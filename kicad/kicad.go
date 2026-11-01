// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package kicad writes a SPICE netlist as a KiCad 8 schematic
// (.kicad_sch, file format 20231120, also read by KiCad 9).
//
// Every element becomes a symbol whose pins end in short wires with net
// labels, so the schematic carries exactly the netlist's connectivity; the
// ground node 0 becomes the power symbol GND, which KiCad's simulator maps
// to node 0. Each symbol carries KiCad's raw-SPICE simulation fields
// (Sim.Device "SPICE", Sim.Params type="…" model="…", as in KiCad's own
// simulation demos), so the element's netlist line — values, waveforms,
// AC specification, controlling-source names, kronspice extensions — is
// passed through unchanged. Dot cards (.model, .tran, .ac, .subckt, …)
// become text directives. The symbols are defined inside the file (library
// "kronspice"), so the schematic opens without extra libraries.
//
// Read checks the result: it recovers the netlist's connectivity from the
// schematic geometry alone (pins, wires, labels, power symbols), which the
// tests compare with the original netlist.
package kicad

import (
	"bufio"
	"bytes"
	"crypto/sha1"
	"errors"
	"fmt"
	"io"
	"math"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/netlist"
)

// ErrNetlist reports a netlist the converter cannot place.
var ErrNetlist = errors.New("kicad: netlist")

// Element is one netlist element as placed: its name, nodes, and the text
// after the nodes (the SPICE "model" part).
type Element struct {
	Name  string
	Kind  byte
	Nodes []string
	Model string
}

// Parsed is a netlist split for the schematic.
type Parsed struct {
	Title      string
	Elements   []Element
	Directives []string // dot cards and .subckt blocks, as written
}

// Split reads a netlist: node counts come from kronspice's parser (so every
// device kind it knows is placed with the right pins), the model text from
// the raw line.
func Split(src string) (*Parsed, error) {
	deck, err := netlist.ParseString(src)
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrNetlist, err)
	}
	nodes := map[string][]string{}
	for _, e := range deck.Elements {
		nodes[strings.ToLower(e.Name)] = e.Nodes
	}
	lines := logical(src)
	p := &Parsed{}
	if len(lines) > 0 {
		p.Title = strings.TrimSpace(lines[0])
	}
	var sub []string
	for _, ln := range lines[1:] {
		t := strings.TrimSpace(ln)
		if t == "" || t[0] == '*' {
			continue
		}
		low := strings.ToLower(t)
		if sub != nil {
			sub = append(sub, t)
			if strings.HasPrefix(low, ".ends") {
				p.Directives = append(p.Directives, strings.Join(sub, "\n"))
				sub = nil
			}
			continue
		}
		if strings.HasPrefix(low, ".subckt") {
			sub = []string{t}
			continue
		}
		if t[0] == '.' {
			if low == ".end" {
				break
			}
			p.Directives = append(p.Directives, t)
			continue
		}
		f := strings.Fields(t)
		name := f[0]
		nd, ok := nodes[strings.ToLower(name)]
		if !ok {
			return nil, fmt.Errorf("%w: element %s not in the parsed deck", ErrNetlist, name)
		}
		k := len(nd)
		switch upper(name[0]) {
		case 'K':
			k = 0 // coupling: no nodes, the rest names the inductors
		case 'E', 'G':
			k = 4 // output pair and controlling pair (the parser keeps the control apart)
		case 'F', 'H':
			k = 2 // output pair; the controlling source name is in the model text
		case 'S':
			k = 4 // switch pair and controlling pair
		case 'W':
			k = 2
		}
		if len(f) < 1+k {
			return nil, fmt.Errorf("%w: %s has %d fields, wants %d nodes", ErrNetlist, name, len(f), k)
		}
		var en []string
		if k > 0 {
			en = append(en, f[1:1+k]...)
		}
		p.Elements = append(p.Elements, Element{Name: name, Kind: upper(name[0]), Nodes: en, Model: strings.Join(f[1+k:], " ")})
	}
	if sub != nil {
		return nil, fmt.Errorf("%w: .subckt without .ends", ErrNetlist)
	}
	return p, nil
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// logical joins continuation lines ("+").
func logical(src string) []string {
	var out []string
	sc := bufio.NewScanner(strings.NewReader(src))
	sc.Buffer(make([]byte, 1<<20), 1<<26)
	for sc.Scan() {
		ln := strings.TrimRight(sc.Text(), "\r")
		if t := strings.TrimSpace(ln); len(out) > 1 && strings.HasPrefix(t, "+") {
			out[len(out)-1] += " " + strings.TrimSpace(t[1:])
			continue
		}
		out = append(out, ln)
	}
	return out
}

// Geometry (mm). KiCad's connection grid is 1.27 mm (50 mil); every
// connection point here is on it.
const (
	grid    = 1.27
	pinLen  = 2.54
	stubLen = 2.54
	cellW   = 25.4
	cellH   = 25.4
)

type pinDef struct {
	x, y  float64 // symbol coordinates (y up)
	angle int     // direction from the connection point into the body
}

// shape gives the symbol name and pin layout for an element kind and pin
// count: two-pin parts vertical (pin 1 top), others a box with pins on the
// left (first half) and right.
func shape(kind byte, n int) (string, []pinDef) {
	if n == 0 {
		return "K0", nil
	}
	if n == 2 {
		name := "TWO"
		switch kind {
		case 'R':
			name = "R"
		case 'C':
			name = "C"
		case 'L':
			name = "L"
		case 'V', 'I':
			name = string(kind)
		case 'D':
			name = "D"
		}
		return name, []pinDef{{0, 3.81, 270}, {0, -3.81, 90}}
	}
	left := (n + 1) / 2
	right := n - left
	rows := left
	if right > rows {
		rows = right
	}
	var pins []pinDef
	top := float64(rows-1) * 2.54 / 2
	for i := 0; i < left; i++ {
		pins = append(pins, pinDef{-7.62, snap(top - float64(i)*2.54), 0})
	}
	for i := 0; i < right; i++ {
		pins = append(pins, pinDef{7.62, snap(top - float64(i)*2.54), 180})
	}
	return fmt.Sprintf("BOX%d", n), pins
}

func snap(v float64) float64 { return math.Round(v/grid) * grid }

// uuid derives a stable UUID from text, so the same netlist always gives
// the same file (determinism, CLAUDE.md §1 goal 5).
func uuid(parts ...string) string {
	h := sha1.Sum([]byte(strings.Join(parts, "\x00")))
	h[6] = (h[6] & 0x0f) | 0x50
	h[8] = (h[8] & 0x3f) | 0x80
	return fmt.Sprintf("%x-%x-%x-%x-%x", h[0:4], h[4:6], h[6:8], h[8:10], h[10:16])
}

func q(s string) string {
	return `"` + strings.NewReplacer(`\`, `\\`, `"`, `\"`, "\n", `\n`).Replace(s) + `"`
}

func f(v float64) string {
	s := fmt.Sprintf("%.4f", v)
	s = strings.TrimRight(strings.TrimRight(s, "0"), ".")
	if s == "-0" {
		s = "0"
	}
	return s
}

const font = `(effects (font (size 1.27 1.27)))`
const fontHide = `(effects (font (size 1.27 1.27)) (hide yes))`

// libSymbol writes one symbol definition of the "kronspice" library.
func libSymbol(w io.Writer, name string, pins []pinDef) {
	fmt.Fprintf(w, "\t\t(symbol %s\n", q("kronspice:"+name))
	fmt.Fprintf(w, "\t\t\t(pin_numbers hide) (pin_names (offset 0)) (exclude_from_sim no) (in_bom yes) (on_board yes)\n")
	fmt.Fprintf(w, "\t\t\t(property \"Reference\" \"X\" (at 3.81 1.27 0) %s)\n", font)
	fmt.Fprintf(w, "\t\t\t(property \"Value\" %s (at 3.81 -1.27 0) %s)\n", q(name), font)
	fmt.Fprintf(w, "\t\t\t(property \"Footprint\" \"\" (at 0 0 0) %s)\n", fontHide)
	fmt.Fprintf(w, "\t\t\t(property \"Datasheet\" \"\" (at 0 0 0) %s)\n", fontHide)
	fmt.Fprintf(w, "\t\t\t(symbol %s\n", q(name+"_0_1"))
	stroke := `(stroke (width 0.254) (type default)) (fill (type none))`
	switch {
	case name == "R":
		fmt.Fprintf(w, "\t\t\t\t(rectangle (start -1.016 -2.54) (end 1.016 2.54) %s)\n", stroke)
	case name == "C":
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -2.032 0.762) (xy 2.032 0.762)) %s)\n", stroke)
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -2.032 -0.762) (xy 2.032 -0.762)) %s)\n", stroke)
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy 0 2.54) (xy 0 0.762)) %s)\n", stroke)
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy 0 -2.54) (xy 0 -0.762)) %s)\n", stroke)
	case name == "L":
		for i := 0; i < 4; i++ {
			y0 := 2.54 - float64(i)*1.27
			fmt.Fprintf(w, "\t\t\t\t(arc (start 0 %s) (mid 0.635 %s) (end 0 %s) %s)\n", f(y0), f(y0-0.635), f(y0-1.27), stroke)
		}
	case name == "V" || name == "I":
		fmt.Fprintf(w, "\t\t\t\t(circle (center 0 0) (radius 2.54) %s)\n", stroke)
		if name == "V" {
			fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -0.635 1.27) (xy 0.635 1.27)) %s)\n", stroke)
			fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy 0 1.905) (xy 0 0.635)) %s)\n", stroke)
			fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -0.635 -1.27) (xy 0.635 -1.27)) %s)\n", stroke)
		} else {
			fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy 0 -1.524) (xy 0 1.524) (xy -0.635 0.762)) %s)\n", stroke)
			fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy 0 1.524) (xy 0.635 0.762)) %s)\n", stroke)
		}
	case name == "D":
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -1.27 1.27) (xy 1.27 1.27) (xy 0 -1.27) (xy -1.27 1.27)) %s)\n", stroke)
		fmt.Fprintf(w, "\t\t\t\t(polyline (pts (xy -1.27 -1.27) (xy 1.27 -1.27)) %s)\n", stroke)
	case name == "TWO":
		fmt.Fprintf(w, "\t\t\t\t(rectangle (start -1.905 -2.54) (end 1.905 2.54) %s)\n", stroke)
	case name == "K0":
		fmt.Fprintf(w, "\t\t\t\t(rectangle (start -2.54 -1.27) (end 2.54 1.27) %s)\n", stroke)
	default:
		top, bot := 0.0, 0.0
		for _, p := range pins {
			top, bot = math.Max(top, p.y), math.Min(bot, p.y)
		}
		fmt.Fprintf(w, "\t\t\t\t(rectangle (start -5.08 %s) (end 5.08 %s) %s)\n", f(bot-1.27), f(top+1.27), stroke)
	}
	fmt.Fprintf(w, "\t\t\t)\n\t\t\t(symbol %s\n", q(name+"_1_1"))
	// Pins reach from the connection point to the body: 1.27 mm on the
	// two-pin parts (body ±2.54), 2.54 mm on the boxes (body ±5.08).
	length := pinLen
	if len(pins) == 2 {
		length = 1.27
	}
	for i, p := range pins {
		fmt.Fprintf(w, "\t\t\t\t(pin passive line (at %s %s %d) (length %s) (name %s %s) (number %s %s))\n",
			f(p.x), f(p.y), p.angle, f(length), q("~"), font, q(fmt.Sprint(i+1)), font)
	}
	fmt.Fprintf(w, "\t\t\t)\n\t\t)\n")
}

const gndSymbol = `		(symbol "power:GND"
			(power) (pin_names (offset 0)) (exclude_from_sim no) (in_bom yes) (on_board yes)
			(property "Reference" "#PWR" (at 0 -6.35 0) (effects (font (size 1.27 1.27)) (hide yes)))
			(property "Value" "GND" (at 0 -3.81 0) (effects (font (size 1.27 1.27))))
			(property "Footprint" "" (at 0 0 0) (effects (font (size 1.27 1.27)) (hide yes)))
			(property "Datasheet" "" (at 0 0 0) (effects (font (size 1.27 1.27)) (hide yes)))
			(symbol "GND_0_1"
				(polyline (pts (xy 0 0) (xy 0 -1.27) (xy 1.27 -1.27) (xy 0 -2.54) (xy -1.27 -1.27) (xy 0 -1.27)) (stroke (width 0) (type default)) (fill (type none)))
			)
			(symbol "GND_1_1"
				(pin power_in line (at 0 0 270) (length 0) hide (name "GND" (effects (font (size 1.27 1.27)))) (number "1" (effects (font (size 1.27 1.27)))))
			)
		)
`

// Options for Write.
type Options struct {
	Project string // KiCad project name (instances); default "kronspice"
	Columns int    // symbols per row in label mode; default 8
	// Labels forces label-only connectivity (short stubs with net labels,
	// no drawn wires). It is used anyway above MaxWired elements.
	Labels bool
}

// MaxWired is the largest netlist drawn with wires; larger ones (field
// networks with thousands of junctions) use net labels only.
const MaxWired = 1000

// sheet writes schematic items with stable UUIDs.
type sheet struct {
	bw            *bytes.Buffer
	root, project string
	pwr, n        int
}

func (s *sheet) id(parts ...string) string {
	s.n++
	return uuid(append([]string{s.root, fmt.Sprint(s.n)}, parts...)...)
}

func (s *sheet) wire(x1, y1, x2, y2 float64) {
	if x1 == x2 && y1 == y2 {
		return
	}
	fmt.Fprintf(s.bw, "\t(wire (pts (xy %s %s) (xy %s %s)) (stroke (width 0) (type default)) (uuid %s))\n",
		f(x1), f(y1), f(x2), f(y2), q(s.id("wire")))
}

func (s *sheet) junction(x, y float64) {
	fmt.Fprintf(s.bw, "\t(junction (at %s %s) (diameter 0) (color 0 0 0 0) (uuid %s))\n", f(x), f(y), q(s.id("junction")))
}

func (s *sheet) label(name string, x, y float64, just string) {
	fmt.Fprintf(s.bw, "\t(label %s (at %s %s 0) (effects (font (size 1.27 1.27)) (justify %s)) (uuid %s))\n",
		q(name), f(x), f(y), just, q(s.id("label", name)))
}

func (s *sheet) gnd(x, y float64) {
	s.pwr++
	gid := s.id("gnd")
	ref := fmt.Sprintf("#PWR%04d", s.pwr)
	fmt.Fprintf(s.bw, "\t(symbol (lib_id \"power:GND\") (at %s %s 0) (unit 1) (exclude_from_sim no) (in_bom yes) (on_board yes) (dnp no)\n", f(x), f(y))
	fmt.Fprintf(s.bw, "\t\t(uuid %s)\n\t\t(property \"Reference\" %s (at %s %s 0) %s)\n", q(gid), q(ref), f(x), f(y+6.35), fontHide)
	fmt.Fprintf(s.bw, "\t\t(property \"Value\" \"GND\" (at %s %s 0) %s)\n", f(x), f(y+3.81), fontHide)
	fmt.Fprintf(s.bw, "\t\t(property \"Footprint\" \"\" (at %s %s 0) %s)\n\t\t(property \"Datasheet\" \"\" (at %s %s 0) %s)\n", f(x), f(y), fontHide, f(x), f(y), fontHide)
	fmt.Fprintf(s.bw, "\t\t(pin \"1\" (uuid %s))\n", q(uuid("gndpin", gid)))
	fmt.Fprintf(s.bw, "\t\t(instances (project %s (path %s (reference %s) (unit 1))))\n\t)\n", q(s.project), q("/"+s.root), q(ref))
}

// symbol places element e; ref and value text at the given points.
func (s *sheet) symbol(e Element, sym string, npins int, x, y float64, rx, ry, vx, vy float64, just string) {
	id := uuid("sym", s.root, e.Name)
	fmt.Fprintf(s.bw, "\t(symbol (lib_id %s) (at %s %s 0) (unit 1) (exclude_from_sim no) (in_bom yes) (on_board yes) (dnp no)\n", q("kronspice:"+sym), f(x), f(y))
	fmt.Fprintf(s.bw, "\t\t(uuid %s)\n", q(id))
	fmt.Fprintf(s.bw, "\t\t(property \"Reference\" %s (at %s %s 0) (effects (font (size 1.27 1.27)) (justify %s)))\n", q(e.Name), f(rx), f(ry), just)
	fmt.Fprintf(s.bw, "\t\t(property \"Value\" %s (at %s %s 0) (effects (font (size 1.27 1.27)) (justify %s)))\n", q(shortValue(e.Model)), f(vx), f(vy), just)
	fmt.Fprintf(s.bw, "\t\t(property \"Footprint\" \"\" (at %s %s 0) %s)\n", f(x), f(y), fontHide)
	fmt.Fprintf(s.bw, "\t\t(property \"Datasheet\" \"\" (at %s %s 0) %s)\n", f(x), f(y), fontHide)
	fmt.Fprintf(s.bw, "\t\t(property \"Sim.Device\" \"SPICE\" (at %s %s 0) %s)\n", f(x), f(y), fontHide)
	fmt.Fprintf(s.bw, "\t\t(property \"Sim.Params\" %s (at %s %s 0) %s)\n",
		q(fmt.Sprintf("type=%q model=%q lib=%q", string(e.Kind), e.Model, "")), f(x), f(y), fontHide)
	var sp []string
	for k := 0; k < npins; k++ {
		sp = append(sp, fmt.Sprintf("%d=%d", k+1, k+1))
	}
	fmt.Fprintf(s.bw, "\t\t(property \"Sim.Pins\" %s (at %s %s 0) %s)\n", q(strings.Join(sp, " ")), f(x), f(y), fontHide)
	for k := 0; k < npins; k++ {
		fmt.Fprintf(s.bw, "\t\t(pin %s (uuid %s))\n", q(fmt.Sprint(k+1)), q(uuid("pin", id, fmt.Sprint(k))))
	}
	fmt.Fprintf(s.bw, "\t\t(instances (project %s (path %s (reference %s) (unit 1))))\n\t)\n", q(s.project), q("/"+s.root), q(e.Name))
}

func shortValue(v string) string {
	if r := []rune(v); len(r) > 40 {
		return string(r[:37]) + "…"
	}
	return v
}

func isGround(n string) bool { return n == "0" || strings.EqualFold(n, "gnd") }

// textW estimates the width of a 1.27 mm text, rounded up to the grid.
func textW(t string) float64 { return math.Ceil((1.1*float64(len([]rune(t)))+1.27)/grid) * grid }

// Write converts the netlist src to a KiCad schematic. Up to MaxWired
// elements it is drawn with wires: each row of parts has one horizontal
// rail per net above it (the net's label at the rail's left end, GND for
// ground), every pin is wired to its net's rail, with a junction at each
// tap; a two-pin part's grounded lower pin gets a GND symbol directly.
// Rails of the same net in different rows are joined by their labels.
// Larger netlists, or Options.Labels, use short stubs with net labels.
func Write(w io.Writer, src string, opt Options) error {
	p, err := Split(src)
	if err != nil {
		return err
	}
	if opt.Project == "" {
		opt.Project = "kronspice"
	}
	if opt.Columns <= 0 {
		opt.Columns = 8
	}
	root := uuid("root", src)
	s := &sheet{bw: &bytes.Buffer{}, root: root, project: opt.Project}
	var width, height float64
	mode := "drawn wires to one rail per net"
	if opt.Labels || len(p.Elements) > MaxWired {
		width, height = s.labels(p, opt)
		mode = "net labels"
	} else {
		width, height = s.wired(p)
	}
	bw := bufio.NewWriter(w)
	fmt.Fprintf(bw, "(kicad_sch\n\t(version 20231120)\n\t(generator \"kronspice\")\n\t(generator_version \"1.0\")\n\t(uuid %s)\n", q(root))
	fmt.Fprintf(bw, "\t(paper \"User\" %s %s)\n", f(math.Max(297, math.Ceil(width/10)*10)), f(math.Max(210, math.Ceil(height/10)*10)))
	fmt.Fprintf(bw, "\t(title_block (title %s) (comment 1 %s))\n", q(p.Title), q("generated by kronspice kicad.Write from the netlist; connectivity by "+mode))
	fmt.Fprintf(bw, "\t(lib_symbols\n")
	used := map[string][]pinDef{}
	var names []string
	for _, e := range p.Elements {
		n, pins := shape(e.Kind, len(e.Nodes))
		if _, ok := used[n]; !ok {
			used[n] = pins
			names = append(names, n)
		}
	}
	sort.Strings(names)
	for _, n := range names {
		libSymbol(bw, n, used[n])
	}
	bw.WriteString(gndSymbol)
	fmt.Fprintf(bw, "\t)\n")
	bw.Write(s.bw.Bytes())
	fmt.Fprintf(bw, "\t(sheet_instances (path \"/\" (page \"1\")))\n)\n")
	return bw.Flush()
}

// directives writes the dot cards as text from y down; returns the end.
func (s *sheet) directives(p *Parsed, x, y float64) float64 {
	for i, d := range p.Directives {
		fmt.Fprintf(s.bw, "\t(text %s (exclude_from_sim no) (at %s %s 0) (effects (font (size 1.27 1.27)) (justify left bottom)) (uuid %s))\n",
			q(d), f(x), f(y), q(uuid("text", s.root, fmt.Sprint(i))))
		y = snap(y + 5.08 + 2.54*float64(strings.Count(d, "\n")))
	}
	return y
}

// labels is the label-only layout: a grid of cells, each pin ending in a
// short stub with its net label (or GND).
func (s *sheet) labels(p *Parsed, opt Options) (float64, float64) {
	maxX := 0.0
	for i, e := range p.Elements {
		sym, pins := shape(e.Kind, len(e.Nodes))
		col, row := i%opt.Columns, i/opt.Columns
		x := snap(25.4 + float64(col)*cellW*map[bool]float64{true: 1.6, false: 1}[len(pins) > 2])
		y := snap(38.1 + float64(row)*cellH)
		if len(pins) > 2 {
			x = snap(25.4 + float64(col)*cellW)
		}
		maxX = math.Max(maxX, x+cellW)
		s.symbol(e, sym, len(pins), x, y, x+3.81, y-1.27, x+3.81, y+1.27, "left")
		for k, pd := range pins {
			px, py := x+pd.x, y-pd.y
			dx, dy := 0.0, 0.0
			switch pd.angle {
			case 270:
				dy = -stubLen
			case 90:
				dy = stubLen
			case 0:
				dx = -stubLen
			case 180:
				dx = stubLen
			}
			ex, ey := snap(px+dx), snap(py+dy)
			s.wire(px, py, ex, ey)
			if isGround(e.Nodes[k]) {
				s.gnd(ex, ey)
				continue
			}
			just := "left bottom"
			if dx < 0 {
				just = "right bottom"
			}
			s.label(e.Nodes[k], ex, ey, just)
		}
	}
	rows := (len(p.Elements) + opt.Columns - 1) / opt.Columns
	y := s.directives(p, 25.4, snap(38.1+float64(rows)*cellH+12.7))
	return math.Max(maxX+25.4, 300), y + 25.4
}

// wired is the drawn layout.
func (s *sheet) wired(p *Parsed) (float64, float64) {
	const rowW = 380.0 // mm of parts per row before wrapping
	type place struct {
		e        Element
		sym      string
		pins     []pinDef
		w, off   float64 // slot width, symbol x offset in the slot
		top, bot float64 // extent above and below the symbol origin
		tw       float64
		nl, nr   int // box pins left and right
	}
	var places []place
	for _, e := range p.Elements {
		sym, pins := shape(e.Kind, len(e.Nodes))
		pl := place{e: e, sym: sym, pins: pins, tw: math.Max(textW(e.Name), textW(shortValue(e.Model)))}
		switch {
		case len(pins) == 0:
			pl.off, pl.w, pl.top, pl.bot = 5.08, pl.tw+10.16, 2.54, 2.54
		case len(pins) == 2:
			pl.off = pl.tw + 2.54
			pl.w = pl.off + 10.16
			pl.top, pl.bot = 3.81, 3.81+2.54+2.54
		default:
			pl.nl = (len(pins) + 1) / 2
			pl.nr = len(pins) - pl.nl
			left := 7.62 + 2.54*float64(pl.nl)
			right := math.Max(7.62+2.54*float64(pl.nr), pl.tw-5.08)
			pl.off = left + 2.54
			pl.w = pl.off + right + 5.08
			for _, pd := range pins {
				pl.top = math.Max(pl.top, pd.y+1.27)
				pl.bot = math.Max(pl.bot, -pd.y+1.27)
			}
			pl.bot += 6.35
		}
		pl.off, pl.w = snap(pl.off), snap(pl.w)
		places = append(places, pl)
	}
	// Rows of parts.
	var rows [][]place
	acc := 0.0
	for _, pl := range places {
		if len(rows) == 0 || acc+pl.w > rowW {
			rows = append(rows, nil)
			acc = 0
		}
		rows[len(rows)-1] = append(rows[len(rows)-1], pl)
		acc += pl.w
	}
	y0, maxX := 25.4, 0.0
	for _, row := range rows {
		// Nets on this row's rails, in order of first use; ground last
		// (only for grounded pins that cannot take a GND symbol directly).
		var nets []string
		seen := map[string]bool{}
		groundRail := false
		labelW := 0.0
		for _, pl := range row {
			for k, n := range pl.e.Nodes {
				if isGround(n) {
					if !(len(pl.pins) == 2 && k == 1) {
						groundRail = true
					}
					continue
				}
				if !seen[n] {
					seen[n] = true
					nets = append(nets, n)
					labelW = math.Max(labelW, textW(n))
				}
			}
		}
		if groundRail {
			nets = append(nets, "0")
		}
		railY := map[string]float64{}
		for r, n := range nets {
			railY[n] = snap(y0 + 2.54*float64(r))
		}
		railX := snap(12.7 + labelW)
		top, bot := 0.0, 0.0
		for _, pl := range row {
			top, bot = math.Max(top, pl.top), math.Max(bot, pl.bot)
		}
		railBottom := y0 + 2.54*float64(max(len(nets)-1, 0))
		yc := snap(railBottom + 7.62 + top)
		taps := map[string][]float64{}
		xs := railX + 5.08
		for _, pl := range row {
			x := snap(xs + pl.off)
			xs += pl.w
			e := pl.e
			switch {
			case len(pl.pins) == 0:
				s.symbol(e, pl.sym, 0, x, yc, x+3.81, yc-1.27, x+3.81, yc+1.27, "left")
			case len(pl.pins) == 2:
				s.symbol(e, pl.sym, 2, x, yc, x-2.54, yc-1.27, x-2.54, yc+1.27, "right")
				// Pin 1 (top) straight up to its rail.
				n1 := e.Nodes[0]
				if isGround(n1) {
					n1 = "0"
				}
				s.wire(x, yc-3.81, x, railY[n1])
				taps[n1] = append(taps[n1], x)
				// Pin 2 (bottom): GND below, or round the part to its rail.
				py := yc + 3.81
				if isGround(e.Nodes[1]) {
					s.wire(x, py, x, py+2.54)
					s.gnd(x, py+2.54)
					break
				}
				n2 := e.Nodes[1]
				lx := x + 5.08
				s.wire(x, py, x, py+2.54)
				s.wire(x, py+2.54, lx, py+2.54)
				s.wire(lx, py+2.54, lx, railY[n2])
				taps[n2] = append(taps[n2], lx)
			default:
				bottom := yc
				for _, pd := range pl.pins {
					bottom = math.Max(bottom, yc-pd.y+1.27)
				}
				s.symbol(e, pl.sym, len(pl.pins), x, yc, x-5.08, bottom+3.81, x-5.08, bottom+6.35, "left")
				for k, pd := range pl.pins {
					n := e.Nodes[k]
					if isGround(n) {
						n = "0"
					}
					px, py := x+pd.x, yc-pd.y
					// Lower pins take lanes further out, so no wire crosses
					// another pin's lead.
					var lx float64
					if k < pl.nl {
						lx = px - 2.54*float64(k+1)
					} else {
						lx = px + 2.54*float64(k-pl.nl+1)
					}
					s.wire(px, py, lx, py)
					s.wire(lx, py, lx, railY[n])
					taps[n] = append(taps[n], lx)
				}
			}
			maxX = math.Max(maxX, xs)
		}
		for _, n := range nets {
			ts := taps[n]
			end := railX
			for _, t := range ts {
				end = math.Max(end, t)
			}
			s.wire(railX, railY[n], end, railY[n])
			for _, t := range ts {
				if t < end {
					s.junction(t, railY[n])
				}
			}
			if n == "0" {
				s.gnd(railX, railY[n])
			} else {
				s.label(n, railX, railY[n], "right bottom")
			}
		}
		y0 = snap(yc + bot + 12.7)
	}
	y := s.directives(p, 25.4, y0)
	return maxX + 25.4, y + 25.4
}
