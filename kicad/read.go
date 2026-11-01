// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kicad

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// node is an S-expression: an atom (Atom set) or a list.
type node struct {
	Atom string
	List []*node
	Str  bool // atom was a quoted string
}

func (n *node) head() string {
	if len(n.List) > 0 && n.List[0].List == nil {
		return n.List[0].Atom
	}
	return ""
}

func (n *node) find(h string) *node {
	for _, c := range n.List {
		if c.head() == h {
			return c
		}
	}
	return nil
}

func (n *node) all(h string) []*node {
	var out []*node
	for _, c := range n.List {
		if c.head() == h {
			out = append(out, c)
		}
	}
	return out
}

func (n *node) num(i int) float64 {
	if i >= len(n.List) {
		return 0
	}
	v, _ := strconv.ParseFloat(n.List[i].Atom, 64)
	return v
}

// parse reads one S-expression.
func parse(s string) (*node, error) {
	p := 0
	var rd func() (*node, error)
	skip := func() {
		for p < len(s) && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' || s[p] == '\r') {
			p++
		}
	}
	rd = func() (*node, error) {
		skip()
		if p >= len(s) {
			return nil, fmt.Errorf("unexpected end")
		}
		switch s[p] {
		case '(':
			p++
			n := &node{List: []*node{}}
			for {
				skip()
				if p >= len(s) {
					return nil, fmt.Errorf("unclosed list")
				}
				if s[p] == ')' {
					p++
					return n, nil
				}
				c, err := rd()
				if err != nil {
					return nil, err
				}
				n.List = append(n.List, c)
			}
		case ')':
			return nil, fmt.Errorf("unexpected ) at %d", p)
		case '"':
			p++
			var b strings.Builder
			for p < len(s) && s[p] != '"' {
				if s[p] == '\\' && p+1 < len(s) {
					p++
					switch s[p] {
					case 'n':
						b.WriteByte('\n')
					default:
						b.WriteByte(s[p])
					}
				} else {
					b.WriteByte(s[p])
				}
				p++
			}
			if p >= len(s) {
				return nil, fmt.Errorf("unclosed string")
			}
			p++
			return &node{Atom: b.String(), Str: true}, nil
		default:
			st := p
			for p < len(s) && !strings.ContainsRune(" \t\n\r()", rune(s[p])) {
				p++
			}
			return &node{Atom: s[st:p]}, nil
		}
	}
	n, err := rd()
	if err != nil {
		return nil, err
	}
	skip()
	if p != len(s) {
		return nil, fmt.Errorf("trailing text at %d", p)
	}
	return n, nil
}

type pt struct{ x, y int64 } // 1 µm units

func key(x, y float64) pt { return pt{int64(math.Round(x * 1000)), int64(math.Round(y * 1000))} }

// Read recovers, from a schematic written by Write, each element's name
// and the net of each of its pins (ground as "0"), using only geometry:
// pin positions from the library symbols, wires, labels and power
// symbols. It also reports dangling pins (connected to nothing) and pins
// or labels off the 1.27 mm grid.
func Read(sch string) (map[string][]string, error) {
	root, err := parse(sch)
	if err != nil {
		return nil, fmt.Errorf("kicad: s-expression: %v", err)
	}
	if root.head() != "kicad_sch" {
		return nil, fmt.Errorf("kicad: not a schematic")
	}
	lib := map[string][]pt{} // lib_id → pin offsets in symbol coordinates (µm), by pin index
	libPins := map[string][][2]float64{}
	if ls := root.find("lib_symbols"); ls != nil {
		for _, s := range ls.all("symbol") {
			id := s.List[1].Atom
			var pins [][2]float64
			for _, u := range s.all("symbol") {
				for _, p := range u.all("pin") {
					at := p.find("at")
					num := p.find("number")
					k, _ := strconv.Atoi(num.List[1].Atom)
					for len(pins) < k {
						pins = append(pins, [2]float64{})
					}
					pins[k-1] = [2]float64{at.num(1), at.num(2)}
				}
			}
			libPins[id] = pins
			lib[id] = nil
		}
	}
	// Union-find over points.
	parent := map[pt]pt{}
	var findp func(pt) pt
	findp = func(a pt) pt {
		if _, ok := parent[a]; !ok {
			parent[a] = a
		}
		for parent[a] != a {
			parent[a] = parent[parent[a]]
			a = parent[a]
		}
		return a
	}
	union := func(a, b pt) { parent[findp(a)] = findp(b) }
	onGrid := func(x, y float64) bool {
		r := func(v float64) bool { return math.Abs(v/grid-math.Round(v/grid)) < 1e-6 }
		return r(x) && r(y)
	}
	// Wires join their two ends. A connection point that lies inside
	// another wire (a T, with or without a junction dot) joins that wire;
	// wires that merely cross do not connect, as in KiCad.
	type seg struct{ a, b pt }
	vert := map[int64][]seg{} // by x
	horz := map[int64][]seg{} // by y
	var points []pt
	for _, w := range root.all("wire") {
		pts := w.find("pts").all("xy")
		if len(pts) != 2 {
			return nil, fmt.Errorf("kicad: wire with %d points", len(pts))
		}
		a, b := key(pts[0].num(1), pts[0].num(2)), key(pts[1].num(1), pts[1].num(2))
		union(a, b)
		points = append(points, a, b)
		switch {
		case a.x == b.x:
			if a.y > b.y {
				a, b = b, a
			}
			vert[a.x] = append(vert[a.x], seg{a, b})
		case a.y == b.y:
			if a.x > b.x {
				a, b = b, a
			}
			horz[a.y] = append(horz[a.y], seg{a, b})
		default:
			return nil, fmt.Errorf("kicad: diagonal wire")
		}
	}
	tee := func(k pt) {
		for _, sg := range vert[k.x] {
			if sg.a.y < k.y && k.y < sg.b.y {
				union(k, sg.a)
			}
		}
		for _, sg := range horz[k.y] {
			if sg.a.x < k.x && k.x < sg.b.x {
				union(k, sg.a)
			}
		}
	}
	for _, j := range root.all("junction") {
		at := j.find("at")
		points = append(points, key(at.num(1), at.num(2)))
	}
	name := map[pt]string{}
	for _, l := range root.all("label") {
		at := l.find("at")
		if !onGrid(at.num(1), at.num(2)) {
			return nil, fmt.Errorf("kicad: label %s off grid", l.List[1].Atom)
		}
		k := key(at.num(1), at.num(2))
		findp(k)
		points = append(points, k)
		name[k] = l.List[1].Atom
	}
	type inst struct {
		ref  string
		pins []pt
	}
	var parts []inst
	for _, s := range root.all("symbol") {
		id := s.find("lib_id").List[1].Atom
		at := s.find("at")
		x, y := at.num(1), at.num(2)
		if at.num(3) != 0 {
			return nil, fmt.Errorf("kicad: rotated symbols are not read")
		}
		ref := ""
		for _, p := range s.all("property") {
			if p.List[1].Atom == "Reference" {
				ref = p.List[2].Atom
			}
		}
		if id == "power:GND" {
			k := key(x, y)
			findp(k)
			points = append(points, k)
			name[k] = "0"
			continue
		}
		pins, ok := libPins[id]
		if !ok {
			return nil, fmt.Errorf("kicad: %s uses undefined symbol %s", ref, id)
		}
		var ps []pt
		for _, pd := range pins {
			px, py := x+pd[0], y-pd[1]
			if !onGrid(px, py) {
				return nil, fmt.Errorf("kicad: pin of %s off grid at (%g, %g)", ref, px, py)
			}
			k := key(px, py)
			findp(k)
			points = append(points, k)
			ps = append(ps, k)
		}
		parts = append(parts, inst{ref, ps})
	}
	for _, k := range points {
		tee(k)
	}
	// Net name per connected group: its label or ground; two different
	// names in one group are a short.
	group := map[pt]string{}
	for k, n := range name {
		r := findp(k)
		if old, ok := group[r]; ok && old != n {
			return nil, fmt.Errorf("kicad: nets %q and %q are joined", old, n)
		}
		group[r] = n
	}
	out := map[string][]string{}
	for _, p := range parts {
		var nets []string
		for _, k := range p.pins {
			n, ok := group[findp(k)]
			if !ok {
				return nil, fmt.Errorf("kicad: %s has a pin connected to no label", p.ref)
			}
			nets = append(nets, n)
		}
		out[p.ref] = nets
	}
	return out, nil
}
