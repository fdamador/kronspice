// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package netlist reads SPICE netlists.
//
// Supported so far: the title line, comments (* and ;), continuation lines
// (+), R, C and L elements (C and L with IC=), independent V and I sources
// with DC, AC and PULSE/SIN/EXP/PWL specifications, and the .op, .tran,
// .title and .end control lines. Unsupported elements and dot commands are
// reported as errors with their line number rather than silently ignored.
package netlist

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"regexp"
	"strings"
)

// ErrSyntax is returned for malformed or unsupported netlist lines.
var ErrSyntax = errors.New("netlist: syntax error")

// maxLineLength bounds a single logical line so hostile input cannot exhaust memory.
const maxLineLength = 1 << 20

// ElementLine is one parsed element card.
type ElementLine struct {
	Name  string   // full name as written, e.g. "R1"
	Kind  byte     // upper-case first letter: 'R', 'C', 'L', 'V', 'I'
	Nodes []string // node names, lower-cased
	Value float64  // SI value: ohms, farads, henries, or a source's DC value
	HasDC bool     // a source's DC value was given

	WaveKind string    // "pulse", "sin", "exp", "pwl" or "" for none
	WaveArgs []float64 // waveform arguments in SI units

	ACMag, ACPhase float64 // AC specification (magnitude, phase in degrees)

	IC    float64 // IC= value for capacitors (volts) and inductors (amperes)
	HasIC bool

	// Params are instance parameters given as name=value (MOSFET l, w, ad,
	// as, pd, ps, nrd, nrs, m), lower-cased.
	Params map[string]float64

	// Couples are the two inductors a K element couples (lower-cased).
	Couples []string
	// Ext is the exterior channel named by EXT= on a resistor line (where
	// its dissipated power goes; "" = thermal), or on a capacitor or
	// inductor line (where the power of its SCALE variation goes;
	// "" = mechanical), lower-cased.
	Ext string
	// ScaleKind and ScaleArgs are a capacitor's or inductor's SCALE=
	// waveform ("pulse", "sin", "exp", "pwl"; "" for none): the time
	// variation s(t) of its charge or flux (kronspice extension).
	ScaleKind string
	ScaleArgs []float64

	Model string  // model name for semiconductor elements (lower-cased)
	Area  float64 // area factor (diodes, BJTs); 0 means the default of 1
	// Nodes are in SPICE order: two for most elements; c b e [s] for Q.
	// E and G use CtrlNodes, F and H CtrlSource; Value is their gain.

	// Behavioral (B) elements, ngspice syntax restricted to a piecewise-
	// linear function of one controlling quantity:
	//   Bname n+ n- I=pwl(V(a[,b]), x1, y1, x2, y2, ...)
	//   Bname n+ n- V=pwl(I(Vname), x1, y1, ...)
	// (and the other two combinations).
	BOut       byte     // 'I' (current output) or 'V' (voltage output)
	BCtrl      byte     // 'V' (node voltage) or 'I' (voltage-source current)
	CtrlNodes  []string // controlling nodes for V(a,b); b = "0" for V(a)
	CtrlSource string   // controlling voltage source for I(Vname), lower-cased
	Curve      []float64
	// A capacitor written with Q=pwl(v1 q1 ...) or an inductor with
	// FLUX=pwl(i1 λ1 ...) keeps its charge or flux curve in Curve.

	Line int // first physical line number (1-based)
}

// Analysis is one analysis request.
type Analysis struct {
	Kind string // "op", "tran", "ac" or "dc"
	// Transient: print step, stop time, start of output, maximum internal step.
	Step, Stop, Start, MaxStep float64
	UIC                        bool // skip the operating point; use IC= values
	// AC: sweep type ("dec", "oct", "lin"), points (per decade/octave, or
	// total for lin), start and stop frequency [Hz].
	Variation     string
	Points        int
	FStart, FStop float64
	// DC: up to two nested source sweeps (Sweeps[0] is the inner one).
	Sweeps []DCSweep
	Line   int
}

// DCSweep is one source swept by .dc.
type DCSweep struct {
	Source            string
	Start, Stop, Step float64
}

// ModelCard is one .model line: ".model name type (param=value …)".
type ModelCard struct {
	Name   string             // lower-cased
	Type   string             // lower-cased device type, e.g. "d"
	Params map[string]float64 // lower-cased parameter names, SI values
	Line   int
}

// Deck is a parsed netlist.
type Deck struct {
	Title    string
	Elements []ElementLine
	Analyses []Analysis
	Models   []ModelCard
	// Domains maps node names (lower-cased) to the energy domain they
	// belong to, from ".domain <channel> node …" lines (kronspice extension);
	// nodes not listed are electrical.
	Domains map[string]string
}

type logicalLine struct {
	text string
	line int
}

// Parse reads a netlist. The first line is always the title, as in SPICE.
func Parse(r io.Reader) (*Deck, error) {
	lines, err := readLogical(r)
	if err != nil {
		return nil, err
	}
	d := &Deck{}
	if len(lines) == 0 {
		return d, nil
	}
	d.Title = strings.TrimSpace(lines[0].text)
	for _, ll := range lines[1:] {
		text := stripComment(ll.text)
		fields := strings.Fields(text)
		if len(fields) == 0 {
			continue
		}
		head := strings.ToLower(fields[0])
		if strings.HasPrefix(head, ".") {
			done, err := d.control(head, fields, ll)
			if err != nil {
				return nil, err
			}
			if done {
				break
			}
			continue
		}
		el, err := parseElement(text, ll.line)
		if err != nil {
			return nil, err
		}
		d.Elements = append(d.Elements, el)
	}
	return d, nil
}

// ParseString is a convenience wrapper around Parse.
func ParseString(s string) (*Deck, error) { return Parse(strings.NewReader(s)) }

func (d *Deck) control(head string, fields []string, ll logicalLine) (end bool, err error) {
	switch head {
	case ".end":
		return true, nil
	case ".op":
		d.Analyses = append(d.Analyses, Analysis{Kind: "op", Line: ll.line})
	case ".tran":
		a, err := parseTran(fields[1:], ll.line)
		if err != nil {
			return false, err
		}
		d.Analyses = append(d.Analyses, a)
	case ".ac":
		a, err := parseAC(fields[1:], ll.line)
		if err != nil {
			return false, err
		}
		d.Analyses = append(d.Analyses, a)
	case ".dc":
		a, err := parseDC(fields[1:], ll.line)
		if err != nil {
			return false, err
		}
		d.Analyses = append(d.Analyses, a)
	case ".model":
		m, err := parseModel(stripComment(ll.text), ll.line)
		if err != nil {
			return false, err
		}
		d.Models = append(d.Models, m)
	case ".title":
		d.Title = strings.TrimSpace(strings.Join(fields[1:], " "))
	case ".domain":
		// .domain mechanical|acoustic|… node node …   (kronspice extension)
		if len(fields) < 3 {
			return false, fmt.Errorf("%w: line %d: .domain needs: domain node [node …]", ErrSyntax, ll.line)
		}
		dom := strings.ToLower(fields[1])
		if d.Domains == nil {
			d.Domains = map[string]string{}
		}
		for _, n := range fields[2:] {
			n = strings.ToLower(n)
			if prev, ok := d.Domains[n]; ok && prev != dom {
				return false, fmt.Errorf("%w: line %d: node %s declared %s and %s", ErrSyntax, ll.line, n, prev, dom)
			}
			d.Domains[n] = dom
		}
	default:
		return false, fmt.Errorf("%w: line %d: unsupported control line %q", ErrSyntax, ll.line, fields[0])
	}
	return false, nil
}

// parseTran reads ".tran tstep tstop [tstart [tmax]] [uic]".
func parseTran(args []string, line int) (Analysis, error) {
	a := Analysis{Kind: "tran", Line: line}
	var nums []float64
	for _, f := range args {
		if strings.EqualFold(f, "uic") {
			a.UIC = true
			continue
		}
		v, err := ParseValue(f)
		if err != nil {
			return a, fmt.Errorf("line %d: .tran: %w", line, err)
		}
		nums = append(nums, v)
	}
	if len(nums) < 2 || len(nums) > 4 {
		return a, fmt.Errorf("%w: line %d: .tran needs: tstep tstop [tstart [tmax]] [uic]", ErrSyntax, line)
	}
	a.Step, a.Stop = nums[0], nums[1]
	if len(nums) > 2 {
		a.Start = nums[2]
	}
	if len(nums) > 3 {
		a.MaxStep = nums[3]
	}
	if a.Step <= 0 || a.Stop <= 0 || a.Start < 0 || a.Start >= a.Stop || a.MaxStep < 0 {
		return a, fmt.Errorf("%w: line %d: .tran needs 0 < tstep, 0 ≤ tstart < tstop, tmax ≥ 0", ErrSyntax, line)
	}
	return a, nil
}

// MaxSweepPoints bounds AC and DC sweeps so hostile input cannot exhaust memory.
const MaxSweepPoints = 1_000_000

// parseAC reads ".ac dec|oct|lin n fstart fstop".
func parseAC(args []string, line int) (Analysis, error) {
	a := Analysis{Kind: "ac", Line: line}
	if len(args) != 4 {
		return a, fmt.Errorf("%w: line %d: .ac needs: dec|oct|lin points fstart fstop", ErrSyntax, line)
	}
	a.Variation = strings.ToLower(args[0])
	if a.Variation != "dec" && a.Variation != "oct" && a.Variation != "lin" {
		return a, fmt.Errorf("%w: line %d: .ac sweep must be dec, oct or lin", ErrSyntax, line)
	}
	n, err := ParseValue(args[1])
	if err != nil {
		return a, fmt.Errorf("line %d: .ac points: %w", line, err)
	}
	if n < 1 || n > MaxSweepPoints || n != float64(int(n)) {
		return a, fmt.Errorf("%w: line %d: .ac points must be a whole number ≥ 1", ErrSyntax, line)
	}
	a.Points = int(n)
	if a.FStart, err = ParseValue(args[2]); err != nil {
		return a, fmt.Errorf("line %d: .ac fstart: %w", line, err)
	}
	if a.FStop, err = ParseValue(args[3]); err != nil {
		return a, fmt.Errorf("line %d: .ac fstop: %w", line, err)
	}
	if a.FStart <= 0 && a.Variation != "lin" || a.FStart < 0 || a.FStop < a.FStart {
		return a, fmt.Errorf("%w: line %d: .ac needs 0 < fstart ≤ fstop (fstart may be 0 for lin)", ErrSyntax, line)
	}
	return a, nil
}

// parseDC reads ".dc src start stop step [src2 start2 stop2 step2]".
func parseDC(args []string, line int) (Analysis, error) {
	a := Analysis{Kind: "dc", Line: line}
	if len(args) != 4 && len(args) != 8 {
		return a, fmt.Errorf("%w: line %d: .dc needs: source start stop step [source2 start2 stop2 step2]", ErrSyntax, line)
	}
	for i := 0; i < len(args); i += 4 {
		sw := DCSweep{Source: args[i]}
		vals := make([]float64, 3)
		for j := range vals {
			v, err := ParseValue(args[i+1+j])
			if err != nil {
				return a, fmt.Errorf("line %d: .dc: %w", line, err)
			}
			vals[j] = v
		}
		sw.Start, sw.Stop, sw.Step = vals[0], vals[1], vals[2]
		if sw.Step == 0 || (sw.Stop-sw.Start)/sw.Step < 0 || (sw.Stop-sw.Start)/sw.Step > MaxSweepPoints {
			return a, fmt.Errorf("%w: line %d: .dc step must be non-zero and move from start toward stop", ErrSyntax, line)
		}
		a.Sweeps = append(a.Sweeps, sw)
	}
	return a, nil
}

// MaxModelParams bounds the parameters on one .model line.
const MaxModelParams = 200

// parseModel reads ".model name type [(] p1=v1 p2=v2 … [)]".
func parseModel(text string, line int) (ModelCard, error) {
	f := tokenize(text)
	m := ModelCard{Line: line, Params: map[string]float64{}}
	if len(f) < 3 {
		return m, fmt.Errorf("%w: line %d: .model needs: name type (params)", ErrSyntax, line)
	}
	m.Name, m.Type = strings.ToLower(f[1]), strings.ToLower(f[2])
	rest := f[3:]
	if len(rest)%3 != 0 || len(rest)/3 > MaxModelParams {
		return m, fmt.Errorf("%w: line %d: .model parameters must be name=value", ErrSyntax, line)
	}
	for i := 0; i < len(rest); i += 3 {
		if rest[i+1] != "=" {
			return m, fmt.Errorf("%w: line %d: .model parameter %q needs =value", ErrSyntax, line, rest[i])
		}
		v, err := ParseValue(rest[i+2])
		if err != nil {
			return m, fmt.Errorf("line %d: .model %s: %w", line, rest[i], err)
		}
		key := strings.ToLower(rest[i])
		if _, dup := m.Params[key]; dup {
			return m, fmt.Errorf("%w: line %d: .model parameter %s given twice", ErrSyntax, line, rest[i])
		}
		m.Params[key] = v
	}
	return m, nil
}

// tokenize splits an element line into fields, treating parentheses, commas
// and "=" as separators (so "PULSE(0 1 0)", "PULSE 0 1 0" and "IC=2" all work).
func tokenize(s string) []string {
	r := strings.NewReplacer("(", " ", ")", " ", ",", " ", "=", " = ")
	return strings.Fields(r.Replace(s))
}

func parseElement(text string, line int) (ElementLine, error) {
	f := tokenize(text)
	if len(f) == 0 {
		return ElementLine{Line: line}, fmt.Errorf("%w: line %d: no element name", ErrSyntax, line)
	}
	name := f[0]
	kind := upper(name[0])
	el := ElementLine{Name: name, Kind: kind, Line: line}
	switch kind {
	case 'R', 'C', 'L':
		if len(f) < 4 {
			return el, fmt.Errorf("%w: line %d: %c element needs: name n1 n2 value", ErrSyntax, line, kind)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		if err := parsePassiveSpec(&el, f[3:], line); err != nil {
			return el, err
		}
	case 'D':
		// Dname n+ n- model [area]
		if len(f) < 4 || len(f) > 5 {
			return el, fmt.Errorf("%w: line %d: diode needs: Dname anode cathode model [area]", ErrSyntax, line)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		el.Model = strings.ToLower(f[3])
		if len(f) == 5 {
			a, err := ParseValue(f[4])
			if err != nil {
				return el, fmt.Errorf("line %d: diode area: %w", line, err)
			}
			el.Area = a
		}
	case 'V', 'I':
		if len(f) < 4 {
			return el, fmt.Errorf("%w: line %d: source needs: %cname n+ n- [DC] value | waveform", ErrSyntax, line, kind)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		if err := parseSourceSpec(&el, f[3:], line); err != nil {
			return el, err
		}
	case 'B':
		if err := parseBehavioral(&el, text, line); err != nil {
			return el, err
		}
	case 'E', 'G':
		// Ename n+ n- nc+ nc- gain
		if len(f) != 6 {
			return el, fmt.Errorf("%w: line %d: %c element needs: name n+ n- nc+ nc- gain", ErrSyntax, line, kind)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		el.CtrlNodes = []string{strings.ToLower(f[3]), strings.ToLower(f[4])}
		v, err := ParseValue(f[5])
		if err != nil {
			return el, fmt.Errorf("line %d: gain: %w", line, err)
		}
		el.Value = v
	case 'Y':
		// Yname p1+ p1- p2+ p2- g   (ideal gyrator / transducer, kronspice
		// extension): i1 = g·V2, i2 = −g·V1.
		if len(f) != 6 {
			return el, fmt.Errorf("%w: line %d: Y element needs: name p1+ p1- p2+ p2- gyration", ErrSyntax, line)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2]), strings.ToLower(f[3]), strings.ToLower(f[4])}
		v, err := ParseValue(f[5])
		if err != nil {
			return el, fmt.Errorf("line %d: gyration: %w", line, err)
		}
		el.Value = v
	case 'P':
		// Pname a b m+ m- R KAPPA=κ W0=ω0 [VF=v] [THERMAL=1 ETAMAX=η]   (ball-bearing
		// motor, kronspice extension; coefficients from a measurement fit)
		if len(f) < 6 {
			return el, fmt.Errorf("%w: line %d: P element needs: name a b m+ m- R KAPPA=.. W0=.. [VF=..] [THERMAL=1 ETAMAX=..]", ErrSyntax, line)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2]), strings.ToLower(f[3]), strings.ToLower(f[4])}
		v, err := ParseValue(f[5])
		if err != nil {
			return el, fmt.Errorf("line %d: resistance: %w", line, err)
		}
		el.Value = v
		rest := f[6:]
		if len(rest)%3 != 0 || len(rest)/3 > 5 {
			return el, fmt.Errorf("%w: line %d: bearing parameters must be name=value", ErrSyntax, line)
		}
		for i := 0; i < len(rest); i += 3 {
			if rest[i+1] != "=" {
				return el, fmt.Errorf("%w: line %d: bearing parameter %q needs =value", ErrSyntax, line, rest[i])
			}
			pv, err := ParseValue(rest[i+2])
			if err != nil {
				return el, fmt.Errorf("line %d: %s: %w", line, rest[i], err)
			}
			if el.Params == nil {
				el.Params = map[string]float64{}
			}
			k := strings.ToLower(rest[i])
			if _, dup := el.Params[k]; dup {
				return el, fmt.Errorf("%w: line %d: %s given twice", ErrSyntax, line, rest[i])
			}
			el.Params[k] = pv
		}
	case 'N':
		// Nname w1+ w1- [w2+ w2- …] m+ m- model   (machine: windings with
		// θ-dependent inductances and a mechanical port, kronspice extension)
		if len(f) < 6 || len(f)%2 != 0 {
			return el, fmt.Errorf("%w: line %d: N element needs: name w1+ w1- [w2+ w2- …] port+ port- model", ErrSyntax, line)
		}
		for _, n := range f[1 : len(f)-1] {
			el.Nodes = append(el.Nodes, strings.ToLower(n))
		}
		el.Model = strings.ToLower(f[len(f)-1])
	case 'K':
		// Kname Lname1 Lname2 k   (mutual inductance, coupling coefficient)
		if len(f) != 4 {
			return el, fmt.Errorf("%w: line %d: K element needs: Kname L1 L2 coupling", ErrSyntax, line)
		}
		el.Couples = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		v, err := ParseValue(f[3])
		if err != nil {
			return el, fmt.Errorf("line %d: coupling: %w", line, err)
		}
		el.Value = v
	case 'F', 'H':
		// Fname n+ n- vname gain
		if len(f) != 5 {
			return el, fmt.Errorf("%w: line %d: %c element needs: name n+ n- vname gain", ErrSyntax, line, kind)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2])}
		el.CtrlSource = strings.ToLower(f[3])
		v, err := ParseValue(f[4])
		if err != nil {
			return el, fmt.Errorf("line %d: gain: %w", line, err)
		}
		el.Value = v
	case 'M':
		// Mname d g s b model [L=val] [W=val] [AD=…] [AS=…] [PD=…] [PS=…]
		// [NRD=…] [NRS=…] [M=…]   (SPICE3f3 User's Manual §2.4.10)
		if len(f) < 6 {
			return el, fmt.Errorf("%w: line %d: MOSFET needs: Mname d g s b model [L=..] [W=..] ...", ErrSyntax, line)
		}
		el.Nodes = []string{strings.ToLower(f[1]), strings.ToLower(f[2]), strings.ToLower(f[3]), strings.ToLower(f[4])}
		el.Model = strings.ToLower(f[5])
		rest := f[6:]
		if len(rest)%3 != 0 || len(rest)/3 > 16 {
			return el, fmt.Errorf("%w: line %d: MOSFET parameters must be name=value", ErrSyntax, line)
		}
		for i := 0; i < len(rest); i += 3 {
			if rest[i+1] != "=" {
				return el, fmt.Errorf("%w: line %d: MOSFET parameter %q needs =value", ErrSyntax, line, rest[i])
			}
			v, err := ParseValue(rest[i+2])
			if err != nil {
				return el, fmt.Errorf("line %d: %s: %w", line, rest[i], err)
			}
			if el.Params == nil {
				el.Params = map[string]float64{}
			}
			k := strings.ToLower(rest[i])
			if _, dup := el.Params[k]; dup {
				return el, fmt.Errorf("%w: line %d: %s given twice", ErrSyntax, line, rest[i])
			}
			el.Params[k] = v
		}
	case 'Q':
		// Qname c b e [s] model [area]. With five fields after the name the
		// last is an area when it is a number, else a model after a
		// substrate node.
		t := f[1:]
		if len(t) < 4 || len(t) > 6 {
			return el, fmt.Errorf("%w: line %d: BJT needs: Qname c b e [s] model [area]", ErrSyntax, line)
		}
		var area string
		switch len(t) {
		case 4:
			el.Nodes, el.Model = t[:3], t[3]
		case 5:
			if _, err := ParseValue(t[4]); err == nil {
				el.Nodes, el.Model, area = t[:3], t[3], t[4]
			} else {
				el.Nodes, el.Model = t[:4], t[4]
			}
		case 6:
			el.Nodes, el.Model, area = t[:4], t[4], t[5]
		}
		for i := range el.Nodes {
			el.Nodes[i] = strings.ToLower(el.Nodes[i])
		}
		el.Nodes = append([]string(nil), el.Nodes...)
		el.Model = strings.ToLower(el.Model)
		if area != "" {
			a, err := ParseValue(area)
			if err != nil {
				return el, fmt.Errorf("line %d: BJT area: %w", line, err)
			}
			el.Area = a
		}
	default:
		return el, fmt.Errorf("%w: line %d: unsupported element %q", ErrSyntax, line, name)
	}
	return el, nil
}

// parsePassiveSpec reads what follows the nodes of an R, C or L line:
//
//	R: value [EXT=channel]
//	C: value | Q=pwl(v1 q1 v2 q2 ...)     [SCALE=wave(...)] [EXT=channel] [IC=v]
//	L: value | FLUX=pwl(i1 λ1 i2 λ2 ...)  [SCALE=wave(...)] [EXT=channel] [IC=i]
//
// EXT= on a resistor names where its dissipated power goes; Q=, FLUX=,
// SCALE= and EXT= on a capacitor or inductor (kronspice extensions) make it a
// nonlinear, time-varying element whose charge or flux is s(t)·curve.
func parsePassiveSpec(el *ElementLine, f []string, line int) error {
	kind := el.Kind
	keyAt := func(i int) string {
		if i+1 < len(f) && f[i+1] == "=" {
			return strings.ToLower(f[i])
		}
		return ""
	}
	// numbers reads values from f[i:] up to the next "key =" or the end.
	numbers := func(i int, what string) ([]float64, int, error) {
		var out []float64
		for ; i < len(f) && keyAt(i) == ""; i++ {
			v, err := ParseValue(f[i])
			if err != nil {
				return nil, i, fmt.Errorf("line %d: %s: %w", line, what, err)
			}
			if len(out) >= 2*MaxCurvePoints {
				return nil, i, fmt.Errorf("%w: line %d: more than %d %s points", ErrSyntax, line, MaxCurvePoints, what)
			}
			out = append(out, v)
		}
		return out, i, nil
	}
	curveKey := map[byte]string{'C': "q", 'L': "flux"}[kind]
	i := 0
	if k := keyAt(0); k != "" && k == curveKey {
		if len(f) < 3 || !strings.EqualFold(f[2], "pwl") {
			return fmt.Errorf("%w: line %d: %s= needs pwl(x1 y1 x2 y2 ...)", ErrSyntax, line, strings.ToUpper(k))
		}
		pts, next, err := numbers(3, "pwl")
		if err != nil {
			return err
		}
		if len(pts) < 4 || len(pts)%2 != 0 {
			return fmt.Errorf("%w: line %d: pwl needs at least two (x, y) pairs", ErrSyntax, line)
		}
		el.Curve, i = pts, next
	} else {
		v, err := ParseValue(f[0])
		if err != nil {
			return fmt.Errorf("line %d: %w", line, err)
		}
		el.Value, i = v, 1
	}
	seen := map[string]bool{}
	for i < len(f) {
		k := keyAt(i)
		ok := k == "ext" || (kind != 'R' && (k == "ic" || k == "scale"))
		if !ok || seen[k] || i+2 >= len(f) {
			return fmt.Errorf("%w: line %d: unexpected %q after value", ErrSyntax, line, strings.Join(f[i:], " "))
		}
		seen[k] = true
		arg := f[i+2]
		switch k {
		case "ext":
			el.Ext = strings.ToLower(arg)
			i += 3
		case "ic":
			ic, err := ParseValue(arg)
			if err != nil {
				return fmt.Errorf("line %d: IC: %w", line, err)
			}
			el.IC, el.HasIC = ic, true
			i += 3
		case "scale":
			if !isWaveName(arg) {
				return fmt.Errorf("%w: line %d: SCALE needs PULSE, SIN, EXP or PWL", ErrSyntax, line)
			}
			args, next, err := numbers(i+3, "SCALE")
			if err != nil {
				return err
			}
			el.ScaleKind, el.ScaleArgs, i = strings.ToLower(arg), args, next
		}
	}
	return nil
}

// MaxCurvePoints bounds the points of a behavioral PWL curve.
const MaxCurvePoints = 100000

var behavioralRE = regexp.MustCompile(`(?is)^\s*(\S+)\s+(\S+)\s+(\S+)\s+([IV])\s*=\s*pwl\s*\(\s*([VI])\s*\(([^()]*)\)\s*(.*)\)\s*$`)

// parseBehavioral reads the supported subset of ngspice B sources: a
// piecewise-linear function of one node voltage or one voltage-source
// current (ngspice manual, "Bxxxx: Nonlinear dependent source", pwl).
func parseBehavioral(el *ElementLine, text string, line int) error {
	m := behavioralRE.FindStringSubmatch(text)
	if m == nil {
		return fmt.Errorf("%w: line %d: B element needs: Bname n+ n- I=pwl(V(a,b), x1,y1, ...) or V=pwl(I(Vname), ...)", ErrSyntax, line)
	}
	el.Nodes = []string{strings.ToLower(m[2]), strings.ToLower(m[3])}
	el.BOut = upper(m[4][0])
	el.BCtrl = upper(m[5][0])
	args := strings.FieldsFunc(m[6], func(r rune) bool { return r == ',' || r == ' ' || r == '\t' })
	switch el.BCtrl {
	case 'V':
		if len(args) < 1 || len(args) > 2 {
			return fmt.Errorf("%w: line %d: V() needs one or two nodes", ErrSyntax, line)
		}
		el.CtrlNodes = []string{strings.ToLower(args[0]), "0"}
		if len(args) == 2 {
			el.CtrlNodes[1] = strings.ToLower(args[1])
		}
	case 'I':
		if len(args) != 1 || upper(args[0][0]) != 'V' {
			return fmt.Errorf("%w: line %d: I() needs one voltage source name", ErrSyntax, line)
		}
		el.CtrlSource = strings.ToLower(args[0])
	}
	pts := strings.FieldsFunc(m[7], func(r rune) bool { return r == ',' || r == ' ' || r == '\t' })
	if len(pts) > 2*MaxCurvePoints {
		return fmt.Errorf("%w: line %d: more than %d curve points", ErrSyntax, line, MaxCurvePoints)
	}
	for _, p := range pts {
		v, err := ParseValue(p)
		if err != nil {
			return fmt.Errorf("line %d: pwl: %w", line, err)
		}
		el.Curve = append(el.Curve, v)
	}
	if len(el.Curve) < 4 || len(el.Curve)%2 != 0 {
		return fmt.Errorf("%w: line %d: pwl needs at least two (x, y) pairs", ErrSyntax, line)
	}
	return nil
}

func isWaveName(s string) bool {
	switch strings.ToLower(s) {
	case "pulse", "sin", "exp", "pwl":
		return true
	}
	return false
}

// parseSourceSpec reads "[DC] v [AC mag [phase]] [PULSE|SIN|EXP|PWL args...]"
// in any order, as SPICE3 does.
func parseSourceSpec(el *ElementLine, f []string, line int) error {
	i := 0
	num := func() (float64, bool, error) {
		if i >= len(f) {
			return 0, false, nil
		}
		low := strings.ToLower(f[i])
		if low == "dc" || low == "ac" || isWaveName(low) {
			return 0, false, nil
		}
		v, err := ParseValue(f[i])
		if err != nil {
			return 0, false, fmt.Errorf("line %d: %w", line, err)
		}
		i++
		return v, true, nil
	}
	for i < len(f) {
		low := strings.ToLower(f[i])
		switch {
		case low == "dc":
			i++
			v, ok, err := num()
			if err != nil {
				return err
			}
			if !ok {
				return fmt.Errorf("%w: line %d: DC needs a value", ErrSyntax, line)
			}
			el.Value, el.HasDC = v, true
		case low == "ac":
			i++
			mag, ok, err := num()
			if err != nil {
				return err
			}
			if !ok {
				mag = 1
			}
			el.ACMag = mag
			if ph, ok, err := num(); err != nil {
				return err
			} else if ok {
				el.ACPhase = ph
			}
		case isWaveName(low):
			if el.WaveKind != "" {
				return fmt.Errorf("%w: line %d: more than one waveform", ErrSyntax, line)
			}
			el.WaveKind = low
			i++
			for {
				v, ok, err := num()
				if err != nil {
					return err
				}
				if !ok {
					break
				}
				el.WaveArgs = append(el.WaveArgs, v)
			}
		default:
			// A bare leading number is the DC value.
			if el.HasDC {
				return fmt.Errorf("%w: line %d: unexpected %q", ErrSyntax, line, f[i])
			}
			v, err := ParseValue(f[i])
			if err != nil {
				return fmt.Errorf("line %d: %w", line, err)
			}
			el.Value, el.HasDC = v, true
			i++
		}
	}
	if !el.HasDC && el.WaveKind == "" && el.ACMag == 0 {
		return fmt.Errorf("%w: line %d: source %s has no value", ErrSyntax, line, el.Name)
	}
	return nil
}

func readLogical(r io.Reader) ([]logicalLine, error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 64*1024), maxLineLength)
	var out []logicalLine
	n := 0
	for sc.Scan() {
		n++
		text := strings.TrimRight(sc.Text(), "\r")
		trimmed := strings.TrimLeft(text, " \t")
		if n > 1 && strings.HasPrefix(trimmed, "+") && len(out) > 1 {
			last := &out[len(out)-1]
			if len(last.text)+len(trimmed) > maxLineLength {
				return nil, fmt.Errorf("%w: line %d: logical line too long", ErrSyntax, n)
			}
			last.text += " " + trimmed[1:]
			continue
		}
		if n > 1 && strings.HasPrefix(trimmed, "*") {
			continue // full-line comment
		}
		out = append(out, logicalLine{text: text, line: n})
	}
	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrSyntax, err)
	}
	return out, nil
}

func stripComment(s string) string {
	if i := strings.IndexByte(s, ';'); i >= 0 {
		return s[:i]
	}
	return s
}

func upper(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 'a' + 'A'
	}
	return b
}
