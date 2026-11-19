// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"errors"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// ampEmbryo: a 12 V rail on terminal vcc, a 1 V AC source with 1 kΩ
// source resistance on in, a 10 kΩ load on out, the 2N3904 model of
// examples/amplifier_gain, a diode and an NMOS model.
const ampEmbryoSrc = `amplifier embryo
VCC vcc 0 DC 12
VIN vs 0 DC 0 AC 1
RS vs in 1k
RL out 0 10k
.model Q2N3904 NPN(IS=6.734f BF=416.4 VAF=74.03 RB=10 CJC=3.638p CJE=4.493p TF=301.2p)
.model D1N4148 D(IS=2.52n N=1.752 RS=0.568)
.model NM NMOS(LEVEL=1 VTO=1 KP=50u)
.ac dec 10 10 1meg
`

func ampEmbryo(t testing.TB) *circuit.Circuit {
	d, err := netlist.ParseString(ampEmbryoSrc)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// ampTopology: node indices 0 ground, 1 in, 2 out, 3 vcc, 4.. internal.
func ampTopology() *Topology {
	t := &Topology{Terminals: []string{"in", "out", "vcc"}, Internal: 6, MaxParts: 14, InitParts: 4, Kinds: []PartKind{
		{Kind: circuit.Resistor, Min: 10, Max: 1e6, Series: "E24"},
		{Kind: circuit.Capacitor, Min: 1e-9, Max: 1e-4, Series: "E12"},
		{Kind: circuit.BJT, Min: 1, Max: 1, Model: "Q2N3904"},
		{Kind: circuit.Diode, Min: 1, Max: 1, Model: "D1N4148"},
		{Kind: circuit.MOSFET, Min: 1e-6, Max: 1e-3, Model: "NM", Length: 2e-6},
		{Kind: circuit.VoltageSource, Min: 0.1, Max: 15},
		{Kind: circuit.CurrentSource, Min: 1e-6, Max: 0.1},
	}}
	t.defaults()
	return t
}

func TestActivePartsDecode(t *testing.T) {
	top, emb := ampTopology(), ampEmbryo(t)
	if err := top.Validate(emb); err != nil {
		t.Fatal(err)
	}
	if !top.polarized() {
		t.Fatal("sources and transistors must make the topology polarized")
	}
	g := Genome{Parts: []Part{
		{Kind: 'Q', A: 2, B: 0, C: 1, X: 0},                 // collector out, emitter ground, base in
		{Kind: 'R', A: 3, B: 2, X: math.Log10(4700)},        // collector load from vcc
		{Kind: 'M', A: 3, B: 4, C: 1, X: math.Log10(20e-6)}, // drain vcc, source n1, gate in
		{Kind: 'R', A: 4, B: 0, X: 3},
		{Kind: 'V', A: 4, B: 0, X: math.Log10(2)},
		{Kind: 'D', A: 1, B: 0, X: 0},
		{Kind: 'I', A: 3, B: 1, X: -4},
	}}
	c, vals := top.Decode(emb, g)
	if err := c.Validate(); err != nil {
		t.Fatal(err)
	}
	n := len(emb.Elements)
	names := func(ns ...int) []string {
		var out []string
		for _, k := range ns {
			out = append(out, c.NodeNames[k])
		}
		return out
	}
	q := c.Elements[n]
	if q.Kind != circuit.BJT || q.Model != "q2n3904" || q.Area != 1 || len(q.Extra) != 1 ||
		!eqs(names(q.Nodes[0], q.Nodes[1], q.Extra[0]), "out", "0", "in") {
		t.Errorf("Q part: %+v", q)
	}
	m := c.Elements[n+2]
	if m.Kind != circuit.MOSFET || m.Params["w"] != 2e-5 || m.Params["l"] != 2e-6 ||
		!eqs(names(m.Nodes[0], m.Nodes[1], m.Extra[0], m.Extra[1]), "vcc", "g_n1", "in", "g_n1") {
		t.Errorf("M part: %+v (bulk on the source)", m)
	}
	if v := c.Elements[n+4]; v.Kind != circuit.VoltageSource || v.Value != 2 || !v.HasDC {
		t.Errorf("V part: %+v", v)
	}
	if d := c.Elements[n+5]; d.Kind != circuit.Diode || d.Model != "d1n4148" || d.Area != 1 {
		t.Errorf("D part: %+v", d)
	}
	if i := c.Elements[n+6]; i.Kind != circuit.CurrentSource || math.Abs(i.Value-1e-4) > 1e-18 {
		t.Errorf("I part: %+v", i)
	}
	if len(vals) != 7 {
		t.Errorf("values %v", vals)
	}
}

func eqs(got []string, want ...string) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if got[i] != want[i] {
			return false
		}
	}
	return true
}

// A transistor whose base node nothing else touches is pruned: its base
// would float.
func TestActivePrune(t *testing.T) {
	top := ampTopology()
	parts := []Part{{Kind: 'Q', A: 2, B: 0, C: 5}, {Kind: 'R', A: 3, B: 2, X: 3}}
	if got := top.prune(parts); len(got) != 1 || got[0].Kind != 'R' {
		t.Errorf("pruned to %+v", got)
	}
	parts[0].C = 1
	if got := top.prune(parts); len(got) != 2 {
		t.Errorf("a based transistor was pruned: %+v", got)
	}
}

// A hand-built common-emitter stage decoded from parts amplifies: the
// decode wires collector, base and emitter as the part says.
func TestActiveCommonEmitter(t *testing.T) {
	top, emb := ampTopology(), ampEmbryo(t)
	// n1 base, n2 collector, n3 emitter.
	g := Genome{Parts: []Part{
		{Kind: 'C', A: 1, B: 4, X: -5},                // input coupling 10 µF
		{Kind: 'R', A: 3, B: 4, X: math.Log10(100e3)}, // bias divider
		{Kind: 'R', A: 4, B: 0, X: math.Log10(22e3)},
		{Kind: 'R', A: 3, B: 5, X: math.Log10(4.7e3)}, // collector load
		{Kind: 'Q', A: 5, B: 6, C: 4, X: 0},
		{Kind: 'R', A: 6, B: 0, X: math.Log10(1e3)}, // emitter degeneration
		{Kind: 'C', A: 6, B: 0, X: -4},              // bypassed: 100 µF
		{Kind: 'C', A: 5, B: 2, X: -5},              // output coupling
	}}
	c, _ := top.Decode(emb, g)
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 10, FStart: 1e3, FStop: 1e3}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(r, "out", "0")
	gain := cmplx.Abs(h[0])
	if gain < 10 {
		t.Errorf("common-emitter gain %g at 1 kHz, want > 10", gain)
	}
	t.Logf("hand-built CE stage: |V(out)/V(vs)| = %.1f at 1 kHz", gain)
}

// Developmental: (Q value cc_c cc_e cc_b) replaces the component under the
// head by a transistor and three leads.
func TestActiveDevelop(t *testing.T) {
	top := ampTopology()
	top.Encoding = Developmental
	top.Wires = [][2]string{{"in", "out"}}
	emb := ampEmbryo(t)
	if err := top.Validate(emb); err != nil {
		t.Fatal(err)
	}
	q := func(b []Op) []Op { return dJoin([]Op{{F: dvQ}}, dK(0), dEnd(), dEnd(), b) }
	cases := []struct {
		name string
		prog []Op
		want [][]string
	}{
		// All leads stay wires: collector and base on in, emitter on out.
		{"bare", q(dEnd()), [][]string{{"Q_g1", "in", "out", "in"}}},
		// Base lead cut: the base floats and the transistor is pruned.
		{"base cut", q(dCut()), nil},
		// Base lead grown into a resistor: base biased from the collector end.
		{"base resistor", q(dPart(dvR, 0, dEnd())), [][]string{{"Q_g1", "in", "out", "g_n1"}, {"R_g2", "g_n1", "in"}}},
		// FLIP before Q swaps collector and emitter.
		{"flipped", dFn(dvFlip, 0, q(dEnd())), [][]string{{"Q_g1", "out", "in", "out"}}},
		// A source part: V from in (+) to out.
		{"source", dPart(dvV, 0, dEnd()), [][]string{{"V_g1", "in", "out"}}},
	}
	for _, tc := range cases {
		if !top.validProgram(tc.prog) {
			t.Fatalf("%s: invalid %s", tc.name, FormatProgram(tc.prog))
		}
		c, _ := top.Decode(emb, Genome{Program: [][]Op{tc.prog}})
		var got [][]string
		for _, e := range c.Elements[len(emb.Elements):] {
			row := []string{e.Name, c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]}
			if e.Kind == circuit.BJT {
				row = append(row, c.NodeNames[e.Extra[0]])
			}
			got = append(got, row)
		}
		if len(got) != len(tc.want) {
			t.Errorf("%s: %s grew %v, want %v", tc.name, FormatProgram(tc.prog), got, tc.want)
			continue
		}
		for i := range got {
			if !eqs(got[i], tc.want[i]...) {
				t.Errorf("%s: %s grew %v, want %v", tc.name, FormatProgram(tc.prog), got, tc.want)
			}
		}
	}
}

func TestActiveOperatorInvariants(t *testing.T) {
	emb := ampEmbryo(t)
	for _, enc := range []Encoding{Direct, Developmental} {
		top := ampTopology()
		top.Encoding, top.Wires = enc, [][2]string{{"in", "out"}, {"vcc", "out"}}
		r := rngFor(21, 1)
		pop := []Genome{top.randomGenome(r), top.randomGenome(r)}
		seen := map[circuit.Kind]bool{}
		for i := 0; i < 2000; i++ {
			a, b := pop[r.IntN(len(pop))], pop[r.IntN(len(pop))]
			c1, c2 := top.crossover(r, a, b)
			top.mutate(r, &c1, 0.3)
			top.mutate(r, &c2, 0.3)
			top.repair(&c1)
			top.repair(&c2)
			pop = append(pop, c1, c2)
			if len(pop) > 60 {
				pop = pop[2:]
			}
			for _, g := range []Genome{c1, c2} {
				for _, p := range g.Parts {
					if !threeTerminal(p.Kind) && p.C != 0 {
						t.Fatalf("%v: two-terminal part with C = %d", enc, p.C)
					}
				}
				c, _ := top.Decode(emb, g)
				for _, e := range c.Elements[len(emb.Elements):] {
					seen[e.Kind] = true
				}
				_ = c.Validate() // may fail (penalized), must not panic
			}
		}
		for _, k := range []circuit.Kind{circuit.BJT, circuit.MOSFET, circuit.VoltageSource, circuit.CurrentSource, circuit.Diode} {
			if !seen[k] {
				t.Errorf("%v: no %v part ever grown", enc, k)
			}
		}
	}
}

func TestActiveValidate(t *testing.T) {
	emb := ampEmbryo(t)
	bad := map[string]func(*Topology){
		"Q without model":   func(t *Topology) { t.Kinds[2].Model = "" },
		"Q with diode":      func(t *Topology) { t.Kinds[2].Model = "D1N4148" },
		"M with BJT":        func(t *Topology) { t.Kinds[4].Model = "Q2N3904" },
		"R with model":      func(t *Topology) { t.Kinds[0].Model = "Q2N3904" },
		"length on R":       func(t *Topology) { t.Kinds[0].Length = 1e-6 },
		"kappa on a source": func(t *Topology) { t.Kinds[5].KappaA = 0.1 },
		"unknown kind":      func(t *Topology) { t.Kinds[0].Kind = circuit.VCVS },
	}
	for name, f := range bad {
		top := ampTopology()
		f(top)
		if err := top.Validate(emb); !errors.Is(err, ErrBadTopology) {
			t.Errorf("%s: %v", name, err)
		}
	}
}
