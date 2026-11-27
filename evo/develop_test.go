// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"encoding/json"
	"errors"
	"math"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

func testDev() *Topology {
	t := testTopology()
	t.Encoding = Developmental
	t.Wires = [][2]string{{"in", "out"}}
	t.MaxSize = 120
	t.InitDepth = 5
	t.defaults()
	return t
}

// Program builders for the tests.
func dEnd() []Op        { return []Op{{F: dvEnd}} }
func dCut() []Op        { return []Op{{F: dvCut}} }
func dK(v float64) []Op { return []Op{{F: dvConst, V: v}} }
func dJoin(p ...[]Op) []Op {
	var o []Op
	for _, x := range p {
		o = append(o, x...)
	}
	return o
}
func dPart(f uint8, v float64, cc []Op) []Op { return dJoin([]Op{{F: f}}, dK(v), cc) }
func dFn(f uint8, k int, sub ...[]Op) []Op   { return dJoin(append([][]Op{{{F: f, K: k}}}, sub...)...) }

// elements lists the grown elements as "name a b".
func grown(t *testing.T, top *Topology, prog []Op) (*circuit.Circuit, []circuit.Element) {
	t.Helper()
	emb := testEmbryo(t)
	g := Genome{Program: [][]Op{prog}}
	if !top.validProgram(prog) {
		t.Fatalf("invalid program %s", FormatProgram(prog))
	}
	c, _ := top.Decode(emb, g)
	return c, c.Elements[len(emb.Elements):]
}

func ends(c *circuit.Circuit, e circuit.Element) [2]string {
	return [2]string{c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]}
}

func TestDevelopKnownCircuits(t *testing.T) {
	top := testDev()
	type el struct {
		name string
		a, b string
	}
	cases := []struct {
		name string
		prog []Op
		want []el
	}{
		// The wire stays a wire: in and out are joined by a 0 V source.
		{"end", dEnd(), []el{{"V_w1", "in", "out"}}},
		// A cut wire leaves in and out unconnected.
		{"cut", dCut(), nil},
		{"capacitor", dPart(dvC, 0, dEnd()), []el{{"C_g1", "in", "out"}}},
		// SERIES: Z, wire, Z′ — the wire merges the two new nodes.
		{"series", dFn(dvSeries, 0, dPart(dvL, 0, dEnd()), dEnd(), dPart(dvC, 0, dEnd())),
			[]el{{"L_g1", "in", "g_n1"}, {"C_g2", "g_n1", "out"}}},
		// VIA to ground: a T-section, L–L series with C to ground between.
		{"via ground", dFn(dvVia, 0, dPart(dvL, 0, dEnd()), dPart(dvL, 0, dEnd()), dPart(dvC, 0, dEnd())),
			[]el{{"L_g1", "in", "g_n1"}, {"L_g2", "g_n1", "out"}, {"C_g3", "g_n1", "0"}}},
		// PARALLEL with the copy cut: only the original R remains; the two
		// new wires dangle and vanish.
		{"parallel cut", dFn(dvParallel, 0, dPart(dvR, 0, dEnd()), dEnd(), dCut(), dEnd()),
			[]el{{"R_g1", "in", "out"}}},
		// PARALLEL with both kept: R ∥ C between in and out.
		{"parallel", dFn(dvParallel, 0, dPart(dvR, 0, dEnd()), dEnd(), dPart(dvC, 0, dEnd()), dEnd()),
			[]el{{"R_g1", "in", "out"}, {"C_g2", "in", "out"}}},
		// Three wires meeting ground: in and out are shorted to ground.
		{"all to ground", dFn(dvVia, 0, dEnd(), dEnd(), dEnd()),
			[]el{{"V_w1", "0", "in"}, {"V_w2", "0", "out"}}},
		// NOP passes the head on.
		{"nop", dFn(dvNop, 0, dPart(dvR, 0, dEnd())), []el{{"R_g1", "in", "out"}}},
	}
	for _, tc := range cases {
		c, els := grown(t, top, tc.prog)
		var got []el
		for _, e := range els {
			n := ends(c, e)
			got = append(got, el{e.Name, n[0], n[1]})
		}
		if !reflect.DeepEqual(got, tc.want) {
			t.Errorf("%s: %s grew %v, want %v", tc.name, FormatProgram(tc.prog), got, tc.want)
		}
	}
}

// The value subtree sets decades from the middle of the kind's range,
// clamped to the range and snapped to its series.
func TestDevelopValues(t *testing.T) {
	top := testDev() // C range 1e-10..1e-4 (middle 1e-7), E24
	val := func(v []Op) float64 {
		_, els := grown(t, top, dJoin([]Op{{F: dvC}}, v, dEnd()))
		return els[0].Value
	}
	if v := val(dK(0)); v != 1e-7 {
		t.Errorf("middle: %g", v)
	}
	// (+ 1 (− 0.5 0.25)) = 1.25 decades: 1e-7 · 10^1.25 = 1.778e-6 → E24 1.8 µF.
	if v := val(dJoin([]Op{{F: dvAdd}}, dK(1), []Op{{F: dvSub}}, dK(0.5), dK(0.25))); math.Abs(v-1.8e-6) > 1e-18 {
		t.Errorf("(+ 1 (- 0.5 0.25)): %g", v)
	}
	if v := val(dK(40)); v != 1e-4 {
		t.Errorf("clamped high: %g", v)
	}
}

// Heads run breadth-first, and MaxParts caps the parts: of SERIES's three
// components the first two become resistors, the third (last in the
// queue) stays a wire.
func TestDevelopBudget(t *testing.T) {
	top := testDev()
	top.MaxParts = 2
	c, els := grown(t, top, dFn(dvSeries, 0, dPart(dvR, 0, dEnd()), dPart(dvR, 0, dEnd()), dPart(dvR, 0, dEnd())))
	if len(els) != 2 || ends(c, els[0]) != [2]string{"in", "g_n1"} || ends(c, els[1]) != [2]string{"g_n1", "out"} {
		var got [][2]string
		for _, e := range els {
			got = append(got, ends(c, e))
		}
		t.Fatalf("grew %v", got)
	}
}

func TestFormatProgram(t *testing.T) {
	p := dFn(dvVia, 0, dPart(dvL, 0.5, dEnd()), dCut(), dEnd())
	if s := FormatProgram(p); s != "(VIA0 (L 0.5 END) SAFE_CUT END)" {
		t.Errorf("%s", s)
	}
}

// Every operator keeps every genome decodable and well-formed, within
// MaxSize, and the same seed gives the same genomes.
func TestDevelopOperatorInvariants(t *testing.T) {
	top, emb := testDev(), testEmbryo(t)
	top.Wires = [][2]string{{"in", "out"}, {"out", "0"}}
	run := func() []Genome {
		r := rngFor(11, 1)
		pop := []Genome{top.randomGenome(r), top.randomGenome(r)}
		for i := 0; i < 3000; i++ {
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
		}
		return pop
	}
	pop := run()
	grewParts := false
	for _, g := range pop {
		if len(g.Program) != 2 || len(g.Parts) != 0 {
			t.Fatalf("genome shape %d programs, %d parts", len(g.Program), len(g.Parts))
		}
		for _, p := range g.Program {
			if !top.validProgram(p) {
				t.Fatalf("invalid program %s", FormatProgram(p))
			}
		}
		c, vals := top.Decode(emb, g)
		if len(vals) > top.MaxParts {
			t.Fatalf("%d parts > %d", len(vals), top.MaxParts)
		}
		grewParts = grewParts || len(vals) > 0
		_ = c.Validate() // may fail (penalized), must not panic
	}
	if !grewParts {
		t.Error("no genome grew any part")
	}
	if !reflect.DeepEqual(pop, run()) {
		t.Fatal("operators are not deterministic for a fixed seed")
	}
	// The original parents are never modified by the operators.
	r := rngFor(3, 0)
	a := top.randomGenome(r)
	keep := a.Clone()
	c1, _ := top.crossover(r, a, a)
	top.mutate(r, &c1, 0.3)
	if !reflect.DeepEqual(a, keep) {
		t.Fatal("an operator modified its parent")
	}
}

func TestDevelopRepair(t *testing.T) {
	top := testDev()
	g := Genome{Program: [][]Op{
		{{F: dvR}, {F: dvEnd}, {F: dvEnd}}, // construction op in a value position
		dEnd(),                             // one program too many
	}, Parts: []Part{{Kind: 'R', A: 1, B: 2}}}
	top.repair(&g)
	if !reflect.DeepEqual(g.Program, [][]Op{dEnd()}) || g.Parts != nil {
		t.Errorf("repaired to %+v", g)
	}
	for name, p := range map[string][]Op{
		"empty":          nil,
		"unknown op":     {{F: 99}},
		"truncated":      {{F: dvSeries}, {F: dvEnd}},
		"trailing":       {{F: dvEnd}, {F: dvEnd}},
		"NaN constant":   dPart(dvC, math.NaN(), dEnd()),
		"via out of set": dFn(dvVia, 7, dEnd(), dEnd(), dEnd()),
		"too long":       longNops(top.MaxSize),
	} {
		if top.validProgram(p) {
			t.Errorf("%s accepted", name)
		}
	}
	// A kind that is not configured is refused.
	top.Kinds = top.Kinds[:1] // L only
	if top.validProgram(dPart(dvC, 0, dEnd())) {
		t.Error("C accepted with only L configured")
	}
}

// longNops is a valid tree of n NOPs and an END: one op over MaxSize.
func longNops(n int) []Op {
	o := make([]Op, n, n+1)
	for i := range o {
		o[i] = Op{F: dvNop}
	}
	return append(o, Op{F: dvEnd})
}

func TestDevelopValidate(t *testing.T) {
	emb := testEmbryo(t)
	bad := map[string]func(*Topology){
		"no wires":      func(t *Topology) { t.Wires = nil },
		"unknown end":   func(t *Topology) { t.Wires = [][2]string{{"in", "vs"}} },
		"self wire":     func(t *Topology) { t.Wires = [][2]string{{"out", "OUT"}} },
		"bad max_size":  func(t *Topology) { t.MaxSize = -1 },
		"bad encoding":  func(t *Topology) { t.Encoding = 7 },
		"bad initdepth": func(t *Topology) { t.InitDepth = 99 },
	}
	for name, f := range bad {
		top := testDev()
		f(top)
		if err := top.Validate(emb); !errors.Is(err, ErrBadTopology) {
			t.Errorf("%s: %v", name, err)
		}
	}
	if err := testDev().Validate(emb); err != nil {
		t.Fatal(err)
	}
	d, _ := netlist.ParseString("e\nV1 in 0 AC 1\nR1 in g_n42 1k\nR2 g_n42 out 1k\nR3 out 0 1k\n")
	c, _ := circuit.FromDeck(d)
	if err := testDev().Validate(c); !errors.Is(err, ErrBadTopology) {
		t.Errorf("reserved node name: %v", err)
	}
	if _, err := ParseEncoding("lisp"); !errors.Is(err, ErrBadTopology) {
		t.Errorf("ParseEncoding: %v", err)
	}
}

// A developmental genome survives a JSON checkpoint unchanged, and a run
// resumed from generation 2 ends bit-identical to an uninterrupted one.
func TestDevelopResume(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "filter_topology", "config_developmental.json"))
	if err != nil {
		t.Fatal(err)
	}
	g := ex.Problem.Topology.randomGenome(rngFor(1, 0))
	raw, _ := json.Marshal(g)
	var back Genome
	if err := json.Unmarshal(raw, &back); err != nil || !reflect.DeepEqual(back.Program, g.Program) {
		t.Fatalf("JSON round trip: %v", err)
	}
	ex.Options.Population, ex.Options.Generations = 40, 6
	full, err := Run(context.Background(), ex.Problem, ex.Options, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	stop := errors.New("stop")
	path := filepath.Join(t.TempDir(), "cp.json")
	_, err = Run(context.Background(), ex.Problem, ex.Options, nil, func(st *State) error {
		if st.Generation == 2 {
			if err := SaveCheckpoint(path, st); err != nil {
				return err
			}
			return stop
		}
		return nil
	})
	if !errors.Is(err, stop) {
		t.Fatal(err)
	}
	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	resumed, err := Run(context.Background(), ex.Problem, ex.Options, cp, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(full.Best, resumed.Best) || !reflect.DeepEqual(full.Population, resumed.Population) {
		t.Fatal("resumed developmental run differs")
	}
}

// A grown circuit with a short (a wire left between two embryo nodes)
// writes a netlist that parses back to the same elements.
func TestDevelopNetlistRoundTrip(t *testing.T) {
	top := testDev()
	// (VIA0 (L 0 END) END END): L in→m; m joined to out and to ground.
	c, _ := grown(t, top, dFn(dvVia, 0, dPart(dvL, 0, dEnd()), dEnd(), dEnd()))
	var b strings.Builder
	if err := c.WriteNetlist(&b); err != nil {
		t.Fatal(err)
	}
	d, err := netlist.ParseString(b.String())
	if err != nil {
		t.Fatalf("%v\n%s", err, b.String())
	}
	back, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	if len(back.Elements) != len(c.Elements) || back.Index("V_w1") < 0 || back.Index("L_g1") < 0 {
		t.Fatalf("round trip:\n%s", b.String())
	}
}
