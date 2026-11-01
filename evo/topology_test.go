// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

func testEmbryo(t testing.TB) *circuit.Circuit {
	d, err := netlist.ParseString("embryo\nVIN vs 0 DC 0 AC 1\nRS vs in 1k\nRL out 0 1k\n.ac dec 20 10 100k\n")
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func testTopology() *Topology {
	t := &Topology{Terminals: []string{"in", "out"}, Internal: 4, MaxParts: 12, InitParts: 3, Kinds: []PartKind{
		{Kind: circuit.Inductor, Min: 1e-5, Max: 10, Series: "E24"},
		{Kind: circuit.Capacitor, Min: 1e-10, Max: 1e-4, Series: "E24"},
		{Kind: circuit.Resistor, Min: 1, Max: 1e6, Series: "E24"},
	}}
	t.defaults()
	return t
}

// Node indices: 0 ground, 1 in, 2 out, 3.. internal g_n1...
func TestTopologyPrune(t *testing.T) {
	top := testTopology()
	parts := []Part{
		{Kind: 'R', A: 1, B: 2, X: 3},  // in–out: kept
		{Kind: 'C', A: 2, B: 2, X: -7}, // self-loop: dropped
		{Kind: 'C', A: 1, B: 3, X: -7}, // in–n1 ...
		{Kind: 'L', A: 3, B: 4, X: -1}, // ... n1–n2, n2 dangles: the chain goes
		{Kind: 'R', A: 2, B: 5, X: 3},  // out–n3 ...
		{Kind: 'C', A: 5, B: 0, X: -7}, // ... n3–ground: a real branch, kept
	}
	got := top.prune(parts)
	want := []Part{parts[0], parts[4], parts[5]}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("pruned to %+v, want %+v", got, want)
	}
}

func TestTopologyDecode(t *testing.T) {
	top, emb := testTopology(), testEmbryo(t)
	before := emb.Clone()
	g := Genome{Parts: []Part{{Kind: 'L', A: 1, B: 2, X: math.Log10(0.21)}, {Kind: 'C', A: 2, B: 3, X: -7}, {Kind: 'R', A: 3, B: 0, X: 3}}}
	c, vals := top.Decode(emb, g)
	if !reflect.DeepEqual(emb, before) {
		t.Fatal("Decode modified the embryo")
	}
	if len(c.Elements) != len(emb.Elements)+3 {
		t.Fatalf("%d elements", len(c.Elements))
	}
	e := c.Elements[len(emb.Elements)]
	if e.Name != "L_g1" || c.NodeNames[e.Nodes[0]] != "in" || c.NodeNames[e.Nodes[1]] != "out" || e.Value != 0.22 {
		t.Fatalf("first part %+v (E24 snap of 0.21 H is 0.22 H)", e)
	}
	if c.NodeNames[c.Elements[len(emb.Elements)+1].Nodes[1]] != internalName(1) || vals[2] != 1000 {
		t.Fatalf("internal node or value: %v %v", c.NodeNames, vals)
	}
	if err := c.Validate(); err != nil {
		t.Fatal(err)
	}
}

// Every operator keeps a genome decodable: values in range, nodes in
// range, kinds from the configuration, size within MaxParts; and the same
// seed gives the same genomes.
func TestTopologyOperatorInvariants(t *testing.T) {
	top, emb := testTopology(), testEmbryo(t)
	run := func() []Genome {
		r := rngFor(7, 1)
		pop := []Genome{top.randomGenome(r), top.randomGenome(r)}
		for i := 0; i < 3000; i++ {
			a, b := pop[r.IntN(len(pop))], pop[r.IntN(len(pop))]
			c1, c2 := top.crossover(r, a, b)
			top.mutate(r, &c1, 0.2)
			top.mutate(r, &c2, 0.2)
			top.repair(&c1)
			top.repair(&c2)
			pop = append(pop, c1, c2)
			if len(pop) > 50 {
				pop = pop[2:]
			}
		}
		return pop
	}
	pop := run()
	for _, g := range pop {
		if len(g.Parts) > top.MaxParts {
			t.Fatalf("%d parts > %d", len(g.Parts), top.MaxParts)
		}
		for _, p := range g.Parts {
			pk, ok := top.kind(p.Kind)
			if !ok || p.A < 0 || p.B < 0 || p.A >= top.nodes() || p.B >= top.nodes() {
				t.Fatalf("bad part %+v", p)
			}
			if v := top.value(p); v < pk.Min || v > pk.Max {
				t.Fatalf("value %g outside [%g, %g]", v, pk.Min, pk.Max)
			}
		}
		c, _ := top.Decode(emb, g)
		_ = c.Validate() // may fail (floating islands are penalized), must not panic
	}
	if !reflect.DeepEqual(pop, run()) {
		t.Fatal("operators are not deterministic for a fixed seed")
	}
}

func TestTopologyValidate(t *testing.T) {
	emb := testEmbryo(t)
	bad := map[string]func(*Topology){
		"unknown terminal": func(t *Topology) { t.Terminals = []string{"nowhere"} },
		"ground terminal":  func(t *Topology) { t.Terminals = []string{"0"} },
		"diode part":       func(t *Topology) { t.Kinds[0].Kind = circuit.Diode },
		"duplicate kind":   func(t *Topology) { t.Kinds[1].Kind = circuit.Inductor },
		"bad range":        func(t *Topology) { t.Kinds[0].Min = 0 },
		"bad series":       func(t *Topology) { t.Kinds[0].Series = "E7" },
	}
	for name, f := range bad {
		top := testTopology()
		f(top)
		if err := top.Validate(emb); !errors.Is(err, ErrBadTopology) {
			t.Errorf("%s: %v", name, err)
		}
	}
	// An embryo that already uses an internal node name.
	d, _ := netlist.ParseString("e\nV1 in 0 AC 1\nR1 in g_n1 1k\nR2 g_n1 out 1k\nR3 out 0 1k\n")
	c, _ := circuit.FromDeck(d)
	if err := testTopology().Validate(c); !errors.Is(err, ErrBadTopology) {
		t.Errorf("name collision: %v", err)
	}
}

// A first-order RC low-pass inside its mask scores 0; outside, at least
// the violation floor.
func TestFilterMask(t *testing.T) {
	c, err := circuit.NewBuilder("rc").Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, HasDC: true, ACMag: 1}, "in", "0").
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1.59e-7).Build() // fc ≈ 1 kHz
	if err != nil {
		t.Fatal(err)
	}
	spec := analysis.ACSpec{Variation: "dec", Points: 20, FStart: 10, FStop: 1e5}
	eval := func(b []Band) float64 {
		v, err := FilterMask("out", "0", spec, b).Eval(&Eval{Ctx: context.Background(), Circuit: c})
		if err != nil {
			t.Fatal(err)
		}
		return v
	}
	if v := eval([]Band{{From: 10, To: 300, MinDB: -1, MaxDB: 0.1}, {From: 1e4, To: 1e5, MinDB: math.Inf(-1), MaxDB: -19}}); v != 0 {
		t.Errorf("RC inside its mask: error %g", v)
	}
	if v := eval([]Band{{From: 10, To: 300, MinDB: -1, MaxDB: 0.1}, {From: 1e4, To: 1e5, MinDB: math.Inf(-1), MaxDB: -40}}); v < MaskViolationFloor {
		t.Errorf("RC cannot reach −40 dB at 10 kHz: error %g", v)
	}
}

func TestTopologyConfigErrors(t *testing.T) {
	emb := "e\nVIN vs 0 AC 1\nRS vs in 1k\nRL out 0 1k\n.ac dec 10 10 100k\n"
	top := &TopologyConfig{Terminals: []string{"in", "out"}, Parts: []PartKindConfig{{Kind: "C", Min: 1e-9, Max: 1e-6}}}
	mask := ObjectiveConfig{Kind: "mask", Node: "out", Bands: []BandConfig{{From: 10, To: 1000}}}
	cases := map[string]Config{
		"genes and topology": {Topology: top, Genes: []GeneConfig{{Element: "RS", Min: 1, Max: 2}}, Objectives: []ObjectiveConfig{mask}},
		"parts without topology": {Genes: []GeneConfig{{Element: "RS", Min: 1, Max: 2}},
			Objectives: []ObjectiveConfig{{Kind: "parts"}}},
		"band outside sweep": {Topology: top, Objectives: []ObjectiveConfig{{Kind: "mask", Node: "out",
			Bands: []BandConfig{{From: 1, To: 10}}}}},
		"part kind Q": {Topology: &TopologyConfig{Terminals: []string{"in"}, Parts: []PartKindConfig{{Kind: "Q", Min: 1, Max: 2}}},
			Objectives: []ObjectiveConfig{mask}},
	}
	for name, cfg := range cases {
		if _, err := NewExperiment(cfg, emb); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

// A topology run stopped part-way and resumed from its checkpoint ends
// bit-identical to an uninterrupted run.
func TestTopologyResume(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "filter_topology", "config.json"))
	if err != nil {
		t.Fatal(err)
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
		t.Fatal("resumed topology run differs")
	}
}

// The M8 acceptance test, run as `cmd/evolve` runs it: from an embryo with
// no parts between source and load, evolution grows a low-pass filter that
// meets the mask — passband 10 Hz–1 kHz within −7.5…−5.5 dB, stopband
// 3–100 kHz at or below −46 dB — checked on the written best netlist at
// ten times the sweep density evolution saw.
func TestFilterTopologyExample(t *testing.T) {
	if testing.Short() {
		t.Skip("about 40 s: the full 200-generation run")
	}
	ex, err := LoadExperiment(filepath.Join("..", "examples", "filter_topology", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 10, 25, 6, 0, 0, 0, time.UTC))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, nil, rd.Observe)
	if err != nil {
		t.Fatal(err)
	}
	rd.Close()
	src, err := os.ReadFile(filepath.Join(rd.Path, "best.cir"))
	if err != nil {
		t.Fatal(err)
	}
	d, err := netlist.ParseString(string(src))
	if err != nil {
		t.Fatalf("%v\n%s", err, src)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 200, FStart: 10, FStop: 1e5}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(r, "out", "0")
	worstPass, worstStop := 0.0, math.Inf(-1)
	minPass := 0.0
	for k, f := range r.Freq {
		db := 20 * math.Log10(cmplx.Abs(h[k]))
		switch {
		case f <= 1000*(1+1e-9): // the sweep's 1 kHz point carries rounding error
			if db < -7.5 || db > -5.5 {
				t.Errorf("passband %g Hz: %.3f dB", f, db)
			}
			minPass = math.Min(minPass, db)
			worstPass = math.Max(worstPass, math.Abs(db+6.0206))
		case f >= 3000*(1-1e-9):
			if db > -46 {
				t.Errorf("stopband %g Hz: %.3f dB", f, db)
			}
			worstStop = math.Max(worstStop, db)
		}
	}
	parts := len(c.Elements) - len(ex.Problem.Template.Elements)
	var names []string
	for _, e := range c.Elements[len(ex.Problem.Template.Elements):] {
		names = append(names, e.Name)
	}
	t.Logf("best after %d generations: %d parts (%s), passband down to %.3f dB, stopband at most %.2f dB",
		st.Generation, parts, strings.Join(names, " "), minPass, worstStop)
	if st.Best.Terms[0] != 0 {
		t.Errorf("mask error %g", st.Best.Terms[0])
	}
}
