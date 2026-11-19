// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// eceEmbryo is the filter fixture with the ECE2 spin-connection terms
// declared (hypothetical, H48) and a 1 V DC drive.
func eceEmbryo(t testing.TB) *circuit.Circuit {
	d, err := netlist.ParseString("embryo\n.ece form=ece2\nVIN vs 0 DC 1 AC 1\nRS vs in 1k\nRL out 0 1k\n")
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// polarDev is a developmental topology whose resistors carry κ_a = 0.75 at
// terminal A and κ_b = 0 at B: with ECE2 (s = −2) the body sees
// −0.5·φ_A − φ_B.
func polarDev() *Topology {
	t := testDev()
	t.Kinds = []PartKind{{Kind: circuit.Resistor, Min: 1, Max: 1e6, KappaA: 0.75}}
	return t
}

func opVolts(t *testing.T, c *circuit.Circuit) (vin, vout, model float64) {
	t.Helper()
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	vin, _ = r.Voltage("in")
	vout, _ = r.Voltage("out")
	return vin, vout, r.Ledger.ModelSupplied
}

// FLIP reverses the component under the head.
func TestFlipDevelop(t *testing.T) {
	top := testDev()
	c, els := grown(t, top, dFn(dvFlip, 0, dPart(dvR, 0, dEnd())))
	if len(els) != 1 || ends(c, els[0]) != [2]string{"out", "in"} {
		t.Fatalf("(FLIP (R 0 END)) grew %+v", els)
	}
	if s := FormatProgram(dFn(dvFlip, 0, dEnd())); s != "(FLIP END)" {
		t.Errorf("%s", s)
	}
}

// For plain R, L, C (Maxwell) a flipped part is the same circuit: every
// node voltage is identical, to the bit.
func TestFlipMaxwellNoOp(t *testing.T) {
	top, emb := testDev(), eceEmbryo(t)
	emb.ECE = nil
	prog := dFn(dvVia, 0, dPart(dvR, 0, dEnd()), dPart(dvL, 0, dEnd()), dPart(dvR, 1, dEnd()))
	flipped := dFn(dvVia, 0, dFn(dvFlip, 0, dPart(dvR, 0, dEnd())), dPart(dvL, 0, dEnd()), dFn(dvFlip, 0, dPart(dvR, 1, dEnd())))
	a, _ := top.Decode(emb, Genome{Program: [][]Op{prog}})
	b, _ := top.Decode(emb, Genome{Program: [][]Op{flipped}})
	ai, ao, _ := opVolts(t, a)
	bi, bo, _ := opVolts(t, b)
	if ai != bi || ao != bo {
		t.Errorf("flip changed a Maxwell circuit: in %g → %g, out %g → %g", ai, bi, ao, bo)
	}
}

// For a polarized part the orientation is the circuit. One resistor R =
// 1 kΩ from in to out between RS = RL = 1 kΩ, body current
// (−0.5·φ_A − φ_B)/R from A to B:
//
//	A at in:  V(in) = 4/3, V(out) = −1/3 (inverting)
//	A at out: V(in) = 1/3, V(out) = 2/3 (above the passive limit ½)
//
// Plain, it is 2/3 and 1/3. Both polarized cases draw model-supplied power
// (booked to the spin-connection term, never as inflow).
func TestFlipPolarized(t *testing.T) {
	top, emb := polarDev(), eceEmbryo(t)
	r1k := dPart(dvR, 0, dEnd()) // middle of 1 Ω..1 MΩ: 1 kΩ
	cases := []struct {
		name      string
		prog      []Op
		vin, vout float64
	}{
		{"A at in", r1k, 4.0 / 3, -1.0 / 3},
		{"A at out", dFn(dvFlip, 0, r1k), 1.0 / 3, 2.0 / 3},
	}
	for _, tc := range cases {
		c, vals := top.Decode(emb, Genome{Program: [][]Op{tc.prog}})
		if len(vals) != 1 || vals[0] != 1000 {
			t.Fatalf("%s: values %v", tc.name, vals)
		}
		if e := c.Elements[len(c.Elements)-1]; e.Params["kappa_a"] != 0.75 || e.Params["kappa_b"] != 0 {
			t.Fatalf("%s: params %v", tc.name, e.Params)
		}
		vin, vout, model := opVolts(t, c)
		if math.Abs(vin-tc.vin) > 1e-9 || math.Abs(vout-tc.vout) > 1e-9 {
			t.Errorf("%s: V(in) %.9g V(out) %.9g, want %.9g %.9g", tc.name, vin, vout, tc.vin, tc.vout)
		}
		if !(math.Abs(model) > 1e-6) {
			t.Errorf("%s: model-supplied %g W, want the spin-connection term's power booked", tc.name, model)
		}
	}
	// κ_a = κ_b: not polarized, orientation irrelevant again. (κ = 0.6 makes
	// a −5 kΩ part; 0.75 would make −2 kΩ, which with RS and RL is exactly
	// singular.)
	sym := polarDev()
	sym.Kinds[0].KappaA, sym.Kinds[0].KappaB = 0.6, 0.6
	a, _ := sym.Decode(emb, Genome{Program: [][]Op{r1k}})
	b, _ := sym.Decode(emb, Genome{Program: [][]Op{dFn(dvFlip, 0, r1k)}})
	ai, ao, _ := opVolts(t, a)
	bi, bo, _ := opVolts(t, b)
	if math.Abs(ai-bi) > 1e-12 || math.Abs(ao-bo) > 1e-12 || sym.polarized() {
		t.Errorf("κ_a = κ_b: flip changed in %g → %g, out %g → %g", ai, bi, ao, bo)
	}
}

// FLIP enters random programs (and the direct encoding's flip mutation
// fires) only when a part kind is polarized, so runs without ECE parts
// draw exactly the random numbers they drew before FLIP existed.
func TestFlipOnlyWhenPolarized(t *testing.T) {
	count := func(top *Topology) int {
		r := rngFor(5, 0)
		n := 0
		for i := 0; i < 300; i++ {
			for _, op := range top.randomProgram(r) {
				if op.F == dvFlip {
					n++
				}
			}
		}
		return n
	}
	if n := count(testDev()); n != 0 {
		t.Errorf("plain R, L, C: %d FLIPs drawn", n)
	}
	if n := count(polarDev()); n == 0 {
		t.Error("polarized parts: no FLIP drawn")
	}
	// Direct encoding: the flip mutation reverses parts only when polarized.
	direct := testTopology()
	direct.Kinds = []PartKind{{Kind: circuit.Resistor, Min: 1, Max: 1e6, KappaA: 0.75}}
	r := rngFor(9, 0)
	g := Genome{Parts: []Part{{Kind: 'R', A: 1, B: 2, X: 3}}}
	flipped := false
	for i := 0; i < 500 && !flipped; i++ {
		h := g.Clone()
		direct.apply(r, &h, opFlip, 0.2)
		flipped = h.Parts[0].A == 2 && h.Parts[0].B == 1
	}
	if !flipped {
		t.Error("opFlip did not reverse the part")
	}
}

func TestFlipValidate(t *testing.T) {
	top := polarDev()
	if err := top.Validate(eceEmbryo(t)); err != nil {
		t.Fatal(err)
	}
	if err := polarDev().Validate(testEmbryo(t)); err == nil {
		t.Error("kappa without .ece accepted")
	}
	bad := polarDev()
	bad.Kinds[0].KappaB = math.NaN()
	if err := bad.Validate(eceEmbryo(t)); err == nil {
		t.Error("NaN kappa accepted")
	}
}

// Evolution uses FLIP: the target V(out) = 1.5 V from a 1 V drive lies
// beyond every passive circuit (at most ½ V with RS = RL = 1 kΩ). Plain
// resistors cannot reach it; resistors with κ_a = 0.75 (ECE2) can, and the
// best design's power comes from the spin-connection term, booked as
// model-supplied.
func TestFlipEvolvesBeyondPassive(t *testing.T) {
	if testing.Short() {
		t.Skip("two short GA runs")
	}
	run := func(top *Topology, emb *circuit.Circuit) *State {
		p := Problem{Topology: top, Template: emb, Objectives: []Objective{OPVoltage("out", 1.5)}}
		opt := Options{Population: 100, Generations: 40, Seed: 1, Target: 1e-6, HasTarget: true}
		st, err := Run(context.Background(), p, opt, nil, nil)
		if err != nil {
			t.Fatal(err)
		}
		return st
	}
	plain := testDev()
	plain.Kinds = []PartKind{{Kind: circuit.Resistor, Min: 1, Max: 1e6}}
	ps := run(plain, eceEmbryo(t))
	if ps.Best.Failure == "" && ps.Best.Fitness < 0.1 {
		t.Errorf("plain resistors reached V(out) = 1.5 V (error %g): impossible passively", ps.Best.Fitness)
	}
	pol := polarDev()
	st := run(pol, eceEmbryo(t))
	b := st.Best
	if b.Failure != "" || b.Fitness > 1e-3 {
		t.Fatalf("polarized: best error %g (%s %s)", b.Fitness, b.Failure, b.Error)
	}
	if !(b.ModelSupplied > 0) {
		t.Errorf("polarized best: model-supplied %g W, want > 0", b.ModelSupplied)
	}
	flips := 0
	for _, op := range b.Genome.Program[0] {
		if op.F == dvFlip {
			flips++
		}
	}
	c, _ := pol.Decode(eceEmbryo(t), b.Genome)
	t.Logf("plain best error %g; polarized best error %g after generation %d, %d FLIPs in %d ops, model-supplied %.3g W",
		ps.Best.Fitness, b.Fitness, st.Generation, flips, len(b.Genome.Program[0]), b.ModelSupplied)
	for _, e := range c.Elements[3:] {
		t.Logf("  %s %s→%s %g", e.Name, c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]], e.Value)
	}
}
