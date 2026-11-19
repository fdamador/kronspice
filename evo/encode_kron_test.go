// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"math/rand/v2"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func testKron() *Topology {
	t := testTopology()
	t.Encoding = Connection
	t.Library = []circuit.Kind{circuit.Inductor, circuit.Inductor, circuit.Capacitor, circuit.Capacitor, circuit.Resistor}
	return t
}

// An open slot stays in the genome and out of the circuit; closing it
// brings the same branch back.
func TestKronDecodeOpenSlots(t *testing.T) {
	top, emb := testKron(), testEmbryo(t)
	if err := top.Validate(emb); err != nil {
		t.Fatal(err)
	}
	// in –L– n1 –L– out, C from n1 to ground, second C and R open.
	g := Genome{Parts: []Part{
		{Kind: 'L', A: 1, B: 3, X: -3}, {Kind: 'L', A: 3, B: 2, X: -3},
		{Kind: 'C', A: 3, B: 0, X: -7}, {Kind: 'C', A: 2, B: 0, X: -7, Off: true},
		{Kind: 'R', A: 1, B: 2, X: 3, Off: true},
	}}
	c, vals := top.Decode(emb, g)
	if n := len(c.Elements) - len(emb.Elements); n != 3 || len(vals) != 3 {
		t.Fatalf("%d parts, want 3", n)
	}
	if g.size() != 3 {
		t.Errorf("size %d, want 3 connected slots", g.size())
	}
	g.Parts[4].Off = false
	c, _ = top.Decode(emb, g)
	if n := len(c.Elements) - len(emb.Elements); n != 4 {
		t.Fatalf("after closing the R slot: %d parts, want 4", n)
	}
}

// identify joins two junctions and tear splits one: on the T network
// above, identifying n1 with out shorts the second inductor (pruned as a
// self-loop); tearing then restores a junction count of the original.
func TestKronIdentifyTear(t *testing.T) {
	top, emb := testKron(), testEmbryo(t)
	_ = top.Validate(emb)
	base := Genome{Parts: []Part{
		{Kind: 'L', A: 1, B: 3, X: -3}, {Kind: 'L', A: 3, B: 2, X: -3},
		{Kind: 'C', A: 3, B: 0, X: -7}, {Kind: 'C', A: 2, B: 0, X: -7, Off: true},
		{Kind: 'R', A: 1, B: 2, X: 3, Off: true},
	}}
	junctions := func(g Genome) int {
		n := 0
		for _, u := range usedNodes(activeParts(g.Parts), top.nodes()) {
			if u {
				n++
			}
		}
		return n
	}
	sawIdentify, sawTear := false, false
	for seed := uint64(1); seed < 200 && !(sawIdentify && sawTear); seed++ {
		r := rand.New(rand.NewPCG(seed, 7))
		g := base.Clone()
		top.kronApply(r, &g, kmIdentify, 0)
		top.kronCanonical(&g)
		if junctions(g) == junctions(base)-1 {
			sawIdentify = true
			// The terminal incidences are kept: every slot still has its ends.
			for i := range g.Parts {
				if g.Parts[i].Off != base.Parts[i].Off || g.Parts[i].Kind != base.Parts[i].Kind {
					t.Fatal("identify changed a slot's state or kind")
				}
			}
		}
		h := base.Clone()
		top.kronApply(r, &h, kmTear, 0)
		top.kronCanonical(&h)
		if junctions(h) == junctions(base)+1 {
			sawTear = true
		}
	}
	if !sawIdentify || !sawTear {
		t.Errorf("identify seen %v, tear seen %v", sawIdentify, sawTear)
	}
}

// Every operator keeps the slots aligned with the library, the labels
// canonical and in range, and every genome decodable.
func TestKronOperatorInvariants(t *testing.T) {
	top, emb := testKron(), testEmbryo(t)
	if err := top.Validate(emb); err != nil {
		t.Fatal(err)
	}
	lib := top.slotKinds()
	check := func(g Genome, what string) {
		t.Helper()
		if len(g.Parts) != len(lib) {
			t.Fatalf("%s: %d slots, want %d", what, len(g.Parts), len(lib))
		}
		term := 1 + len(top.Terminals)
		next := term
		for _, off := range []bool{false, true} {
			for i, p := range g.Parts {
				if p.Kind != byte(lib[i]) {
					t.Fatalf("%s: slot %d kind %c, want %c", what, i, p.Kind, lib[i])
				}
				if p.Off != off {
					continue
				}
				for _, n := range p.terminals() {
					if n < 0 || n >= top.nodes() {
						t.Fatalf("%s: node %d out of range", what, n)
					}
					if n >= term {
						if n > next {
							t.Fatalf("%s: labels not canonical: %d before %d (%+v)", what, n, next, g.Parts)
						}
						if n == next {
							next++
						}
					}
				}
			}
		}
		top.Decode(emb, g)
	}
	r := rand.New(rand.NewPCG(3, 4))
	pop := make([]Genome, 20)
	for i := range pop {
		pop[i] = top.randomGenome(r)
		check(pop[i], "random")
	}
	for step := 0; step < 2000; step++ {
		a, b := pop[r.IntN(len(pop))], pop[r.IntN(len(pop))]
		c1, c2 := top.crossover(r, a, b)
		check(c1, "crossover")
		check(c2, "crossover")
		top.mutate(r, &c1, 0.3)
		top.repair(&c1)
		check(c1, "mutate")
		pop[r.IntN(len(pop))] = c1
	}
}

// Crossover is aligned: each child slot comes whole from one parent.
func TestKronCrossoverAligned(t *testing.T) {
	top := testKron()
	_ = top.Validate(testEmbryo(t))
	r := rand.New(rand.NewPCG(9, 9))
	a, b := top.randomGenome(r), top.randomGenome(r)
	c1, c2 := top.crossover(r, a, b)
	for i := range c1.Parts {
		x, y := c1.Parts[i], c2.Parts[i]
		// Values and open state travel with the slot (labels may be renumbered).
		fromA := x.X == a.Parts[i].X && y.X == b.Parts[i].X && x.Off == a.Parts[i].Off
		fromB := x.X == b.Parts[i].X && y.X == a.Parts[i].X && x.Off == b.Parts[i].Off
		if !fromA && !fromB {
			t.Fatalf("slot %d mixed: %+v %+v from %+v %+v", i, x, y, a.Parts[i], b.Parts[i])
		}
	}
}

func TestKronValidateAndConfig(t *testing.T) {
	emb := testEmbryo(t)
	top := testKron()
	top.Library = append(top.Library, circuit.BJT)
	if err := top.Validate(emb); err == nil {
		t.Error("library slot of an unlisted kind accepted")
	}
	d := testTopology()
	d.Library = []circuit.Kind{circuit.Resistor}
	if err := d.Validate(emb); err == nil {
		t.Error("library accepted for the direct encoding")
	}
	if e, err := ParseEncoding("kron"); err != nil || e != Connection || e.String() != "kron" {
		t.Errorf("ParseEncoding(kron) = %v, %v", e, err)
	}
	// Default library: MaxParts slots round-robin over the kinds.
	top = testKron()
	top.Library = nil
	if got := top.slotKinds(); len(got) != top.MaxParts || got[3] != circuit.Inductor {
		t.Errorf("default library %v", got)
	}
}

// m = b − N + S: the T network (in–L–n1–L–out, C from n1) on the embryo
// (source to ground, RS, RL to ground) has meshes source–RS–L–C and
// C–L–RL; adding a resistor from in to out adds a third.
func TestMeshes(t *testing.T) {
	top, emb := testKron(), testEmbryo(t)
	_ = top.Validate(emb)
	g := Genome{Parts: []Part{
		{Kind: 'L', A: 1, B: 3, X: -3}, {Kind: 'L', A: 3, B: 2, X: -3},
		{Kind: 'C', A: 3, B: 0, X: -7}, {Kind: 'C', A: 2, B: 0, X: -7, Off: true},
		{Kind: 'R', A: 1, B: 2, X: 3, Off: true},
	}}
	c, _ := top.Decode(emb, g)
	if m := Meshes(c); m != 2 {
		t.Errorf("meshes %d, want 2", m)
	}
	g.Parts[4].Off = false
	c, _ = top.Decode(emb, g)
	if m := Meshes(c); m != 3 {
		t.Errorf("meshes %d, want 3", m)
	}
	if m := Meshes(emb); m != 0 {
		t.Errorf("bare embryo meshes %d, want 0", m)
	}
	v, ok := DescMeshes(0, 10, 10).Of(Problem{Topology: top, Template: emb}, Individual{Genome: g})
	if !ok || v != 3 {
		t.Errorf("descriptor %g %v", v, ok)
	}
}
