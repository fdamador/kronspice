// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"math"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

func TestGridStep(t *testing.T) {
	step := func(v float64, ser string, integer bool, dir int, min, max float64) (float64, bool) {
		x, ok := gridStep(math.Log10(v), min, max, ser, integer, 0.02, dir)
		return math.Pow(10, x), ok
	}
	near := func(a, b float64) bool { return math.Abs(a-b) <= 1e-9*b }
	for _, c := range []struct {
		v      float64
		ser    string
		dir    int
		want   float64
		moved  bool
		lo, hi float64
	}{
		{1000, "E24", 1, 1100, true, 1, 1e6},
		{1000, "E24", -1, 910, true, 1, 1e6},
		{9100, "E24", 1, 10000, true, 1, 1e6}, // across a decade
		{1e6, "E24", 1, 1e6, false, 1, 1e6},   // at the top
		{4.7e-9, "E12", -1, 3.9e-9, true, 1e-12, 1e-6},
	} {
		got, ok := step(c.v, c.ser, false, c.dir, c.lo, c.hi)
		if ok != c.moved || (ok && !near(got, c.want)) {
			t.Errorf("%g %s dir %d: %g %v, want %g %v", c.v, c.ser, c.dir, got, ok, c.want, c.moved)
		}
	}
	if got, ok := step(7, "", true, 1, 1, 10); !ok || !near(got, 8) {
		t.Errorf("integer 7+1: %g", got)
	}
	if _, ok := step(10, "", true, 1, 1, 10); ok {
		t.Error("integer past the top moved")
	}
	if got, ok := step(100, "", false, -1, 1, 1e3); !ok || !near(got, 100/math.Pow(10, 0.02)) {
		t.Errorf("continuous step: %g", got)
	}
}

// dividerProblem: 10 V across R1 + R2; the mid-point must stay at 5 V
// (weight 10) while the supply current should reach 1 mA. From R1 = R2 =
// 10 kΩ every single-part E24 step breaks the ratio and costs more than it
// gains; stepping both parts down together keeps 5 V and moves the
// current toward 1 mA, down to 5.1 kΩ each.
func dividerProblem(t testing.TB) Problem {
	d, err := netlist.ParseString("divider\nV1 top 0 DC 10\nR1 top mid 10k\nR2 mid 0 10k\n")
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	mid := OPVoltage("mid", 5)
	mid.Weight = 10
	return Problem{Template: c, Space: Space{
		{Element: "R1", Min: 1e3, Max: 1e5, Series: "E24"},
		{Element: "R2", Min: 1e3, Max: 1e5, Series: "E24"},
	}, Objectives: []Objective{mid, supplyTarget(0.01)}}
}

func TestLocalSearchMovesPartsTogether(t *testing.T) {
	p := dividerProblem(t)
	p.Objectives[1] = supplyTarget(0.01) // 10 mW: 1 mA at 10 V
	opt := Options{}
	if err := opt.defaults(2); err != nil {
		t.Fatal(err)
	}
	ev := &evaluator{p: p, opt: opt}
	start := ev.evaluate(context.Background(), []Genome{{X: []float64{4, 4}, Sign: []int8{1, 1}}})[0]
	if start.Failure != "" {
		t.Fatal(start.Error)
	}
	one, _ := ev.improve(context.Background(), start, LocalSearch{Size: 1, MaxEvals: 1000, Step: 0.02})
	if one.Fitness != start.Fitness {
		t.Fatalf("single-part moves improved %g → %g; the test problem is meant to stall them", start.Fitness, one.Fitness)
	}
	two, evals := ev.improve(context.Background(), start, LocalSearch{Size: 2, MaxEvals: 1000, Step: 0.02})
	if !reflect.DeepEqual(two.Values, []float64{5100, 5100}) {
		t.Fatalf("pair moves reached %v (fitness %g), want 5.1 kΩ each", two.Values, two.Fitness)
	}
	t.Logf("start %g; size 1 stalls; size 2 reaches %v, fitness %g, in %d evaluations", start.Fitness, two.Values, two.Fitness, evals)
}

// supplyTarget scores the power the sources deliver against target [W].
func supplyTarget(target float64) Objective {
	sp := SupplyPower("V1")
	return Objective{Name: "supply", Weight: 1, Eval: func(e *Eval) (float64, error) {
		p, err := sp.Eval(e)
		if err != nil {
			return 0, err
		}
		return relErr(p, target, 1e-3), nil
	}}
}

func TestLocalSearchInRun(t *testing.T) {
	p := dividerProblem(t)
	p.Objectives[1] = supplyTarget(0.01)
	opt := Options{Population: 8, Generations: 6, Seed: 4, Workers: 3,
		LocalSearch: &LocalSearch{Size: 2, Every: 2}}
	full, err := Run(context.Background(), p, opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	plain := opt
	plain.LocalSearch = nil
	ga, err := Run(context.Background(), p, plain, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if full.Best.Fitness > ga.Best.Fitness || full.Evaluations <= ga.Evaluations {
		t.Errorf("with local search %g (%d evals), without %g (%d evals)", full.Best.Fitness, full.Evaluations, ga.Best.Fitness, ga.Evaluations)
	}
	// Checkpointed at generation 3 and resumed: bit-identical.
	stop := errors.New("stop")
	path := filepath.Join(t.TempDir(), "cp.json")
	_, err = Run(context.Background(), p, opt, nil, func(st *State) error {
		if st.Generation == 3 {
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
	resumed, err := Run(context.Background(), p, opt, cp, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(full.Population, resumed.Population) || full.Evaluations != resumed.Evaluations {
		t.Fatal("resumed run with local search differs")
	}
	// Options checked.
	bad := opt
	bad.Algorithm = NSGA2
	if _, err := Run(context.Background(), p, bad, nil, nil); !errors.Is(err, ErrBadOptions) {
		t.Errorf("NSGA-II with local search: %v", err)
	}
	bad = opt
	bad.LocalSearch = &LocalSearch{Size: 9}
	if _, err := Run(context.Background(), p, bad, nil, nil); !errors.Is(err, ErrBadOptions) {
		t.Errorf("size 9: %v", err)
	}
}

// The direct topology genome's part values are searched too; the
// developmental genome is left to the GA.
func TestLocalSearchTopology(t *testing.T) {
	top := testTopology()
	p := Problem{Topology: top, Template: testEmbryo(t)}
	g := Genome{Parts: []Part{{Kind: 'R', A: 1, B: 2, X: 3}, {Kind: 'C', A: 2, B: 0, X: -7}}}
	n, move := p.lsGenes(g, 0.02)
	nb := neighbours(g, n, 2, move)
	if n != 2 || len(nb) != 2*2+1*4 {
		t.Fatalf("%d genes, %d neighbours", n, len(nb))
	}
	if v := top.value(nb[0].Parts[0]); v != 1100 {
		t.Errorf("first neighbour R %g, want 1.1 kΩ", v)
	}
	p.Topology = testDev()
	if _, move := p.lsGenes(Genome{}, 0.02); move != nil {
		t.Error("developmental genome searched")
	}
}
