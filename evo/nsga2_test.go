// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"reflect"
	"sort"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func TestFastNonDominatedSort(t *testing.T) {
	pts := [][]float64{{1, 5}, {2, 2}, {3, 3}, {5, 1}, {2, 2}, {4, 4}, {6, 6}}
	got := FastNonDominatedSort(pts)
	want := [][]int{{0, 1, 3, 4}, {2}, {5}, {6}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("fronts = %v, want %v", got, want)
	}
}

// Front 0 of the fast sort is exactly the brute-force non-dominated set
// (with duplicates kept), and every point of front r+1 is dominated by a
// point of front r.
func TestFastNonDominatedSortProperties(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	dom := func(p, q []float64) bool {
		lt := false
		for k := range p {
			if p[k] > q[k] {
				return false
			}
			lt = lt || p[k] < q[k]
		}
		return lt
	}
	for trial := 0; trial < 200; trial++ {
		n, m := 1+r.IntN(40), 1+r.IntN(4)
		pts := make([][]float64, n)
		for i := range pts {
			pts[i] = make([]float64, m)
			for k := range pts[i] {
				pts[i][k] = float64(r.IntN(6)) // many ties
			}
		}
		fronts := FastNonDominatedSort(pts)
		seen := 0
		for fi, f := range fronts {
			seen += len(f)
			if !sort.IntsAreSorted(f) {
				t.Fatalf("front not sorted: %v", f)
			}
			for _, a := range f {
				for _, b := range f {
					if dom(pts[a], pts[b]) {
						t.Fatalf("trial %d: %v dominates %v inside front %d", trial, pts[a], pts[b], fi)
					}
				}
				if fi > 0 {
					ok := false
					for _, b := range fronts[fi-1] {
						ok = ok || dom(pts[b], pts[a])
					}
					if !ok {
						t.Fatalf("trial %d: %v in front %d not dominated by front %d", trial, pts[a], fi, fi-1)
					}
				}
			}
		}
		if seen != n {
			t.Fatalf("trial %d: %d of %d points placed", trial, seen, n)
		}
		var brute []int
		for i := range pts {
			d := false
			for j := range pts {
				d = d || dom(pts[j], pts[i])
			}
			if !d {
				brute = append(brute, i)
			}
		}
		if !reflect.DeepEqual(fronts[0], brute) {
			t.Fatalf("trial %d: front 0 %v, brute force %v", trial, fronts[0], brute)
		}
	}
}

func TestCrowdingDistance(t *testing.T) {
	// Points on a line f2 = 4 − f1; interior distance = (gap₁/4) + (gap₂/4).
	pts := [][]float64{{0, 4}, {1, 3}, {3, 1}, {4, 0}}
	d := CrowdingDistance(pts, []int{0, 1, 2, 3})
	if !math.IsInf(d[0], 1) || !math.IsInf(d[3], 1) {
		t.Fatalf("extremes = %v", d)
	}
	if math.Abs(d[1]-1.5) > 1e-15 || math.Abs(d[2]-1.5) > 1e-15 {
		t.Fatalf("interior = %v, want 1.5, 1.5", d)
	}
	if d := CrowdingDistance(pts, []int{1, 2}); !math.IsInf(d[0], 1) || !math.IsInf(d[1], 1) {
		t.Fatalf("two-point front = %v", d)
	}
}

// Constrained dominance: a failed candidate ranks after every feasible
// one, and failed candidates order by penalty.
func TestRankPopulationFailures(t *testing.T) {
	obs := []Objective{{Name: "a", Weight: 1}, {Name: "b", Weight: 1}}
	pop := []Individual{
		{Failure: "singular", Fitness: 1e6},
		{Terms: []float64{9, 9}, Fitness: 18},
		{Failure: "no_convergence", Fitness: 1e5},
		{Terms: []float64{1, 2}, Fitness: 3},
		{Terms: []float64{2, 1}, Fitness: 3},
	}
	rankPopulation(pop, obs)
	ranks := []int{pop[0].Rank, pop[1].Rank, pop[2].Rank, pop[3].Rank, pop[4].Rank}
	if !reflect.DeepEqual(ranks, []int{3, 1, 2, 0, 0}) {
		t.Fatalf("ranks = %v", ranks)
	}
	if f := ParetoFront(pop); len(f) != 2 || f[0].Fitness != 3 {
		t.Fatalf("front = %+v", f)
	}
	all := []Individual{{Failure: "x", Fitness: 5}, {Failure: "y", Fitness: 5}}
	rankPopulation(all, obs)
	if all[0].Rank != 0 || all[1].Rank != 0 || len(ParetoFront(all)) != 0 {
		t.Fatalf("all-failed ranks = %+v", all)
	}
}

// zdtGenes maps 30 log-scale genes onto x ∈ [0, 1]: x = 1 + log10(v)/3 on
// v ∈ [1e−3, 1].
func zdtX(vals []float64) []float64 {
	x := make([]float64, len(vals))
	for i, v := range vals {
		x[i] = 1 + math.Log10(v)/3
	}
	return x
}

// zdtProblem is a ZDT test problem (E. Zitzler, K. Deb, L. Thiele,
// "Comparison of multiobjective evolutionary algorithms: empirical
// results," Evol. Comput. 8(2), 2000) on a trivial circuit, so the whole
// engine path (decode, validate, worker pool) runs.
func zdtProblem(n int, h func(f1, g float64) float64) Problem {
	space := make(Space, n)
	for i := range space {
		space[i] = Gene{Name: fmt.Sprintf("x%d", i+1), Min: 1e-3, Max: 1}
	}
	build := func([]float64) (*circuit.Circuit, error) {
		b := circuit.NewBuilder("zdt")
		b.VoltageSource("V1", "a", "0", 1)
		b.Resistor("R1", "a", "0", 1)
		return b.Build()
	}
	g := func(x []float64) float64 {
		s := 0.0
		for _, v := range x[1:] {
			s += v
		}
		return 1 + 9*s/float64(len(x)-1)
	}
	return Problem{Build: build, Space: space, Objectives: []Objective{
		{Name: "f1", Weight: 1, Eval: func(e *Eval) (float64, error) { return zdtX(e.Values)[0], nil }},
		{Name: "f2", Weight: 1, Eval: func(e *Eval) (float64, error) {
			x := zdtX(e.Values)
			gg := g(x)
			return gg * h(x[0], gg), nil
		}},
	}}
}

func zdt1(f1, g float64) float64 { return 1 - math.Sqrt(f1/g) }
func zdt2(f1, g float64) float64 { return 1 - (f1/g)*(f1/g) }

// The evolved front approaches the true ZDT front (g = 1) and spreads
// along it: convex (ZDT1) and non-convex (ZDT2), where a weighted sum can
// only ever find the two ends.
func TestNSGA2ZDT(t *testing.T) {
	for _, tc := range []struct {
		name string
		h    func(f1, g float64) float64
	}{{"ZDT1", zdt1}, {"ZDT2", zdt2}} {
		t.Run(tc.name, func(t *testing.T) {
			p := zdtProblem(10, tc.h)
			st, err := Run(context.Background(), p, Options{Algorithm: NSGA2, Population: 60, Generations: 150,
				Seed: 7, Mutation: 0.1, Sigma: 0.3}, nil, nil)
			if err != nil {
				t.Fatal(err)
			}
			if len(st.Front) < 30 {
				t.Fatalf("front has only %d members", len(st.Front))
			}
			// Distance to the true front: f2 − h(f1, 1) over the front.
			worst, lo, hi := 0.0, 1.0, 0.0
			interior := 0
			for _, ind := range st.Front {
				f1, f2 := ind.Terms[0], ind.Terms[1]
				worst = math.Max(worst, f2-tc.h(f1, 1))
				lo, hi = math.Min(lo, f1), math.Max(hi, f1)
				if f1 > 0.2 && f1 < 0.8 {
					interior++
				}
			}
			t.Logf("%s: %d front members, f1 ∈ [%.3f, %.3f], %d interior, worst gap to true front %.4f, %d evaluations",
				tc.name, len(st.Front), lo, hi, interior, worst, st.Evaluations)
			if worst > 0.05 {
				t.Errorf("front is %.3f above the true front", worst)
			}
			if lo > 0.02 || hi < 0.98 || interior < 10 {
				t.Errorf("front does not spread: f1 ∈ [%g, %g], %d interior points", lo, hi, interior)
			}
		})
	}
}

// NSGA-II runs are deterministic across worker counts and resume exactly
// from a checkpoint.
func TestNSGA2DeterministicAndResumable(t *testing.T) {
	p := zdtProblem(6, zdt1)
	opt := Options{Algorithm: NSGA2, Population: 24, Generations: 10, Seed: 11, Workers: 1}
	a, err := Run(context.Background(), p, opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	opt.Workers = 8
	b, err := Run(context.Background(), p, opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(a, b) {
		t.Fatal("result depends on the worker count")
	}
	var buf bytes.Buffer
	stop := errors.New("stop")
	_, err = Run(context.Background(), p, opt, nil, func(st *State) error {
		if st.Generation == 4 {
			if err := WriteCheckpoint(&buf, st); err != nil {
				return err
			}
			return stop
		}
		return nil
	})
	if !errors.Is(err, stop) {
		t.Fatal(err)
	}
	cp, err := ReadCheckpoint(&buf)
	if err != nil {
		t.Fatal(err)
	}
	c, err := Run(context.Background(), p, opt, cp, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(a, c) {
		t.Fatal("resumed NSGA-II run differs")
	}
}

func TestParseAlgorithm(t *testing.T) {
	for s, want := range map[string]Algorithm{"": WeightedSum, "weighted": WeightedSum, "NSGA2": NSGA2, "nsga-ii": NSGA2} {
		if got, err := ParseAlgorithm(s); err != nil || got != want {
			t.Errorf("%q: %v, %v", s, got, err)
		}
	}
	if _, err := ParseAlgorithm("spea2"); !errors.Is(err, ErrBadOptions) {
		t.Errorf("spea2: %v", err)
	}
}

// Constraint violators rank after feasible candidates and before failed
// ones, ordered by weighted violation; constraints are not objectives.
func TestRankPopulationConstraints(t *testing.T) {
	obs := []Objective{{Name: "a", Weight: 1}, {Name: "c", Weight: 2, Constraint: true}}
	pop := []Individual{
		{Failure: "singular", Fitness: 1e6},
		{Terms: []float64{0, 0.5}},  // best objective but violates: v = 1
		{Terms: []float64{5, 0}},    // feasible
		{Terms: []float64{0, 0.25}}, // v = 0.5
		{Terms: []float64{3, 0}},    // feasible, dominates [2]
	}
	rankPopulation(pop, obs)
	ranks := []int{pop[0].Rank, pop[1].Rank, pop[2].Rank, pop[3].Rank, pop[4].Rank}
	if !reflect.DeepEqual(ranks, []int{4, 3, 1, 2, 0}) {
		t.Fatalf("ranks = %v", ranks)
	}
	if pop[1].Violation != 1 || pop[3].Violation != 0.5 {
		t.Fatalf("violations %g, %g", pop[1].Violation, pop[3].Violation)
	}
	if f := ParetoFront(pop); len(f) != 1 || f[0].Terms[0] != 3 {
		t.Fatalf("front = %+v", f)
	}
}

// ZDT1 with a constraint f1 ≥ 0.5: the front lies on the feasible half.
func TestNSGA2Constrained(t *testing.T) {
	p := zdtProblem(10, zdt1)
	p.Objectives = append(p.Objectives, Objective{Name: "f1 ≥ 0.5", Weight: 1, Constraint: true,
		Eval: func(e *Eval) (float64, error) { return math.Max(0, 0.5-zdtX(e.Values)[0]), nil }})
	st, err := Run(context.Background(), p, Options{Algorithm: NSGA2, Population: 60, Generations: 100,
		Seed: 3, Mutation: 0.1, Sigma: 0.3}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	lo := 1.0
	for _, ind := range st.Front {
		lo = math.Min(lo, ind.Terms[0])
	}
	if len(st.Front) < 20 || lo < 0.5 || lo > 0.51 {
		t.Fatalf("front of %d, smallest f1 %g (want just ≥ 0.5)", len(st.Front), lo)
	}
	only := Problem{Build: p.Build, Space: p.Space, Objectives: p.Objectives[2:]}
	if _, err := Run(context.Background(), only, Options{Algorithm: NSGA2}, nil, nil); !errors.Is(err, ErrBadOptions) {
		t.Fatalf("constraint-only problem: %v", err)
	}
}
