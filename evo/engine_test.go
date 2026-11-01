// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"bytes"
	"context"
	"errors"
	"math"
	"reflect"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

func deck(t testing.TB, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// RC low-pass, the M5 acceptance problem: R in E24, C in E12, target
// cutoff 1 kHz.
func rcProblem(t testing.TB) Problem {
	tpl := deck(t, "rc\nV1 in 0 AC 1\nR1 in out 10k\nC1 out 0 10n\n.ac dec 100 10 100k\n")
	spec := analysis.ACSpec{Variation: "dec", Points: 100, FStart: 10, FStop: 1e5}
	return Problem{
		Template: tpl,
		Space: Space{
			{Element: "R1", Min: 100, Max: 1e6, Series: "E24"},
			{Element: "C1", Min: 1e-10, Max: 1e-5, Series: "E12"},
		},
		Objectives: []Objective{Cutoff("out", "0", spec, 1000, true)},
	}
}

func TestRCLowpassHitsCutoff(t *testing.T) {
	st, err := Run(context.Background(), rcProblem(t),
		Options{Population: 40, Generations: 40, Seed: 1, Target: 0.002, HasTarget: true}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	r, c := st.Best.Values[0], st.Best.Values[1]
	fc := 1 / (2 * math.Pi * r * c) // exact cutoff of the decoded values
	if math.Abs(fc-1000) > 10 {
		t.Fatalf("best R = %g, C = %g: fc = %.2f Hz, not within 1%% of 1 kHz (fitness %g)", r, c, fc, st.Best.Fitness)
	}
	t.Logf("R = %g Ω, C = %g F, fc = %.2f Hz after %d generations, %d evaluations",
		r, c, fc, st.Generation, st.Evaluations)
}

// Same seed ⇒ identical run, whatever the number of workers.
func TestDeterministicAcrossWorkers(t *testing.T) {
	p := rcProblem(t)
	run := func(w int) *State {
		st, err := Run(context.Background(), p, Options{Population: 24, Generations: 6, Seed: 7, Workers: w}, nil, nil)
		if err != nil {
			t.Fatal(err)
		}
		return st
	}
	a, b := run(1), run(8)
	if !reflect.DeepEqual(a, b) {
		t.Fatal("1 worker and 8 workers gave different runs")
	}
	c, _ := Run(context.Background(), p, Options{Population: 24, Generations: 6, Seed: 8}, nil, nil)
	if reflect.DeepEqual(a.Best, c.Best) {
		t.Fatal("different seeds gave the same best")
	}
}

// A run checkpointed after generation 3 and resumed ends exactly where an
// uninterrupted run does.
func TestCheckpointResumeIsExact(t *testing.T) {
	p := rcProblem(t)
	opt := Options{Population: 20, Generations: 8, Seed: 3}
	full, err := Run(context.Background(), p, opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	var buf bytes.Buffer
	stop := errors.New("stop")
	_, err = Run(context.Background(), p, opt, nil, func(st *State) error {
		if st.Generation == 3 {
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
	resumed, err := Run(context.Background(), p, opt, cp, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(full, resumed) {
		t.Fatalf("resumed run differs: best %+v vs %+v", full.Best, resumed.Best)
	}
}

// Broken candidates are penalized by class, never crash the run.
func TestFailuresArePenalized(t *testing.T) {
	// R2 may take either sign; at R2 = −R1 the divider node has no
	// conductance (singular). A 1-ns timeout on a long transient times out.
	tpl := deck(t, "f\nV1 in 0 1\nR1 in out 1k\nR2 out 0 1k\n.op\n")
	ev := &evaluator{p: Problem{Template: tpl, Space: Space{{Element: "R2", Min: 1, Max: 1e4, Sign: Either}},
		Objectives: []Objective{OPVoltage("out", 0.5)}}, opt: Options{Timeout: time.Second}}
	ind := ev.one(context.Background(), Genome{X: []float64{3}, Sign: []int8{-1}})
	if ind.Failure != "singular" || ind.Fitness != PenaltySingular {
		t.Fatalf("singular candidate: %+v", ind)
	}
	ind = ev.one(context.Background(), Genome{X: []float64{3}, Sign: []int8{1}})
	if ind.Failure != "" || ind.Fitness > 1e-9 {
		t.Fatalf("good candidate: %+v", ind)
	}
	slow := deck(t, "s\nV1 in 0 SIN(0 1 1meg)\nR1 in out 1k\nC1 out 0 1n\n.tran 1n 10m\n")
	ev = &evaluator{p: Problem{Template: slow, Space: Space{{Element: "R1", Min: 1, Max: 1e4}},
		Objectives: []Objective{TranAverage("out", analysis.TranSpec{Step: 1e-9, Stop: 10e-3}, 0, 10e-3, 0)}},
		opt: Options{Timeout: time.Millisecond}}
	ind = ev.one(context.Background(), Genome{X: []float64{3}, Sign: []int8{1}})
	if ind.Failure != "timeout" || ind.Fitness != PenaltyTimeout {
		t.Fatalf("slow candidate: %+v", ind)
	}
}

// The gain target needs |H| = 2 from a divider, which only a negative lower
// leg can give: evolution must find the sign. The model-supplied power is
// reported with the winner.
func TestEvolvesNegativeResistance(t *testing.T) {
	tpl := deck(t, "n\nV1 in 0 DC 1 AC 1\nR1 in out 1k\nR2 out 0 1k\n.ac lin 1 1k 1k\n")
	spec := analysis.ACSpec{Variation: "lin", Points: 1, FStart: 1e3, FStop: 1e3}
	p := Problem{Template: tpl,
		Space:      Space{{Element: "R2", Min: 10, Max: 1e5, Sign: Either}},
		Objectives: []Objective{GainAt("out", "0", spec, 1e3, 20*math.Log10(2))}}
	st, err := Run(context.Background(), p, Options{Population: 30, Generations: 30, Seed: 5, Sigma: 0.3}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	r2 := st.Best.Values[0]
	// |R2/(1k + R2)| = 2 ⇒ R2 = −2 kΩ or −667 Ω.
	if r2 >= 0 || st.Best.Fitness > 1e-3 {
		t.Fatalf("best R2 = %g, fitness %g", r2, st.Best.Fitness)
	}
	if st.Best.ModelSupplied <= 0 {
		t.Fatal("negative resistance's power not reported as model-supplied")
	}
}

// Nonlinear: evolve the filter capacitor of a half-wave rectifier for at
// most 0.5 V ripple with the smallest capacitor that does it (a second,
// lightly weighted objective pulls C down).
func TestEvolvesRectifierFilter(t *testing.T) {
	tpl := deck(t, "r\nV1 in 0 SIN(0 10 1k)\nD1 in out dr\nR1 out 0 1k\nC1 out 0 1u\n.model dr D(IS=1e-14)\n.tran 10u 6m\n")
	spec := analysis.TranSpec{Step: 10e-6, Stop: 6e-3}
	small := Objective{Name: "small C", Weight: 1e-3, Eval: func(e *Eval) (float64, error) {
		return e.Circuit.Elements[e.Circuit.Index("C1")].Value / 1e-6, nil
	}}
	p := Problem{Template: tpl,
		Space:      Space{{Element: "C1", Min: 1e-7, Max: 1e-3, Series: "E12"}},
		Objectives: []Objective{TranRipple("out", spec, 3e-3, 6e-3, 0.5), small}}
	st, err := Run(context.Background(), p, Options{Population: 16, Generations: 12, Seed: 2}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	// Ripple ≈ Vpeak/(f·R·C) ⇒ C ≈ 9.3 V/(1 kHz · 1 kΩ · 0.5 V) ≈ 18.6 µF;
	// the smallest E12 value that meets it is 22 µF (18 µF just misses).
	if c := st.Best.Values[0]; c != 22e-6 && c != 18e-6 {
		t.Fatalf("best C = %g (terms %v)", c, st.Best.Terms)
	}
}

func BenchmarkGenerationRC(b *testing.B) {
	p := rcProblem(b)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := Run(context.Background(), p, Options{Population: 50, Generations: 1, Seed: uint64(i)}, nil, nil); err != nil {
			b.Fatal(err)
		}
	}
}
