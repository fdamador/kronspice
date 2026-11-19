// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// The M7 acceptance test, run as `cmd/evolve` runs it: the 2N3904
// common-emitter example evolves, from E-series values only, to the spec
// 30 ± 0.5 dB mid-band gain, −3 dB edges at 100 Hz and 200 kHz within 5 %
// and V_CE ≥ 3 V, checked on the written best netlist re-simulated from
// scratch. The tolerances are those of the 5 % parts the design is built
// from: with E24/E12 values in ~4 % steps and five coupled targets the
// weighted error settles near 0.02, the floor of the grid (single-part
// moves from the GA's best do not improve it), so the spec, not the
// fitness value, is what the test holds. Over seeds 1–12, 10 meet it; the
// other two miss the gain by under 0.1 dB beyond 0.5 dB (doc/M7.md).
func TestAmplifierExample(t *testing.T) {
	if testing.Short() {
		t.Skip("about 25 s: the full 300-generation run")
	}
	ex, err := LoadExperiment(filepath.Join("..", "examples", "amplifier_gain", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	checkAmplifier(t, ex, true)
}

// Seeds 3 and 10 miss the gain spec by under 0.1 dB with the GA alone
// (doc/M7.md). M7 guessed a multi-part local search would close the gap;
// measured, it does not: seed 3's best is a strict local optimum for every
// move of up to three parts by one E-series step, and of any pair by up to
// four steps — it sits in another bias basin (I_C 4.2 mA against seed 1's
// 1.1 mA). Run every 10 or 25 generations (Lamarckian) the search lowers
// the error (seed 3 0.058 → 0.044, seed 10 0.048 → 0.037) but not to seed
// 1's 0.020. The test records those numbers; it holds only that the search
// never makes the result worse.
func TestAmplifierLocalSearch(t *testing.T) {
	if testing.Short() {
		t.Skip("four full 300-generation runs")
	}
	for _, c := range []struct {
		seed  uint64
		every int
	}{{3, 10}, {10, 25}} {
		base, err := LoadExperiment(filepath.Join("..", "examples", "amplifier_gain", "config.json"))
		if err != nil {
			t.Fatal(err)
		}
		base.Options.Seed = c.seed
		ga, err := Run(context.Background(), base.Problem, base.Options, nil, nil)
		if err != nil {
			t.Fatal(err)
		}
		ex, _ := LoadExperiment(filepath.Join("..", "examples", "amplifier_gain", "config.json"))
		ex.Options.Seed = c.seed
		ex.Options.LocalSearch = &LocalSearch{Size: 2, Every: c.every}
		t.Run(fmt.Sprint("seed", c.seed), func(t *testing.T) {
			ls := checkAmplifier(t, ex, false)
			t.Logf("GA alone %.5g, with local search every %d: %.5g", ga.Best.Fitness, c.every, ls.Best.Fitness)
			if ls.Best.Fitness > ga.Best.Fitness {
				t.Errorf("local search made the result worse")
			}
		})
	}
}

// checkAmplifier runs the amplifier experiment as cmd/evolve does and
// checks the written best netlist against the engineering spec (errors
// reported when strict, logged otherwise).
func checkAmplifier(t *testing.T, ex *Experiment, strict bool) *State {
	t.Helper()
	fail := t.Errorf
	if !strict {
		fail = t.Logf
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 10, 24, 22, 0, 0, 0, time.UTC))
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
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 20, FStart: 10, FStop: 1e8}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(r, "o", "0")
	g, _ := analysis.GainAt(r.Freq, h, 1e4)
	lo, hi, err := analysis.Cutoff3dB(r.Freq, h)
	if err != nil {
		t.Fatal(err)
	}
	if db := 20 * math.Log10(g); math.Abs(db-30) > 0.5 {
		fail("gain at 10 kHz %.3f dB, want 30 ± 0.5", db)
	}
	if math.Abs(hi-2e5)/2e5 > 0.05 || math.Abs(lo-100)/100 > 0.05 {
		fail("band edges %.4g Hz and %.4g Hz, want 100 Hz and 200 kHz within 5 %%", lo, hi)
	}
	if r.Ledger.MaxRelativeResidual > 1e-9 {
		t.Errorf("AC ledger residual %g", r.Ledger.MaxRelativeResidual)
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	q, _ := op.Element("Q1")
	if q.Voltage < 3 {
		fail("V_CE %g V, want ≥ 3 V", q.Voltage)
	}
	t.Logf("best after %d generations (%d evaluations): fitness %.5g, gain %.3f dB, edges %.4g Hz / %.4g Hz, V_CE %s, I_C %s",
		st.Generation, st.Evaluations, st.Best.Fitness, 20*math.Log10(g), lo, hi, output.Eng(q.Voltage, "V"), output.Eng(q.Current, "A"))
	return st
}
