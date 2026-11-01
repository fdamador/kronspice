// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// The M9 acceptance test: the multi-objective MOSFET amplifier
// (examples/mos_amplifier_pareto) evolves, under NSGA-II, a Pareto front
// trading mid-band gain, bandwidth and supply power, subject to gain ≥ 20 dB.
// Every front netlist written to front/*.cir is re-simulated from scratch:
// it must meet the constraint, reproduce the errors the GA recorded, balance
// its energy ledger, and no member may dominate another. The front must span
// the trade-off: a design above 34 dB, one above 30 MHz and one below 50 µW
// (seeds 1–6 reach 34.5–35.2 dB, 34–36 MHz and 20–34 µW; doc/M9.md).
func TestMOSParetoExample(t *testing.T) {
	if testing.Short() {
		t.Skip("about 6 s: the full 100-generation NSGA-II run")
	}
	ex, err := LoadExperiment(filepath.Join("..", "examples", "mos_amplifier_pareto", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 10, 25, 7, 0, 0, 0, time.UTC))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, nil, rd.Observe)
	if err != nil {
		t.Fatal(err)
	}
	rd.Close()
	if len(st.Front) < 20 {
		t.Fatalf("front has %d designs", len(st.Front))
	}
	spec := analysis.ACSpec{Variation: "dec", Points: 10, FStart: 100, FStop: 1e9}
	type pt struct{ gain, bw, p float64 }
	var pts []pt
	var errs [][]float64
	maxGain, maxBW, minP := 0.0, 0.0, math.Inf(1)
	files, _ := filepath.Glob(filepath.Join(rd.Path, "front", "*.cir"))
	if len(files) != len(st.Front) {
		t.Fatalf("%d front netlists for %d designs", len(files), len(st.Front))
	}
	for i, f := range files {
		src, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		d, err := netlist.ParseString(string(src))
		if err != nil {
			t.Fatalf("%s: %v", f, err)
		}
		c, err := circuit.FromDeck(d)
		if err != nil {
			t.Fatal(err)
		}
		r, err := analysis.AC(c, spec, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		h, _ := analysis.Transfer(r, "o", "0")
		g, _ := analysis.GainAt(r.Freq, h, 1e4)
		_, hi, err := analysis.Cutoff3dB(r.Freq, h)
		if err != nil {
			t.Fatal(err)
		}
		hi = math.Min(hi, 1e9)
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		vdd, _ := op.Element("VDD")
		p := pt{20 * math.Log10(g), hi, -vdd.PowerAbsorbed}
		pts = append(pts, p)
		if p.gain < 20-1e-9 {
			t.Errorf("%s: gain %.3f dB violates the 20 dB constraint", f, p.gain)
		}
		if r.Ledger.MaxRelativeResidual > 1e-9 || op.Ledger.RelativeResidual() > 1e-9 {
			t.Errorf("%s: ledger residual AC %g, OP %g", f, r.Ledger.MaxRelativeResidual, op.Ledger.RelativeResidual())
		}
		// The written netlist reproduces the GA's errors (front files are
		// sorted by the first objective, like front.csv).
		want := sortedFront(st.Front)[i].Terms
		got := []float64{math.Max(0, 40-p.gain) / 20, math.Max(0, math.Log10(1e9/p.bw)), math.Max(0, p.p)}
		for k := range got {
			if math.Abs(got[k]-want[k]) > 1e-9*math.Max(1, math.Abs(want[k])) {
				t.Errorf("%s: objective %d re-simulates to %g, GA recorded %g", f, k, got[k], want[k])
			}
		}
		errs = append(errs, got)
		maxGain, maxBW, minP = math.Max(maxGain, p.gain), math.Max(maxBW, p.bw), math.Min(minP, p.p)
	}
	// No member dominates another in the objectives the GA ranked (the
	// capped shortfalls: gains above the 40 dB ceiling count as equal).
	for a := range errs {
		for b := range errs {
			if a != b && len(FastNonDominatedSort([][]float64{errs[a], errs[b]})) > 1 {
				t.Errorf("front member %d dominates %d: %+v vs %+v", a, b, pts[a], pts[b])
			}
		}
	}
	if maxGain < 34 || maxBW < 30e6 || minP > 50e-6 {
		t.Errorf("front does not span the trade-off: max gain %.2f dB, max bandwidth %.3g Hz, min power %.3g W",
			maxGain, maxBW, minP)
	}
	t.Logf("%d designs on the front after %d generations (%d evaluations): gain up to %.2f dB, bandwidth up to %.3g Hz, power down to %.3g W",
		len(st.Front), st.Generation, st.Evaluations, maxGain, maxBW, minP)
}
