// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// Every case: resonance and Q match the closed forms (to the sweep's
// 0.0275 Hz grid), the engine's terms equal the ordinary-element twin, and
// the AC and transient ledgers close with the terms' power booked as
// model-supplied (README, Results).
func TestCases(t *testing.T) {
	for _, k := range append(append([]Case{}, Cases...), UFTCases...) {
		r, err := Run(k)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.FPeak-r.FWant) > 0.03 || math.Abs(r.Q-r.QWant) > 2e-4*r.QWant {
			t.Errorf("%s: f_r %.3f Hz (want %.3f), Q %.4f (want %.4f)", k.Name, r.FPeak, r.FWant, r.Q, r.QWant)
		}
		if r.TwinDiff > 1e-12 || r.ACResidual > 1e-12 {
			t.Errorf("%s: twin difference %.2g, AC residual %.2g", k.Name, r.TwinDiff, r.ACResidual)
		}
		// The trapezoidal mismatch scales with the energy moved through the
		// circuit. Where the terms supply more than the source (kl0: 1.1 J
		// model-supplied against 0.44 J in, 1.47 J heat), that is the heat,
		// so the tolerance is relative to the larger of source input and heat.
		scale := math.Max(r.EIn, r.EHeat)
		if bal := r.EIn + r.EModel - r.EOut - r.EHeat - r.EStored; math.Abs(bal) > 1e-5*scale || math.Abs(r.Mismatch) > 1e-5*scale {
			t.Errorf("%s: transient books open by %.3g J (mismatch %.3g J)", k.Name, bal, r.Mismatch)
		}
		if k.Form == "" && (r.PModel != 0 || r.EModel != 0) {
			t.Errorf("plain circuit books model-supplied power %.3g W, %.3g J", r.PModel, r.EModel)
		}
	}
}

// ECE and ECE2 give opposite torsion effects for the same ω₀: ECE adds
// 0.25 Ω (Q 8.70, the term absorbs 0.08 W), ECE2 adds −0.5 Ω (Q 21.74,
// the term supplies 1.0 W, model-supplied, while the source gives 1.0 W).
func TestOppositeSigns(t *testing.T) {
	ece, err := Run(Cases[1])
	if err != nil {
		t.Fatal(err)
	}
	ece2, err := Run(Cases[2])
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ece.PModel+0.08) > 1e-3 || math.Abs(ece2.PModel-1.0) > 2e-3 || math.Abs(ece2.PSource-1.0) > 2e-3 {
		t.Errorf("model-supplied ECE %.4f W, ECE2 %.4f W (source %.4f W)", ece.PModel, ece2.PModel, ece2.PSource)
	}
}

// UFT 456's conductivities on the ECE2 torsion case (README, "UFT 456
// conductivities"): kl0 takes the loop from 0.5 Ω to 0.25 Ω (4.000 A),
// k0's 1 mS leak adds G/(ωC)² ≈ 0.118 Ω (1.618 A), both together 2.716 A.
// The model-supplied power is the series terms' ½|I|²·|r| less the leak's
// ½|V_C|²·G, so the leak's loss is booked in the same column, not as heat.
func TestUFT456(t *testing.T) {
	want := []float64{1.6178, 4.0000, 2.7164}
	for n, k := range UFTCases {
		r, err := Run(k)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.IPeak-want[n]) > 1e-4 {
			t.Errorf("%s: peak |I| %.5f A, want %.4f", k.Name, r.IPeak, want[n])
		}
		w := 2 * math.Pi * r.FPeak
		g := k.shunt()
		vc := r.IPeak / math.Hypot(g, w*C)
		pm := 0.5*r.IPeak*r.IPeak*(-k.series()) - 0.5*vc*vc*g
		if math.Abs(r.PModel-pm) > 1e-9*math.Abs(pm) {
			t.Errorf("%s: model-supplied %.6f W, series terms less leak %.6f W", k.Name, r.PModel, pm)
		}
		if heat := 0.5 * r.IPeak * r.IPeak * R; math.Abs(r.PHeat-heat) > 1e-9*heat {
			t.Errorf("%s: heat %.6f W, want ½|I|²R = %.6f W", k.Name, r.PHeat, heat)
		}
	}
}

// The netlists, their KiCad schematics and results sheets go to export/
// (KRONSPICE_EXPORT=1); the sheet is the AC solution at f0.
func TestExport(t *testing.T) {
	var items []export.Item
	for _, k := range append(append([]Case{}, Cases...), UFTCases...) {
		name := "ece_rlc_" + strings.NewReplacer(" (maxwell)", "", " + ", "_", " ", "_").Replace(strings.ToLower(k.Name))
		items = append(items, export.Item{Name: name, Src: k.Netlist(), Opt: sheet.Options{Freq: F0}})
	}
	export.Test(t, "export", items...)
}
