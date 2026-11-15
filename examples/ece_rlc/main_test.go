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
	for _, k := range Cases {
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
		if bal := r.EIn + r.EModel - r.EOut - r.EHeat - r.EStored; math.Abs(bal) > 1e-5*r.EIn || math.Abs(r.Mismatch) > 1e-5*r.EIn {
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

// The five netlists, their KiCad schematics and results sheets go to
// export/ (KRONSPICE_EXPORT=1); the sheet is the AC solution at f0.
func TestExport(t *testing.T) {
	var items []export.Item
	for _, k := range Cases {
		name := "ece_rlc_" + strings.NewReplacer(" (maxwell)", "", " ", "_").Replace(strings.ToLower(k.Name))
		items = append(items, export.Item{Name: name, Src: k.Netlist(), Opt: sheet.Options{Freq: F0}})
	}
	export.Test(t, "export", items...)
}
