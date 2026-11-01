// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"os"
	"testing"
)

// Eq. (9) with the paper's values gives about 30 rad/s.
func TestResonance(t *testing.T) {
	if w := Paper().OmegaRes(); math.Abs(w-29.84) > 0.01 {
		t.Errorf("ω_res = %.3f", w)
	}
}

// The paper's equations: the energy that appears in R and the coil beyond
// the drive's work is the work of the omitted back torque.
func TestPaperMissingEnergy(t *testing.T) {
	p := Paper()
	tr := RunPaper(p, 0, 0, p.Alpha, 0.6, 1e-6, 1e4)
	if !tr.Blowup {
		t.Fatal("expected the printed equations to run away")
	}
	missing := tr.EOhm + tr.EL - (tr.EDrive - tr.ERot)
	if math.Abs(missing-tr.EBackTorque) > 1e-3*tr.EBackTorque {
		t.Errorf("missing %.6g J, back-torque work %.6g J", missing, tr.EBackTorque)
	}
}

// With the back torque (kronspice), the current stays bounded, the speed
// stays near ω_res, and the books close.
func TestCircuitBounded(t *testing.T) {
	p := Paper()
	tr, mm, err := p.RunCircuit(p.L, p.Alpha, 0, 3)
	if err != nil {
		t.Fatal(err)
	}
	for k := range tr.I {
		if math.Abs(tr.I[k]) > 100 {
			t.Fatalf("I = %g A at %g s", tr.I[k], tr.T[k])
		}
	}
	if math.Abs(tr.W[len(tr.W)-1]-p.OmegaRes()) > 3 {
		t.Errorf("ω_end %g, ω_res %g", tr.W[len(tr.W)-1], p.OmegaRes())
	}
	other := tr.EDrive - tr.EOhm - tr.EL - tr.ERot
	if math.Abs(other) > 1e-3*tr.EDrive || math.Abs(mm) > 1e-4*tr.EDrive {
		t.Errorf("books: unaccounted %g J, ledger mismatch %g J of %g J", other, mm, tr.EDrive)
	}
}

// Second run: the resistor receives what the flywheel and drive gave.
func TestFlywheelDischarge(t *testing.T) {
	p := Paper()
	for _, a := range []float64{0, 1} {
		tr, _, err := p.RunCircuit(p.L, a, 40, 0.5)
		if err != nil {
			t.Fatal(err)
		}
		give := -tr.ERot + tr.EDrive - tr.EL
		if math.Abs(tr.EOhm-give) > 1e-3*give {
			t.Errorf("α=%g: R %g J, shaft gave %g J", a, tr.EOhm, give)
		}
	}
}

// TestDumpChart writes the comparison series when ECKARDT_DUMP is set.
func TestDumpChart(t *testing.T) {
	if os.Getenv("ECKARDT_DUMP") == "" {
		t.Skip()
	}
	p := Paper()
	pp := RunPaper(p, 0, 0, p.Alpha, 0.6, 1e-6, 1e6)
	tr, _, _ := p.RunCircuit(p.L, p.Alpha, 0, 3)
	for _, s := range []float64{0.30, 0.35, 0.40, 0.42, 0.44, 0.45, 0.46} {
		i1, w1 := at(pp, s)
		i2, w2 := at(tr, s)
		fmt.Printf("%.2f %g %g %g %g\n", s, i1, w1, i2, w2)
	}
	fmt.Println("--")
	for s := 0.0; s <= 3.0001; s += 0.1 {
		i2, w2 := at(tr, s)
		fmt.Printf("%.1f %.4f %.4f\n", s, i2, w2)
	}
}
