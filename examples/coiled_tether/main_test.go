// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestMotionalEMF: in a uniform field the EMF between a wire's ends is
// (v×B)·(r_end − r_start), whatever the path: vBL for the straight tether,
// vB × the coil's axial length for either coil, zero out and back.
func TestMotionalEMF(t *testing.T) {
	s, c, b, fo := MotionalEMFs()
	ord, _ := Coils()
	e := OrbitSpeed * FieldB
	if math.Abs(s-e*FiveMiles)/(e*FiveMiles) > 1e-12 {
		t.Errorf("straight %.6g V", s)
	}
	// Coil ends: z from 0 to N·pitch (ordinary); bifilar from 0 to N·pitch too.
	if math.Abs(c-e*ord.Length())/(e*ord.Length()) > 1e-9 {
		t.Errorf("coiled %.6g V, want %.6g", c, e*ord.Length())
	}
	if math.Abs(b-e*ord.Length())/(e*ord.Length()) > 2e-3 {
		t.Errorf("bifilar %.6g V, want %.6g", b, e*ord.Length())
	}
	if math.Abs(fo) > 1e-9 {
		t.Errorf("folded %.3g V, want 0", fo)
	}
}

// TestResonance: at resonance the coil voltage is Q × the drive, the Sun's
// input equals heat (+ load) at steady state, the whole-run ledger closes,
// and the bifilar coil holds Tesla's N²/4 times the ordinary coil's
// turn-to-turn energy (US 512,340, p. 2) while the total stored energy is
// the same.
func TestResonance(t *testing.T) {
	ord, bif := Coils()
	cb := bif.SelfCapacitance()
	ro, err := RunResonant("ord", ord, cb, 0)
	if err != nil {
		t.Fatal(err)
	}
	rb, err := RunResonant("bif", bif, cb, 0)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range []Resonant{ro, rb} {
		if g := r.VCoil / r.VDrive; math.Abs(g-r.Q)/r.Q > 0.01 {
			t.Errorf("%s: step-up %.2f, Q %.2f", r.Name, g, r.Q)
		}
		if math.Abs(r.PIn-SolarPower())/SolarPower() > 0.01 || math.Abs(r.PIn-r.PHeat)/r.PIn > 0.01 {
			t.Errorf("%s: in %.1f W, heat %.1f W, array %.1f W", r.Name, r.PIn, r.PHeat, SolarPower())
		}
		if math.Abs(r.EIn-r.EDiss-r.EStored) > 1e-3*r.EIn {
			t.Errorf("%s: ledger %.4g − %.4g − %.4g", r.Name, r.EIn, r.EDiss, r.EStored)
		}
	}
	n := float64(ord.Turns)
	if ratio := rb.WTurns / ro.WTurns; math.Abs(ratio-n*n/4)/(n*n/4) > 1e-9 {
		t.Errorf("turn-energy ratio %.6g, want N²/4 = %.6g", ratio, n*n/4)
	}
	if math.Abs(rb.WTotal-ro.WTotal)/ro.WTotal > 1e-6 {
		t.Errorf("total stored %.6g vs %.6g", rb.WTotal, ro.WTotal)
	}
	// A matched load takes half the Sun's input, never more.
	rl, err := RunResonant("load", bif, cb, bif.Resistance(rb.F))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(rl.PLoad-0.5*rl.PIn)/rl.PIn > 0.01 {
		t.Errorf("load %.1f W of %.1f W", rl.PLoad, rl.PIn)
	}
}

func TestRuns(t *testing.T) { main() }
