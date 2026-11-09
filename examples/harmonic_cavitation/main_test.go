// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The calibration brackets inception correctly and the best speed lies
// inside the cavitating range.
func TestCalibrationOrder(t *testing.T) {
	c := Calibrate(Xue2024)
	if !(c.SigmaNone > c.SigmaOn && c.SigmaOn > c.SigmaBest && c.SigmaBest > 0) {
		t.Fatalf("σ order wrong: none %.3f, on %.3f, best %.3f", c.SigmaNone, c.SigmaOn, c.SigmaBest)
	}
	// Surface speeds: 3000 rpm at R = 40 mm is 12.57 m/s.
	if u := omega(3000) * 0.040; math.Abs(u-12.566) > 1e-3 {
		t.Fatalf("surface speed %.4f", u)
	}
}

// The resonance rule applied to Xue's own rotor recovers its groove count
// to within 10 % (the dominant frequency sits 6.5 % below groove passing).
func TestResonanceRuleSelfConsistent(t *testing.T) {
	c := Calibrate(Xue2024)
	n := 2 * math.Pi * Xue2024.R * c.StLock / Xue2024.W
	if math.Abs(n/float64(Xue2024.Cavities)-1) > 0.10 {
		t.Fatalf("rule gives N = %.2f against %d", n, Xue2024.Cavities)
	}
}

// Wood: pure liquid gives c_l; the bubbly minimum near α = 0.5 is
// √(p/(ρ α(1−α))) ≈ √(4p/ρ) when the liquid term is negligible.
func TestWood(t *testing.T) {
	if c := WoodSpeed(0, pAtm); math.Abs(c/cLiq-1) > 1e-12 {
		t.Fatalf("α = 0: %.3f", c)
	}
	want := math.Sqrt(4 * pAtm / rho)
	if c := WoodSpeed(0.5, pAtm); math.Abs(c/want-1) > 0.01 {
		t.Fatalf("α = 0.5: %.3f against %.3f", c, want)
	}
	a := SonicVoid(30, pAtm)
	if math.Abs(WoodSpeed(a, pAtm)-30) > 1e-6 {
		t.Fatalf("SonicVoid inverse failed: α %.5f", a)
	}
	if !math.IsNaN(SonicVoid(5, pAtm)) {
		t.Fatal("5 m/s is below the mixture minimum and must give NaN")
	}
}
