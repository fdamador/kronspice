// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Tewari's void radius from eq. 2.8 (his 4 × 10⁻¹¹ cm, eq. 2.9), and the
// hydrogen contact formula reproduces the Darwin term α⁴mc²/2n³.
func TestTewariAndDarwin(t *testing.T) {
	tw := TewariVoid()
	if math.Abs(tw.ReCm/4.029e-11-1) > 1e-3 {
		t.Errorf("r_e = %g cm", tw.ReCm)
	}
	darwin := math.Pow(Alpha, 4) * MC2Hz / (2 * 8)
	if got := ContactShiftHz(DarwinR2, 2); math.Abs(got/darwin-1) > 1e-12 {
		t.Errorf("contact shift with ¾ƛ² = %g Hz, Darwin %g Hz", got, darwin)
	}
	// LEP2 transfer: 200 GeV / ħc in units of 1/ƛ.
	if q := LEPq(); math.Abs(q/3.914e5-1) > 1e-3 {
		t.Errorf("q ƛ = %g", q)
	}
}

// The measurements on a known shape: a uniform ball of radius R has
// ⟨r²⟩ = 3R²/5 and form factor 3(sin x − x cos x)/x³.
func TestBallMeasures(t *testing.T) {
	g := NewGrid(0.1, 6, 6)
	s := &State{G: g, Psi: make([]complex128, g.NR*g.NZ), dens: make([]float64, g.NR*g.NZ)}
	R := 3.0
	for i, r := range g.Rho {
		for j, z := range g.Z {
			if math.Hypot(r, z) < R {
				s.Psi[g.at(i, j)] = 0
			} else {
				s.Psi[g.at(i, j)] = 1
			}
		}
	}
	rho := s.Charge(Deficit, 0)
	if r2 := s.R2(rho); math.Abs(r2/(0.6*R*R)-1) > 0.01 {
		t.Errorf("⟨r²⟩ = %g, want %g", r2, 0.6*R*R)
	}
	x := 2.0 * R
	want := 3 * (math.Sin(x) - x*math.Cos(x)) / (x * x * x)
	if f := s.FormFactor(rho, 2); math.Abs(f-want) > 0.01 {
		t.Errorf("F = %g, want %g", f, want)
	}
	enc := s.Enclosed(rho, []float64{R / 2, 2 * R})
	if math.Abs(enc[0]-0.125) > 0.01 || math.Abs(enc[1]-1) > 1e-9 {
		t.Errorf("enclosed %v", enc)
	}
}

// RingFrequency reads a damped ringing.
func TestRingFrequency(t *testing.T) {
	var x []float64
	for i := 0; i < 600; i++ {
		tt := float64(i) * 0.1
		x = append(x, math.Exp(-tt/25)*math.Cos(0.38*tt))
	}
	if w, n := RingFrequency(x, 0.1, 30, 21); math.Abs(w/0.38-1) > 0.02 || n < 4 {
		t.Errorf("ω = %g from %d turning points", w, n)
	}
}

// The bubble's charge is many ƛ across, and two voids' interaction dies
// off far faster than Coulomb's.
func TestBubbleCharge(t *testing.T) {
	s := relaxedBall(EpsRef, []int{0}, 5)
	ee, el := s.Energies()
	lbar := EpsRef / (ee + el)
	for _, p := range []Placement{Deficit, Wall, Light} {
		if r2 := s.R2(s.Charge(p, 0)) / (lbar * lbar); r2 < 5 {
			t.Errorf("%v: ⟨r²⟩ = %g ƛ²", p, r2)
		}
	}
	f := s.Profile()
	far := NewState(grid(), EpsRef, V0, []int{0, 1}, Pair(f, 8))
	far.refreshLight(40)
	e2, _ := far.Energies()
	inter := e2 + 2*EpsRef*math.Sqrt(far.Lam[0]) - 2*(ee+el)
	if math.Abs(inter) > 0.01*Alpha*EpsRef/16 {
		t.Errorf("interaction at 2d = 16: %g against Coulomb %g", inter, Alpha*EpsRef/16)
	}
}
