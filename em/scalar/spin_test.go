// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// tm110 is a box with the TM110 mode set in E_z.
func tm110(t *testing.T, spin *SpinConnection) (*Region, float64) {
	t.Helper()
	N, D := [3]int{17, 15, 9}, [3]float64{0.01, 0.01, 0.01}
	r, err := New(Config{N: N, D: D, Mode: Standard, Spin: spin})
	if err != nil {
		t.Fatal(err)
	}
	lx, ly := float64(N[0]-1)*D[0], float64(N[1]-1)*D[1]
	r.SetE(func(x, y, z float64) [3]float64 {
		return [3]float64{0, 0, 100 * math.Sin(math.Pi*x/lx) * math.Sin(math.Pi*y/ly)}
	})
	w := em.SpeedOfLight * math.Hypot(math.Pi/lx, math.Pi/ly)
	return r, w
}

// Gate 1: a spin connection declared zero is Maxwell's (em.FDTD's) scheme
// bit for bit, with no spin work.
func TestSpinZeroIsMaxwell(t *testing.T) {
	a, _ := tm110(t, nil)
	b, _ := tm110(t, &SpinConnection{
		Kappa0: func(x, y, z float64) float64 { return 0 },
		Kappa:  func(x, y, z float64) [3]float64 { return [3]float64{} },
	})
	cur := func(t float64) float64 { return 0.3 * math.Sin(5e9*t) }
	_ = a.AddSource(EdgeSource{Comp: 1, I: 5, J: 6, K: 4, Current: cur})
	_ = b.AddSource(EdgeSource{Comp: 1, I: 5, J: 6, K: 4, Current: cur})
	for s := 0; s < 300; s++ {
		a.Step()
		b.Step()
	}
	for c := 0; c < 3; c++ {
		for n := range a.E[c] {
			if a.E[c][n] != b.E[c][n] || a.B[c][n] != b.B[c][n] {
				t.Fatalf("κ = 0 differs from Maxwell at component %d, index %d", c, n)
			}
		}
	}
	if w := b.Ledger().SpinWork; w != 0 {
		t.Errorf("κ = 0 but spin work %.3g J", w)
	}
}

// The scalar part κ₀ is a conductivity σ = 2κ₀/(μ₀c): the mode's energy
// changes as exp(−2κ₀c t), falling for κ₀ > 0 and growing for κ₀ < 0 (the
// "energy from spacetime" reading), and every joule of the change is the
// spin connection's work, booked in SpinWork.
func TestSpinScalarConductivity(t *testing.T) {
	for _, sgn := range []float64{1, -1} {
		_, w := tm110(t, nil)
		k0 := sgn * 0.01 * w / em.SpeedOfLight // σ/ε₀ = 2κ₀c = 0.02 ω
		r, _ := tm110(t, &SpinConnection{Kappa0: func(x, y, z float64) float64 { return k0 }})
		period := 2 * math.Pi / w
		r.Run(0) // start the ledger
		r.Step()
		w0 := r.Ledger().W
		t0 := r.Time()
		// sample over a whole period at the end to average the E/B exchange
		r.Run(t0 + 15*period)
		var sum float64
		var n int
		var worst float64
		for r.Time() < t0+16*period {
			r.Step()
			l := r.Ledger()
			sum += l.W * math.Exp(2*k0*em.SpeedOfLight*(r.Time()-t0))
			n++
			worst = math.Max(worst, math.Abs(l.Residual())/l.W)
		}
		got := sum / float64(n) / w0 // W·e^{+2κ₀ct}/W0 ≈ 1 if the rate is right
		l := r.Ledger()
		t.Logf("κ₀ = %+.3g 1/m: W/W0 = %.4g after %.0f periods, rate-corrected %.4f; spin work %+.4g J of ΔW %+.4g J; residual %.2g",
			k0, l.W/l.W0, (r.Time()-t0)/period, got, l.SpinWork, l.W-l.W0, worst)
		if math.Abs(got-1) > 0.02 {
			t.Errorf("κ₀ = %g: energy rate off, corrected ratio %.4f", k0, got)
		}
		if worst > 1e-12 {
			t.Errorf("κ₀ = %g: ledger residual %.3g", k0, worst)
		}
		if math.Signbit(l.SpinWork) != (sgn > 0) {
			t.Errorf("κ₀ = %g: spin work has the wrong sign (%.3g)", k0, l.SpinWork)
		}
		if d := r.CoulombDefect(); d != 0 {
			t.Errorf("κ = 0: Coulomb defect %.3g, want 0", d)
		}
	}
}

// The vector part κ: the Hall-like current (2/μ₀)κ×B. With κ ∥ E (the
// slides' free-field case, κ ⊥ B) it starts doing no work; the field it
// drives across the mode later exchanges energy with it, all booked. ECE2's
// Coulomb law ∇·E = −2κ·E does not follow from its Ampère law: for the TM110
// mode (∇·E = 0) with κ ∥ z it fails at once by O(1).
func TestSpinVectorPart(t *testing.T) {
	_, w := tm110(t, nil)
	kz := 0.02 * w / em.SpeedOfLight
	r, _ := tm110(t, &SpinConnection{Kappa: func(x, y, z float64) [3]float64 { return [3]float64{0, 0, kz} }})
	worst := 0.0
	for s := 0; s < 2000; s++ {
		r.Step()
		l := r.Ledger()
		worst = math.Max(worst, math.Abs(l.Residual())/l.W0)
	}
	l := r.Ledger()
	d := r.CoulombDefect()
	t.Logf("κ_z = %.3g 1/m: W/W0 = %.4f, spin work %+.3g of W0 %.3g J, residual %.2g; Coulomb defect %.3g",
		kz, l.W/l.W0, l.SpinWork, l.W0, worst, d)
	if worst > 1e-12 {
		t.Errorf("ledger residual %.3g", worst)
	}
	if d < 0.5 {
		t.Errorf("Coulomb defect %.3g: ECE2's ∇·E = −2κ·E should fail by O(1) here", d)
	}
}
