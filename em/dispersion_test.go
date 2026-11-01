// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"math/cmplx"
	"testing"
)

// TestWaterDebye: limits, the loss peak at ωτ = 1 with ε″ = Δε/2, and
// causality (Kramers–Kronig): ε′(0) − ε∞ = (2/π)∫ε″(ω)/ω dω.
func TestWaterDebye(t *testing.T) {
	w := WaterDebye25
	if err := w.Validate(); err != nil {
		t.Fatal(err)
	}
	if math.Abs(w.Static()-78.6) > 1e-12 || math.Abs(real(w.Eps(1e18))-5.4) > 1e-3 {
		t.Errorf("limits %.6g, %.6g", w.Static(), real(w.Eps(1e18)))
	}
	tau := w.Debye[0].Tau
	e := w.Eps(1 / tau)
	if math.Abs(-imag(e)-73.2/2) > 1e-9 || math.Abs(real(e)-(5.4+73.2/2)) > 1e-9 {
		t.Errorf("at ωτ = 1: %v", e)
	}
	if f := 1 / (2 * math.Pi * tau); math.Abs(f-19.175e9)/19.175e9 > 1e-3 {
		t.Errorf("relaxation frequency %.4g Hz", f)
	}
	// Kramers–Kronig on a log grid (ω from 1e6 to 1e18 rad/s).
	var kk float64
	n := 200000
	lo, hi := math.Log(1e6), math.Log(1e18)
	for i := 0; i < n; i++ {
		x := lo + (hi-lo)*(float64(i)+0.5)/float64(n)
		kk += -imag(w.Eps(math.Exp(x))) * (hi - lo) / float64(n) // ∫ε″/ω dω = ∫ε″ d(ln ω)
	}
	kk *= 2 / math.Pi
	if math.Abs(kk-73.2)/73.2 > 1e-4 {
		t.Errorf("Kramers–Kronig %.6g, want 73.2", kk)
	}
}

// TestDispersiveTimeTwin: driven by a sine, the time-domain polarization
// equations settle to the amplitude and phase ε(ω) gives — the twin of the
// frequency-domain form — for Debye (exact integrator) and Lorentz (RK4).
func TestDispersiveTimeTwin(t *testing.T) {
	m := Dispersive{Name: "test", EpsInf: 2, Debye: []DebyeTerm{{10, 1e-9}},
		Lorentz: []LorentzTerm{{3, 2 * math.Pi * 2e9, 2e8}}}
	if err := m.Validate(); err != nil {
		t.Fatal(err)
	}
	for _, f := range []float64{5e7, 3e8, 1.9e9, 6e9} {
		w := 2 * math.Pi * f
		per := 1 / f
		h := per / 2000
		n := int(200e-9/h) + 4000 // settle 20 Lorentz decay times (2/γ = 10 ns), then two cycles
		d := m.Response(func(t float64) float64 { return math.Sin(w * t) }, h, n)
		// Project the last whole cycles onto sin and cos.
		var s, c float64
		cyc := 2000
		for k := n - 2*cyc; k < n; k++ {
			tk := float64(k) * h
			s += d[k] * math.Sin(w*tk)
			c += d[k] * math.Cos(w*tk)
		}
		s, c = 2*s/float64(2*cyc), 2*c/float64(2*cyc)
		// E = sin ωt = Im e^{jωt} gives D = Im[ε e^{jωt}] = Re ε·sin ωt + Im ε·cos ωt.
		want := m.Eps(w)
		got := complex(s, c)
		if cmplx.Abs(got-want)/cmplx.Abs(want) > 2e-4 {
			t.Errorf("f %.3g: time domain %v, ε(ω) %v", f, got, want)
		}
	}
}

// TestPodolsky: the static potential solves (1 − a²∇²)∇²φ = 0 away from
// the charge, its flux is q/ε0 far away, a → 0 is Coulomb, and the field
// is minus the gradient of the potential.
func TestPodolsky(t *testing.T) {
	const q, a = 1e-9, 1e-3
	k := q / (4 * math.Pi * Eps0)
	lap := func(f func(r float64) float64, r, h float64) float64 { // radial Laplacian (1/r)(rf)″
		g := func(x float64) float64 { return x * f(x) }
		return (g(r+h) - 2*g(r) + g(r-h)) / (h * h) / r
	}
	phi := func(r float64) float64 { return PodolskyPotential(q, r, a) }
	u := func(r float64) float64 { return lap(phi, r, 1e-6) }
	for _, r := range []float64{0.5e-3, 1e-3, 3e-3} {
		// ∇²φ analytic: −k e^(−r/a)/(a² r)
		if want := -k * math.Exp(-r/a) / (a * a * r); math.Abs(u(r)-want)/math.Abs(want) > 1e-4 {
			t.Errorf("r %g: ∇²φ %.6g, want %.6g", r, u(r), want)
		}
		res := u(r) - a*a*lap(u, r, 2e-5)
		if math.Abs(res)/math.Abs(u(r)) > 2e-3 {
			t.Errorf("r %g: (1 − a²∇²)∇²φ = %.3g of ∇²φ", r, res/u(r))
		}
		e := PodolskyField(q, r, a)
		g := -(phi(r+1e-8) - phi(r-1e-8)) / 2e-8
		if math.Abs(e-g)/e > 1e-6 {
			t.Errorf("r %g: E %.8g, −∇φ %.8g", r, e, g)
		}
	}
	if r := 50 * a; math.Abs(PodolskyField(q, r, a)*r*r/k-1) > 1e-15 {
		t.Errorf("far flux not q/ε0")
	}
	if PodolskyPotential(q, 0.01, 0) != k/0.01 {
		t.Errorf("a = 0 is not Coulomb")
	}
	if math.Abs(PodolskyPotential(q, 0, a)-k/a) > 0 || math.Abs(PodolskyField(q, 1e-9, a)-k/(2*a*a))/(k/(2*a*a)) > 1e-5 {
		t.Errorf("centre values")
	}
	if u := PodolskySelfEnergy(q, a); math.Abs(u-0.5*q*PodolskyPotential(q, 0, a))/u > 1e-14 {
		t.Errorf("self-energy ≠ ½qφ(0)")
	}
	mx, ms := PodolskyModes(1, a)
	if mx != SpeedOfLight || math.Abs(ms-SpeedOfLight*math.Sqrt(1+1e6)) > 1e-3 {
		t.Errorf("modes %g %g", mx, ms)
	}
	if b := PodolskyBound(1, 1e-6); math.Abs(PodolskyRelativeDeviation(1, b)-1e-6) > 1e-18 {
		t.Errorf("bound inverse")
	}
	if l := PodolskyLength(1e9); math.Abs(l-1.973269804e-16) > 1e-25 {
		t.Errorf("1 GeV length %g", l)
	}
}
