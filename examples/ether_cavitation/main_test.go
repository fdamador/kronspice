// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The potential's derivatives against central differences.
func TestEOSDerivatives(t *testing.T) {
	for _, e := range []EOS{Q(1.7), S(0.8, 0.13), S(1.2, -0.07).Twist(0.3), Q(1).Twist(0.2)} {
		for _, a := range []float64{-0.9, 0.05, 0.4, 0.77, 1.1} {
			const d = 1e-5
			d1 := (e.V(a+d) - e.V(a-d)) / (2 * d)
			d2 := (e.D1(a+d) - e.D1(a-d)) / (2 * d)
			if math.Abs(d1-e.D1(a)) > 1e-8 || math.Abs(d2-e.D2(a)) > 1e-7 {
				t.Errorf("%c a=%.2f: D1 %.9f vs %.9f, D2 %.9f vs %.9f", e.Kind, a, e.D1(a), d1, e.D2(a), d2)
			}
		}
	}
	if v := S(1, 0.3).V(0); math.Abs(v-0.3) > 1e-15 {
		t.Errorf("S: V(0) = %g, want ε", v)
	}
}

// The ladder's residual is minus the gradient of its energy, so its
// operating points are the energy's stationary points.
func TestLadderGradient(t *testing.T) {
	e := S(1, 0.05).Twist(0.2)
	for _, l := range []Ladder{Planar(12, 0.3), Radial(12, 0.3)} {
		a := make([]float64, 12)
		for i := range a {
			a[i] = 0.1 + 0.07*float64(i) + 0.03*math.Sin(float64(i))
		}
		f := l.Residual(e, a)
		for i := 1; i < 11; i++ {
			const d = 1e-6
			a[i] += d
			ep := l.Energy(e, a, 0)
			a[i] -= 2 * d
			em := l.Energy(e, a, 0)
			a[i] += d
			g := (ep - em) / (2 * d) / l.Scale
			if math.Abs(g+f[i]) > 1e-7 {
				t.Errorf("node %d: dE/da %.9f, −F %.9f", i, g, -f[i])
			}
		}
	}
}

// Wall tensions and profiles against the closed forms.
func TestWalls(t *testing.T) {
	w, err := QCavityWall(1, 20, 401)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(w.Sigma/w.SigmaExact-1) > 5e-5 || w.ProfileErr > 1e-4 {
		t.Errorf("Q cavity wall: σ %.7f vs %.7f, profile %.1e", w.Sigma, w.SigmaExact, w.ProfileErr)
	}
	s, _, err := SCoexistenceWall(1, 12, 481)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(s.Sigma/s.SigmaExact-1) > 1e-4 || s.ProfileErr > 2e-4 {
		t.Errorf("S wall: σ %.7f vs %.7f, profile %.1e", s.Sigma, s.SigmaExact, s.ProfileErr)
	}
	if q := SigmaQuad(Q(2), 0, 1, 0); math.Abs(q-2.0/3*math.Sqrt(2)) > 1e-10 {
		t.Errorf("quadrature σ_Q(λ=2) = %.12f", q)
	}
}

// Tensile strength: Eckhaus criterion = maximum of the twist stress;
// Q: k_c = 1/(√6 ξ), j_c = 4√(κ2λ)/(3√6).
func TestTensileStrength(t *testing.T) {
	kc, kj, jc := Spinodal(Q(1))
	if math.Abs(kc-1/math.Sqrt(6)) > 1e-9 || math.Abs(kj-kc) > 1e-4 || math.Abs(jc-4/(3*math.Sqrt(6))) > 1e-8 {
		t.Errorf("Q: k_c %.9f, k_j %.6f, j_c %.9f", kc, kj, jc)
	}
	for _, eps := range []float64{0.05, 0.2} {
		e := S(1, eps)
		kc, kj, _ := Spinodal(e)
		if math.Abs(kj/kc-1) > 1e-3 {
			t.Errorf("S ε=%.2f: Eckhaus k_c %.5f vs max-j %.5f", eps, kc, kj)
		}
		if kb := Binodal(e); kb >= kc || math.Abs(kb/math.Sqrt(eps)-1) > 0.03 {
			t.Errorf("S ε=%.2f: binodal %.4f, spinodal %.4f", eps, kb, kc)
		}
	}
}

// Relaxation of a twisted ring: below k_c it holds; above, the amplitude
// reaches zero (a phase slip, the 1-D cavity) and twist is released. An S²
// ether unwinds over the pole without a cavity.
func TestRing(t *testing.T) {
	if testing.Short() {
		t.Skip("ring relaxation")
	}
	hold := Ring(Q(1), 400, 0.5, 12, 2, 1e-3, 800, 1)
	slip := Ring(Q(1), 400, 0.5, 14, 2, 1e-3, 800, 1)
	if hold.W1 != 12 || hold.MinA < 0.8 {
		t.Errorf("k = 0.92 k_c: winding %d, min a %.3f", hold.W1, hold.MinA)
	}
	if slip.W1 >= 14 || slip.MinA > 0.05 {
		t.Errorf("k = 1.08 k_c: winding %d, min a %.3f", slip.W1, slip.MinA)
	}
	esc := Ring(Q(1), 400, 0.5, 6, 3, 1e-3, 3000, 1)
	if esc.W1 != 0 || esc.MinA < 0.9 {
		t.Errorf("S² escape: winding %d, min a %.3f", esc.W1, esc.MinA)
	}
	for _, r := range []RingRun{hold, slip, esc} {
		if r.LedgerRel > 1e-3 {
			t.Errorf("ledger %.1e", r.LedgerRel)
		}
	}
}

// Reversal slab: Q tears continuously (L_t = L_c ≈ 4.93 ξ); S tears with
// hysteresis (L_t < L_c).
func TestReversal(t *testing.T) {
	q, err := CriticalLengths(Q(1), 0.02)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(q.Lt/q.Lc-1) > 1e-3 || math.Abs(q.Lc-4.934) > 0.01 {
		t.Errorf("Q: L_c %.4f, L_t %.4f", q.Lc, q.Lt)
	}
	s, err := CriticalLengths(S(1, 0.2), 0.02)
	if err != nil {
		t.Fatal(err)
	}
	if s.Lt > 0.5*s.Lc {
		t.Errorf("S ε=0.2: L_t %.4f, L_c %.4f, want hysteresis", s.Lt, s.Lc)
	}
}

// Critical bubble: a saddle with exactly one unstable mode, tending to the
// thin-wall R* = 2σ/Δ and E* = 16πσ³/3Δ² as Δ → 0.
func TestCriticalBubble(t *testing.T) {
	prevR, prevE := 0.0, 0.0
	for _, eps := range []float64{-0.1, -0.05, -0.02} {
		b, err := CriticalBubble(S(1, eps), 0.05)
		if err != nil {
			t.Fatal(err)
		}
		rr, re := b.R/b.RTW, b.E/b.ETW
		if b.Neg != 1 || rr <= prevR || re <= prevE || rr > 1 || re > 1 {
			t.Errorf("ε=%.2f: neg %d, R/R_tw %.4f, E/E_tw %.4f", eps, b.Neg, rr, re)
		}
		prevR, prevE = rr, re
	}
	if prevR < 0.98 || prevE < 0.92 {
		t.Errorf("ε=−0.02: R/R_tw %.4f, E/E_tw %.4f", prevR, prevE)
	}
}

// Kron twin: the kronspice circuit's operating point equals the ladder
// solution (a minimum and a saddle), its ledger balances.
func TestKronTwin(t *testing.T) {
	e := S(1, -0.1)
	b, err := CriticalBubble(e, 0.1)
	if err != nil {
		t.Fatal(err)
	}
	k, err := KronTwin(Radial(len(b.Profile), b.H), e, b.Profile, e.V(b.A0))
	if err != nil {
		t.Fatal(err)
	}
	if k.MaxDiff > 2e-5 || math.Abs(k.Energy/b.E-1) > 1e-6 || k.LedgerResidual > 1e-9 {
		t.Errorf("bubble twin: max diff %.1e, E %.6f vs %.6f, ledger %.1e", k.MaxDiff, k.Energy, b.E, k.LedgerResidual)
	}
}
