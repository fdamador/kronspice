// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Step 1d: the critical cavity bubble in 3-D (sextic ether).
//
// A spherical cavity a(r) → 0 inside, a → a0 (the ether) outside, is the
// saddle between "no cavity" and "a growing cavity": the analogue of the
// critical bubble of classical nucleation, with energy barrier E*. It
// exists when the cavity phase is favoured, Δ = V(a0) − V(0) > 0:
//
//	ether metastable, no twist:   Δ = −ε (ε < 0)
//	ether under twist k (local):  Δ = V(a0) + κ2k²a0² − ε  (κ2k² beyond the binodal)
//
// Thin wall (Δ → 0), classical nucleation (Gibbs; Blake threshold in
// liquids) and the static bubble of Coleman's problem:
//
//	R* = 2σ/Δ,      E* = 16πσ³/(3Δ²),
//
// with σ the coexistence wall tension (wall.go). Twisted case: the twist of
// n round a real bubble is not uniform (n must pass round it); here it is
// taken as a uniform local term κ2k²a², labelled as that approximation;
// the full 3-D twisted bubble is step 3.
//
// Method: shooting on the radial equation 2κ2(a″ + 2a′/r) = V′(a) from the
// linear interior a = A·sinh(mr)/(mr), m² = V_u(0)/κ2, bisection on ln A
// (overshoot past a0 against turning back), then Newton on the radial
// ladder from that profile for the discrete saddle.

import "math"

// Bubble is a critical cavity.
type Bubble struct {
	R, E     float64 // radius where a² = a0²/2; barrier
	RTW, ETW float64 // thin-wall values
	Delta    float64
	Sigma    float64
	A0       float64
	Neg      int // negative eigenvalues of the radial Hessian (saddle: 1)
	Profile  []float64
	H        float64
	Iter     int
	ShootR   float64
}

// shoot integrates from the centre with amplitude A; returns +1 overshoot,
// −1 undershoot, and the radius reached.
func shoot(e EOS, a0, lnA, h, rmax float64) (int, float64) {
	_, vu0, _ := e.Vu(0)
	m := math.Sqrt(vu0 / kappa2)
	A := math.Exp(lnA)
	r := h
	sh := math.Sinh(m * r)
	a := A * sh / (m * r)
	da := A * (m*r*math.Cosh(m*r) - sh) / (m * r * r)
	f := func(r, a, da float64) (float64, float64) {
		return da, e.D1(a)/(2*kappa2) - 2*da/r
	}
	for r < rmax {
		k1a, k1v := f(r, a, da)
		k2a, k2v := f(r+h/2, a+h/2*k1a, da+h/2*k1v)
		k3a, k3v := f(r+h/2, a+h/2*k2a, da+h/2*k2v)
		k4a, k4v := f(r+h, a+h*k3a, da+h*k3v)
		a += h / 6 * (k1a + 2*k2a + 2*k3a + k4a)
		da += h / 6 * (k1v + 2*k2v + 2*k3v + k4v)
		r += h
		if a > a0 {
			return 1, r
		}
		if da < 0 && a > 1e-3 {
			return -1, r
		}
	}
	return 1, r // still in the cavity at rmax: left too late, as an overshoot
}

// CriticalBubble finds the saddle for e (twist included in e.K).
func CriticalBubble(e EOS, h float64) (Bubble, error) {
	var b Bubble
	a0, ok := e.Ether()
	if !ok {
		return b, errNoConverge
	}
	b.A0 = a0
	b.Delta = e.V(a0) - e.V(0)
	if b.Delta <= 0 {
		return b, errNoConverge
	}
	// Coexistence tension of this EOS: shift ε so that V(0) = V(a0), and
	// integrate 2√(κ2(V − V(a0))) from 0 to a0 (exact for no twist; for a
	// twist the shift is the local approximation's).
	ec := e
	ec.Eps += b.Delta
	b.Sigma = SigmaQuad(ec, 0, a0, ec.V(a0))
	b.RTW = 2 * b.Sigma / b.Delta
	b.ETW = 16 * math.Pi * b.Sigma * b.Sigma * b.Sigma / (3 * b.Delta * b.Delta)
	rmax := 3*b.RTW + 30*e.Xi()
	lo, hi := -700.0, 0.0
	var reach float64
	for it := 0; it < 200; it++ {
		m := (lo + hi) / 2
		s, r := shoot(e, a0, m, h/4, rmax)
		if s > 0 { // left the cavity too late: too little friction lost
			lo = m
		} else {
			hi = m
		}
		reach = r
		if hi-lo < 1e-15 {
			break
		}
	}
	// Rebuild the profile on the ladder grid up to the departure radius.
	b.ShootR = reach
	box := reach + 12*e.Xi()
	n := int(math.Round(box/h)) + 1
	lad := Radial(n, h)
	prof := make([]float64, n)
	{
		_, vu0, _ := e.Vu(0)
		m := math.Sqrt(vu0 / kappa2)
		A := math.Exp(lo)
		r := h / 4
		sh := math.Sinh(m * r)
		a := A * sh / (m * r)
		da := A * (m*r*math.Cosh(m*r) - sh) / (m * r * r)
		f := func(r, a, da float64) (float64, float64) { return da, e.D1(a)/(2*kappa2) - 2*da/r }
		hs := h / 4
		prof[0] = A
		done := false
		for i := 1; i < n; i++ {
			target := float64(i) * h
			for !done && r < target-1e-12 {
				k1a, k1v := f(r, a, da)
				k2a, k2v := f(r+hs/2, a+hs/2*k1a, da+hs/2*k1v)
				k3a, k3v := f(r+hs/2, a+hs/2*k2a, da+hs/2*k2v)
				k4a, k4v := f(r+hs, a+hs*k3a, da+hs*k3v)
				a += hs / 6 * (k1a + 2*k2a + 2*k3a + k4a)
				da += hs / 6 * (k1v + 2*k2v + 2*k3v + k4v)
				r += hs
				if a > a0 || da < 0 && a > 1e-3 || r > reach-2*e.Xi() {
					done = true
				}
			}
			if done {
				prof[i] = a0
			} else {
				prof[i] = a
			}
		}
	}
	prof[n-1] = a0
	it, err := lad.Solve(e, prof, 1e-10, 300)
	if err != nil {
		return b, err
	}
	b.Iter = it
	b.Profile, b.H = prof, h
	b.E = lad.Energy(e, prof, e.V(a0))
	// Radius at a² = a0²/2.
	for i := 0; i+1 < n; i++ {
		u0, u1 := prof[i]*prof[i], prof[i+1]*prof[i+1]
		t := a0 * a0 / 2
		if u0 < t && u1 >= t {
			b.R = (float64(i) + (t-u0)/(u1-u0)) * h
			break
		}
	}
	// Hessian inertia (radial, interior free nodes 0..n−2).
	m := n - 1
	di, off := make([]float64, m), make([]float64, m-1)
	for i := 0; i < m; i++ {
		d := lad.W[i] * e.D2(prof[i])
		if i > 0 {
			d += 2 * kappa2 * lad.A[i-1] / h
		}
		d += 2 * kappa2 * lad.A[i] / h
		di[i] = d
		if i < m-1 {
			off[i] = -2 * kappa2 * lad.A[i] / h
		}
	}
	b.Neg = negatives(di, off)
	return b, nil
}
