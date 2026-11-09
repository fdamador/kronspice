// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Step 1a: the ether's surface tension, and step 1c: the reversal slab.
//
// Walls. With n uniform the amplitude obeys 2κ2 a″ = V′(a), first integral
// κ2 a′² = V(a) − V(a∞) for a profile reaching a stationary point, so a wall
// between a = 0 and a = 1 has tension σ = ∫₀¹ 2√(κ2 V) da (Bogomolny):
//
//	Q, cavity wall (a forced to 0 on a sheet, e.g. a defect plane):
//	    a = tanh(x/2ξ),            σ_Q = (2/3)√(κ2λ)
//	S, coexistence (ε = 0), ether | cavity phase:
//	    a² = 1/(1 + e^(−2x/ξ)),    σ_S = (1/2)√(κ2λ)
//
// In Q the cavity is not a phase (V(0) = λ/4 above the ether), so a cavity
// exists only where something forces a = 0; in S at ε = 0 the two phases
// coexist and the wall is free.
//
// Reversal slab. n is held at +n̂ at x = 0 and at −n̂ at x = L. The ether can
// reverse by turning n through half a great circle (bending energy
// ≈ κ2π²/L with a ≈ 1) or by passing φ straight through zero (a cavity
// sheet, energy ≈ 2σ for L ≫ ξ). The straight (cavity) solution φ = f(x)n̂
// is a stationary point for every L; it is stable against turning out of
// the line (perturbation η ⊥ n̂, quadratic form κ2η′² + V_u(f²)η²) only for
// L < L_c. So a reversal squeezed into less than L_c tears the ether: the
// field must pass through zero. Same problem for an S² ether (n ∈ S²): with
// antipodal ends every half great circle is equivalent, so turning out of
// the plane does not help, and L_c is the same.

import "math"

// WallResult compares a numerical wall with its closed form.
type WallResult struct {
	Sigma, SigmaExact float64
	ProfileErr        float64 // max |a − a_exact|
	Iter              int
}

// QCavityWall solves the Q half-kink a(0) = 0, a(L) = 1 on n nodes.
func QCavityWall(lambda, L float64, n int) (WallResult, error) {
	e := Q(lambda)
	xi := e.Xi()
	h := L / float64(n-1)
	lad := Planar(n, h)
	a := make([]float64, n)
	for i := range a {
		a[i] = math.Tanh(float64(i) * h / (2 * xi) * 1.3) // deliberately off
	}
	a[0], a[n-1] = 0, 1
	it, err := lad.Solve(e, a, 1e-11, 200)
	if err != nil {
		return WallResult{}, err
	}
	var pe float64
	for i := range a {
		pe = math.Max(pe, math.Abs(a[i]-math.Tanh(float64(i)*h/(2*xi))))
	}
	return WallResult{Sigma: lad.Energy(e, a, 0), SigmaExact: 2.0 / 3 * math.Sqrt(kappa2*lambda), ProfileErr: pe, Iter: it}, nil
}

// SCoexistenceWall solves the S wall at ε = 0 on [−L, L].
func SCoexistenceWall(lambda, L float64, n int) (WallResult, float64, error) {
	e := S(lambda, 0)
	xi := e.Xi()
	h := 2 * L / float64(n-1)
	lad := Planar(n, h)
	a := make([]float64, n)
	exact := func(x float64) float64 { return math.Sqrt(1 / (1 + math.Exp(-2*x/xi))) }
	for i := range a {
		a[i] = exact((-L + float64(i)*h) * 0.8)
	}
	a[0], a[n-1] = 0, 1
	it, err := lad.Solve(e, a, 1e-11, 200)
	if err != nil {
		return WallResult{}, 0, err
	}
	// Centre: where a² = 1/2; compare the profile about it (the wall is
	// free to slide, a zero mode).
	var x0 float64
	for i := 0; i+1 < n; i++ {
		u0, u1 := a[i]*a[i], a[i+1]*a[i+1]
		if u0 < 0.5 && u1 >= 0.5 {
			x0 = -L + (float64(i)+(0.5-u0)/(u1-u0))*h
		}
	}
	var pe float64
	for i := range a {
		pe = math.Max(pe, math.Abs(a[i]-exact(-L+float64(i)*h-x0)))
	}
	return WallResult{Sigma: lad.Energy(e, a, 0), SigmaExact: 0.5 * math.Sqrt(kappa2*lambda), ProfileErr: pe, Iter: it}, x0, nil
}

// SigmaQuad is σ = ∫ 2√(κ2(V − Vmin)) da between a1 and a2 by Gauss–
// Legendre on 200 panels: the closed-form check for any equation of state.
func SigmaQuad(e EOS, a1, a2, vref float64) float64 {
	gx := [4]float64{-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526}
	gw := [4]float64{0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538}
	const np = 200
	h := (a2 - a1) / np
	var s float64
	for p := 0; p < np; p++ {
		c := a1 + (float64(p)+0.5)*h
		for q := range gx {
			v := e.V(c+gx[q]*h/2) - vref
			if v > 0 {
				s += gw[q] * h / 2 * 2 * math.Sqrt(kappa2*v)
			}
		}
	}
	return s
}
