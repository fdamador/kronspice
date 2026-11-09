// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The soft ether: φ = a·n with an amplitude a ≥ 0 and an equation of state
// V(a). Units: κ2 = 1 (the ether's stiffness, as in examples/ether_knot), so
// lengths are in ℓ = √(κ4/κ2) and energies in E0 = √(κ2κ4) once κ4 = 1 too;
// the amplitude stiffness is λ = β·E0/ℓ³ with β = λκ4/κ2², and the healing
// length is ξ = √(κ2/λ) = ℓ/√β. In 1-D (n uniform, or n twisted in a plane
// at a fixed rate k) κ4 does not enter.
//
// Both potentials are written in u = a²:
//
//	Q (quartic):  V = (λ/4)(u − 1)²                 a = 0 is a maximum
//	S (sextic):   V = λ u(1 − u)² + ε (1 − u)²(1 + 2u)
//
// S has stationary points at u = 0 and u = 1 for every ε; h(u) = (1−u)²(1+2u)
// is the cubic with h(0) = 1, h(1) = 0, h′(0) = h′(1) = 0, so ε moves V(0) =
// ε without moving either minimum. ε > 0: the cavity phase a = 0 costs ε per
// volume (the ether is the stable phase, the cavity metastable, as vapour
// in a liquid below its boiling point); ε < 0: the ether itself is
// metastable (false vacuum, Coleman 1977). a = 1 stays a minimum for
// ε > −λ/3.
//
// A twist of n at rate k in a plane adds κ2 k² a² to the energy density
// (|∂φ|² = (∂a)² + a²(∂n)², (∂n)² = k²): the bending of the ether is paid
// only where the ether is (a > 0). This is the ether's tension: it favours
// the cavity exactly as negative pressure favours vapour.

import "math"

const kappa2 = 1.0

// EOS is an equation of state with an optional uniform twist K (rate of
// turning of n, 1/ℓ).
type EOS struct {
	Kind   byte // 'Q' or 'S'
	Lambda float64
	Eps    float64 // S only
	K      float64 // twist rate
}

// Q returns the quartic ether with amplitude stiffness λ.
func Q(lambda float64) EOS { return EOS{Kind: 'Q', Lambda: lambda} }

// S returns the sextic ether with stiffness λ and cavity energy ε.
func S(lambda, eps float64) EOS { return EOS{Kind: 'S', Lambda: lambda, Eps: eps} }

// Twist returns e with the twist rate k.
func (e EOS) Twist(k float64) EOS { e.K = k; return e }

// Xi is the healing length √(κ2/λ).
func (e EOS) Xi() float64 { return math.Sqrt(kappa2 / e.Lambda) }

// Vu returns V, ∂V/∂u, ∂²V/∂u² at u = a² (twist included).
func (e EOS) Vu(u float64) (v, vu, vuu float64) {
	l := e.Lambda
	switch e.Kind {
	case 'Q':
		v, vu, vuu = l/4*(u-1)*(u-1), l/2*(u-1), l/2
	default:
		w := 1 - u
		v = l*u*w*w + e.Eps*w*w*(1+2*u)
		vu = l*w*(1-3*u) - 6*e.Eps*u*w
		vuu = l*(6*u-4) - 6*e.Eps*(1-2*u)
	}
	k2 := kappa2 * e.K * e.K
	return v + k2*u, vu + k2, vuu
}

// V, D1, D2 are the potential and its first two derivatives in a (valid for
// a < 0 too: the field φ passing through zero along a line).
func (e EOS) V(a float64) float64 { v, _, _ := e.Vu(a * a); return v }
func (e EOS) D1(a float64) float64 {
	_, vu, _ := e.Vu(a * a)
	return 2 * a * vu
}
func (e EOS) D2(a float64) float64 {
	u := a * a
	_, vu, vuu := e.Vu(u)
	return 2*vu + 4*u*vuu
}

// Ether returns the metastable or stable ether amplitude a0 near 1 (the
// local minimum of V + κ2k²a² closest to a = 1) and false when the twist
// has removed it (past the spinodal).
func (e EOS) Ether() (float64, bool) {
	if e.Kind == 'Q' {
		u := 1 - 2*kappa2*e.K*e.K/e.Lambda
		if u <= 0 {
			return 0, false
		}
		return math.Sqrt(u), true
	}
	// S: root of Vu(u) = 0 on the branch through u = 1, found by walking
	// down from u = 1 while ∂Vu/∂u > 0 (the minimum branch).
	f := func(u float64) float64 { _, vu, _ := e.Vu(u); return vu }
	if e.K == 0 {
		return 1, true
	}
	// Vu(1) = κ2k² > 0; on the minimum branch Vu increases with u, so the
	// root lies below 1 where Vu < 0. Scan down for the first sign change
	// before Vuu turns negative.
	hi := 1.0
	for u := 1.0; u > 0; u -= 1e-4 {
		_, vu, vuu := e.Vu(u)
		if vuu <= 0 {
			return 0, false
		}
		if vu < 0 {
			lo := u
			for it := 0; it < 100; it++ {
				m := (lo + hi) / 2
				if f(m) < 0 {
					lo = m
				} else {
					hi = m
				}
			}
			return math.Sqrt((lo + hi) / 2), true
		}
		hi = u
	}
	return 0, false
}
