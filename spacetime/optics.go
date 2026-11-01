// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import "math"

// Light in a static gravitational field moves as in a medium of index n
// with ε_r = μ_r = n (Plebanski 1960; Landau & Lifshitz §90). The
// impedance √(μ/ε) stays η0, so gravity bends and delays waves without
// reflecting them.

// IndexWeak is the weak-field index n = 1 + (1+γ)·GM/(r c²) at radius r
// (PPN; γ = 1 in general relativity; C. M. Will, Living Rev. Relativ. 17,
// 4, 2014, §3.4).
func IndexWeak(gm, r, gamma float64) float64 { return 1 + (1+gamma)*gm/(r*C*C) }

// IndexIsotropic is the exact Schwarzschild index in isotropic radius ρ,
// n = (1 + r_s/4ρ)³/(1 − r_s/4ρ) (Eddington, The Mathematical Theory of
// Relativity, 1923, §43). It reduces to IndexWeak(γ = 1) for ρ ≫ r_s.
func IndexIsotropic(gm, rho float64) float64 {
	u := SchwarzschildRadius(gm) / (4 * rho)
	return (1 + u) * (1 + u) * (1 + u) / (1 - u)
}

// ShapiroStraight is the one-way excess delay [s] of a signal along a
// straight path passing a mass gm at closest distance b, from a point a
// distance x1 before closest approach to one x2 after it (both ≥ 0 along
// the path), weak field: ((1+γ)GM/c³)·ln[(r1 + x1)(r2 + x2)/b²],
// r = √(x² + b²) (I. I. Shapiro, Phys. Rev. Lett. 13, 789, 1964).
func ShapiroStraight(gm, b, x1, x2, gamma float64) float64 {
	r1, r2 := math.Hypot(x1, b), math.Hypot(x2, b)
	return (1 + gamma) * gm / (C * C * C) * math.Log((r1+x1)*(r2+x2)/(b*b))
}

// PathDelay integrates the excess delay ∫(n − 1) dl / c [s] along the
// straight path of ShapiroStraight for any index profile n(r), by
// Simpson's rule in the substitution x = b·sinh(s) (n samples).
func PathDelay(n func(r float64) float64, b, x1, x2 float64, steps int) float64 {
	if steps%2 == 1 {
		steps++
	}
	s1, s2 := -math.Asinh(x1/b), math.Asinh(x2/b)
	h := (s2 - s1) / float64(steps)
	f := func(s float64) float64 {
		x := b * math.Sinh(s)
		return (n(math.Hypot(x, b)) - 1) * b * math.Cosh(s)
	}
	sum := f(s1) + f(s2)
	for i := 1; i < steps; i++ {
		w := 2.0
		if i%2 == 1 {
			w = 4
		}
		sum += w * f(s1+float64(i)*h)
	}
	return sum * h / 3 / C
}

// Deflection is the total bending angle [rad] of a ray with impact
// parameter b passing a spherical index profile n(r), integrated along the
// straight path from the transverse gradient of n (eikonal, weak
// deflection): α = −∫ (dn/dr)·(b/r) dl. For IndexWeak it equals
// 2(1+γ)GM/(b c²) — 1.75″ at the solar limb for γ = 1.
func Deflection(n func(r float64) float64, b float64, steps int) float64 {
	if steps%2 == 1 {
		steps++
	}
	smax := math.Asinh(1e6) // path from −10⁶ b to +10⁶ b
	h := 2 * smax / float64(steps)
	f := func(s float64) float64 {
		x := b * math.Sinh(s)
		r := math.Hypot(x, b)
		dr := 1e-4 * r
		dndr := (n(r+dr) - n(r-dr)) / (2 * dr)
		return -dndr * b / r * b * math.Cosh(s)
	}
	sum := f(-smax) + f(smax)
	for i := 1; i < steps; i++ {
		w := 2.0
		if i%2 == 1 {
			w = 4
		}
		sum += w * f(-smax+float64(i)*h)
	}
	return sum * h / 3
}
