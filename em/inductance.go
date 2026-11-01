// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// Inductance of filamentary coils (M9i): Neumann's formula and its
// closed forms.

// EllipticKE returns the complete elliptic integrals K(k) and E(k) of the
// first and second kind (modulus k, 0 ≤ k < 1) by the arithmetic–geometric
// mean (M. Abramowitz, I. A. Stegun, Handbook of Mathematical Functions,
// §17.6): a₀ = 1, b₀ = k′, c₀ = k; K = π/(2a_N), E = K·(1 − Σ 2ⁿ⁻¹cₙ²).
func EllipticKE(k float64) (K, E float64) {
	return ellipticKE(k, math.Sqrt((1-k)*(1+k)))
}

// ellipticKE takes the complementary modulus k′ = √(1 − k²) directly, so
// k → 1 keeps its precision when k′ is known accurately.
func ellipticKE(k, kp float64) (K, E float64) {
	a, b, c := 1.0, kp, k
	sum := 0.5 * c * c
	pow := 0.5
	for i := 0; i < 60; i++ {
		an, bn, cn := 0.5*(a+b), math.Sqrt(a*b), 0.5*(a-b)
		a, b, c = an, bn, cn
		pow *= 2
		sum += pow * c * c
		if math.Abs(c) < 1e-17*a {
			break
		}
	}
	K = math.Pi / (2 * a)
	return K, K * (1 - sum)
}

// MutualCircles is the mutual inductance of two coaxial circles of radii
// a and b at axial distance d (Maxwell, Treatise Vol. II art. 701):
// M = μ₀√(ab)·[(2/k − k)·K(k) − (2/k)·E(k)], k² = 4ab/((a + b)² + d²).
func MutualCircles(a, b, d float64) float64 {
	den := (a+b)*(a+b) + d*d
	k := math.Sqrt(4 * a * b / den)
	kp := math.Sqrt(((a-b)*(a-b) + d*d) / den)
	K, E := ellipticKE(k, kp)
	return Mu0 * math.Sqrt(a*b) * ((2/k-k)*K - 2/k*E)
}

// Neumann is the mutual inductance of two closed polygonal filaments,
// M = μ₀/4π ∮∮ dl₁·dl₂/|r₁ − r₂|, by Gauss–Legendre quadrature on every
// pair of segments (the filaments must not touch).
func Neumann(l1, l2 Loop) float64 {
	s := 0.0
	n1, n2 := len(l1.P), len(l2.P)
	for i := 0; i < n1; i++ {
		a, b := l1.P[i], l1.P[(i+1)%n1]
		d1 := sub3(b, a)
		for j := 0; j < n2; j++ {
			c, d := l2.P[j], l2.P[(j+1)%n2]
			d2 := sub3(d, c)
			dd := dot(d1, d2)
			if dd == 0 {
				continue
			}
			for p := range gl8x {
				x := add3(a, scale3(d1, 0.5*(1+gl8x[p])))
				for q := range gl8x {
					y := add3(c, scale3(d2, 0.5*(1+gl8x[q])))
					r := sub3(x, y)
					s += 0.25 * gl8w[p] * gl8w[q] * dd / math.Sqrt(dot(r, r))
				}
			}
		}
	}
	return Mu0 / (4 * math.Pi) * s
}

// SolenoidSheet is the inductance of a cylindrical current sheet of radius
// R and length l carrying N turns' worth of current uniformly: the double
// integral of MutualCircles over its length,
// L = (N/l)²·2∫₀ˡ (l − u)·M(R, R, u) du, with the logarithmic singularity of
// M at u = 0 handled by graded Gauss–Legendre panels.
func SolenoidSheet(R, l float64, N float64) float64 {
	f := func(u float64) float64 { return (l - u) * MutualCircles(R, R, u) }
	// Panels [0, h0], [h0, 2h0], [2h0, 4h0], … up to l, then Gauss on each;
	// the first panel uses the substitution u = h0·t² to tame the log.
	h0 := l * 1e-9
	s := 0.0
	for p := range gl8x {
		t := 0.5 * (1 + gl8x[p])
		u := h0 * t * t
		s += 0.5 * gl8w[p] * f(u) * 2 * h0 * t
	}
	for a := h0; a < l; {
		b := math.Min(2*a, l)
		for p := range gl8x {
			u := a + (b-a)*0.5*(1+gl8x[p])
			s += 0.5 * gl8w[p] * (b - a) * f(u)
		}
		a = b
	}
	return (N / l) * (N / l) * 2 * s
}

// Nagaoka is Lorenz's closed form for the same current sheet,
// L = μ₀πR²N²/l · K_N, with Nagaoka's coefficient
// K_N = 4/(3πk′)·[(k′²/k²)·K(k) − ((k′² − k²)/k²)·E(k) − k],
// k = 2R/√(4R² + l²), k′ = √(1 − k²) (H. Nagaoka, J. Coll. Sci. Tokyo 27,
// 1909; F. W. Grover, Inductance Calculations, Van Nostrand 1946, Ch. 13).
func Nagaoka(R, l float64, N float64) float64 {
	k := 2 * R / math.Sqrt(4*R*R+l*l)
	kp := math.Sqrt(1 - k*k)
	K, E := EllipticKE(k)
	kn := 4 / (3 * math.Pi * kp) * ((kp*kp/(k*k))*K - ((kp*kp-k*k)/(k*k))*E - k)
	return Mu0 * math.Pi * R * R * N * N / l * kn
}
