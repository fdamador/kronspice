// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Part 11: g − 2 on the photon ring.
//
// The ring's light carries energy E_c round a loop at c. Its angular
// momentum is (E_c/c)·⟨ρ⟩, its charge e riding with it at c gives a moment
// μ = e·c·⟨ρ⟩/2, and g = 2mμ/(e·S). With S = ħ/2 and the same ⟨ρ⟩ for
// charge and light this is
//
//	g = M c² / E_c,
//
// the rest energy over the part of it that circulates. If all of it
// circulates (H23 as written) g = 1, as claude_electron step 0 found.
// The rest of the mass has to sit in something that does not go round:
// the confinement that holds the light, and the charge's own static field.

// Alpha is the fine structure constant (CODATA 2022).
const Alpha = 7.2973525643e-3

// Measured g − 2 and its QED terms.
const (
	AeMeasured  = 1.15965218059e-3 // (g − 2)/2, Fan et al. 2023
	AeSchwinger = Alpha / (2 * math.Pi)
)

// SplitG is g for light held by a squeeze k·R^n plus a static,
// non-circulating self-energy B/R next to the light's A/R:
//
//	E(R) = (A + B)/R + k·R^n,  equilibrium (A + B)/R = n·k·R^n,
//	M c² = (A + B)/R·(1 + 1/n),  E_c = A/R,
//	g = (1 + 1/n)·(1 + B/A).
//
// g = 2 needs n = 1 (a constant-force string) and B = 0.
func SplitG(n, bOverA float64) float64 { return (1 + 1/n) * (1 + bOverA) }

// StringRing is the n = 1 ring of mass m: the light carries half the rest
// energy, its angular momentum (mc/2)·R = ħ/2 puts it on R = ħ/(mc), and
// a photon of energy mc²/2 has wavelength 2λ_C — exactly twice the
// circumference: the double loop. Returned: radius [m], string tension
// [N], photon wavelength / circumference.
func StringRing(m float64) (R, k, ratio float64) {
	R = Hbar / (m * C0)
	k = m * C0 * C0 / 2 / R // (A + B)/R² with B = 0
	lam := 2 * math.Pi * Hbar * C0 / (m * C0 * C0 / 2)
	return R, k, lam / (2 * math.Pi * R)
}

// WvdMMoment is H23's loop (radius λ_C/4π, all energy circulating):
// μ/μ_B and g.
func WvdMMoment() (muOverMuB, g float64) {
	r := Compton() / (4 * math.Pi)
	mu := 1.602176634e-19 * C0 * r / 2
	muB := 1.602176634e-19 * Hbar / (2 * Me)
	return mu / muB, 2 * Me * mu / (1.602176634e-19 * Hbar / 2)
}

// ---------------------------------------------------------------- anomaly

// The anomaly comes from B, the charge's self-energy, which sits in its
// static field and does not circulate: U = α·ħc·β/R, so B/A = 2αβ (A =
// ħc/2 from S = ħ/2). Schwinger's α/2π needs β = 1/(4π).
const BetaSchwinger = 1 / (4 * math.Pi)

// ellipK is the complete elliptic integral K(m), m = k², by the AGM.
func ellipK(m float64) float64 {
	a, b := 1.0, math.Sqrt(1-m)
	for i := 0; i < 40 && math.Abs(a-b) > 1e-15*a; i++ {
		a, b = (a+b)/2, math.Sqrt(a*b)
	}
	return math.Pi / (2 * a)
}

// TorusInvD is ⟨1/d⟩ for charge spread uniformly over the surface of a
// torus of major radius 1 and tube radius a (≤ 1): two coaxial circles at
// tube angles v, v' interact through 4K(m)/√((ρ+ρ')²+Δz²) averaged over
// 2π, and the circles are weighted by their length. Staggered midpoint
// grids keep the log singularity off the nodes.
func TorusInvD(a float64, n int) float64 {
	w := func(v float64) float64 { return (1 + a*math.Cos(v)) / (2 * math.Pi) }
	h := 2 * math.Pi / float64(n)
	s := 0.0
	for i := 0; i < n; i++ {
		v := (float64(i) + 0.5) * h
		r1, z1 := 1+a*math.Cos(v), a*math.Sin(v)
		for j := 0; j < n; j++ {
			vp := (float64(j) + 0.25) * h // offset quarter cell from the field grid
			r2, z2 := 1+a*math.Cos(vp), a*math.Sin(vp)
			den := (r1+r2)*(r1+r2) + (z1-z2)*(z1-z2)
			m := 4 * r1 * r2 / den
			if m >= 1 {
				m = 1 - 1e-16
			}
			s += w(v) * w(vp) * (2 / math.Pi) * ellipK(m) / math.Sqrt(den) * h * h
		}
	}
	return s
}

// ThinTorusInvD is the thin-torus limit ln(8/a)/π (from the capacitance
// 4π²ε₀R/ln(8R/a)).
func ThinTorusInvD(a float64) float64 { return math.Log(8/a) / math.Pi }

// Beta is β = ½⟨1/d⟩·R for a charge distribution with mean inverse
// separation invD (in units of 1/R).
func Beta(invD float64) float64 { return invD / 2 }

// Anomaly is (g − 2)/2 for the string ring with self-energy factor β:
// g/2 = 1 + 2αβ.
func Anomaly(beta float64) float64 { return 2 * Alpha * beta }

// ShellForSchwinger is the radius (in units of R) of a uniform charged
// shell whose self-energy gives α/2π: β = 1/(2s) ⇒ s = 2π. With the
// magnetic self-energy of the moving charge counted too (equal to the
// electric for v = c), β doubles and s = 4π.
func ShellForSchwinger(magnetic bool) float64 {
	if magnetic {
		return 4 * math.Pi
	}
	return 2 * math.Pi
}

// NForAnomaly is the squeeze exponent that alone would give the measured
// g with B = 0: 1 + 1/n = 2(1 + a_e).
func NForAnomaly() float64 { return 1 / (1 + 2*AeMeasured) }
