// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// SI constants (CODATA 2018).
const (
	Hbar = 1.054571817e-34
	C    = 299792458.0
	Eps0 = 8.8541878128e-12
	Qe   = 1.602176634e-19
)

// PerfectEnergy is Casimir's energy per area between perfect mirrors at
// gap a (Forward eq. 1, πhc/1440a³ = π²ħc/720a³) [J/m²], negative.
func PerfectEnergy(a float64) float64 { return -math.Pi * math.Pi * Hbar * C / (720 * a * a * a) }

// PerfectPressure is the attractive pressure (Forward eq. 2) [Pa], negative.
func PerfectPressure(a float64) float64 {
	return -math.Pi * math.Pi * Hbar * C / (240 * a * a * a * a)
}

// LifshitzEnergy is the zero-temperature Lifshitz energy per area between
// two identical half-spaces with the plasma-model permittivity
// ε(iξ) = 1 + ωp²/ξ² (plasma wavelength lp = 2πc/ωp) [J/m²]:
//
//	E = ħc/(4π²a³) ∫∫ y dy dx Σ_{TE,TM} ln(1 − r² e^{−2ρ}),  ρ = √(x²+y²)
//
// in the scaled variables x = ξa/c, y = k·a (E. M. Lifshitz, JETP 2, 73,
// 1956; the form of Lambrecht & Reynaud, Eur. Phys. J. D 8, 309, 2000).
// With lp → 0 it is Casimir's result.
func LifshitzEnergy(a, lp float64) float64 {
	wp := 2 * math.Pi * a / lp // ωp·a/c
	// polar coordinates: x = ρ cos φ, y = ρ sin φ; dx dy y = ρ² sin φ dρ dφ
	const nr, nf = 400, 64
	rmax := 25.0
	sum := 0.0
	for i := 0; i < nr; i++ {
		// Gauss-like midpoint in u with ρ = u²·rmax for resolution near 0
		u := (float64(i) + 0.5) / nr
		rho := u * u * rmax
		drho := 2 * u * rmax / nr
		e := math.Exp(-2 * rho)
		inner := 0.0
		for j := 0; j < nf; j++ {
			phi := (float64(j) + 0.5) / nf * math.Pi / 2
			x, y := rho*math.Cos(phi), rho*math.Sin(phi)
			r2TE, r2TM := 1.0, 1.0
			if lp > 0 {
				eps := 1 + wp*wp/(x*x)
				kap := math.Sqrt(y*y + eps*x*x)
				rte := (rho - kap) / (rho + kap)
				rtm := (eps*rho - kap) / (eps*rho + kap)
				r2TE, r2TM = rte*rte, rtm*rtm
			}
			inner += math.Sin(phi) * (math.Log1p(-r2TE*e) + math.Log1p(-r2TM*e))
		}
		sum += rho * rho * inner * (math.Pi / 2 / nf) * drho
	}
	return Hbar * C / (4 * math.Pi * math.Pi * a * a * a) * sum
}

// Material is a metal for the leaves.
type Material struct {
	Name    string
	Lp      float64 // plasma wavelength [m] (0: perfect mirror)
	Density float64 // [kg/m³]
}

// Aluminium: plasma wavelength 107 nm (Lambrecht & Reynaud 2000),
// density 2,700 kg/m³.
var (
	Aluminium = Material{"aluminium", 107e-9, 2700}
	Perfect   = Material{"perfect mirror", 0, 2700}
)

// Table is a log-spaced table of the energy per area against gap, for fast
// interpolation of energy and pressure.
type Table struct {
	la, e []float64 // ln a, energy
}

// NewTable tabulates the energy from amin to amax.
func NewTable(m Material, amin, amax float64, n int) Table {
	var t Table
	for i := 0; i < n; i++ {
		a := amin * math.Pow(amax/amin, float64(i)/float64(n-1))
		t.la = append(t.la, math.Log(a))
		t.e = append(t.e, LifshitzEnergy(a, m.Lp)/PerfectEnergy(a)) // η_E
	}
	return t
}

// eta interpolates η_E(a) linearly in ln a (clamped at the ends).
func (t Table) eta(a float64) (eta, deta float64) {
	x := math.Log(a)
	n := len(t.la)
	if x <= t.la[0] {
		return t.e[0], 0
	}
	if x >= t.la[n-1] {
		return t.e[n-1], 0
	}
	i := int((x - t.la[0]) / (t.la[1] - t.la[0]))
	if i >= n-1 {
		i = n - 2
	}
	f := (x - t.la[i]) / (t.la[i+1] - t.la[i])
	d := (t.e[i+1] - t.e[i]) / (t.la[i+1] - t.la[i]) / a
	return t.e[i] + f*(t.e[i+1]-t.e[i]), d
}

// Energy is the energy per area at gap a [J/m²].
func (t Table) Energy(a float64) float64 {
	h, _ := t.eta(a)
	return h * PerfectEnergy(a)
}

// Pressure is −dE/da [Pa] (negative: attraction).
func (t Table) Pressure(a float64) float64 {
	h, dh := t.eta(a)
	// E = η·E0, E0 ∝ a⁻³: dE/da = dη·E0 + η·(−3E0/a)
	e0 := PerfectEnergy(a)
	return -(dh*e0 - 3*h*e0/a)
}
