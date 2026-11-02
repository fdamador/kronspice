// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// SI constants (CODATA 2018).
const (
	C0   = 299792458.0
	Eps0 = 8.8541878128e-12
	Mu0  = 1.25663706212e-6
	Qe   = 1.602176634e-19
	Me   = 9.1093837015e-31
	EoM  = Qe / Me // |e/m|
)

// Nachamkin's relations (§III, §V), and what they are in standard plasma
// physics. Inside the sphere the fields have spatial wavenumber α and time
// frequency ω (vacuum wavenumber k = ω/c). His eq. 14 gives the current
// J = (α − k²/α)·H, and electrons oscillating along E (his eq. 37) carry it
// only if their charge density is (his eq. 49)
//
//	ρ_e = (m/e)·ε0·c²·(k² − α²)          (m/e < 0 for electrons)
//
// i.e. negative only above the critical frequency ω0 = αc (his eq. 15).
// With ρ_e = −n·e this is n·e²/(m·ε0) = ω² − α²c²: the electrons' plasma
// frequency ωp obeys
//
//	ω² = ωp² + α²c²
//
// — the dispersion relation of an electromagnetic wave in a cold plasma
// (e.g. F. F. Chen, Introduction to Plasma Physics, §4.16). Nachamkin's
// "critical frequency" is the cutoff of a wave whose wavenumber is α, and
// his plasmoid field is a standing wave of that plasma.

// Example is Nachamkin's Appendix B: a 10 cm radius, one-megajoule
// plasmoid in air at 0 °C, 1 atm.
type Example struct {
	R, Alpha, K, Omega  float64
	RhoE, N, Electrons  float64 // charge density, electron density, count
	Wp                  float64 // plasma frequency from n
	DispersionMismatch  float64 // (ω² − ωp² − α²c²)/ω²
	C                   float64 // his small-oscillation factor (eq. 50)
	E0                  float64 // field for that factor, eq. 50: E = C·R·ω²/|e/m|
	Quiver, QuiverOverC float64 // electron velocity amplitude ωξ = ω·C·R [m/s], ÷ c
	FieldEnergy         float64 // ½ε0E0²·(sphere volume)·(mean square of F ÷ max)
	AirBreakdownRatio   float64 // E0 ÷ 3 MV/m (DC breakdown of air at STP, order of magnitude)
}

// AppendixB recomputes his numbers.
func AppendixB() Example {
	e := Example{R: 0.1, K: 2686, C: 0.01}
	e.Alpha = J1Zero / e.R
	e.Omega = e.K * C0
	e.RhoE = -(1 / EoM) * Eps0 * C0 * C0 * (e.K*e.K - e.Alpha*e.Alpha)
	e.N = -e.RhoE / Qe
	vol := 4.0 / 3 * math.Pi * e.R * e.R * e.R
	e.Electrons = e.N * vol
	e.Wp = math.Sqrt(e.N * Qe * Qe / (Me * Eps0))
	e.DispersionMismatch = (e.Omega*e.Omega - e.Wp*e.Wp - e.Alpha*e.Alpha*C0*C0) / (e.Omega * e.Omega)
	e.E0 = e.C * e.R * e.Omega * e.Omega / EoM
	e.Quiver = e.Omega * e.C * e.R
	e.QuiverOverC = e.Quiver / C0
	e.FieldEnergy = 0.5 * Eps0 * e.E0 * e.E0 * vol * meanSquare()
	e.AirBreakdownRatio = e.E0 / 3e6
	return e
}

// meanSquare is ⟨|F|²⟩ over the sphere ÷ max |F|² for the l = 1 field.
func meanSquare() float64 {
	const n = 200
	var s, w, mx float64
	for i := 0; i < n; i++ {
		x := (float64(i) + 0.5) / n * J1Zero
		for j := 0; j < n; j++ {
			th := (float64(j) + 0.5) / n * math.Pi
			fr, ft, fp := Beltrami(1, x, th)
			f2 := fr*fr + ft*ft + fp*fp
			dv := x * x * math.Sin(th)
			s += f2 * dv
			w += dv
			mx = math.Max(mx, f2)
		}
	}
	return s / w / mx
}
