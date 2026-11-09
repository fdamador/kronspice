// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"sync"
)

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

// AppB is one pass through Nachamkin's Appendix B (R0 = 10 cm, 1 MJ at
// 0 °C, 1 atm): half the energy in the field (eq. B2), E0 from the field
// integral (B3), ω from his small-oscillation factor ξ = 0.01 (eq. 58,
// B5), k = ω/c moved to a zero of j1 (B6), the charge density from eq. 38
// (B7). IntegralM is the value used for ∫|M|²dV: his printed B3 is
// 35.87/α³; his B4 follows only from 35.87/α².
type AppB struct {
	Label                 string
	R, Alpha, IntegralM   float64
	E0                    float64 // B4 [V/m]
	Omega, K0             float64 // B5 [rad/s], k = ω/c [1/m]
	XNearest              float64 // nearest zero of j1 to k·R
	NNearest              int
	K                     float64 // the k used for B6/B7 [1/m]
	RhoE, Electrons       float64 // B7 [C/m³], count
	IonFrac               float64 // of the 1.126e23 molecules
	Beta                  float64 // electron speed ω·ξ·R ÷ c (eq. 58 at the final ω)
	PeakOverBreakdown     float64 // |E0·M| at the centre ÷ 3 MV/m
	Wp, DispersionResidue float64
}

// AppendixB recomputes the example. integralM is the ∫|M|²dV to use;
// kFixed, if > 0, is the k used for B6 (the paper's 2686 /m), else the
// nearest zero of j1 to the B5 value.
func AppendixB(label string, integralM, kFixed float64) AppB {
	e := AppB{Label: label, R: 0.1, IntegralM: integralM}
	e.Alpha = J1Zero / e.R
	const xi = 0.01
	e.E0 = math.Sqrt(4 * 0.5e6 / (Eps0 * integralM))
	e.Omega = math.Sqrt(EoM * e.E0 / (xi * e.R))
	e.K0 = e.Omega / C0
	e.XNearest, e.NNearest = nearestZero(e.K0 * e.R)
	e.K = e.XNearest / e.R
	if kFixed > 0 {
		e.K = kFixed
	}
	w := e.K * C0
	e.RhoE = -(1 / EoM) * Eps0 * C0 * C0 * (e.K*e.K - e.Alpha*e.Alpha)
	vol := 4.0 / 3 * math.Pi * e.R * e.R * e.R
	e.Electrons = -e.RhoE / Qe * vol
	e.IonFrac = e.Electrons / (Loschmidt * vol)
	e.Beta = w * xi * e.R / C0
	e.PeakOverBreakdown = e.E0 * math.Sqrt(MMax2) / AirBreakdown
	e.Wp = math.Sqrt(-e.RhoE / Qe * Qe * Qe / (Me * Eps0))
	e.DispersionResidue = (w*w - e.Wp*e.Wp - e.Alpha*e.Alpha*C0*C0) / (w * w)
	return e
}

// meanSquare is ⟨|F|²⟩ over the sphere ÷ max |F|² for the l = 1 field
// (computed once).
func meanSquare() float64 {
	msOnce.Do(func() { msValue = computeMeanSquare() })
	return msValue
}

var (
	msOnce  sync.Once
	msValue float64
)

func computeMeanSquare() float64 {
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
