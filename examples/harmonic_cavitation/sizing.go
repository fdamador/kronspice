// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/fluid"
)

// Water at 300 K. ρ and c from IAPWS-95 (fluid.Water); saturation pressure
// from IAPWS-IF97 (Wagner et al., J. Eng. Gas Turbines Power 122, 150
// (2000), Table 35: p_s(300 K) = 0.353658941e-2 MPa).
var (
	rho  = fluid.Water.Density
	cLiq = fluid.Water.SoundSpeed
)

const (
	pVap = 3536.58941 // [Pa] saturation pressure at 300 K, IAPWS-IF97
	pAtm = 101325.0   // [Pa]
)

// Sigma is the cavitation number (p − p_v)/(½ρU²).
func Sigma(p, U float64) float64 { return (p - pVap) / (0.5 * rho * U * U) }

// omega converts rpm to rad/s.
func omega(rpm float64) float64 { return 2 * math.Pi * rpm / 60 }

// Reactor is a published rotational cavitation reactor used to calibrate
// the sizing rules (LITERATURE.md §1).
type Reactor struct {
	Name     string
	R        float64 // rotor radius at the cavities [m]
	W        float64 // cavity mouth width along the motion [m]
	Cavities int     // cavities round the rotor
	RPMNone  float64 // highest speed seen without cavitation
	RPMOn    float64 // lowest speed seen with cavitation
	RPMBest  float64 // speed of the strongest pulsation
	FDom     float64 // dominant pulsation frequency [Hz] ...
	RPMAtF   float64 // ... measured at this speed
	PRef     float64 // pressure at the cavities [Pa] (declared: not reported)
}

// Xue2024 is the cylindrical reactor of L. Xue et al., Ultrason. Sonochem.
// 109, 106999 (2024), Sec. 2.1 and 3.1–3.3: rotor Ø 80 mm, 20 grooves 6 mm
// wide × 3 mm deep, stator gap 5 mm; no bubbles at 3000 rpm, bubbles from
// 3400 rpm, best pulsation at 3800 rpm, dominant 1122 Hz at 3600 rpm.
// The paper does not give the pressure at the rotor: 1 atm is declared.
var Xue2024 = Reactor{
	Name: "Xue et al. 2024", R: 0.040, W: 0.006, Cavities: 20,
	RPMNone: 3000, RPMOn: 3400, RPMBest: 3800, FDom: 1122, RPMAtF: 3600,
	PRef: pAtm,
}

// Calibration holds what one reactor says about a cavity row.
type Calibration struct {
	SigmaNone float64 // σ at the fastest non-cavitating speed (inception is below this)
	SigmaOn   float64 // σ at the slowest cavitating speed (inception is above this)
	SigmaBest float64 // σ at the best speed
	StLock    float64 // f·W/U of the dominant pulsation
	PassRatio float64 // dominant f / cavity-passing frequency
}

// Calibrate reads the σ band and the pulsation Strouhal number from a
// reactor, with U the rotor surface speed ωR (the reactor papers' velocity).
func Calibrate(r Reactor) Calibration {
	U := func(rpm float64) float64 { return omega(rpm) * r.R }
	return Calibration{
		SigmaNone: Sigma(r.PRef, U(r.RPMNone)),
		SigmaOn:   Sigma(r.PRef, U(r.RPMOn)),
		SigmaBest: Sigma(r.PRef, U(r.RPMBest)),
		StLock:    r.FDom * r.W / U(r.RPMAtF),
		PassRatio: r.FDom / (float64(r.Cavities) * r.RPMAtF / 60),
	}
}

// WoodSpeed is the low-frequency sound speed of a bubbly mixture at void
// fraction α and pressure p, with isothermal bubble contents (A. B. Wood,
// A Textbook of Sound, 1930; C. E. Brennen, Cavitation and Bubble Dynamics,
// OUP 1995, ch. 6, homogeneous mixture, polytropic index 1):
// 1/(ρ_m c²) = α/p + (1 − α)/(ρ_l c_l²), ρ_m = (1 − α)ρ_l.
func WoodSpeed(alpha, p float64) float64 {
	rm := (1 - alpha) * rho
	return math.Sqrt(1 / (rm * (alpha/p + (1-alpha)/(rho*cLiq*cLiq))))
}

// SonicVoid is the void fraction (on the bubbly branch, α ≤ 0.5) at which
// the mixture's sound speed equals U. It returns NaN when U is below the
// mixture's minimum speed (≈ √(4p/ρ) near α = 0.5), when no α matches.
func SonicVoid(U, p float64) float64 {
	lo, hi := 1e-9, 0.5
	if WoodSpeed(hi, p) > U || WoodSpeed(lo, p) < U {
		return math.NaN()
	}
	for i := 0; i < 80; i++ { // c falls monotonically with α on this branch
		m := 0.5 * (lo + hi)
		if WoodSpeed(m, p) > U {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi)
}
