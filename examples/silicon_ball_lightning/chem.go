// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The chemical store and the time it takes to empty.
//
// Energy. A gram of an element burned in air releases its oxidation heat;
// making that gram into vapor (the arc's or the lightning's job) costs the
// heat to warm, melt, warm and boil it. The vapor condenses into
// nanoparticles near the source, giving the boiling heat back there; what
// the ball then carries away is the oxidation heat, released as the
// particles burn. The ratio "oxidation heat ÷ vaporization cost" says how
// much chemical energy the source buys per joule spent making the vapor.
//
// Time. Silicon oxidizes through a growing SiO2 shell. Deal & Grove
// (J. Appl. Phys. 36, 3770, 1965) give the shell thickness x after time t
// for dry O2 on silicon as x² + A·x = B·(t + τ0), with
//
//	B   = 772 µm²/h · exp(−1.23 eV/kT)       (parabolic, diffusion-limited)
//	B/A = 6.23×10⁶ µm/h · exp(−2.00 eV/kT)  (linear, reaction-limited)
//
// fitted on planar wafers at 700–1200 °C in 1 atm of O2. Both constants
// are close to proportional to the oxygen partial pressure (Deal & Grove
// 1965), so in air they are scaled by pO2 (0.21 atm at sea level, 0.16 at
// the Qinghai site, 2530 m). One µm of oxide consumes
// 0.44 µm of silicon, so a particle of radius r is burned through when
// x ≈ r/0.44 (planar approximation, thin shells), at t = x²/B + x/(B/A).

import "math"

const kB = 8.617333262e-5 // eV/K

// OxidationHeat is the heat released per gram burned [J/g].
func (e Element) OxidationHeat() float64 { return e.HOxPerMol * 1e3 / e.M }

// VaporizationCost is the heat to take one gram from 298 K to vapor at its
// boiling point [J/g].
func (e Element) VaporizationCost() float64 {
	h := e.CpS*(e.Tm-298) + e.HFus*1e3 + e.CpL*(e.Tb-e.Tm) + e.HVap*1e3
	return h / e.M
}

// Return is oxidation heat ÷ vaporization cost.
func (e Element) Return() float64 { return e.OxidationHeat() / e.VaporizationCost() }

// MassFor is the grams of element whose burning releases energy E [J].
func (e Element) MassFor(E float64) float64 { return E / e.OxidationHeat() }

// DealGroveB is the parabolic rate constant [m²/s] at temperature T [K].
func DealGroveB(T float64) float64 {
	return 772 * 1e-12 / 3600 * math.Exp(-1.23/(kB*T))
}

// DealGroveBA is the linear rate constant B/A [m/s].
func DealGroveBA(T float64) float64 {
	return 6.23e6 * 1e-6 / 3600 * math.Exp(-2.00/(kB*T))
}

// BurnTime is the time [s] for a silicon particle of radius r [m] to oxidize
// through at temperature T [K] in oxygen of partial pressure pO2 [atm]
// (Deal–Grove, planar approximation).
func BurnTime(r, T, pO2 float64) float64 {
	x := r / 0.44
	return x*x/(pO2*DealGroveB(T)) + x/(pO2*DealGroveBA(T))
}
