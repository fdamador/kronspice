// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Magnetism of a gas, neutral and ionized, in SI volume susceptibility χ
// (M = χH, dimensionless).

// ChiNeutral is the volume susceptibility of the neutral gas at p [Pa] and
// T [K]. The molar CGS value converts as χ_m(SI) = 4π·10⁻⁶·χ_m(CGS)
// [m³/mol] (CRC, "Magnetic susceptibility", units note), times the molar
// density p/RT. Paramagnetic O₂ follows Curie's law (∝ 1/T); its table
// value is at 293 K, which is where it is used. The diamagnetic gases are
// temperature independent.
func (g Gas) ChiNeutral(p, T float64) float64 {
	chi := g.ChiMolar * 1e-6 // cm³/mol, CGS
	if g.Unpaired > 0 {
		chi *= Troom / T
	}
	return 4 * math.Pi * 1e-6 * chi * p / (Rgas * T)
}

// ChiCurie is the Curie susceptibility of n [1/m³] moments of μ_eff [μ_B]
// at temperature T: χ = μ₀·n·μ_eff²·μ_B² / (3·k·T) (Kittel, Introduction
// to Solid State Physics, 8th ed., eq. 14.21).
func ChiCurie(n, muEff, T float64) float64 {
	return mu0 * n * muEff * muEff * muB * muB / (3 * kB * T)
}

// ChiElectrons is the susceptibility of n free, non-degenerate electrons
// at temperature T: spin (Curie, μ_eff² = 3) plus Landau diamagnetism,
// which is −1/3 of the spin part in this limit (Landau 1930; Kittel,
// Introduction to Solid State Physics, ch. 14). Net χ = ⅔·μ₀·n·μ_B²/(kT).
func ChiElectrons(n, T float64) float64 {
	return 2.0 / 3 * mu0 * n * muB * muB / (kB * T)
}

// Ionized is the gas made singly ionized at a given fraction and
// temperature: what removing one electron per molecule does to the
// magnetism, before any confinement is considered.
type Ionized struct {
	Gas
	P, Fraction, T  float64 // fill pressure [Pa] at room temperature, ion fraction, temperature [K]
	N, Ni           float64 // molecules, ions [1/m³] (density fixed at the fill)
	Neutral         float64 // χ of the remaining neutrals at T
	Ions, Electrons float64 // χ of the ions and of the free electrons
	Total           float64
}

// Ionized returns the magnetism of the gas filled to p at room temperature,
// a fraction f of its molecules then singly ionized at temperature T
// (sealed volume: the density stays, the pressure rises).
func (g Gas) Ionized(p, f, T float64) Ionized {
	n := Density(p, Troom)
	ni := f * n
	nn := n - ni
	r := Ionized{Gas: g, P: p, Fraction: f, T: T, N: n, Ni: ni}
	r.Neutral = g.ChiNeutral(nn*kB*T, T)
	r.Ions = ChiCurie(ni, g.IonMuEff, T)
	r.Electrons = ChiElectrons(ni, T)
	r.Total = r.Neutral + r.Ions + r.Electrons
	return r
}

// MagneticPressure is B²/2μ₀ [Pa] (Jackson, Classical Electrodynamics,
// 3rd ed., §6.7: energy density of the field = pressure on a conductor).
func MagneticPressure(B float64) float64 { return B * B / (2 * mu0) }

// FieldToHold is the field whose magnetic pressure equals p [T], plasma
// β = 2μ₀p/B² = 1. A field holds a plasma only where β < 1; the plasma's
// diamagnetic current lowers the field inside by just enough that
// B_in²/2μ₀ + p = B_out²/2μ₀ (pressure balance; Chen, Introduction to
// Plasma Physics, 3rd ed., §3.5).
func FieldToHold(p float64) float64 { return math.Sqrt(2 * mu0 * p) }

// PlasmaPressure is the kinetic pressure of a sealed fill (density from p
// at room temperature) heated to T with a fraction f singly ionized: the
// neutrals, ions and free electrons each count, p = (1 + f)·n·k·T
// (electrons and ions at the same T, declared).
func PlasmaPressure(p, f, T float64) float64 {
	return (1 + f) * Density(p, Troom) * kB * T
}

// AdiabaticT is the temperature after compressing a gas from p1 to p2 in
// one fast stroke starting at T1, T2 = T1·(p2/p1)^((γ−1)/γ) (reversible
// adiabatic, ideal gas).
func (g Gas) AdiabaticT(T1, p1, p2 float64) float64 {
	return T1 * math.Pow(p2/p1, (g.Gamma-1)/g.Gamma)
}

// IonizeEnergy is the energy to singly ionize every molecule in volume V
// [m³] filled to p at room temperature [J]: n·V·E_i. A spark ionizes only a
// thin channel; this is the bound for the whole fill.
func (g Gas) IonizeEnergy(p, V float64) float64 {
	return Density(p, Troom) * V * g.Ionize * qe
}

// ForceOverWeight is the magnetic force on the gas in a field with
// B·dB/dz [T²/m], f = (χ/μ₀)·B·dB/dz per unit volume (Jackson §5.16, for
// |χ| ≪ 1), divided by its weight ρg at the same p, T. Positive pulls into
// the field. 100 T²/m is a strong laboratory magnet (2 T falling over
// 4 cm).
func (g Gas) ForceOverWeight(p, T, BdB float64) float64 {
	f := g.ChiNeutral(p, T) / mu0 * BdB
	rho := p * g.Molar / (Rgas * T)
	return f / (rho * 9.80665)
}

// HHOChemical is the heat of combustion [J] stored in volume V [m³] of
// oxyhydrogen filled to p at room temperature, at the lower heating
// value of hydrogen, 241.8 kJ per mole of H₂ burned to water vapour (NIST
// Chemistry WebBook, ΔfH° of H₂O(g) = −241.826 kJ/mol).
func HHOChemical(p, V float64) float64 {
	molH2 := 2.0 / 3 * p * V / (Rgas * Troom)
	return molH2 * 241.826e3
}
