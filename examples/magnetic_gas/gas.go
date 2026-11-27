// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// Physical constants (CODATA 2018, exact or to the digits used).
const (
	kB    = 1.380649e-23     // Boltzmann [J/K]
	qe    = 1.602176634e-19  // elementary charge [C]
	mu0   = 1.25663706212e-6 // vacuum permeability [H/m]
	muB   = 9.2740100783e-24 // Bohr magneton [J/T]
	NA    = 6.02214076e23    // Avogadro [1/mol]
	Rgas  = NA * kB          // molar gas constant [J/(mol·K)]
	atm   = 101325.0         // [Pa]
	torr  = atm / 760        // [Pa]
	psi   = 6894.757293168   // [Pa]
	Troom = 293.15           // [K]
)

// Gas holds the numbers this example needs for one gas. Every value is
// measured and cited unless its comment says declared or estimate.
type Gas struct {
	Name, Formula string
	Molar         float64 // molar mass [kg/mol]
	// Paschen minimum: the lowest sparking voltage [V] and the pressure ×
	// gap where it occurs [Torr·cm], uniform field, from M. S. Naidu and
	// V. Kamaraju, High Voltage Engineering, p. 27 (table reproduced at
	// pe2bz.philpem.me.uk, "Paschen's Law").
	Vmin, PDmin float64
	// Townsend coefficients A [1/(cm·Torr)] and B [V/(cm·Torr)] from
	// Yu. P. Raizer, Gas Discharge Physics (Springer 1991), Table 4.1,
	// where listed (0 = not listed). Used only as a cross-check.
	A, B float64
	// ChiMolar is the molar magnetic susceptibility, CGS, in units of
	// 1e-6 cm³/mol at room temperature (CRC Handbook of Chemistry and
	// Physics, "Magnetic susceptibility of the elements and inorganic
	// compounds").
	ChiMolar float64
	Unpaired int // unpaired electrons in the neutral ground state
	// Ion is the singly charged ion: its unpaired electrons and effective
	// moment μ_eff in Bohr magnetons. Ar⁺ ²P₃/₂: g_J = 4/3, J = 3/2, so
	// μ_eff = g_J·√(J(J+1)) = 2.582. H₂⁺ ²Σ: spin only, √3 = 1.732. O₂⁺
	// ²Π_g: spin only, √3 (orbital part not modelled: an approximation).
	Ion         string
	IonUnpaired int
	IonMuEff    float64
	Ionize      float64 // first ionization energy [eV] (NIST Chemistry WebBook)
	EpsMinus1   float64 // ε − 1 at 20 °C, 1 atm (CRC, "Permittivity of gases")
	Gamma       float64 // heat-capacity ratio cp/cv near room temperature
	Fuel        bool    // contains fuel and oxidizer: a spark ignites it
	Estimate    bool    // Paschen values estimated from its components
}

// The gases compared. Argon is the reference case.
var (
	Argon = Gas{Name: "argon", Formula: "Ar", Molar: 39.948e-3,
		Vmin: 137, PDmin: 0.9, A: 12, B: 180,
		ChiMolar: -19.6, Unpaired: 0, Ion: "Ar⁺", IonUnpaired: 1, IonMuEff: 2.582,
		Ionize: 15.7596, EpsMinus1: 5.17e-4, Gamma: 5.0 / 3}
	Hydrogen = Gas{Name: "hydrogen", Formula: "H₂", Molar: 2.016e-3,
		Vmin: 273, PDmin: 1.15, A: 5, B: 130,
		ChiMolar: -3.99, Unpaired: 0, Ion: "H₂⁺", IonUnpaired: 1, IonMuEff: math.Sqrt(3),
		Ionize: 15.4259, EpsMinus1: 2.54e-4, Gamma: 1.41}
	Oxygen = Gas{Name: "oxygen", Formula: "O₂", Molar: 31.998e-3,
		Vmin: 450, PDmin: 0.7,
		ChiMolar: 3449, Unpaired: 2, Ion: "O₂⁺", IonUnpaired: 1, IonMuEff: math.Sqrt(3),
		Ionize: 12.0697, EpsMinus1: 4.94e-4, Gamma: 1.40}
	// HHO is oxyhydrogen, 2 H₂ : 1 O₂ by moles (the stoichiometric mix an
	// electrolyser makes). Its susceptibility, permittivity and molar mass
	// are the mole-weighted mix (exact for a dilute gas). No measured
	// Paschen curve was found: Vmin and pd_min are the mole-weighted
	// values of H₂ and O₂, an ESTIMATE. Its ion is taken as O₂⁺, the lowest
	// ionization energy present.
	HHO = mix("HHO (2 H₂ + O₂)", "HHO", Hydrogen, Oxygen, 2.0/3)
	// Air is used only to check the Paschen method against a well-known
	// uniform-field breakdown (about 30 kV across 1 cm at 1 atm).
	Air = Gas{Name: "air", Formula: "air", Molar: 28.97e-3, Vmin: 327, PDmin: 0.567, A: 15, B: 365}
)

// Gases is the comparison set, argon first.
var Gases = []Gas{Argon, Hydrogen, Oxygen, HHO}

func mix(name, formula string, a, b Gas, xa float64) Gas {
	xb := 1 - xa
	w := func(u, v float64) float64 { return xa*u + xb*v }
	return Gas{Name: name, Formula: formula, Molar: w(a.Molar, b.Molar),
		Vmin: w(a.Vmin, b.Vmin), PDmin: w(a.PDmin, b.PDmin),
		ChiMolar: w(a.ChiMolar, b.ChiMolar), Unpaired: b.Unpaired, Ion: b.Ion,
		IonUnpaired: b.IonUnpaired, IonMuEff: b.IonMuEff, Ionize: math.Min(a.Ionize, b.Ionize),
		EpsMinus1: w(a.EpsMinus1, b.EpsMinus1), Gamma: 1.40, Fuel: true, Estimate: true}
}

// Breakdown is the uniform-field sparking voltage [V] at pressure p [Pa]
// across gap d [m] from the Townsend criterion written through the
// measured Paschen minimum. With x = pd / pd_min the Townsend form
// V = B·pd / (ln(A·pd) − ln ln(1 + 1/γ)) becomes
//
//	V = V_min · x / (1 + ln x)
//
// (at the minimum V/pd = B and A·pd = e·ln(1 + 1/γ); Raizer 1991, §4.4).
// It needs no A, B or γ, only the two measured numbers, so every gas is
// treated by the same rule. Left of the minimum the curve turns up and
// below x = 1/e no breakdown is predicted (+Inf).
func (g Gas) Breakdown(p, d float64) float64 {
	x := p / torr * d * 100 / g.PDmin
	if x <= 1/math.E {
		return math.Inf(1)
	}
	return g.Vmin * x / (1 + math.Log(x))
}

// BreakdownRaizer is the Townsend form with Raizer's A and B and γ = 0.01,
// for gases he lists (NaN otherwise). A cross-check only.
func (g Gas) BreakdownRaizer(p, d float64) float64 {
	if g.A == 0 {
		return math.NaN()
	}
	pd := p / torr * d * 100
	den := math.Log(g.A*pd) - math.Log(math.Log(1+1/0.01))
	if den <= 0 {
		return math.Inf(1)
	}
	return g.B * pd / den
}

// Density is the number density of molecules [1/m³] at p [Pa], T [K]
// (ideal gas; at 34 atm and room temperature the compressibility factor
// of these gases is within a few per cent of 1).
func Density(p, T float64) float64 { return p / (kB * T) }

// Label is the gas name with its estimate flag.
func (g Gas) Label() string {
	if g.Estimate {
		return g.Formula + " (est.)"
	}
	return g.Formula
}

// Describe prints the gas in one line.
func (g Gas) Describe() string {
	return fmt.Sprintf("%-6s %-16s Paschen minimum %3.0f V at %.2f Torr·cm; ionization %.2f eV; %d unpaired e⁻ neutral, %d in %s",
		g.Formula, g.Name, g.Vmin, g.PDmin, g.Ionize, g.Unpaired, g.IonUnpaired, g.Ion)
}
