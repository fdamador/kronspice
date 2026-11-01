// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command space_vortex puts P. Tewari's space-vortex electron (hypothesis
// H13) on the electron test bench of doc/kelvin-references.md. Source:
// P. Tewari, "Structural relation between the vacuum space and the
// electron," Physics Essays 31(1), 108–129 (2018), DOI
// 10.4006/0836-1398-31.1.108 (free copy on tewari.org). The equations are
// transcribed from the paper's text layer and are to be checked against
// the typeset page (doc/tewari-references.md lists which).
//
// The model: the electron is an irrotational vortex of an incompressible,
// massless "space fluid" whose speed reaches c at the surface of a central
// spherical void of radius r_e. From that it states, in its own units
// (charge in cm³/s, mass in cm⁴/s):
//
//	q_e = π² r_e² c               (eq. 4)
//	m_e = (4π/3) r_e³ c           (eq. 6)
//	L   = (4/5) m_e c r_e         (eq. 15, angular momentum)
//	μ   = (3/4) q_e c r_e         (eq. 16, magnetic moment)
//
// (page images checked: eqs. 10–12 give 1 g = 8.4, 8.8, 8.6 × 10⁶ cm⁴/s;
// eq. 30 gives the proton charge q_p = 12 q_e, p. 118)
// and, for hydrogen, an electron orbit of radius ≈ 24 r_e at 0.69 c
// (eq. 35). Each prediction is compared with the measured value (CODATA
// 2018, NIST). Standard counterpart: the Dirac electron with QED
// corrections, which gives the measured values.
//
//	go run ./examples/space_vortex
package main

import (
	"fmt"
	"math"
)

// Measured constants, Gaussian units (CODATA 2018 recommended values,
// NIST; Lyman-α from the NIST Atomic Spectra Database).
const (
	C       = 2.99792458e10    // cm/s
	E       = 4.80320471e-10   // statC
	Me      = 9.1093837015e-28 // g
	Hbar    = 1.054571817e-27  // erg s
	MuB     = 9.2740100783e-21 // erg/G, eħ/2m_e c
	GE      = 2.00231930436256 // |g| of the free electron
	A0      = 5.29177210903e-9 // cm, Bohr radius
	LymanA  = 121.567e-7       // cm, hydrogen 1s–2p
	RClass  = 2.8179403262e-13 // cm, classical electron radius
	AlphaFS = 7.2973525693e-3
)

// Model holds Tewari's electron as he states it.
type Model struct {
	Re float64 // void radius [cm]
}

// FromCharge is the model with r_e fixed from the measured charge by
// eq. 4 read numerically (q in "cm³/s" = statC): r_e = √(q/π²c).
func FromCharge() Model { return Model{Re: math.Sqrt(E / (math.Pi * math.Pi * C))} }

// Charge is eq. 4 — an input, not a prediction, once r_e is fitted to it.
func (m Model) Charge() float64 { return math.Pi * math.Pi * m.Re * m.Re * C }

// MassUnits is eq. 6 in the model's own unit (cm⁴/s); the paper defines
// its gram by dividing this by the measured mass, so mass is an input too.
func (m Model) MassUnits() float64 { return 4 * math.Pi / 3 * m.Re * m.Re * m.Re * C }

// AngularMomentum is eq. 15 with the measured mass in grams [erg s].
func (m Model) AngularMomentum() float64 { return 0.8 * Me * C * m.Re }

// Moment is eq. 16 read in Gaussian units. As printed, μ = (3/4)q c r_e
// has the dimensions of an SI moment (current × area); the Gaussian moment
// of the same current loop is that divided by c, μ = (3/4) q r_e.
func (m Model) Moment() float64 { return 0.75 * E * m.Re }

// GFactor is g = 2 m_e c μ/(e L), from eqs. 15–16: (15/8), independent of
// r_e and of the unit reading, since both use the same q, m and r_e.
func (m Model) GFactor() float64 { return 2 * Me * C * m.Moment() / (E * m.AngularMomentum()) }

// Alpha is the fine-structure constant e²/(2ε₀hc) evaluated with the
// model's own constants (q = π²r_e²c, ε₀ = π/2c, h = (4/5)m_e c r_e with
// m_e = (4π/3)r_e³c): 15π²/16, independent of r_e (Tewari, Universal
// Principles of Space and Matter, 2017, eq. 5.2, p. 130).
func (m Model) Alpha() float64 {
	q := m.Charge()
	eps0 := math.Pi / (2 * C)
	h := 0.8 * m.MassUnits() * C * m.Re
	return q * q / (2 * eps0 * h * C)
}

// HydrogenOrbit returns the orbit radius 24 r_e, speed 0.69c and the
// frequency v/2πr a classical orbiting charge radiates at (eq. 35).
func (m Model) HydrogenOrbit() (r, v, f float64) {
	r, v = 24*m.Re, 0.69*C
	return r, v, v / (2 * math.Pi * r)
}

// Check is one comparison.
type Check struct {
	Name           string
	Model, Measure float64
	Unit           string
	Note           string
}

// Ratio is model / measured.
func (c Check) Ratio() float64 { return c.Model / c.Measure }

// Checks runs the electron test bench on the model.
func (m Model) Checks() []Check {
	_, _, fOrb := m.HydrogenOrbit()
	fLy := C / LymanA
	return []Check{
		{"charge (eq. 4)", m.Charge(), E, "statC", "input: r_e is fitted to e"},
		{"void radius vs classical radius", m.Re, RClass, "cm", "paper: ≈142× larger; no measured size, scattering bounds ≲1e-18 cm"},
		{"angular momentum (eq. 15) vs spin ħ/2", m.AngularMomentum(), Hbar / 2, "erg s", "spin ½"},
		{"magnetic moment (eq. 16) vs measured", m.Moment(), GE / 2 * MuB, "erg/G", "|μ_e| = (g/2) μ_B"},
		{"g-factor (eqs. 15–16)", m.GFactor(), GE, "", "model gives 15/8 for any r_e"},
		{"proton charge (eq. 30) vs measured e", 12 * m.Charge(), E, "statC", "atoms neutral to ~1e-21 e; model neutralises hydrogen by moments (eqs. 32–35)"},
		{"g-factor with the ring moment ½qvρ (eq. 16 ÷ 2)", m.GFactor() / 2, GE, "", "p. 115 integrates dq·v·ρ; a ring current's moment is ½·dq·v·ρ"},
		{"fine-structure constant (UPSM eq. 5.2: 15π²/16)", m.Alpha(), AlphaFS, "", "book p. 130 attributes 1/137 to experimental error; measured to 1.5e-10"},
		{"density limit (UPSM §2.8, p. 54) vs white dwarf Sirius B", 3.42e3, 2e6, "g/cm³", "Sirius B mean density ~2e6 g/cm³ (M ≈ 1.02 M☉, R ≈ 0.008 R☉; order of magnitude)"},
		{"hydrogen orbit radius (eq. 35) vs Bohr radius", 24 * m.Re, A0, "cm", ""},
		{"hydrogen orbit speed vs ground-state αc", 0.69 * C, AlphaFS * C, "cm/s", ""},
		{"hydrogen radiation: orbit frequency vs Lyman-α", fOrb, fLy, "Hz", "classical orbit radiates at v/2πr"},
	}
}

func main() {
	m := FromCharge()
	fmt.Println("Tewari, \"Structural relation between the vacuum space and the electron\", Physics Essays 31, 108 (2018)")
	fmt.Printf("void radius from the measured charge (eq. 4): r_e = %.4g cm\n\n", m.Re)
	fmt.Printf("%-48s %12s %12s %9s  %s\n", "quantity", "model", "measured", "ratio", "note")
	for _, c := range m.Checks() {
		fmt.Printf("%-48s %12.4g %12.4g %9.4g  %s\n", c.Name, c.Model, c.Measure, c.Ratio(), c.Note)
	}
	fmt.Println("\nZeeman effect: not treated in the paper. The test bench needs the model's prediction for the splitting")
	fmt.Println("and polarisation of a spectral line in a field B (measured: Δω = eB/2m_e c per unit of m_l, σ± circular).")
}
