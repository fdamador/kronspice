// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package material holds cited physical constants and material properties
// used to turn physical designs into circuit parameters. Every value names
// its source; values are SI and near room temperature unless stated.
//
// Leaf package: it imports nothing from kronspice.
package material

import "math"

// Exact SI constants (CODATA 2018 / SI 2019) and µ0, ε0 (CODATA 2018).
const (
	Mu0     = 1.25663706212e-6 // vacuum permeability [H/m], CODATA 2018
	Eps0    = 8.8541878128e-12 // vacuum permittivity [F/m], CODATA 2018
	Faraday = 96485.33212      // Faraday constant [C/mol] (exact: N_A·e)
)

// Metal is a conductor's properties.
type Metal struct {
	Name string
	// Seebeck is the Seebeck coefficient relative to platinum at room
	// temperature [V/K]. Only differences between two metals are physical
	// for a thermocouple; using the same reference for both makes the
	// difference exact to the table's accuracy.
	Seebeck float64
	// Resistivity at 20 °C [Ω·m].
	Resistivity float64
	// ThermalConductivity near 300 K [W/(m·K)].
	ThermalConductivity float64
	Density             float64 // [kg/m³]
	MolarMass           float64 // [kg/mol]
	// Valence is the number of electrons per atom in the electrode reaction
	// M → M^(n+) + n e⁻, and StandardPotential the standard reduction
	// potential E° of M^(n+)/M versus the standard hydrogen electrode [V].
	Valence           int
	StandardPotential float64
	Source            string
}

// Sources.
const (
	srcSeebeck = "Seebeck coefficients relative to Pt: The Seebeck Coefficient, Electronics Cooling (Nov 2006), " +
		"as tabulated in Wikipedia 'Seebeck coefficient' (Bi −72, Al +3.5, Cu +6.5 µV/K)"
	srcResist = "Resistivity at 20 °C: CRC Handbook of Chemistry and Physics, 'Electrical resistivity of pure metals'; " +
		"thermal conductivity at 300 K: CRC Handbook, 'Thermal conductivity of metals and semiconductors'"
	srcDensity   = "Density and molar mass: CRC Handbook of Chemistry and Physics; IUPAC standard atomic weights"
	srcPotential = "Standard reduction potentials: CRC Handbook, 'Electrochemical series' (P. Vanýsek): " +
		"Al³⁺/Al −1.662 V, Bi³⁺/Bi +0.308 V"
)

// Metals used by the examples.
var (
	Bismuth = Metal{Name: "bismuth", Seebeck: -72e-6, Resistivity: 1.29e-6, ThermalConductivity: 7.97, Density: 9780, MolarMass: 208.98e-3,
		Valence: 3, StandardPotential: 0.308, Source: srcSeebeck + "; " + srcResist + "; " + srcDensity + "; " + srcPotential}
	Aluminium = Metal{Name: "aluminium", Seebeck: 3.5e-6, Resistivity: 2.65e-8, ThermalConductivity: 237, Density: 2700, MolarMass: 26.98e-3,
		Valence: 3, StandardPotential: -1.662, Source: srcSeebeck + "; " + srcResist + "; " + srcDensity + "; " + srcPotential}
	Copper = Metal{Name: "copper", Seebeck: 6.5e-6, Resistivity: 1.68e-8, ThermalConductivity: 401, Density: 8960, MolarMass: 63.546e-3,
		Valence: 2, StandardPotential: 0.342, Source: srcSeebeck + "; " + srcResist + "; " + srcDensity +
			"; CRC electrochemical series Cu²⁺/Cu +0.342 V"}
)

// ThermocoupleEMF is the open-circuit EMF of a junction of metals a and b
// with a temperature difference dT between the junction and the far ends
// [V]: (S_b − S_a)·dT. It is positive when b's free end is positive with
// the junction hotter (dT > 0).
func ThermocoupleEMF(a, b Metal, dT float64) float64 { return (b.Seebeck - a.Seebeck) * dT }

// CellEMF is the standard open-circuit EMF of a galvanic cell with anode a
// (oxidized) and cathode c [V]: E°(c) − E°(a). It is the thermodynamic
// upper bound; a real cell is lower (concentration, passivation — aluminium
// in particular forms an oxide film — and overpotentials).
func CellEMF(anode, cathode Metal) float64 {
	return cathode.StandardPotential - anode.StandardPotential
}

// ChargeCapacity is the charge a mass m of metal delivers when fully
// oxidized [C] (Faraday's law): m·n·F/M.
func ChargeCapacity(m Metal, mass float64) float64 {
	return mass * float64(m.Valence) * Faraday / m.MolarMass
}

// SpectralLine is an atomic emission line.
type SpectralLine struct {
	Element       string
	AirWavelength float64 // [m], standard air
	Intensity     float64 // relative, as tabulated
	Source        string
}

const srcNIST = "NIST Handbook of Basic Atomic Spectroscopic Data (Sansonetti & Martin), strong lines, air wavelengths"

// Strongest neutral-atom lines of the metals used here (NIST handbook,
// table 2 of each element, intensity 1000).
var (
	LineCopper    = SpectralLine{Element: "Cu I", AirWavelength: 324.754e-9, Intensity: 1000, Source: srcNIST}
	LineAluminium = SpectralLine{Element: "Al I", AirWavelength: 396.152e-9, Intensity: 1000, Source: srcNIST}
	LineBismuth   = SpectralLine{Element: "Bi I", AirWavelength: 306.7700e-9, Intensity: 1000, Source: srcNIST}
)

// SpeedOfLight in vacuum [m/s] (exact).
const SpeedOfLight = 299792458.0

// AirIndex is the refractive index of standard dry air at vacuum
// wavenumber σ = 1/λ_vac [µm⁻¹]: Edlén's dispersion formula as revised by
// Birch & Downs (Metrologia 30, 1993; 31, 1994):
// (n − 1)·10⁸ = 8342.54 + 2406147/(130 − σ²) + 15998/(38.9 − σ²).
func AirIndex(sigma float64) float64 {
	s2 := sigma * sigma
	return 1 + 1e-8*(8342.54+2406147/(130-s2)+15998/(38.9-s2))
}

// VacuumWavelength converts the tabulated air wavelength [m] (two fixed-point
// iterations of λ_vac = n(λ_vac)·λ_air, converged to < 1e-12 relative).
func (l SpectralLine) VacuumWavelength() float64 {
	v := l.AirWavelength
	for i := 0; i < 4; i++ {
		v = AirIndex(1e-6/v) * l.AirWavelength
	}
	return v
}

// Frequency is c/λ_vac [Hz].
func (l SpectralLine) Frequency() float64 { return SpeedOfLight / l.VacuumWavelength() }

// Subharmonic returns the whole number n for which f_line/n is nearest to
// f, and the offset (f − f_line/n)/f. For f far below the line, adjacent
// subharmonics are only a fraction 1/n apart, so the offset is at most
// 1/(2n) in magnitude.
func (l SpectralLine) Subharmonic(f float64) (n float64, offset float64) {
	fl := l.Frequency()
	n = math.Max(1, math.Round(fl/f))
	return n, (f - fl/n) / f
}
