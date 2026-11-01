// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package fluid holds liquid media and their acoustics (CLAUDE.md §3h):
// property values with cited sources, and the acoustic loads a transducer
// sees in a liquid. It is a leaf package: it returns numbers, and callers
// wire them into circuits through circuit.Builder.
package fluid

// Medium is a liquid's properties at one state point.
type Medium struct {
	Name        string
	Temperature float64 // [K]
	Pressure    float64 // [Pa]
	Density     float64 // ρ [kg/m³]
	SoundSpeed  float64 // c [m/s]
	Source      string  // where the values come from
}

// Impedance is the specific acoustic impedance ρ·c [Pa·s/m] (rayl).
func (m Medium) Impedance() float64 { return m.Density * m.SoundSpeed }

// Water is ordinary water (H₂O) at 300 K, from the IAPWS-95 formulation's
// computer-program verification values: ρ = 996.556 kg/m³ (p = 0.0992 MPa),
// w = 1501.51914 m/s (W. Wagner, A. Pruß, "The IAPWS formulation 1995 for
// the thermodynamic properties of ordinary water substance for general and
// scientific use," J. Phys. Chem. Ref. Data 31, 387 (2002), Table 7).
var Water = Medium{Name: "H2O", Temperature: 300, Pressure: 0.0992418352e6, Density: 996.556,
	SoundSpeed: 1501.51914, Source: "IAPWS-95, Wagner & Pruß, JPCRD 31, 387 (2002), Table 7, T = 300 K"}

// HeavyWater is D₂O at 300 K from the IAPWS reference equation of state
// for heavy water's verification values: 55.126 mol/dm³ (× 20.0276 g/mol
// = 1104.04 kg/m³), w = 1403.75 m/s at p = 0.0529 MPa (S. Herrig, M. Thol,
// A. H. Harvey, E. W. Lemmon, "A reference equation of state for heavy
// water," J. Phys. Chem. Ref. Data 47, 043102 (2018), Table 7). D₂O is
// denser and slower than H₂O, which shifts a resonator's frequency and the
// bubble's dynamics (§3h).
var HeavyWater = Medium{Name: "D2O", Temperature: 300, Pressure: 0.0529e6, Density: 55.126 * 20.0276,
	SoundSpeed: 1403.75, Source: "Herrig et al., JPCRD 47, 043102 (2018), Table 7, T = 300 K"}
