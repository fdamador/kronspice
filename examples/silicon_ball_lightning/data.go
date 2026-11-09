// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Observations and material data.

// Observation is one measured or reported ball (or ball-like glow).
type Observation struct {
	Source    string
	Kind      string
	DiameterM [2]float64 // range [m]
	LifetimeS [2]float64 // range [s]
	Note      string
}

// Observations on record (Cen, Yuan & Xue, PRL 112, 035001, 2014, and the
// works it cites; Abrahamson & Dinniss, Nature 403, 519, 2000).
var Observations = []Observation{
	{"Abrahamson & Dinniss 2000 (the 'average' ball)", "natural, eyewitness", [2]float64{0.3, 0.3}, [2]float64{10, 10},
		"luminosity similar to a 100 W lamp; ends by explosion or fading"},
	{"Cen, Yuan & Xue 2014 (Qinghai, 23 July 2012)", "natural, spectrographs at 0.9 km", [2]float64{5, 5}, [2]float64{1.64, 1.64},
		"apparent (luminous) diameter ≈ 5 m in the stable stage; 8.6 m/s; Si I, Fe I, Ca I lines throughout; N I, O I periodic; 99.4 Hz light oscillation; 35 kV line 20 m away"},
	{"Paiva et al. 2007 (PRL 98, 048501)", "laboratory: arc on pure silicon", [2]float64{0.01, 0.04}, [2]float64{2, 5},
		"arc 20–25 V, 100–140 A, 60 Hz; bluish white or orange white"},
	{"Stephan & Massey 2008 (J. Atmos. Sol.-Terr. Phys. 70, 1589)", "laboratory: silicon", [2]float64{1e-4, 1e-3}, [2]float64{0, 0},
		"about 3140 K; luminous balls of 1–4 cm lit by a liquid core of about 1 mm"},
	{"Dikhtyar & Jerby 2006 (PRL 96, 045002)", "laboratory: microwave hot spot", [2]float64{0.03, 0.03}, [2]float64{0.03, 0.04},
		"fireballs last 30–40 ms after the microwave power is turned off"},
}

// Element is a metal (or reduced species) that can be vaporized by an arc
// or a lightning strike and burn back in air.
//
// Enthalpies of formation at 298 K (kJ/mol of oxide) from the NIST-JANAF
// tables / CRC Handbook; melting and boiling points, heats of fusion and
// vaporization from the CRC Handbook. Cp values are rounded means over the
// solid and liquid ranges (approximate: they make the vaporization cost
// uncertain by roughly ±10 %).
type Element struct {
	Name       string
	M          float64 // molar mass [g/mol]
	Oxide      string
	HOxPerMol  float64 // heat released per mol of element burned to the oxide [kJ/mol] (positive)
	Tm, Tb     float64 // melting, boiling point [K]
	HFus, HVap float64 // [kJ/mol]
	CpS, CpL   float64 // mean heat capacities, solid and liquid [J/(mol·K)]
}

// Elements in soil and in welding electrodes (steel core; ferro-silicon,
// ferro-manganese deoxidizers; rutile flux holds Ti only as TiO2).
var Elements = []Element{
	{"Si", 28.086, "SiO2", 910.7, 1687, 3538, 50.21, 359, 25, 27.2},
	{"Fe", 55.845, "Fe2O3", 824.2 / 2, 1811, 3134, 13.81, 340, 33, 46},
	{"Mn", 54.938, "MnO", 385.2, 1519, 2334, 12.91, 221, 32, 46},
	{"Ca", 40.078, "CaO", 634.9, 1115, 1757, 8.54, 154.7, 28, 31},
	{"Al", 26.982, "Al2O3", 1675.7 / 2, 933.5, 2792, 10.71, 294, 28, 31.8},
	{"Ti", 47.867, "TiO2", 944.0, 1941, 3560, 14.15, 425, 30, 47},
}

// Reduced silicon species Abrahamson & Dinniss name besides Si: SiO(g)
// (ΔHf −100.4 kJ/mol) and SiC (ΔHf −73.2 kJ/mol, CO2 −393.5): heat on
// burning to SiO2 (+ CO2).
var Species = []struct {
	Name    string
	M       float64
	HPerMol float64
}{
	{"Si → SiO2", 28.086, 910.7},
	{"SiO(g) → SiO2", 44.085, 910.7 - 100.4},
	{"SiC → SiO2 + CO2", 40.096, 910.7 + 393.5 - 73.2},
}
