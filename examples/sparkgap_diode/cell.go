// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/material"
)

// Electrode is one electrode material in water.
type Electrode struct {
	Name string
	// Valve metals grow an insulating anodic oxide when made the anode
	// (positive): current is blocked that way until the oxide sparks.
	Valve bool
	// E° of the metal's ion/metal couple versus the standard hydrogen
	// electrode [V] (CRC Handbook, electrochemical series).
	E0 float64
	// Fraction of the plate that touches the water with the working metal
	// (silver tips on a copper plate); 0 means the whole plate.
	Fraction float64
	// Silver layer thickness [m] when the working surface is a silver
	// coating or tip on another metal; 0 for solid electrodes.
	Silver float64
	Base   string // the metal under the silver
}

// active is the working fraction of the plate.
func (e Electrode) active() float64 {
	if e.Fraction <= 0 {
		return 1
	}
	return e.Fraction
}

// Electrodes compared.
var (
	Al = Electrode{Name: "aluminium", Valve: true, E0: material.Aluminium.StandardPotential}
	Pb = Electrode{Name: "lead", E0: -0.126}
	Fe = Electrode{Name: "iron", E0: -0.447}
	Cu = Electrode{Name: "copper", E0: material.Copper.StandardPotential}
	// Silver surfaces: E°(Ag⁺/Ag) = +0.7996 V (CRC electrochemical series).
	// Declared geometry: silver tips covering 10 % of a copper plate, 1 mm
	// thick (solid contact tips; the rest of the copper insulated), and a
	// 10 µm silver plating over a whole iron plate.
	CuAgTips = Electrode{Name: "copper with silver tips", E0: 0.7996, Fraction: 0.1, Silver: 1e-3, Base: "copper"}
	FeAgCoat = Electrode{Name: "silver-plated iron", E0: 0.7996, Silver: 10e-6, Base: "iron"}
)

// Silver (CRC Handbook): density 10,490 kg/m³, molar mass 107.8682 g/mol;
// one electron per atom (Ag → Ag⁺).
const (
	silverDensity = 10490.0
	silverMolar   = 107.8682e-3
)

// SilverLife is how long a silver working surface lasts when it is the
// positive electrode (forward conduction: current leaves it into the
// water) at mean current i [A], if all that charge dissolves silver
// (Faraday's law, worst case: in practice part of it evolves oxygen or
// forms an oxide/phosphate film). Returns the silver mass [kg], the charge
// it can pass [C] and the time [s]; zero without silver.
func (c Cell) SilverLife(e Electrode, i float64) (mass, charge, seconds float64) {
	if e.Silver == 0 || i <= 0 {
		return 0, 0, 0
	}
	mass = silverDensity * e.Silver * c.Area * e.active()
	charge = mass / silverMolar * material.Faraday
	return mass, charge, charge / i
}

// Cell is a water cell between two electrode plates: electrode A and
// electrode K, the current counted from A to K through the water.
//
// Model (each number declared or cited in DefaultCell):
//
//   - water between the plates: resistance R = ρ·gap/area;
//   - electrolysis needs a few volts before current flows (decomposition
//     onset Vd, the 1.23 V reversible voltage plus overpotentials);
//   - two different metals also form a battery, with open-circuit EMF up to
//     E°(K) − E°(A) (reported by Galvanic; a weak chemical source, not a
//     rectifier, so it is left out of the I–V curve — declared);
//   - a valve metal as the positive electrode blocks with its anodic oxide
//     (leakage Rleak) until the oxide breaks down at the sparking voltage
//     Vspark; above it, sparks (plasma at the oxide) carry the current with
//     resistance Rspark. This is the aluminium electrolytic rectifier of
//     Pollak and Nodon (Encyclopaedia Britannica 1911, "Transformers",
//     vol. 27: about 140 V per cell in phosphate solution, efficiency about
//     75–80 %).
//   - the oxide film is also a capacitor, ε_r·ε0·area/(k·Vform), which is
//     why these rectifiers fail at high frequency (the 1911 article notes
//     it): a parallel capacitance per valve electrode.
type Cell struct {
	A, K        Electrode
	Area, Gap   float64 // plate area [m²], gap [m]
	Rho         float64 // electrolyte resistivity [Ω·m]
	Vd          float64 // decomposition onset [V]
	Vspark      float64 // anodic-oxide sparking voltage [V]
	LeakPerArea float64 // oxide leakage current density at Vspark [A/m²]
	Rspark      float64 // resistance while sparking [Ω]
	EpsOx, KOx  float64 // oxide relative permittivity, growth constant [m/V]
}

// DefaultCell: 10 cm × 10 cm plates 1 cm apart in a phosphate solution.
// Declared: ρ = 0.1 Ω·m (a strong salt solution), Vd = 2 V, leakage
// 10 A/m² (1 mA/cm²) at the sparking voltage, Rspark = 2 Ω; Vspark = 140 V
// (Pollak's cells, EB 1911); oxide ε_r = 9 and 1.3 nm per formation volt
// (typical barrier anodic alumina, declared).
func DefaultCell(a, k Electrode) Cell {
	return Cell{A: a, K: k, Area: 0.01, Gap: 0.01, Rho: 0.1, Vd: 2, Vspark: 140, LeakPerArea: 10, Rspark: 2, EpsOx: 9, KOx: 1.3e-9}
}

// R is the water's resistance [Ω], set by the smaller working area (a
// declared simplification: the current spreads from the smaller electrode).
func (c Cell) R() float64 { return c.Rho * c.Gap / (c.Area * math.Min(c.A.active(), c.K.active())) }

// Galvanic is the pair's open-circuit battery EMF, E°(K) − E°(A) with A
// as the oxidized metal (zero when a valve metal is involved: its oxide
// blocks the reaction — declared). Reported, not part of the I–V curve.
func (c Cell) Galvanic() float64 {
	if c.A.Valve || c.K.Valve {
		return 0
	}
	// A is the anode (oxidized) when E°(A) < E°(K): current flows A → K
	// inside the cell, i.e. the cell drives current out of K externally.
	return material.CellEMF(material.Metal{StandardPotential: c.A.E0}, material.Metal{StandardPotential: c.K.E0})
}

// Current is the static current from A to K at voltage v = V(A) − V(K).
func (c Cell) Current(v float64) float64 {
	e := v
	r := c.R()
	valveArea := c.Area
	if c.A.Valve && !c.K.Valve {
		valveArea *= c.A.active()
	} else if c.K.Valve && !c.A.Valve {
		valveArea *= c.K.active()
	}
	leak := c.LeakPerArea * valveArea / c.Vspark // leakage conductance [S]
	block := func(v float64) float64 {           // the valve electrode positive by v > 0
		if v <= c.Vspark {
			return v * leak
		}
		return c.Vspark*leak + (v-c.Vspark)/(c.Rspark+r)
	}
	conduct := func(v float64) float64 { // ordinary electrolysis, v > 0
		if v <= c.Vd {
			return 0
		}
		return (v - c.Vd) / r
	}
	switch {
	case c.A.Valve && c.K.Valve:
		if e >= 0 {
			return block(e)
		}
		return -block(-e)
	case c.A.Valve: // A positive → blocked
		if e >= 0 {
			return block(e)
		}
		return -conduct(-e)
	case c.K.Valve: // K positive (e < 0) → blocked
		if e >= 0 {
			return conduct(e)
		}
		return -block(-e)
	}
	if e >= 0 {
		return conduct(e)
	}
	return -conduct(-e)
}

// FilmC is the anodic oxide's capacitance [F] (0 without a valve metal).
func (c Cell) FilmC() float64 {
	if !c.A.Valve && !c.K.Valve {
		return 0
	}
	return c.EpsOx * material.Eps0 * c.Area / (c.KOx * c.Vspark)
}

// Name describes the pair.
func (c Cell) Name() string { return c.A.Name + " | water | " + c.K.Name }

// curve is the piecewise-linear i(v) for kronspice, sampled densely
// around the breaks.
func (c Cell) curve(vmax float64) []float64 {
	pts := map[float64]bool{-vmax: true, vmax: true, 0: true}
	for _, b := range []float64{c.Vd, c.Vspark} {
		for _, s := range []float64{-1, 1} {
			if x := s * b; math.Abs(x) < vmax {
				pts[x] = true
			}
		}
	}
	var vs []float64
	for v := range pts {
		vs = append(vs, v)
	}
	sortFloats(vs)
	var out []float64
	for _, v := range vs {
		out = append(out, v, c.Current(v))
	}
	return out
}

func sortFloats(x []float64) {
	for i := 1; i < len(x); i++ {
		for j := i; j > 0 && x[j] < x[j-1]; j-- {
			x[j], x[j-1] = x[j-1], x[j]
		}
	}
}

// Add puts the cell between nodes a and k of a builder (as a nonlinear
// resistor and, for valve metals, the oxide's capacitance).
func (c Cell) Add(b *circuit.Builder, name, a, k string) {
	b.NonlinearVC(name, a, k, c.curve(20*c.Vspark))
	if cf := c.FilmC(); cf > 0 {
		b.Capacitor("C"+name, a, k, cf)
	}
}

// Describe prints the cell's parameters.
func (c Cell) Describe() string {
	return fmt.Sprintf("%s: plates %.0f cm² (working area %.0f cm²), %.0f cm apart, water %.3g Ω, decomposition onset %.1f V, galvanic %.3f V, oxide sparks at %.0f V, film %.3g µF",
		c.Name(), c.Area*1e4, c.Area*1e4*math.Min(c.A.active(), c.K.active()), c.Gap*100, c.R(), c.Vd, c.Galvanic(), c.Vspark, c.FilmC()*1e6)
}
