// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/material"
)

// ErrOutsideModelValidity means a geometry is outside the range where an
// analytic formula is accurate. It is returned instead of extrapolating
// (CLAUDE.md §3d principle 4); the GA treats it as an invalid design.
var ErrOutsideModelValidity = errors.New("em: outside model validity")

// Coil is a single-layer solenoid, wound ordinarily or bifilar as in Tesla's
// U.S. Patent 512,340 (two wires side by side, the end of the first joined
// to the start of the second, so the turns alternate between the halves).
//
// These are analytic, quasi-static formulas (CLAUDE.md §3c: inductor row,
// flux linkage of 𝔄; capacitor row, field confined between adjacent
// turns), a first approximation until em/extract solves the fields (M9i).
type Coil struct {
	Turns         int     // total turns (both halves for a bifilar coil)
	Diameter      float64 // mean turn diameter [m]
	WireDiameter  float64 // conductor diameter [m]
	Pitch         float64 // centre-to-centre spacing of adjacent turns [m], > WireDiameter
	InsulationEps float64 // relative permittivity between adjacent turns
	Bifilar       bool
	Conductor     material.Metal
}

// Length is the winding length [m].
func (c Coil) Length() float64 { return float64(c.Turns) * c.Pitch }

// Validate checks the geometry against the formulas' ranges: Wheeler's
// formula is within about 1 % for length > 0.4·diameter (Wheeler 1928),
// and the adjacent-wire capacitance needs pitch > wire diameter.
func (c Coil) Validate() error {
	switch {
	case c.Turns < 2:
		return fmt.Errorf("%w: coil needs at least 2 turns", ErrOutsideModelValidity)
	case !(c.Diameter > 0) || !(c.WireDiameter > 0) || !(c.InsulationEps >= 1):
		return fmt.Errorf("%w: coil dimensions must be positive", ErrOutsideModelValidity)
	case !(c.Pitch > c.WireDiameter):
		return fmt.Errorf("%w: pitch %.3g m must exceed the wire diameter %.3g m", ErrOutsideModelValidity, c.Pitch, c.WireDiameter)
	case c.Length() < 0.4*c.Diameter:
		return fmt.Errorf("%w: winding length %.3g m < 0.4 × diameter %.3g m (Wheeler's formula)", ErrOutsideModelValidity,
			c.Length(), c.Diameter)
	case c.WireDiameter > 0.2*c.Diameter:
		return fmt.Errorf("%w: wire diameter is more than 1/5 of the coil diameter", ErrOutsideModelValidity)
	}
	return nil
}

// Inductance [H] by Wheeler's formula, L[µH] = r²N²/(9r + 10ℓ) with r and
// ℓ in inches (H. A. Wheeler, "Simple inductance formulas for radio coils",
// Proc. IRE 16(10), 1928). A bifilar coil whose halves are joined in series
// carries the same current in the same direction in every turn, so its
// inductance is that of an ordinary coil with the same turns and pitch.
func (c Coil) Inductance() float64 {
	const inch = 0.0254
	r, l := c.Diameter/2/inch, c.Length()/inch
	n := float64(c.Turns)
	return r * r * n * n / (9*r + 10*l) * 1e-6
}

// TurnCapacitance is the capacitance between two adjacent turns [F]:
// two parallel round wires of diameter d at centre spacing p, length πD,
// C = π·ε·πD / acosh(p/d) (the two-wire line formula; e.g. Ramo, Whinnery
// & Van Duzer, Fields and Waves in Communication Electronics, ch. 4).
func (c Coil) TurnCapacitance() float64 {
	return math.Pi * material.Eps0 * c.InsulationEps * math.Pi * c.Diameter / math.Acosh(c.Pitch/c.WireDiameter)
}

// TurnVoltage is the voltage between adjacent turns for a voltage v across
// the coil, in the idealized uniform model of Tesla's patent: v/N for an
// ordinary winding, v/2 for a bifilar one (each turn lies next to the turn
// half a coil away). Tesla's example (US 512,340, p. 1): 1,000 turns at
// 100 V give 0.1 V and 50 V.
func (c Coil) TurnVoltage(v float64) float64 {
	if c.Bifilar {
		return v / 2
	}
	return v / float64(c.Turns)
}

// SelfCapacitance is the terminal capacitance equivalent to the energy
// stored between adjacent turns [F]: (N − 1) gaps, each holding
// ½·C_tt·(TurnVoltage)², referred to the terminal voltage. Ordinary:
// (N−1)·C_tt/N²; bifilar: (N−1)·C_tt/4. The ratio N²/4 is Tesla's
// 250,000 for N = 1,000. (For ordinary coils this idealized model
// underestimates the measured self-capacitance, which Medhurst's empirical
// formula gives; Medhurst, Wireless Engineer 24, 1947.)
func (c Coil) SelfCapacitance() float64 {
	n := float64(c.Turns)
	ct := c.TurnCapacitance()
	if c.Bifilar {
		return (n - 1) * ct / 4
	}
	return (n - 1) * ct / (n * n)
}

// TurnToTurnEnergy is the electrostatic energy between adjacent turns for a
// voltage v across the coil [J].
func (c Coil) TurnToTurnEnergy(v float64) float64 {
	return 0.5 * c.SelfCapacitance() * v * v
}

// WireLength is the conductor length [m].
func (c Coil) WireLength() float64 { return float64(c.Turns) * math.Pi * c.Diameter }

// Resistance is the conductor resistance at frequency f [Ω]: the DC value,
// raised by the skin effect when the skin depth δ = √(ρ/(π f µ0)) is small
// against the wire, R ≈ R_dc·d/(4δ) (Terman, Radio Engineers' Handbook,
// 1943, sec. 2). Proximity effect between turns is not included.
func (c Coil) Resistance(f float64) float64 {
	area := math.Pi * c.WireDiameter * c.WireDiameter / 4
	rdc := c.Conductor.Resistivity * c.WireLength() / area
	if f <= 0 {
		return rdc
	}
	delta := math.Sqrt(c.Conductor.Resistivity / (math.Pi * f * material.Mu0))
	return rdc * math.Max(1, c.WireDiameter/(4*delta))
}

// SelfResonance is 1/(2π√(L·C_self)) [Hz].
func (c Coil) SelfResonance() float64 {
	return 1 / (2 * math.Pi * math.Sqrt(c.Inductance()*c.SelfCapacitance()))
}

// Volume is the volume of the cylinder enclosing the winding [m³].
func (c Coil) Volume() float64 {
	r := (c.Diameter + c.WireDiameter) / 2
	return math.Pi * r * r * (c.Length() + c.WireDiameter)
}

// ParallelPlate is the capacitance of plates of area a at gap g with
// relative permittivity eps, C = ε0·ε·a/g, without fringing [F]. Fringing is
// negligible when the gap is small against the plate size; a gap above a
// tenth of √a is rejected.
func ParallelPlate(a, g, eps float64) (float64, error) {
	if !(a > 0) || !(g > 0) || !(eps >= 1) {
		return 0, fmt.Errorf("%w: plate area, gap and permittivity must be positive", ErrOutsideModelValidity)
	}
	if g > 0.1*math.Sqrt(a) {
		return 0, fmt.Errorf("%w: gap %.3g m is not small against the plate size %.3g m", ErrOutsideModelValidity, g, math.Sqrt(a))
	}
	return material.Eps0 * eps * a / g, nil
}
