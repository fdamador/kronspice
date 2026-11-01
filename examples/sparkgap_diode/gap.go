// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// PolarGap is a point–plane spark gap in pure water. Above the oxide's
// volts, a spark in water is the conduction path, and what makes it one-way
// is the electrode shape, not the metal: a sharp positive point starts
// breakdown at a lower average field than a sharp negative point.
//
//   - Positive point (the point is the anode): the Sandia water relation
//     E·τ^(1/3) = 0.135 MV/cm·µs^(1/3) for a significantly field-enhanced
//     anode (W. A. Stygar et al., Phys. Rev. ST Accel. Beams 9, 070401,
//     2006), with τ the time the voltage stays above 63 % of its peak. It
//     was measured for τ = 0.011–0.6 µs; at 60 Hz τ is milliseconds, so the
//     value used here is an EXTRAPOLATION and is labelled so.
//   - Negative point: NegRatio times the positive value. 2 is declared, from
//     the factor of about two between the negative and positive constants
//     in J. C. Martin's water relations as widely quoted (to be checked
//     against J. C. Martin on Pulsed Power, 1996, before it is relied on).
//   - Before breakdown pure water (0.055 µS/cm) leaks through the spreading
//     resistance of the point, 1/(2π·σ·a) for a hemispherical tip of
//     radius a; after it, an arc of declared voltage and resistance.
type PolarGap struct {
	Gap, Tip, Sigma float64 // gap [m], tip radius [m], water conductivity [S/m]
	Tau             float64 // effective time [s]
	NegRatio        float64
	Varc, Rarc      float64 // arc voltage [V] and resistance [Ω]
	Iext            float64 // the arc goes out when its current falls below this [A]
	Hold            float64 // a fresh arc lasts at least this long (channel de-ionization) [s]
	// Recovery: after an arc goes out the gap's breakdown voltages come back
	// as 1 − exp(−t/Recovery) (hot channel cooling; 0 = instant). Declared.
	Recovery float64
	// Fixed, when > 0, replaces the water relation with a fixed forward
	// breakdown voltage (an air or mist gap).
	Fixed float64
	// Floor: a gap whose channel is still hot re-strikes at k·Vb but never
	// below this. The channel is vapour (or air), and no gas breaks down
	// below its Paschen minimum: water vapour 480 V (Škoro et al., Phys.
	// Rev. E 84, 055401, 2011), air about 330 V (Raizer 1991).
	Floor float64
}

// restrike is the breakdown threshold k·vb of a recovering gap.
func (g PolarGap) restrike(k, vb float64) float64 {
	return math.Max(k*vb, math.Max(g.Floor, 2*g.Varc))
}

// AtFrequency returns the gap with its effective time set for a sine of
// frequency f (time above 63 % of the peak).
func (g PolarGap) AtFrequency(f float64) PolarGap {
	g.Tau = 2 * (math.Pi/2 - math.Asin(0.63)) / (2 * math.Pi * f)
	return g
}

// DefaultGap: 3 mm, 0.5 mm tip, pure water, 60 Hz (τ = time above 63 %
// of a sine's peak), negative ×2, arc 50 V + 0.5 Ω that goes out below
// 0.1 A after lasting at least 10 µs (declared).
func DefaultGap() PolarGap {
	w := 2 * math.Pi * 60
	tau := 2 * (math.Pi/2 - math.Asin(0.63)) / w
	return PolarGap{Gap: 3e-3, Tip: 0.5e-3, Sigma: 5.5e-6, Tau: tau, NegRatio: 2, Varc: 50, Rarc: 0.5, Iext: 0.1, Hold: 10e-6, Floor: 480}
}

// Vpos is the breakdown voltage with the point positive [V].
func (g PolarGap) Vpos() float64 {
	if g.Fixed > 0 {
		return g.Fixed
	}
	e := 0.135 / math.Cbrt(g.Tau*1e6) * 1e8 // V/m
	return e * g.Gap
}

// Vneg is the breakdown voltage with the point negative [V].
func (g PolarGap) Vneg() float64 { return g.NegRatio * g.Vpos() }

// Rleak is the water's spreading resistance from the tip [Ω].
func (g PolarGap) Rleak() float64 { return 1 / (2 * math.Pi * g.Sigma * g.Tip) }

// Describe prints the gap.
func (g PolarGap) Describe() string {
	return fmt.Sprintf("point–plane gap %.0f mm in pure water: breaks down at +%.0f V (point positive, EXTRAPOLATED from Stygar 2006) and −%.0f V (point negative, ×%.0f declared); leaks through %.3g MΩ; arc %.0f V + %.1f Ω",
		g.Gap*1e3, g.Vpos(), g.Vneg(), g.NegRatio, g.Rleak()/1e6, g.Varc, g.Rarc)
}
