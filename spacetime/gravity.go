// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import (
	"errors"
	"fmt"
	"math"
)

// Constants (CODATA 2018 for c and G; IAU 2015 nominal GM values; IERS
// conventions 2010 for Earth).
const (
	C          = 299792458.0    // m/s, exact
	G          = 6.67430e-11    // m³/(kg·s²)
	GMEarth    = 3.986004418e14 // m³/s² (IERS 2010)
	REarthEq   = 6378137.0      // m, equatorial (WGS-84)
	OmegaEarth = 7.2921150e-5   // rad/s (IERS 2010)
	GMSun      = 1.3271244e20   // m³/s² (IAU 2015 nominal)
	RSun       = 6.957e8        // m (IAU 2015 nominal)
)

// Errors.
var (
	ErrInsideHorizon = errors.New("spacetime: at or inside the event horizon")
	ErrNoOrbit       = errors.New("spacetime: no circular orbit (r ≤ 1.5 r_s, the photon sphere)")
	ErrSuperluminal  = errors.New("spacetime: speed at or above c")
	ErrBadInput      = errors.New("spacetime: invalid input")
)

// SchwarzschildRadius r_s = 2GM/c² [m].
func SchwarzschildRadius(gm float64) float64 { return 2 * gm / (C * C) }

// StaticRate is dτ/dt for a clock at rest at radius r outside mass gm:
// √(1 − r_s/r) (Misner, Thorne, Wheeler, Gravitation, §25.2).
func StaticRate(gm, r float64) (float64, error) {
	return MovingRate(gm, r, 0)
}

// MovingRate is dτ/dt for a clock at radius r moving with coordinate
// speed v tangentially (circular motion, orbit or rotating support):
// √(1 − r_s/r − v²/c²) in the Schwarzschild exterior (MTW §25.4; for a
// free circular orbit v² = GM/r and the rate is √(1 − 3r_s/2r)).
func MovingRate(gm, r, v float64) (float64, error) {
	if gm < 0 || r <= 0 || v < 0 || math.IsNaN(gm+r+v) {
		return 0, fmt.Errorf("%w: GM %g, r %g, v %g", ErrBadInput, gm, r, v)
	}
	rs := SchwarzschildRadius(gm)
	if r <= rs {
		return 0, fmt.Errorf("%w: r %g m ≤ r_s %g m", ErrInsideHorizon, r, rs)
	}
	x := 1 - rs/r - v*v/(C*C)
	if x <= 0 {
		return 0, fmt.Errorf("%w: v %g m/s at r %g m", ErrSuperluminal, v, r)
	}
	return math.Sqrt(x), nil
}

// OrbitRate is dτ/dt on a free circular orbit at radius r.
func OrbitRate(gm, r float64) (float64, error) {
	if r <= 1.5*SchwarzschildRadius(gm) {
		return 0, fmt.Errorf("%w: r %g m", ErrNoOrbit, r)
	}
	return MovingRate(gm, r, math.Sqrt(gm/r))
}

// Shift is the fractional frequency shift of a signal from an emitter
// clock running at rateE to a receiver clock running at rateR, both rates
// dτ/dt against the same coordinate time in a static field (light
// conserves its coordinate frequency): f_R/f_E − 1 = rateE/rateR − 1.
// Doppler shifts from relative motion along the line of sight are not
// included. Rates within ~10⁻¹² of 1 lose precision in float64; use
// StaticShift for small shifts between static clocks.
func Shift(rateE, rateR float64) float64 { return rateE/rateR - 1 }

// StaticShift is the fractional frequency shift of a signal sent from a
// static clock at radius rE to one at rR, f_R/f_E − 1 =
// √((1 − r_s/rE)/(1 − r_s/rR)) − 1, evaluated without cancellation so it
// stays accurate for shifts of 10⁻¹⁵ and below (Pound–Rebka scale). It is
// positive (blue shift) when the signal falls inward (rR < rE).
func StaticShift(gm, rE, rR float64) (float64, error) {
	if _, err := StaticRate(gm, rE); err != nil {
		return 0, err
	}
	if _, err := StaticRate(gm, rR); err != nil {
		return 0, err
	}
	rs := SchwarzschildRadius(gm)
	x := (rs/rR - rs/rE) / (1 - rs/rR) // ratio − 1
	return math.Expm1(0.5 * math.Log1p(x)), nil
}

// Tidal is the tidal (stretching) acceleration across a length l aligned
// radially at radius r: 2GM·l/r³ [m/s²] (Newtonian; MTW §1.6).
func Tidal(gm, r, l float64) float64 { return 2 * gm * l / (r * r * r) }
