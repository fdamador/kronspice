// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import (
	"fmt"
	"math"
)

// Fast-rotating machines (homopolar disks, flywheels, rotors, spinning
// hulls). Standard physics: a rotor is a flat-spacetime system seen from
// rotating coordinates (the Langevin metric, ds² = −(1 − ω²r²/c²)c²dt² +
// 2ωr²dφdt + dr² + r²dφ² + dz²). Clocks on it run slow by the speed of the
// rim, light around it shows the Sagnac delay, and in its own frame the
// centrifugal potential acts like gravity (Einstein's equivalence
// principle). The rim speed is bounded by the material long before
// (v/c)² becomes large; BurstSpeed gives that bound.

// RimRate is dτ/dt for a clock riding a rotor at radius r with angular
// speed ω, √(1 − ω²r²/c²) (special relativity; the clock hypothesis that
// the rate depends on speed, not acceleration — tested to 10⁻³ at
// 10¹⁸ g in the CERN muon storage ring, J. Bailey et al., Nature 268,
// 301, 1977).
func RimRate(omega, r float64) (float64, error) {
	v := math.Abs(omega * r)
	if v >= C {
		return 0, fmt.Errorf("%w: rim speed %g m/s", ErrSuperluminal, v)
	}
	return math.Sqrt(1 - v*v/(C*C)), nil
}

// RotorShift is the fractional energy shift ΔE/E of a gamma ray (or any
// signal) emitted at radius rE on a rotor and absorbed at radius rA on the
// same rotor, ΔE/E = k·(ω²/c²)·(rA² − rE²) to order (v/c)². Standard
// prediction k = 1/2 (time dilation of the two clocks; the first-order
// Doppler shifts cancel). The coefficient k is a parameter so the
// reported anomaly (hypothesis H17) can be computed alongside: rotor
// Mössbauer experiments with source at the centre report k = 0.68 ± 0.03
// (A. L. Kholmetskii et al., Phys. Scr. 79, 2009) and
// k = 0.69 ± 0.02 (T. Yarman et al., Can. J. Phys.,
// doi:10.1139/cjp-2015-0063, arXiv:1503.05853); their reanalysis of W. Kündig, Phys. Rev. 129, 2371
// (1963) gives 0.596 ± 0.006 (Phys. Scr. 77, 2008), against
// Kündig's own agreement with 1/2 to 1.1 %. Proposed explanations
// (disputed): clock desynchronization adding 1/6, k = 2/3 (C. Corda,
// arXiv:2203.02282); line broadening from finite slits invalidating the
// count-rate analyses (Y. Friedman, I. Nowik, Phys. Scr. 85, 065702,
// 2012). A positive shift means the absorber sees the line higher
// (blue-shifted) than it emits it.
func RotorShift(omega, rE, rA, k float64) float64 {
	return k * omega * omega * (rA*rA - rE*rE) / (C * C)
}

// SagnacDelay is the arrival-time difference [s] between light sent
// co- and counter-rotating around a closed loop of enclosed area a on a
// platform turning at ω: Δt = 4aω/c² to first order (G. Sagnac, C. R.
// Acad. Sci. 157, 708, 1913; E. J. Post, Rev. Mod. Phys. 39, 475, 1967).
func SagnacDelay(area, omega float64) float64 { return 4 * area * omega / (C * C) }

// SagnacRingExact is the delay for a circular loop of radius r measured on
// the rotor (exact in the lab frame, t± = 2πr/(c ∓ ωr), converted to the
// rim clock by RimRate).
func SagnacRingExact(r, omega float64) (float64, error) {
	v := omega * r
	rate, err := RimRate(omega, r)
	if err != nil {
		return 0, err
	}
	l := 2 * math.Pi * r
	return (l/(C-v) - l/(C+v)) * rate, nil
}

// SagnacSplitting is the frequency difference [Hz] between the co- and
// counter-rotating resonances of a ring resonator of enclosed area a and
// perimeter p, at free-space wavelength λ, turning at ω:
// Δf = 4aω/(λ·p) (ring-laser gyroscope; W. W. Chow et al., Rev. Mod.
// Phys. 57, 61, 1985, eq. 2.8). For a guide whose wave speed differs from
// c the splitting is unchanged to first order (Post 1967 §IV.B).
func SagnacSplitting(area, perimeter, lambda, omega float64) float64 {
	return 4 * area * omega / (lambda * perimeter)
}

// Material is a rotor material for the burst limit.
type Material struct {
	Name     string
	Strength float64 // tensile strength [Pa]
	Density  float64 // [kg/m³]
	Source   string
}

// T1000G is Toray T1000G carbon fibre (fibre only; a composite rotor is
// lower): 6,370 MPa, 1.80 g/cm³ (Toray Composite Materials America, T1000G
// data sheet, 16 Feb 2024).
var T1000G = Material{"Toray T1000G carbon fibre (fibre)", 6.37e9, 1800, "Toray CMA T1000G data sheet, 2024-02-16"}

// BurstSpeed is the rim speed at which a thin rotating hoop reaches its
// tensile strength, v = √(σ/ρ) (hoop stress ρv²; S. Timoshenko,
// Strength of Materials II, 1956, §8). A solid disk bursts lower
// (its peak stress is (3+ν)/8·ρω²R²).
func BurstSpeed(m Material) float64 { return math.Sqrt(m.Strength / m.Density) }

// Rotor summarizes the relativistic size of a spinning machine.
type Rotor struct {
	RimSpeed     float64 // ωr [m/s]
	Beta2        float64 // (ωr/c)²
	RimRate      float64 // dτ/dt on the rim
	LossPerDay   float64 // rim clock lag against the axis [s/day]
	SagnacPerRev float64 // Sagnac delay for one loop at the rim [s]
	BurstMargin  float64 // rim speed / hoop burst speed (> 1: the rotor flies apart)
}

// AnalyzeRotor reports a rotor of radius r at ω made of m.
func AnalyzeRotor(omega, r float64, m Material) (Rotor, error) {
	if !(r > 0) || math.IsNaN(omega) {
		return Rotor{}, fmt.Errorf("%w: r %g, ω %g", ErrBadInput, r, omega)
	}
	rate, err := RimRate(omega, r)
	if err != nil {
		return Rotor{}, err
	}
	v := math.Abs(omega * r)
	b2 := v * v / (C * C)
	lag := -math.Expm1(0.5 * math.Log1p(-b2)) // 1 − rate without cancellation
	out := Rotor{RimSpeed: v, Beta2: b2, RimRate: rate, LossPerDay: lag * 86400,
		SagnacPerRev: SagnacDelay(math.Pi*r*r, math.Abs(omega))}
	if m.Strength > 0 && m.Density > 0 {
		out.BurstMargin = v / BurstSpeed(m)
	}
	return out, nil
}
