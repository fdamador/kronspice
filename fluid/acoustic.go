// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package fluid

import (
	"errors"
	"fmt"
	"math"
)

// ErrOutsideModelValidity is returned for a model used outside its range.
var ErrOutsideModelValidity = errors.New("fluid: outside model validity")

// PistonRadiation is the acoustic radiation load of a rigid circular piston
// of radius a [m] in an infinite baffle, radiating into medium m, as a
// lumped acoustic network for the acoustic domain (across = pressure,
// through = volume velocity):
//
//	inertance Ma [kg/m⁴] in parallel with a resistance Rpar [Pa·s/m³],
//
// Ma = 8ρ/(3π²a) and Rpar = 128ρc/(9π³a²). For ka ≪ 1 its impedance,
// jωMa·Rpar/(Rpar + jωMa) ≈ jωMa + (ωMa)²/Rpar, equals the low-frequency
// radiation impedance of the piston, Z = (ρc/S)·[(ka)²/2 + j·8ka/(3π)]
// with S = πa² (L. E. Kinsler, A. R. Frey, A. B. Coppens, J. V. Sanders,
// Fundamentals of Acoustics, 4th ed., Wiley 2000, §7.5, eqs 7.5.11–7.5.13),
// so the power through Rpar is the sound radiated into the medium. The
// error of both parts is of order (ka)²; fmax is the frequency at which
// ka = 0.5, the end of the stated range.
func PistonRadiation(a float64, m Medium) (Ma, Rpar, fmax float64, err error) {
	if !(a > 0) || !(m.Density > 0) || !(m.SoundSpeed > 0) {
		return 0, 0, 0, fmt.Errorf("%w: radius %g, medium %+v", ErrOutsideModelValidity, a, m)
	}
	rho, c := m.Density, m.SoundSpeed
	Ma = 8 * rho / (3 * math.Pi * math.Pi * a)
	Rpar = 128 * rho * c / (9 * math.Pi * math.Pi * math.Pi * a * a)
	fmax = 0.5 * c / (2 * math.Pi * a)
	return Ma, Rpar, fmax, nil
}

// PistonImpedance is the exact radiation impedance of the baffled piston
// in acoustic units, Z = (ρc/S)·[1 − 2J₁(2ka)/(2ka)] + j·(ρc/S)·X₁(2ka),
// X₁(x) = 2·H₁(x)/x and H₁ the Struve function, from its power series
// H₁(x) = Σ (−1)^k (x/2)^(2k+2) / (Γ(k+3/2)·Γ(k+5/2)) (Kinsler et al.
// §7.5, eq 7.5.12; M. Abramowitz, I. A. Stegun, Handbook of Mathematical
// Functions, 1964, eq 12.1.3), summed to convergence (accurate for
// 2ka ≲ 20).
func PistonImpedance(a float64, m Medium, f float64) complex128 {
	S := math.Pi * a * a
	k := 2 * math.Pi * f / m.SoundSpeed
	x := 2 * k * a
	r1 := 1 - 2*math.J1(x)/x
	h1, term := 0.0, 0.0
	for kk := 0; kk < 200; kk++ {
		lg1, _ := math.Lgamma(float64(kk) + 1.5)
		lg2, _ := math.Lgamma(float64(kk) + 2.5)
		term = math.Exp(float64(2*kk+2)*math.Log(x/2) - lg1 - lg2)
		if kk%2 == 1 {
			term = -term
		}
		h1 += term
		if math.Abs(term) < 1e-17*math.Abs(h1) {
			break
		}
	}
	x1 := 2 * h1 / x
	z := m.Impedance() / S
	return complex(z*r1, z*x1)
}

// SphericalResonator is the first radial (breathing) mode of a liquid
// sphere of radius R [m] with a pressure-release wall, the standing wave of
// single-bubble sonoluminescence flasks: p(r) = P·sin(kr)/(kr) with
// kR = π, so f0 = c/(2R), and the pressure antinode — where the bubble is
// trapped — at the center (M. P. Brenner, S. Hilgenfeldt, D. Lohse,
// "Single-bubble sonoluminescence," Rev. Mod. Phys. 74, 425 (2002), §II.A;
// Kinsler et al. §9.3). As a lumped acoustic node whose across variable is
// the center pressure amplitude P, the mode is a parallel resonator:
//
//	Ca = 2R³/(π·ρ·c²)   (from the peak potential energy ∫ p²/(2ρc²) dV = ½·Ca·P²)
//	Ma = 1/(ω0²·Ca),     Rpar = Q/(ω0·Ca)
//
// with Q the mode's quality factor (wall, viscous and thermal losses; a
// measured value of the actual flask should be used). How strongly a
// transducer on the wall drives this node depends on the mode shape where
// it sits; that coupling is the caller's (an effective piston area).
func SphericalResonator(R float64, m Medium, Q float64) (f0, Ca, Ma, Rpar float64, err error) {
	if !(R > 0) || !(Q > 0) || !(m.Density > 0) || !(m.SoundSpeed > 0) {
		return 0, 0, 0, 0, fmt.Errorf("%w: radius %g, Q %g", ErrOutsideModelValidity, R, Q)
	}
	rho, c := m.Density, m.SoundSpeed
	f0 = c / (2 * R)
	w0 := 2 * math.Pi * f0
	Ca = 2 * R * R * R / (math.Pi * rho * c * c)
	Ma = 1 / (w0 * w0 * Ca)
	Rpar = Q / (w0 * Ca)
	return f0, Ca, Ma, Rpar, nil
}
