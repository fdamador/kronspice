// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// CODATA 2018 values (SI).
const (
	Hbar = 1.054571817e-34   // J·s
	E    = 1.602176634e-19   // C
	Me   = 9.1093837015e-31  // kg
	Eps0 = 8.8541878128e-12  // F/m
	C    = 299792458.0       // m/s
	A0   = 5.29177210903e-11 // Bohr radius [m]
	Ry   = 13.605693122994   // Rydberg energy [eV]
)

// Puthoff's balance (H. E. Puthoff, Phys. Rev. D 35, 3266, 1987).
//
// An electron on a circular orbit of radius r at angular frequency ω is
// two one-dimensional oscillators in quadrature (eq. 17). From the random
// zero-point field of spectral energy density ħω³/(2π²c³) (eq. 1) it
// absorbs on average (eq. 18)
//
//	P_abs = e²ħω³ / (6πε0·m·c³)
//
// and, being accelerated, it radiates (Larmor, eq. 19)
//
//	P_rad = e²r²ω⁴ / (6πε0·c³)
//
// Equal when m·ω·r² = ħ (eq. 20) — Bohr's quantization of angular momentum.

// Pabs is eq. 18 [W].
func Pabs(omega float64) float64 {
	return E * E * Hbar * omega * omega * omega / (6 * math.Pi * Eps0 * Me * C * C * C)
}

// Prad is eq. 19 [W].
func Prad(r, omega float64) float64 {
	return E * E * r * r * math.Pow(omega, 4) / (6 * math.Pi * Eps0 * C * C * C)
}

// Kepler is the Coulomb orbit's angular frequency at radius r:
// m·ω²·r = e²/(4πε0·r²).
func Kepler(r float64) float64 {
	return math.Sqrt(E * E / (4 * math.Pi * Eps0 * Me * r * r * r))
}

// Orbit is the radius where the balance (eq. 20) and the Coulomb force
// hold together: r = 4πε0ħ²/(m e²), the Bohr radius.
func Orbit() (r, omega, energyEV float64) {
	r = 4 * math.Pi * Eps0 * Hbar * Hbar / (Me * E * E)
	omega = Kepler(r)
	energyEV = -E / (8 * math.Pi * Eps0 * r)
	return r, omega, energyEV
}

// Net is P_abs − P_rad on the Coulomb orbit of radius r [W]: positive
// means the orbit gains energy (grows).
func Net(r float64) float64 {
	w := Kepler(r)
	return Pabs(w) - Prad(r, w)
}

// Relax is the time constant [s] of the averaged balance: the orbit energy
// is −e²/(8πε0 r), so dE/dr = e²/(8πε0 r²) and dr/dt = Net(r)/(dE/dr);
// linearized about a0, dr/dt = −(r − a0)/τ.
func Relax() float64 {
	r0, _, _ := Orbit()
	h := 1e-6 * r0
	dEdr := E * E / (8 * math.Pi * Eps0 * r0 * r0)
	slope := (Net(r0+h) - Net(r0-h)) / (2 * h) / dEdr
	return -1 / slope
}

// Collapse is the classical fall time of the orbit without the zero-point
// field (Larmor only, from a0): t = a0³/(4·r_e²·c), r_e the classical
// electron radius.
func Collapse() float64 {
	re := E * E / (4 * math.Pi * Eps0 * Me * C * C)
	return A0 * A0 * A0 / (4 * re * re * C)
}
