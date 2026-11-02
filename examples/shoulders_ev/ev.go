// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// K. R. Shoulders, EV: A Tale of Discovery (Jupiter Technologies, Austin,
// 1987). Numbers used (book chapter-page):
//
//   - preferred EV size about 1 µm diameter, spherical (Preface);
//     witness-plate strikes 1–10 µm (2-7);
//   - "one EV bead carries around 10¹¹ charges" (8-49);
//   - a 3 µm EV leaves 2×10¹⁰ electrons on an electrode; a 10 µm EV gives
//     off 3.5×10¹⁴ electrons over 1 mm; density ≈ 6.6×10²³ electrons/cm³,
//     "approximately that of a solid"; translational speed c/10 (4-9, 4-63);
//   - at most one ion per 100,000 electrons (4-63, property 5);
//   - EVs made from a ~1 pF cathode storing "below 10⁻⁶ joules" (2-7);
//   - the container: zero-point fluctuations (Preface; 9-2).

// SI constants (CODATA 2018).
const (
	C0   = 299792458.0
	Eps0 = 8.8541878128e-12
	Qe   = 1.602176634e-19
	Me   = 9.1093837015e-31
	Hbar = 1.054571817e-34
	Ke   = 1 / (4 * math.Pi * Eps0)
)

// EV is a uniform sphere of N electrons of diameter D moving at speed V.
type EV struct {
	N, D, V float64
}

// Book is Shoulders' standard bead: 10¹¹ electrons, 1 µm, c/10.
func Book() EV { return EV{N: 1e11, D: 1e-6, V: C0 / 10} }

func (e EV) R() float64       { return e.D / 2 }
func (e EV) Charge() float64  { return e.N * Qe }
func (e EV) Volume() float64  { return 4.0 / 3 * math.Pi * math.Pow(e.R(), 3) }
func (e EV) Density() float64 { return e.N / e.Volume() } // electrons/m³

// SelfEnergy is the electrostatic energy of the uniform sphere,
// (3/5)·Q²/(4πε0·R) [J] — the work needed to assemble it.
func (e EV) SelfEnergy() float64 { return 0.6 * Ke * e.Charge() * e.Charge() / e.R() }

// SurfaceField is Q/(4πε0R²) [V/m].
func (e EV) SurfaceField() float64 { return Ke * e.Charge() / (e.R() * e.R()) }

// Pressure is the outward electrostatic stress at the surface, ε0E²/2 [Pa].
func (e EV) Pressure() float64 { E := e.SurfaceField(); return 0.5 * Eps0 * E * E }

// Kinetic is the electrons' total kinetic energy at speed V [J].
func (e EV) Kinetic() float64 {
	g := 1 / math.Sqrt(1-e.V*e.V/(C0*C0))
	return e.N * (g - 1) * Me * C0 * C0
}

// PlasmaFreq is the electron plasma frequency of the cluster [rad/s].
func (e EV) PlasmaFreq() float64 { return math.Sqrt(e.Density() * Qe * Qe / (Me * Eps0)) }

// MagneticFraction is the ratio of magnetic attraction to electric
// repulsion between electrons moving together at speed V: β².
func (e EV) MagneticFraction() float64 { return e.V * e.V / (C0 * C0) }

// Explode integrates the Coulomb explosion of the bare cluster: the outer
// electron of a uniformly expanding sphere always has all the charge
// inside it, so dp/dt = k·Q·e/r² (relativistic momentum), from rest.
// Returns the times for the radius to grow ×2 and ×10 and the outer
// electron's kinetic energy at ×10 [eV].
func (e EV) Explode() (t2, t10, keV float64) {
	f := Ke * e.Charge() * Qe
	mc := Me * C0
	r, p, t := e.R(), 0.0, 0.0
	dt := 1e-3 / e.PlasmaFreq()
	vel := func(p float64) float64 { return p * C0 / math.Hypot(p, mc) }
	for r < 10*e.R() {
		// midpoint (RK2) step
		pm := p + 0.5*dt*f/(r*r)
		rm := r + 0.5*dt*vel(p)
		p += dt * f / (rm * rm)
		r += dt * vel(pm)
		t += dt
		if t2 == 0 && r >= 2*e.R() {
			t2 = t
		}
	}
	gamma := math.Hypot(p, mc) / mc
	return t2, t, (gamma - 1) * Me * C0 * C0 / Qe
}

// Lifetime is Shoulders' EV lifetime estimate (4-9): 1 mm of track at c/10.
func Lifetime() float64 { return 1e-3 / (C0 / 10) }

// EscapeEnergy is the energy an electron gains leaving the surface of the
// bare cluster for far away, k·Q·e/R [eV].
func (e EV) EscapeEnergy() float64 { return Ke * e.Charge() / e.R() }

// GaussRadius is the smallest radius a net charge q can sit inside if
// electrons leaving it gain only v volts: by Gauss's law the potential
// outside any closed surface around q is k·q/r whatever is inside, so
// r ≥ k·q/v.
func GaussRadius(q, v float64) float64 { return Ke * q / v }

// BrillouinField is the magnetic field needed to confine electrons at
// density n against their own repulsion, the Brillouin limit
// n = ε0B²/(2m) [T].
func BrillouinField(n float64) float64 { return math.Sqrt(2 * Me * n / Eps0) }

// CasimirPlates is the ideal Casimir pressure between plates at gap a [Pa]
// (attractive), π²ħc/(240a⁴).
func CasimirPlates(a float64) float64 {
	return math.Pi * math.Pi * Hbar * C0 / (240 * math.Pow(a, 4))
}

// BoyerSphere is the zero-point stress on a perfectly conducting
// spherical shell of radius a: energy +0.04618·ħc/a (T. H. Boyer, Phys. Rev.
// 174, 1764, 1968), so it pushes outward with pressure 0.04618·ħc/(4πa⁴).
func BoyerSphere(a float64) float64 {
	return 0.04618 * Hbar * C0 / (4 * math.Pi * math.Pow(a, 4))
}

// SourceCharge is the most charge a cathode capacitance c can give when its
// stored energy is at most w: Q = √(2wc).
func SourceCharge(c, w float64) float64 { return math.Sqrt(2 * w * c) }
