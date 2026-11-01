// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "github.com/fdamador/kronspice/hyper/ga"

// Lorentz force and power by grade selection (CLAUDE.md §3f, M9n), after
// P. Joot, Geometric Algebra for Electrical Engineers (GAelectrodynamics
// v0.3.6.4), Theorem 3.16 "Lorentz force and power" — implemented from the
// equations; no text or code is taken from that work.
//
// In the algebra of physical space G3 the field is the multivector
// F = E + I·c·B (I = e₁₂₃) and a particle is a charge/velocity multivector
// Q. The product F·Q† (Q† the reverse) holds both the power and the force:
//
//	⟨F·Q†⟩₁ = dp/dt (force),   c·⟨F·Q†⟩₀ = dℰ/dt (power).
//
// Electric charge q at velocity v:   Q = q·(1 + v/c)
//
//	F·Q† = q(E + IcB)(1 + v/c): grade 1 = q(E + v×B), grade 0 = q E·v/c.
//
// Magnetic charge q_m [A·m, SI: force q_m·B on a charge at rest] at
// velocity v:   Q = I·(q_m/c)·(1 + v/c)
//
//	grade 1 = q_m(B − v×E/c²),   c·grade 0 = q_m B·v.
//
// Joot normalizes the magnetic charge so the force reads q_m(cB − v×E/c);
// his q_m is ours divided by c. Magnetic charge has never been observed
// (MoEDAL, IceCube searches); every use of it is hypothetical (duality.go).

// FieldMV is the field multivector F = E + I·c·B.
func FieldMV(e, b [3]float64) ga.G3 {
	c := SpeedOfLight
	return ga.G3{ga.G3E1: e[0], ga.G3E2: e[1], ga.G3E3: e[2],
		// I·(b₁e₁ + b₂e₂ + b₃e₃) = b₁e₂₃ + b₂e₃₁ + b₃e₁₂, e₃₁ = −e₁₃
		ga.G3E23: c * b[0], ga.G3E13: -c * b[1], ga.G3E12: c * b[2]}
}

// ElectricQ is the charge/velocity multivector q·(1 + v/c).
func ElectricQ(q float64, v [3]float64) ga.G3 {
	c := SpeedOfLight
	return ga.G3{ga.G3S: q, ga.G3E1: q * v[0] / c, ga.G3E2: q * v[1] / c, ga.G3E3: q * v[2] / c}
}

// MagneticQ is I·(q_m/c)·(1 + v/c) for a magnetic charge q_m [A·m].
func MagneticQ(qm float64, v [3]float64) ga.G3 {
	I := ga.G3{ga.G3I: 1}
	return I.Mul(ElectricQ(qm, v)).Scale(1 / SpeedOfLight)
}

// LorentzGA returns the force ⟨F·Q†⟩₁ [N] and power c·⟨F·Q†⟩₀ [W] on a
// particle Q in the field F.
func LorentzGA(f, q ga.G3) (force [3]float64, power float64) {
	m := f.Mul(q.Reverse())
	return [3]float64{m[ga.G3E1], m[ga.G3E2], m[ga.G3E3]}, SpeedOfLight * m[ga.G3S]
}

// LorentzVector is q(E + v×B) and q E·v.
func LorentzVector(q float64, v, e, b [3]float64) ([3]float64, float64) {
	vb := cross(v, b)
	return [3]float64{q * (e[0] + vb[0]), q * (e[1] + vb[1]), q * (e[2] + vb[2])}, q * dot(e, v)
}

// MagneticLorentzVector is q_m(B − v×E/c²) and q_m B·v (hypothetical).
func MagneticLorentzVector(qm float64, v, e, b [3]float64) ([3]float64, float64) {
	c2 := SpeedOfLight * SpeedOfLight
	ve := cross(v, e)
	return [3]float64{qm * (b[0] - ve[0]/c2), qm * (b[1] - ve[1]/c2), qm * (b[2] - ve[2]/c2)}, qm * dot(b, v)
}

// MotionalEMF integrates (v×B)·dl along a path of points with the
// velocity and field at each point (trapezoidal), the EMF a moving
// conductor develops by the Lorentz force on its charges — computed by
// grade selection on a unit test charge.
func MotionalEMF(path [][3]float64, vel, field func(p [3]float64) [3]float64) float64 {
	emf := 0.0
	per := func(p [3]float64) [3]float64 {
		f, _ := LorentzGA(FieldMV([3]float64{}, field(p)), ElectricQ(1, vel(p)))
		return f
	}
	for k := 1; k < len(path); k++ {
		a, b := path[k-1], path[k]
		fa, fb := per(a), per(b)
		dl := [3]float64{b[0] - a[0], b[1] - a[1], b[2] - a[2]}
		emf += 0.5 * (dot(fa, dl) + dot(fb, dl))
	}
	return emf
}
