// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
)

// Symmetrical components (Fortescue) in Kron's power-invariant form.
//
// L. V. Bewley, Tensor Analysis of Electric Circuits and Machines, Ronald
// Press 1961, Ch. 4 ("Bewley 1961"). The phase currents are the old
// variables and the sequence currents (0, 1, 2) the new ones, i = C·i′,
// with (Eq. 4‑17)
//
//	C = (1/√3)·[[1, 1, 1], [1, a², a], [1, a, a²]],   a = e^(j·120°).
//
// The 1/√3 makes C unitary (C*ₜ·C = I, Eqs. 4‑18 to 4‑20), so the complex
// power e*·i is invariant and impedances transform as Z′ = C*ₜ·Z·C
// (Eqs. 4‑10 to 4‑13, 4‑21 to 4‑23). Note the power-invariant scaling: the
// conventional (non-invariant) components are those returned here divided
// by √3, and a phase quantity's positive-sequence part is √3 times the
// conventional one.

// SeqOperator is a = e^(j·120°).
var SeqOperator = cmplx.Rect(1, 2*math.Pi/3)

// FortescueMatrix is Bewley's Eq. 4‑17: rows a, b, c; columns 0, 1, 2.
func FortescueMatrix() [3][3]complex128 {
	a := SeqOperator
	a2 := a * a
	s := complex(1/math.Sqrt(3), 0)
	return [3][3]complex128{
		{s, s, s},
		{s, s * a2, s * a},
		{s, s * a, s * a2},
	}
}

// ToSequence returns the power-invariant (0, 1, 2) components of the phase
// quantities x = (x_a, x_b, x_c): x′ = C*ₜ·x, the inverse of x = C·x′.
func ToSequence(x [3]complex128) [3]complex128 {
	c := FortescueMatrix()
	var y [3]complex128
	for i := 0; i < 3; i++ {
		for k := 0; k < 3; k++ {
			y[i] += cmplx.Conj(c[k][i]) * x[k]
		}
	}
	return y
}

// FromSequence returns the phase quantities x = C·x′.
func FromSequence(s [3]complex128) [3]complex128 {
	c := FortescueMatrix()
	var x [3]complex128
	for i := 0; i < 3; i++ {
		for k := 0; k < 3; k++ {
			x[i] += c[i][k] * s[k]
		}
	}
	return x
}

// SequenceImpedance transforms a phase impedance (or admittance) matrix to
// sequence axes, Z′ = C*ₜ·Z·C (Bewley 1961 Eqs. 4‑21 to 4‑23). A
// cyclic-symmetric Z (a machine, a transposed line) comes out diagonal
// (Eq. 4‑24); an unbalanced one couples the sequences (Eq. 4‑26).
func SequenceImpedance(z [3][3]complex128) [3][3]complex128 {
	c := FortescueMatrix()
	var zc, out [3][3]complex128
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				zc[i][j] += z[i][k] * c[k][j]
			}
		}
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				out[i][j] += cmplx.Conj(c[k][i]) * zc[k][j]
			}
		}
	}
	return out
}

// PhaseImpedance is the inverse law, Z = C·Z′·C*ₜ: the phase matrix of a
// network given its sequence impedance matrix (for example a generator's
// diag(Z₀, Z₁, Z₂)).
func PhaseImpedance(zs [3][3]complex128) [3][3]complex128 {
	c := FortescueMatrix()
	var t, out [3][3]complex128
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				t[i][j] += c[i][k] * zs[k][j]
			}
		}
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				out[i][j] += t[i][k] * cmplx.Conj(c[j][k])
			}
		}
	}
	return out
}
