// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"

	"github.com/fdamador/kronspice/hyper/quat"
)

// Quaternion polyphase networks (CLAUDE.md §3b item 2, M9f).
//
// A set of three-phase quantities (x_a, x_b, x_c) is the pure quaternion
// x = x_a·i + x_b·j + x_c·k. The power-invariant Clarke transformation
// (E. Clarke, Circuit Analysis of A-C Power Systems, Wiley 1943, Vol. I,
// Ch. X) is an orthogonal matrix of determinant +1, so it is a rotation of
// that 3-space: the axes α, β, 0 are the images of i, j, k, with the 0
// (zero-sequence) axis along (1, 1, 1)/√3. Park's transformation (R. H.
// Park, "Two-reaction theory of synchronous machines," AIEE Trans. 48,
// 1929) is a further rotation by −θ about the 0 axis, so abc → dq0 is the
// single unit quaternion r(θ) = r_P(θ)·r_C acting by r·x·r* — Kron's
// rotating reference frame of the machine (Kron 1951 Ch. 4; RotatingFrame)
// as a rotor. Power-invariant scaling is used throughout, so the frame
// changes are orthogonal (Kron's power-invariant transformations).
//
// Instantaneous power in quaternion form: s = v·i* for pure quaternions v
// and i. Because i* = −i and uv = −u·v + u×v (hyper/quat),
//
//	s = v·i − v×i:   S.s = v_a·i_a + v_b·i_b + v_c·i_c = p (active power),
//	                 V.s = −v×i, whose length is Akagi's instantaneous
//	                 imaginary power (H. Akagi, Y. Kanazawa, A. Nabae,
//	                 IEEE Trans. Ind. Appl. IA-20(3), 1984; the vector
//	                 form of F. Z. Peng, J.-S. Lai, IEEE Trans. Instrum.
//	                 Meas. 45(1), 1996).
//
// Under any frame rotation r, s transforms as r·s·r*: p is invariant and
// the non-active vector turns with the frame (the quaternion three-phase
// power of A. S. F. Komeno et al., RE&PQJ).

// ClarkeMatrix is the power-invariant abc → αβ0 matrix
// √(2/3)·[[1, −½, −½], [0, √3/2, −√3/2], [1/√2, 1/√2, 1/√2]].
func ClarkeMatrix() [3][3]float64 {
	s := math.Sqrt(2.0 / 3.0)
	h := math.Sqrt(3) / 2
	r := 1 / math.Sqrt2
	return [3][3]float64{{s, -s / 2, -s / 2}, {0, s * h, -s * h}, {s * r, s * r, s * r}}
}

// ParkMatrix is the power-invariant abc → dq0 matrix at electrical angle θ
// (d axis on phase a at θ = 0): Park's rotation [[cos θ, sin θ, 0],
// [−sin θ, cos θ, 0], [0, 0, 1]] after Clarke.
func ParkMatrix(theta float64) [3][3]float64 {
	s, c := math.Sincos(theta)
	rot := [3][3]float64{{c, s, 0}, {-s, c, 0}, {0, 0, 1}}
	return mul3(rot, ClarkeMatrix())
}

// ClarkeRotor is the unit quaternion of the Clarke transformation:
// x_αβ0 = r·x_abc·r*, components (α, β, 0) on (i, j, k).
func ClarkeRotor() quat.Q { return quat.FromMatrix(ClarkeMatrix()) }

// ParkRotor is the unit quaternion of the abc → dq0 transformation at
// electrical angle θ: the Park rotation about the 0 axis, cos(θ/2) −
// k·sin(θ/2), composed with the Clarke rotor.
func ParkRotor(theta float64) quat.Q {
	s, c := math.Sincos(theta / 2)
	return quat.Q{W: c, Z: -s}.Mul(ClarkeRotor())
}

// Phases is the pure quaternion a·i + b·j + c·k of three phase values.
func Phases(a, b, c float64) quat.Q { return quat.Vec(a, b, c) }

// Rotate3 returns r·x·r* for a three-phase quaternion x (a frame change).
func Rotate3(r, x quat.Q) quat.Q { return r.Mul(x).Mul(r.Conj()) }

// InstantPower is the quaternion power s = v·i* of three-phase voltage and
// current quaternions: S.s is the instantaneous active power and V.s the
// instantaneous non-active power vector (−v×i).
func InstantPower(v, i quat.Q) quat.Q { return v.Mul(i.Conj()) }

// ParkFrame returns Kron's moving frame for three-phase winding triples
// (indices into a machine's windings) at electrical angle p·θ: winding
// currents i_abc = Pᵀ(pθ)·i_dq0, the identity on the other windings. It is
// the three-phase counterpart of RotorFrame, built from the Park rotor, and
// with RotatingFrame turns a salient machine's θ-dependent abc
// inductances into constant d, q, 0 inductances.
func ParkFrame(n int, triples [][3]int, p float64) MovingFrame {
	build := func(th float64, deriv bool) [][]float64 {
		c := square(n)
		in := map[int]bool{}
		for _, t := range triples {
			for _, w := range t {
				in[w] = true
			}
		}
		for k := 0; k < n; k++ {
			if !in[k] && !deriv {
				c[k][k] = 1
			}
		}
		var m [3][3]float64
		if deriv {
			// d/dθ of Pᵀ(pθ): the rotor's matrix differentiated
			// analytically (p·dRᵀ/d(pθ)·K).
			s, co := math.Sincos(p * th)
			drot := [3][3]float64{{-s, co, 0}, {-co, -s, 0}, {0, 0, 0}}
			m = mul3(drot, ClarkeMatrix())
			for i := range m {
				for j := range m[i] {
					m[i][j] *= p
				}
			}
		} else {
			m = ParkRotor(p * th).Matrix()
		}
		for _, t := range triples {
			for r := 0; r < 3; r++ {
				for col := 0; col < 3; col++ {
					// C = Pᵀ: row = abc winding, column = dq0 axis.
					c[t[r]][t[col]] = m[col][r]
				}
			}
		}
		return c
	}
	return MovingFrame{C: func(t float64) [][]float64 { return build(t, false) }, D: func(t float64) [][]float64 { return build(t, true) }}
}

func mul3(a, b [3][3]float64) [3][3]float64 {
	var c [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for k := 0; k < 3; k++ {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}
