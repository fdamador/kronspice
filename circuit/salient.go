// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"errors"
	"math"
)

// NewSalientPoleSpec returns the three armature windings (a, b, c) of a
// salient-pole machine whose direct-, quadrature- and zero-axis
// inductances are ld, lq and l0, written in the phase (holonomic) frame
// with θ the electrical rotor angle (d axis on phase a at θ = 0):
//
//	L(θ) = P(θ)ᵀ·diag(ld, lq, l0)·P(θ)
//
// with P the power-invariant Park matrix (kron.ParkMatrix). Expanded, the
// self-inductances are L₁ + L₂·cos 2θ and the a–b mutual −M₁ + L₂·cos 2(θ −
// 60°) (Bewley's θ + 60° counts phase b the other way round), L₁ = (ld + lq + l0)/3, L₂ = (ld − lq)/3, M₁ = (ld + lq
// − 2l0)/6 (L. V. Bewley, Tensor Analysis of Electric Circuits and
// Machines, Ronald Press 1961, Art. 20‑10, Eqs. 20‑46, 20‑47; Kron 1951
// Ch. 4). So the spec has a constant part and a second harmonic only; Park's
// transformation turns it back into the constant diag(ld, lq, l0).
//
// Field and amortisseur windings, if wanted, are added by the caller
// (indices 3 and up after growing the spec).
func NewSalientPoleSpec(ld, lq, l0 float64) (*MachineSpec, error) {
	if !(ld > 0 && lq > 0 && l0 > 0) {
		return nil, errors.New("circuit: NewSalientPoleSpec: ld, lq and l0 must be positive")
	}
	// Power-invariant Clarke matrix (rows α, β, 0).
	s := math.Sqrt(2.0 / 3.0)
	h := math.Sqrt(3) / 2
	r := 1 / math.Sqrt2
	cl := [3][3]float64{{s, -s / 2, -s / 2}, {0, s * h, -s * h}, {s * r, s * r, s * r}}
	// Rᵀ(θ)·diag(ld, lq, l0)·R(θ) = mean + cos 2θ·A + sin 2θ·B on the αβ0 axes.
	m, d := (ld+lq)/2, (ld-lq)/2
	mean := [3][3]float64{{m, 0, 0}, {0, m, 0}, {0, 0, l0}}
	a := [3][3]float64{{d, 0, 0}, {0, -d, 0}, {0, 0, 0}}
	b := [3][3]float64{{0, d, 0}, {d, 0, 0}, {0, 0, 0}}
	back := func(x [3][3]float64) [][]float64 { // Clarkeᵀ·x·Clarke
		out := squareZero(3)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				var v float64
				for p := 0; p < 3; p++ {
					for q := 0; q < 3; q++ {
						v += cl[p][i] * x[p][q] * cl[q][j]
					}
				}
				out[i][j] = v
			}
		}
		return out
	}
	spec := NewMachineSpec(3)
	spec.L0 = back(mean)
	h2 := spec.Harmonic(2)
	h2.Cos = back(a)
	h2.Sin = back(b)
	return spec, nil
}
