// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
)

// Multiwinding transformers from leakage (short-circuit) impedances.
//
// Bewley 1961 Ch. 5, Art. 5‑1. Winding 0 is the reference winding of N₀
// turns; winding j has the turns ratio n_j = N_j/N₀. Ignoring the
// magnetizing current, N₀·i⁰ + Σ N_k·i^k = 0, i.e. i⁰ = −n_k·i^k
// (Eq. 5‑5), so the (n+1) windings reduce to n axes with the equation
//
//	e_j − n_j·e₀ = Z_jk·i^k                                 (Eq. 5‑9)
//	Z_jk = (n_j·n_k/2)·[Z₀₋ⱼ + Z₀₋ₖ − Z_j₋ₖ]
//
// where Z_a₋b is the leakage impedance between windings a and b referred
// to the reference winding (Eq. 5‑4):
//
//	Z_a₋b = z_aa/n_a² + z_bb/n_b² − 2·z_ab/(n_a·n_b)
//
// The point of Eq. 5‑9 is that it needs only the pairwise short-circuit
// impedances a transformer test report gives, not the self and mutual
// impedances (which differ by the large magnetizing reactance).

// PairLeakage returns Z_a₋b of every pair of windings from the self and
// mutual impedances z (Eq. 5‑4). n are the turns ratios to the reference
// winding (n[0] = 1 for the reference itself).
func PairLeakage(z [][]complex128, n []float64) ([][]complex128, error) {
	m := len(n)
	if len(z) != m {
		return nil, fmt.Errorf("kron: PairLeakage: %d windings, %d rows of z", m, len(z))
	}
	out := make([][]complex128, m)
	for a := range out {
		if len(z[a]) != m {
			return nil, fmt.Errorf("kron: PairLeakage: row %d has %d entries", a, len(z[a]))
		}
		if n[a] == 0 {
			return nil, errors.New("kron: PairLeakage: zero turns ratio")
		}
		out[a] = make([]complex128, m)
	}
	for a := 0; a < m; a++ {
		for b := 0; b < m; b++ {
			if a == b {
				continue
			}
			na, nb := complex(n[a], 0), complex(n[b], 0)
			out[a][b] = z[a][a]/(na*na) + z[b][b]/(nb*nb) - 2*z[a][b]/(na*nb)
		}
	}
	return out, nil
}

// LeakageMatrix returns Bewley's Eq. 5‑9 matrix Z_jk for windings 1..n
// from the pairwise leakage impedances pair[a][b] = Z_a₋b (indices 0..n,
// 0 the reference winding) and the turns ratios n[0..n] (n[0] is unused).
// The result is the impedance of the transformer seen from the voltages
// e_j − n_j·e₀, with the magnetizing current ignored.
func LeakageMatrix(pair [][]complex128, n []float64) ([][]complex128, error) {
	m := len(n)
	if m < 2 || len(pair) != m {
		return nil, fmt.Errorf("kron: LeakageMatrix: need ≥ 2 windings and an (n+1)×(n+1) table, have %d and %d", m, len(pair))
	}
	for a := range pair {
		if len(pair[a]) != m {
			return nil, fmt.Errorf("kron: LeakageMatrix: row %d has %d entries", a, len(pair[a]))
		}
	}
	z := make([][]complex128, m-1)
	for j := 1; j < m; j++ {
		z[j-1] = make([]complex128, m-1)
		for k := 1; k < m; k++ {
			nn := complex(n[j]*n[k]/2, 0)
			zjk := complex(0, 0)
			if j != k {
				zjk = pair[j][k]
			}
			z[j-1][k-1] = nn * (pair[0][j] + pair[0][k] - zjk)
		}
	}
	return z, nil
}
