// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"

	"github.com/fdamador/kronspice/circuit"
)

// Curvature returns the Riemann-Christoffel curvature tensor R_abcd of a
// machine's metric at angle θ, for the coordinates (q_1 … q_n, θ) of
// Christoffel: g = L(θ) on the winding axes and J on the θ axis
// (L. V. Bewley, Tensor Analysis of Electric Circuits and Machines, Ronald
// Press 1961, Art. 19‑6, Eq. 19‑44; Kron, "Non-Riemannian Dynamics of
// Rotating Electrical Machinery," J. Math. Phys. 13, 1934). With the
// symbols of the first kind [bc, a] (Christoffel),
//
//	R_abcd = ½(∂b∂c g_ad + ∂a∂d g_bc − ∂a∂c g_bd − ∂b∂d g_ac)
//	         + g^ef·([bc, e]·[ad, f] − [bd, e]·[ac, f])
//
// (index order out[a][b][c][d]; in two dimensions R_1212 = K·det g, K the
// Gaussian curvature). The metric depends on θ alone, so only θθ second
// derivatives appear. A machine whose inductances do not vary has R = 0;
// a salient-pole or slotted one does not, which is the curvature Bewley
// uses in the frame-independent small-oscillation equation (Art. 19‑6).
//
// Curvature is a tensor: unlike the torsion of a non-holonomic frame
// (Torsion), it does not vanish by a change of reference axes. It is the
// geometric object to set beside the curvature of Cartan geometry in the
// ECE comparison (doc/bewley.md).
func Curvature(s *circuit.MachineSpec, J, theta float64) ([][][][]float64, error) {
	if J <= 0 {
		return nil, errors.New("kron: Curvature: inertia must be positive")
	}
	n := s.N
	N := n + 1
	l, d2 := square(n), square(n)
	s.Inductance(theta, l, nil, d2)
	g := square(N)
	for i := 0; i < n; i++ {
		copy(g[i], l[i])
	}
	g[n][n] = J
	gi, err := invertDense(g)
	if err != nil {
		return nil, errors.New("kron: Curvature: singular metric (L(θ) not positive definite)")
	}
	ch := Christoffel(s, theta) // [b][c][a] = [bc, a]
	dd := func(a, b, c, d int) float64 { // ∂c∂d g_ab
		if c != n || d != n || a == n || b == n {
			return 0
		}
		return d2[a][b]
	}
	out := make([][][][]float64, N)
	for a := 0; a < N; a++ {
		out[a] = make([][][]float64, N)
		for b := 0; b < N; b++ {
			out[a][b] = make([][]float64, N)
			for c := 0; c < N; c++ {
				out[a][b][c] = make([]float64, N)
				for d := 0; d < N; d++ {
					v := 0.5 * (dd(a, d, b, c) + dd(b, c, a, d) - dd(b, d, a, c) - dd(a, c, b, d))
					for e := 0; e < N; e++ {
						for f := 0; f < N; f++ {
							if gi[e][f] == 0 {
								continue
							}
							v += gi[e][f] * (ch[b][c][e]*ch[a][d][f] - ch[b][d][e]*ch[a][c][f])
						}
					}
					out[a][b][c][d] = v
				}
			}
		}
	}
	return out, nil
}

// invertDense inverts a small dense matrix (Gauss–Jordan, partial pivoting).
func invertDense(a [][]float64) ([][]float64, error) {
	n := len(a)
	m := make([][]float64, n)
	for i := range a {
		m[i] = make([]float64, 2*n)
		copy(m[i], a[i])
		m[i][n+i] = 1
	}
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if abs(m[r][c]) > abs(m[p][c]) {
				p = r
			}
		}
		if m[p][c] == 0 {
			return nil, errors.New("singular")
		}
		m[c], m[p] = m[p], m[c]
		piv := m[c][c]
		for k := range m[c] {
			m[c][k] /= piv
		}
		for r := 0; r < n; r++ {
			if r == c || m[r][c] == 0 {
				continue
			}
			f := m[r][c]
			for k := range m[r] {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	out := make([][]float64, n)
	for i := range out {
		out[i] = m[i][n:]
	}
	return out, nil
}

func abs(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
