// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Ladder2: the planar ladder with a two-component field φ_i ∈ ℝ² (a and the
// angle of n in its plane), energy
//
//	E = Σ κ2|φ_{i+1} − φ_i|²/H + Σ W_i V(|φ_i|²),
//
// both ends fixed. Newton by block-tridiagonal elimination; stability by the
// inertia of the block LDLᵀ factorization (Haynsworth: the number of
// negative eigenvalues of the Hessian is the sum over the pivot blocks).

import "math"

type v2 = [2]float64
type m2 = [2][2]float64

// grad2 returns ∂E/∂φ_i.
func grad2(l Ladder, e EOS, p []v2) []v2 {
	n := len(p)
	g := make([]v2, n)
	for i := 1; i < n-1; i++ {
		u := p[i][0]*p[i][0] + p[i][1]*p[i][1]
		_, vu, _ := e.Vu(u)
		for c := 0; c < 2; c++ {
			g[i][c] = 2*l.W[i]*vu*p[i][c] +
				2*kappa2*(l.A[i-1]*(p[i][c]-p[i-1][c])-l.A[i]*(p[i+1][c]-p[i][c]))/l.H
		}
	}
	return g
}

// hess2 returns the diagonal blocks and the scalar off-diagonal couplings
// of the Hessian over the interior nodes.
func hess2(l Ladder, e EOS, p []v2) ([]m2, []float64) {
	n := len(p)
	d := make([]m2, n-2)
	o := make([]float64, n-3)
	for i := 1; i < n-1; i++ {
		u := p[i][0]*p[i][0] + p[i][1]*p[i][1]
		_, vu, vuu := e.Vu(u)
		s := 2 * kappa2 * (l.A[i-1] + l.A[i]) / l.H
		var b m2
		for r := 0; r < 2; r++ {
			for c := 0; c < 2; c++ {
				b[r][c] = l.W[i] * 4 * vuu * p[i][r] * p[i][c]
			}
			b[r][r] += s + 2*l.W[i]*vu
		}
		d[i-1] = b
		if i < n-2 {
			o[i-1] = -2 * kappa2 * l.A[i] / l.H
		}
	}
	return d, o
}

func inv2(a m2) m2 {
	det := a[0][0]*a[1][1] - a[0][1]*a[1][0]
	return m2{{a[1][1] / det, -a[0][1] / det}, {-a[1][0] / det, a[0][0] / det}}
}

// negBlocks counts the negative eigenvalues of a symmetric 2×2 block.
func negBlocks(a m2) int {
	tr, det := a[0][0]+a[1][1], a[0][0]*a[1][1]-a[0][1]*a[1][0]
	switch {
	case det < 0:
		return 1
	case tr < 0:
		return 2
	}
	return 0
}

// Solve2 runs damped Newton; returns iterations and the Hessian's number of
// negative eigenvalues at the solution.
func (l Ladder) Solve2(e EOS, p []v2, tol float64, maxit int) (int, int, error) {
	n := len(p)
	for it := 1; it <= maxit; it++ {
		g := grad2(l, e, p)
		var norm float64
		for i := 1; i < n-1; i++ {
			norm = math.Max(norm, math.Hypot(g[i][0], g[i][1])/l.W[i])
		}
		d, o := hess2(l, e, p)
		// Block LDLᵀ: D_k = d_k − o²·D_{k−1}⁻¹ (off-diagonal blocks are o·I).
		m := len(d)
		D := make([]m2, m)
		neg := 0
		for k := 0; k < m; k++ {
			D[k] = d[k]
			if k > 0 {
				iv := inv2(D[k-1])
				for r := 0; r < 2; r++ {
					for c := 0; c < 2; c++ {
						D[k][r][c] -= o[k-1] * o[k-1] * iv[r][c]
					}
				}
			}
			neg += negBlocks(D[k])
		}
		if norm < tol {
			return it, neg, nil
		}
		// Forward: y_k = −g_k − o·D_{k−1}⁻¹·y_{k−1}; back: x_k = D_k⁻¹(y_k − o·x_{k+1}).
		y := make([]v2, m)
		for k := 0; k < m; k++ {
			y[k] = v2{-g[k+1][0], -g[k+1][1]}
			if k > 0 {
				iv := inv2(D[k-1])
				for r := 0; r < 2; r++ {
					y[k][r] -= o[k-1] * (iv[r][0]*y[k-1][0] + iv[r][1]*y[k-1][1])
				}
			}
		}
		x := make([]v2, m)
		var big float64
		for k := m - 1; k >= 0; k-- {
			r := y[k]
			if k < m-1 {
				r[0] -= o[k] * x[k+1][0]
				r[1] -= o[k] * x[k+1][1]
			}
			iv := inv2(D[k])
			x[k] = v2{iv[0][0]*r[0] + iv[0][1]*r[1], iv[1][0]*r[0] + iv[1][1]*r[1]}
			big = math.Max(big, math.Hypot(x[k][0], x[k][1]))
		}
		s := 1.0
		if big > 0.2 {
			s = 0.2 / big
		}
		for k := 0; k < m; k++ {
			p[k+1][0] += s * x[k][0]
			p[k+1][1] += s * x[k][1]
		}
	}
	return maxit, -1, errNoConverge
}

// Energy2 is the discrete energy.
func (l Ladder) Energy2(e EOS, p []v2) float64 {
	var s float64
	for i := 0; i+1 < len(p); i++ {
		dx, dy := p[i+1][0]-p[i][0], p[i+1][1]-p[i][1]
		s += kappa2 * l.A[i] * (dx*dx + dy*dy) / l.H
	}
	for i := range p {
		s += l.W[i] * e.V(p[i][0]*p[i][0]+p[i][1]*p[i][1])
	}
	return s
}

// minAmp is min |φ| over the interior.
func minAmp(p []v2) float64 {
	m := math.Inf(1)
	for _, q := range p[1 : len(p)-1] {
		m = math.Min(m, math.Hypot(q[0], q[1]))
	}
	return m
}
