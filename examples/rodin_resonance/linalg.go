// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
	"sort"
)

// Small dense linear algebra (the project keeps to the standard library).

// cholesky returns the lower factor G of a symmetric positive-definite
// matrix, a = G·Gᵀ.
func cholesky(a [][]float64) ([][]float64, error) {
	n := len(a)
	g := newMat(n)
	for j := 0; j < n; j++ {
		s := a[j][j]
		for k := 0; k < j; k++ {
			s -= g[j][k] * g[j][k]
		}
		if s <= 0 {
			return nil, errors.New("matrix not positive definite")
		}
		g[j][j] = math.Sqrt(s)
		for i := j + 1; i < n; i++ {
			t := a[i][j]
			for k := 0; k < j; k++ {
				t -= g[i][k] * g[j][k]
			}
			g[i][j] = t / g[j][j]
		}
	}
	return g, nil
}

// lowerInverse inverts a lower-triangular matrix.
func lowerInverse(g [][]float64) [][]float64 {
	n := len(g)
	x := newMat(n)
	for c := 0; c < n; c++ {
		for i := c; i < n; i++ {
			s := 0.0
			if i == c {
				s = 1
			}
			for k := c; k < i; k++ {
				s -= g[i][k] * x[k][c]
			}
			x[i][c] = s / g[i][i]
		}
	}
	return x
}

// spdInverse inverts a symmetric positive-definite matrix.
func spdInverse(a [][]float64) ([][]float64, error) {
	g, err := cholesky(a)
	if err != nil {
		return nil, err
	}
	gi := lowerInverse(g)
	n := len(a)
	out := newMat(n)
	parallelFor(n, func(i int) {
		for j := 0; j < n; j++ {
			s := 0.0
			for k := max(i, j); k < n; k++ {
				s += gi[k][i] * gi[k][j]
			}
			out[i][j] = s
		}
	})
	return out, nil
}

func matMul(a, b [][]float64) [][]float64 {
	n, m, p := len(a), len(b), len(b[0])
	out := make([][]float64, n)
	parallelFor(n, func(i int) {
		row := make([]float64, p)
		for k := 0; k < m; k++ {
			if a[i][k] == 0 {
				continue
			}
			for j := 0; j < p; j++ {
				row[j] += a[i][k] * b[k][j]
			}
		}
		out[i] = row
	})
	return out
}

func transpose(a [][]float64) [][]float64 {
	out := make([][]float64, len(a[0]))
	for j := range out {
		out[j] = make([]float64, len(a))
		for i := range a {
			out[j][i] = a[i][j]
		}
	}
	return out
}

// symEigen returns the eigenvalues (ascending) and eigenvectors (columns)
// of a symmetric matrix by cyclic Jacobi rotations (G. H. Golub,
// C. F. Van Loan, Matrix Computations, 3rd ed., §8.4).
func symEigen(h [][]float64) ([]float64, [][]float64) {
	n := len(h)
	a := newMat(n)
	v := newMat(n)
	scale := 0.0
	for i := range a {
		copy(a[i], h[i])
		v[i][i] = 1
		for j := range a[i] {
			scale += a[i][j] * a[i][j]
		}
	}
	for sweep := 0; sweep < 60; sweep++ {
		off := 0.0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off <= 1e-28*scale {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(a[p][q]) < 1e-300 {
					continue
				}
				th := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := math.Copysign(1, th) / (math.Abs(th) + math.Sqrt(th*th+1))
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p], a[k][q] = c*akp-s*akq, s*akp+c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k], a[q][k] = c*apk-s*aqk, s*apk+c*aqk
				}
				for k := 0; k < n; k++ {
					vkp, vkq := v[k][p], v[k][q]
					v[k][p], v[k][q] = c*vkp-s*vkq, s*vkp+c*vkq
				}
			}
		}
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(x, y int) bool { return a[idx[x]][idx[x]] < a[idx[y]][idx[y]] })
	vals := make([]float64, n)
	vecs := newMat(n)
	for c, i := range idx {
		vals[c] = a[i][i]
		for k := 0; k < n; k++ {
			vecs[k][c] = v[k][i]
		}
	}
	return vals, vecs
}

// genEigen solves K·x = λ·C·x for symmetric K and SPD C: with C = G·Gᵀ,
// G⁻¹·K·G⁻ᵀ·y = λ·y and x = G⁻ᵀ·y.
func genEigen(k, c [][]float64) ([]float64, [][]float64, error) {
	g, err := cholesky(c)
	if err != nil {
		return nil, nil, err
	}
	gi := lowerInverse(g)
	a := matMul(matMul(gi, k), transpose(gi))
	for i := range a { // symmetrize round-off
		for j := 0; j < i; j++ {
			m := 0.5 * (a[i][j] + a[j][i])
			a[i][j], a[j][i] = m, m
		}
	}
	vals, y := symEigen(a)
	return vals, matMul(transpose(gi), y), nil
}
