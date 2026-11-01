// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"fmt"
	"math"
)

// MaxDense bounds dense matrices; larger systems belong to the sparse backend.
const MaxDense = 4096

// Dense is a row-major dense matrix solved by LU factorization with partial
// (row) pivoting. All storage is allocated once in NewDense, so Clear, Add and
// Solve do not allocate (CLAUDE.md §8 rule 7).
//
// Algorithm: Doolittle elimination with partial pivoting, as in Golub & Van
// Loan, Matrix Computations, Alg. 3.4.1, followed by forward and back
// substitution.
type Dense struct {
	n   int
	a   []float64 // stamped matrix
	lu  []float64 // working copy, factored in place
	piv []int
	y   []float64
	col []float64 // largest |entry| of each column of the stamped matrix
}

// NewDense allocates an n×n system.
func NewDense(n int) (*Dense, error) {
	if n <= 0 || n > MaxDense {
		return nil, fmt.Errorf("%w: %d", ErrBadSize, n)
	}
	return &Dense{
		n:   n,
		a:   make([]float64, n*n),
		lu:  make([]float64, n*n),
		piv: make([]int, n),
		y:   make([]float64, n),
		col: make([]float64, n),
	}, nil
}

// Size returns the number of unknowns.
func (d *Dense) Size() int { return d.n }

// Clear zeroes the matrix.
func (d *Dense) Clear() {
	for i := range d.a {
		d.a[i] = 0
	}
}

// Add adds v to entry (i, j). Out-of-range indices (such as ground, −1) are ignored.
func (d *Dense) Add(i, j int, v float64) {
	if i < 0 || j < 0 || i >= d.n || j >= d.n {
		return
	}
	d.a[i*d.n+j] += v
}

// AbsProduct sets y_i = Σ_j |A_ij|·|x_j|.
func (d *Dense) AbsProduct(x, y []float64) {
	n := d.n
	for i := 0; i < n; i++ {
		s := 0.0
		row := d.a[i*n : (i+1)*n]
		for j, v := range row {
			s += math.Abs(v) * math.Abs(x[j])
		}
		y[i] = s
	}
}

// At returns entry (i, j) of the stamped (unfactored) matrix.
func (d *Dense) At(i, j int) float64 { return d.a[i*d.n+j] }

// Solve factors a copy of the matrix and solves A·x = b.
func (d *Dense) Solve(b, x []float64) error {
	n := d.n
	if len(b) != n || len(x) != n {
		return fmt.Errorf("%w: rhs %d, x %d, matrix %d", ErrBadSize, len(b), len(x), n)
	}
	copy(d.lu, d.a)
	scale := 0.0
	for _, v := range d.lu {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return ErrNonFinite
		}
		if a := math.Abs(v); a > scale {
			scale = a
		}
	}
	if scale == 0 {
		return ErrSingularMatrix
	}
	// A pivot smaller than n·1e−14 of its column's largest original entry
	// means the system is singular to working precision. The test is per
	// column, not against the whole matrix, so systems that mix units of
	// very different size — the energy domains of one circuit (volts next
	// to pascals and acoustic impedances of 10⁹) — are not mistaken for
	// singular: scaling a column (changing an unknown's unit) leaves it
	// unchanged.
	for j := range d.col {
		d.col[j] = 0
	}
	for i := 0; i < n; i++ {
		for j, v := range d.lu[i*n : (i+1)*n] {
			if a := math.Abs(v); a > d.col[j] {
				d.col[j] = a
			}
		}
	}

	lu := d.lu
	for k := 0; k < n; k++ {
		// Partial pivoting: pick the largest |entry| in column k at or below row k.
		p := k
		maxv := math.Abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(lu[i*n+k]); v > maxv {
				maxv, p = v, i
			}
		}
		if maxv <= d.col[k]*float64(n)*1e-14 {
			return fmt.Errorf("%w: pivot %d is %.3g", ErrSingularMatrix, k, maxv)
		}
		d.piv[k] = p
		if p != k {
			rk, rp := lu[k*n:(k+1)*n], lu[p*n:(p+1)*n]
			for j := range rk {
				rk[j], rp[j] = rp[j], rk[j]
			}
		}
		pivot := lu[k*n+k]
		for i := k + 1; i < n; i++ {
			f := lu[i*n+k] / pivot
			lu[i*n+k] = f
			if f == 0 {
				continue
			}
			ri, rk := lu[i*n:(i+1)*n], lu[k*n:(k+1)*n]
			for j := k + 1; j < n; j++ {
				ri[j] -= f * rk[j]
			}
		}
	}
	// Apply the row permutation to b, then forward substitution (unit lower).
	y := d.y
	copy(y, b)
	for k := 0; k < n; k++ {
		if p := d.piv[k]; p != k {
			y[k], y[p] = y[p], y[k]
		}
	}
	for i := 0; i < n; i++ {
		s := y[i]
		for j := 0; j < i; j++ {
			s -= lu[i*n+j] * y[j]
		}
		y[i] = s
	}
	// Back substitution (upper).
	for i := n - 1; i >= 0; i-- {
		s := y[i]
		for j := i + 1; j < n; j++ {
			s -= lu[i*n+j] * x[j]
		}
		x[i] = s / lu[i*n+i]
	}
	for _, v := range x {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return ErrNonFinite
		}
	}
	return nil
}
