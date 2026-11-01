// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"fmt"
	"math"
	"math/cmplx"
)

// CStamper receives complex matrix contributions (AC analysis).
type CStamper interface {
	Add(i, j int, v complex128)
}

// CDense is the complex counterpart of Dense for AC small-signal analysis:
// row-major storage, LU with partial pivoting on |entry|, no allocation after
// construction (Golub & Van Loan, Alg. 3.4.1, in complex arithmetic).
type CDense struct {
	n   int
	a   []complex128
	lu  []complex128
	piv []int
	col []float64

	factored bool
}

// NewCDense allocates an n×n complex system.
func NewCDense(n int) (*CDense, error) {
	if n <= 0 || n > MaxDense {
		return nil, fmt.Errorf("%w: %d", ErrBadSize, n)
	}
	return &CDense{n: n, a: make([]complex128, n*n), lu: make([]complex128, n*n), piv: make([]int, n),
		col: make([]float64, n)}, nil
}

// Size returns the number of unknowns.
func (d *CDense) Size() int { return d.n }

// Clear zeroes the matrix.
func (d *CDense) Clear() {
	d.factored = false
	for i := range d.a {
		d.a[i] = 0
	}
}

// Add adds v to entry (i, j); out-of-range indices (ground) are ignored.
func (d *CDense) Add(i, j int, v complex128) {
	if i < 0 || j < 0 || i >= d.n || j >= d.n {
		return
	}
	d.a[i*d.n+j] += v
}

// AbsProduct sets y_i = Σ_j |A_ij|·|x_j|.
func (d *CDense) AbsProduct(x []complex128, y []float64) {
	n := d.n
	for i := 0; i < n; i++ {
		s := 0.0
		for j, v := range d.a[i*n : (i+1)*n] {
			s += cmplx.Abs(v) * cmplx.Abs(x[j])
		}
		y[i] = s
	}
}

// At returns entry (i, j) of the stamped matrix.
func (d *CDense) At(i, j int) complex128 { return d.a[i*d.n+j] }

func cfinite(v complex128) bool {
	return !math.IsNaN(real(v)) && !math.IsNaN(imag(v)) && !math.IsInf(real(v), 0) && !math.IsInf(imag(v), 0)
}

// Solve factors a copy of the matrix and solves A·x = b.
func (d *CDense) Solve(b, x []complex128) error {
	n := d.n
	if len(b) != n || len(x) != n {
		return fmt.Errorf("%w: rhs %d, x %d, matrix %d", ErrBadSize, len(b), len(x), n)
	}
	if err := d.Factor(); err != nil {
		return err
	}
	return d.SolveFactored(b, x)
}

// Factor factors a copy of the stamped matrix in place, so several right-
// hand sides can then be solved with SolveFactored (diakoptics solves each
// torn piece for its source vector and every torn branch at once).
func (d *CDense) Factor() error {
	n := d.n
	d.factored = false
	copy(d.lu, d.a)
	scale := 0.0
	for _, v := range d.lu {
		if !cfinite(v) {
			return ErrNonFinite
		}
		if a := cmplx.Abs(v); a > scale {
			scale = a
		}
	}
	if scale == 0 {
		return ErrSingularMatrix
	}
	// Column-relative singularity test, as in Dense.Solve.
	for j := range d.col {
		d.col[j] = 0
	}
	for i := 0; i < n; i++ {
		for j, v := range d.lu[i*n : (i+1)*n] {
			if a := cmplx.Abs(v); a > d.col[j] {
				d.col[j] = a
			}
		}
	}
	lu := d.lu
	for k := 0; k < n; k++ {
		p := k
		maxv := cmplx.Abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := cmplx.Abs(lu[i*n+k]); v > maxv {
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
	d.factored = true
	return nil
}

// SolveFactored solves A·x = b with the factors of the last Factor. It
// only reads the factors (x is its working storage), so once a matrix is
// factored any number of goroutines may solve with it at the same time.
func (d *CDense) SolveFactored(b, x []complex128) error {
	n := d.n
	if len(b) != n || len(x) != n {
		return fmt.Errorf("%w: rhs %d, x %d, matrix %d", ErrBadSize, len(b), len(x), n)
	}
	if !d.factored {
		return fmt.Errorf("%w: not factored", ErrSingularMatrix)
	}
	lu := d.lu
	copy(x, b)
	for k := 0; k < n; k++ {
		if p := d.piv[k]; p != k {
			x[k], x[p] = x[p], x[k]
		}
	}
	for i := 0; i < n; i++ {
		s := x[i]
		for j := 0; j < i; j++ {
			s -= lu[i*n+j] * x[j]
		}
		x[i] = s
	}
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= lu[i*n+j] * x[j]
		}
		x[i] = s / lu[i*n+i]
	}
	for _, v := range x {
		if !cfinite(v) {
			return ErrNonFinite
		}
	}
	return nil
}
