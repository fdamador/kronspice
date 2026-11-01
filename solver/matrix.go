// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package solver holds the linear-algebra backends for modified nodal
// analysis (CLAUDE.md §3).
//
// Two backends sit behind the Matrix interface: a dense LU with partial
// pivoting (Dense, CDense) and a sparse Gilbert–Peierls LU with a
// minimum-degree ordering (Sparse, CSparse), both written against the
// standard library only. Auto picks one by size; every backend passes the same tests.
package solver

import "errors"

// Errors returned by solvers. Sentinel values so the GA can map them to penalties.
var (
	ErrSingularMatrix = errors.New("solver: singular matrix")
	ErrBadSize        = errors.New("solver: invalid matrix size")
	ErrNonFinite      = errors.New("solver: non-finite value in system")
)

// Matrix is an MNA system matrix that can be cleared, stamped and solved.
type Matrix interface {
	// Size returns the number of unknowns.
	Size() int
	// Clear sets every entry to zero, keeping storage.
	Clear()
	// Add adds v to entry (i, j). Indices are 0-based unknown indices.
	Add(i, j int, v float64)
	// Solve solves A·x = b. It must not allocate after construction.
	Solve(b, x []float64) error
	// AbsProduct sets y_i = Σ_j |A_ij|·|x_j| for the stamped matrix: the
	// size of the terms that cancel in each row, from which the analyses
	// estimate round-off in their energy ledgers.
	AbsProduct(x, y []float64)
}

// CMatrix is the complex counterpart of Matrix, for AC analysis.
type CMatrix interface {
	Size() int
	Clear()
	Add(i, j int, v complex128)
	Solve(b, x []complex128) error
	// AbsProduct sets y_i = Σ_j |A_ij|·|x_j| (see Matrix.AbsProduct).
	AbsProduct(x []complex128, y []float64)
}

// Backend selects a linear solver.
type Backend uint8

const (
	// Auto picks Dense below DenseThreshold unknowns and Sparse above.
	Auto Backend = iota
	// DenseLU forces the dense backend.
	DenseLU
	// SparseLU forces the sparse Gilbert–Peierls backend.
	SparseLU
)

func (b Backend) String() string {
	switch b {
	case DenseLU:
		return "dense"
	case SparseLU:
		return "sparse"
	}
	return "auto"
}

// DenseThreshold is the size from which Auto uses the sparse backend. Below
// it the dense LU is as fast or faster (CLAUDE.md §3, "Linear solver
// backends"; see BenchmarkGridSolve and doc/M6.md), and most evolved
// circuits are this small.
const DenseThreshold = 30

// Resolve returns the backend Auto stands for at size n.
func (b Backend) Resolve(n int) Backend {
	if b != Auto {
		return b
	}
	if n < DenseThreshold {
		return DenseLU
	}
	return SparseLU
}

// NewMatrix returns an n×n real system on backend b.
func NewMatrix(n int, b Backend) (Matrix, error) {
	if b.Resolve(n) == SparseLU {
		return NewSparse(n)
	}
	return NewDense(n)
}

// NewCMatrix returns an n×n complex system on backend b.
func NewCMatrix(n int, b Backend) (CMatrix, error) {
	if b.Resolve(n) == SparseLU {
		return NewCSparse(n)
	}
	return NewCDense(n)
}
