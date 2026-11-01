// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/solver"
)

// Form selects the reference frame the network equations are solved in.
type Form uint8

const (
	// MeshForm solves for mesh currents: i = C·i′, Z′ = Cₜ·Z·C,
	// e′ = Cₜ·e (Kron 1939 Ch. IV–V).
	MeshForm Form = iota
	// NodeForm solves for node-pair potentials: E = A·v, Y′ = Aᵀ·Y·A
	// (Kron 1939 Ch. XIV).
	NodeForm
)

func (f Form) String() string {
	if f == NodeForm {
		return "node"
	}
	return "mesh"
}

// Solution is the network's state in the primitive frame: every branch's
// voltage E and current i (phasors for AC, real values for DC), together
// with the transformed tensors and the frame-invariant powers.
type Solution struct {
	Form Form
	E, I Vec[Primitive]
	// Invariant: the impressed power e·i in the primitive frame and the
	// same contraction e′·i′ in the solving frame; Kron's power invariance
	// says they are equal (Kron 1939 eq 4.24). Tellegen is Σ E·i over all
	// branches, zero for any solution (CLAUDE.md §3c).
	PrimitivePower, FramePower complex128
	Tellegen                   complex128
}

// Solve solves the network at angular frequency omega (dc for the
// operating point) in the chosen form.
func (n *Network) Solve(form Form, omega float64, dc bool) (*Solution, error) {
	rel := n.relations(omega, dc)
	switch form {
	case MeshForm:
		return n.solveMesh(rel)
	case NodeForm:
		return n.solveNode(rel)
	}
	return nil, fmt.Errorf("%w: form %d", ErrUnsupported, form)
}

// solveSystem factors and solves a square complex system (never forming an
// inverse, CLAUDE.md §3a rule 2). The rows and columns are first scaled to
// unit maximum (equilibration): the augmented Kron systems mix impedances
// of very different size (an acoustic radiation resistance of 10⁹ Pa·s/m³
// next to unit coefficients of a transducer row), and without scaling the
// LU's relative pivot test would see a well-determined system as
// singular. Two steps of iterative refinement follow, because branch
// quantities recovered afterwards (an inductor's current from two node
// potentials) magnify the solve's error.
func solveSystem(a [][]complex128, rhs []complex128) ([]complex128, error) {
	size := len(rhs)
	x := make([]complex128, size)
	if size == 0 {
		return x, nil
	}
	rs, cs := make([]float64, size), make([]float64, size)
	for i, row := range a {
		mx := 0.0
		for _, v := range row {
			mx = math.Max(mx, cmplx.Abs(v))
		}
		if mx == 0 {
			return nil, fmt.Errorf("%w: equation %d is empty", ErrSingular, i)
		}
		rs[i] = 1 / mx
	}
	for j := 0; j < size; j++ {
		mx := 0.0
		for i := range a {
			mx = math.Max(mx, cmplx.Abs(a[i][j])*rs[i])
		}
		if mx == 0 {
			return nil, fmt.Errorf("%w: unknown %d appears in no equation", ErrSingular, j)
		}
		cs[j] = 1 / mx
	}
	m, err := solver.NewCMatrix(size, solver.Auto)
	if err != nil {
		return nil, err
	}
	for i, row := range a {
		for j, v := range row {
			if v != 0 {
				m.Add(i, j, v*complex(rs[i]*cs[j], 0))
			}
		}
	}
	b := make([]complex128, size)
	for i := range b {
		b[i] = rhs[i] * complex(rs[i], 0)
	}
	y := make([]complex128, size)
	if err := m.Solve(b, y); err != nil {
		if errors.Is(err, solver.ErrSingularMatrix) {
			return nil, fmt.Errorf("%w: %v", ErrSingular, err)
		}
		return nil, err
	}
	for j := range x {
		x[j] = y[j] * complex(cs[j], 0)
	}
	r := make([]complex128, size)
	d := make([]complex128, size)
	for pass := 0; pass < 2; pass++ {
		for i, row := range a {
			s := rhs[i]
			for j, v := range row {
				s -= v * x[j]
			}
			r[i] = s * complex(rs[i], 0)
		}
		if err := m.Solve(r, d); err != nil {
			return nil, err
		}
		for j := range x {
			x[j] += d[j] * complex(cs[j], 0)
		}
	}
	return x, nil
}

func zeros(r, c int) [][]complex128 {
	a := make([][]complex128, r)
	for i := range a {
		a[i] = make([]complex128, c)
	}
	return a
}
