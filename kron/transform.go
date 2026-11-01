// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import "fmt"

// Frame is a reference frame of the network (CLAUDE.md §3a rule 1). A
// tensor's frames are part of its Go type, so a transformation between
// the wrong frames does not compile.
type Frame interface{ frameName() string }

// Primitive is the frame of the unconnected branches ("coils"): one axis
// per branch (Kron 1939 Ch. III).
type Primitive struct{}

// Mesh is the frame of the closed paths: one axis per independent mesh
// (fundamental loop), m = b − N + S of them (Kron 1939 p. 75, eqs 3.1–3.2).
type Mesh struct{}

// Node is the frame of the node-pair (junction-pair) potentials measured
// from the datum (ground): N − 1 axes for a network connected to ground.
type Node struct{}

// Junction is the frame of the junction-pair voltages kept as unknowns in
// the mesh form: the voltages across branches whose current is impressed
// or controlled (current sources, controlled current sources, open
// capacitors at DC, control coils). Along these axes E is unknown and i is
// known, so Kron treats them as junction pairs, not meshes (Kron 1939
// p. 411).
type Junction struct{}

// Apparent is the frame of the currents kept as unknowns in the node form:
// the currents of branches with no admittance — voltage sources,
// controlled voltage sources, shorted inductors at DC — which Kron calls
// apparent coils (Kron 1939 pp. 407–408).
type Apparent struct{}

func (Primitive) frameName() string { return "primitive" }
func (Mesh) frameName() string      { return "mesh" }
func (Node) frameName() string      { return "node" }
func (Junction) frameName() string  { return "junction" }
func (Apparent) frameName() string  { return "apparent" }

// Vec is a vector of components along the axes of frame F.
type Vec[F Frame] struct{ v []complex128 }

// NewVec returns a zero vector of n components.
func NewVec[F Frame](n int) Vec[F] { return Vec[F]{v: make([]complex128, n)} }

// Len is the number of components.
func (x Vec[F]) Len() int { return len(x.v) }

// At returns component i.
func (x Vec[F]) At(i int) complex128 { return x.v[i] }

// Set sets component i.
func (x Vec[F]) Set(i int, c complex128) { x.v[i] = c }

// Dot is the invariant contraction Σ xᵢ·yᵢ of a covariant and a
// contravariant vector of the same frame (no complex conjugate: Kron's
// e·i; for phasors the complex power is formed by the caller).
func Dot[F Frame](x, y Vec[F]) complex128 {
	s := complex128(0)
	for i := range x.v {
		s += x.v[i] * y.v[i]
	}
	return s
}

// Mat is a rank-2 tensor that maps vectors of frame C into frame R (rows
// along R, columns along C), stored densely by rows.
type Mat[R, C Frame] struct {
	r, c int
	a    []complex128
}

// NewMat returns an r × c zero tensor.
func NewMat[R, C Frame](r, c int) Mat[R, C] {
	return Mat[R, C]{r: r, c: c, a: make([]complex128, r*c)}
}

// Dims returns the number of rows and columns.
func (m Mat[R, C]) Dims() (int, int) { return m.r, m.c }

// At returns entry (i, j).
func (m Mat[R, C]) At(i, j int) complex128 { return m.a[i*m.c+j] }

// Set sets entry (i, j).
func (m Mat[R, C]) Set(i, j int, v complex128) { m.a[i*m.c+j] = v }

// Add adds v to entry (i, j).
func (m Mat[R, C]) Add(i, j int, v complex128) { m.a[i*m.c+j] += v }

func dimErr(op string, a, b int) error {
	return fmt.Errorf("%w: %s: %d against %d", ErrDimension, op, a, b)
}

// Mul is the product a·b: frame K is contracted.
func Mul[R, K, C Frame](a Mat[R, K], b Mat[K, C]) (Mat[R, C], error) {
	if a.c != b.r {
		return Mat[R, C]{}, dimErr("Mul", a.c, b.r)
	}
	out := NewMat[R, C](a.r, b.c)
	for i := 0; i < a.r; i++ {
		for k := 0; k < a.c; k++ {
			x := a.a[i*a.c+k]
			if x == 0 {
				continue
			}
			row := b.a[k*b.c : (k+1)*b.c]
			dst := out.a[i*b.c : (i+1)*b.c]
			for j, y := range row {
				dst[j] += x * y
			}
		}
	}
	return out, nil
}

// TMul is the product aᵀ·b, contracting the row frame K of both: the
// covariant transformation of a covariant tensor's first index.
func TMul[K, R, C Frame](a Mat[K, R], b Mat[K, C]) (Mat[R, C], error) {
	if a.r != b.r {
		return Mat[R, C]{}, dimErr("TMul", a.r, b.r)
	}
	out := NewMat[R, C](a.c, b.c)
	for k := 0; k < a.r; k++ {
		arow := a.a[k*a.c : (k+1)*a.c]
		brow := b.a[k*b.c : (k+1)*b.c]
		for i, x := range arow {
			if x == 0 {
				continue
			}
			dst := out.a[i*b.c : (i+1)*b.c]
			for j, y := range brow {
				dst[j] += x * y
			}
		}
	}
	return out, nil
}

// MulVec is a·x: the contravariant transformation of x from frame C into
// frame R (for example i = C·i′, Kron 1939 eq 4.23).
func MulVec[R, C Frame](a Mat[R, C], x Vec[C]) (Vec[R], error) {
	if a.c != len(x.v) {
		return Vec[R]{}, dimErr("MulVec", a.c, len(x.v))
	}
	out := NewVec[R](a.r)
	for i := 0; i < a.r; i++ {
		s := complex128(0)
		for j, y := range a.a[i*a.c : (i+1)*a.c] {
			s += y * x.v[j]
		}
		out.v[i] = s
	}
	return out, nil
}

// TMulVec is aᵀ·x: the covariant transformation of x from frame K into
// frame C (for example e′ = Cₜ·e, Kron 1939 eqs 4.26–4.27).
func TMulVec[K, C Frame](a Mat[K, C], x Vec[K]) (Vec[C], error) {
	if a.r != len(x.v) {
		return Vec[C]{}, dimErr("TMulVec", a.r, len(x.v))
	}
	out := NewVec[C](a.c)
	for k := 0; k < a.r; k++ {
		y := x.v[k]
		if y == 0 {
			continue
		}
		for j, c := range a.a[k*a.c : (k+1)*a.c] {
			out.v[j] += c * y
		}
	}
	return out, nil
}

// Transform is Kron's transformation law of a primitive impedance or
// admittance tensor into frame F through the connection tensor c:
// z′ = Cₜ·z·C (Kron 1939 Ch. IV pp. 103–104, following eq 4.37; for a
// rectangular, singular C, Ch. V pp. 119–121). The same law gives the
// node-pair admittance Y′ = Aᵀ·Y·A with c = A (A = Cₜ⁻¹ in Kron's
// notation, Ch. XIV p. 360).
func Transform[F Frame](c Mat[Primitive, F], z Mat[Primitive, Primitive]) (Mat[F, F], error) {
	zc, err := Mul(z, c)
	if err != nil {
		return Mat[F, F]{}, err
	}
	return TMul(c, zc)
}
