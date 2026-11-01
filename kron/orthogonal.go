// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math/cmplx"
)

// Orthogonal is the frame of Kron's orthogonal network: closed paths
// (meshes, carrying the currents) and open paths (junction pairs,
// carrying the voltages) in one reference frame (Kron 1939 Ch. XVI). For b
// branches, N junctions and S sub-networks there are m = b − N + S closed
// paths and p = N − S open paths (Kron 1939 p. 75, eqs 3.1–3.2); together
// they span the branch space, so the transformation [C | P] is square and
// non-singular. The two families are orthogonal, Pᵀ·C = 0: a closed path
// carries no junction-pair voltage and an open path no mesh current.
type Orthogonal struct{}

func (Orthogonal) frameName() string { return "orthogonal" }

// OrthogonalFrame is a network's orthogonal reference frame.
type OrthogonalFrame struct {
	Closed Mat[Primitive, Mesh] // C, one column per closed path
	Open   Mat[Primitive, Node] // P, one column per open path (junction pair)
	// Square is [C | P], b × b.
	Square Mat[Primitive, Orthogonal]
}

// Frame builds the orthogonal frame. The open paths are the junction pairs
// from each junction to its sub-network's datum (ground for the
// sub-network that contains it): the incidence tensor with one datum
// column removed per sub-network, so P has N − S columns even for a torn
// or disconnected network.
func (n *Network) Frame() (*OrthogonalFrame, error) {
	g := n.g
	C := g.loops()
	isRoot := make([]bool, g.n)
	for _, r := range g.roots {
		isRoot[r] = true
	}
	col := make([]int, g.n)
	p := 0
	for x := 0; x < g.n; x++ {
		col[x] = -1
		if !isRoot[x] {
			col[x] = p
			p++
		}
	}
	b := len(g.from)
	P := NewMat[Primitive, Node](b, p)
	for k := 0; k < b; k++ {
		if c := col[g.from[k]]; c >= 0 {
			P.Add(k, c, 1)
		}
		if c := col[g.to[k]]; c >= 0 {
			P.Add(k, c, -1)
		}
	}
	_, m := C.Dims()
	if m+p != b {
		return nil, fmt.Errorf("%w: %d closed + %d open paths for %d branches", ErrDimension, m, p, b)
	}
	sq := NewMat[Primitive, Orthogonal](b, b)
	for k := 0; k < b; k++ {
		for j := 0; j < m; j++ {
			sq.Set(k, j, C.At(k, j))
		}
		for j := 0; j < p; j++ {
			sq.Set(k, m+j, P.At(k, j))
		}
	}
	return &OrthogonalFrame{Closed: C, Open: P, Square: sq}, nil
}

// Coordinates expresses a primitive vector x in the orthogonal frame:
// x = C·xc + P·xo, solving the square system [C | P]·(xc, xo) = x. Branch
// currents that satisfy KCL have xo = 0 (they lie in the closed-path
// space); branch voltages that satisfy KVL have xc = 0 (open-path space).
// A singular square frame returns ErrSingular.
func (f *OrthogonalFrame) Coordinates(x Vec[Primitive]) (closed Vec[Mesh], open Vec[Node], err error) {
	b, _ := f.Square.Dims()
	a := zeros(b, b)
	for i := 0; i < b; i++ {
		for j := 0; j < b; j++ {
			a[i][j] = f.Square.At(i, j)
		}
	}
	rhs := append([]complex128(nil), x.v...)
	y, err := solveSystem(a, rhs)
	if err != nil {
		return Vec[Mesh]{}, Vec[Node]{}, err
	}
	_, m := f.Closed.Dims()
	return Vec[Mesh]{v: y[:m]}, Vec[Node]{v: y[m:]}, nil
}

// Invariants are Kron's two power invariants of the orthogonal network
// (Kron 1939 p. 409, pp. 412–413): the power input e·i of the impressed
// emfs, carried by the closed paths, and the power output E·I of the
// impressed currents, carried by the open paths. Each is computed in the
// primitive frame and again from orthogonal coordinates (e′·i′ with
// e′ = Cₜ·e, and E′·I′ with I′ = Pᵀ·I); invariance says the pairs agree.
type Invariants struct {
	InputPrimitive, InputFrame   complex128 // e·i, e′·i′
	OutputPrimitive, OutputFrame complex128 // E·I, E′·I′
}

// SourcePower is the power the impressed sources deliver to the network,
// e·i − E·I: emfs push current around the closed paths, impressed currents
// are withdrawn across the open paths (Kron's outside loads, p. 407).
func (iv Invariants) SourcePower() complex128 { return iv.InputFrame - iv.OutputFrame }

// Invariants computes both invariants of a solution. Impressed emfs are the
// values of the voltage sources (Kron's e = −E on a source branch, so that
// E + e = 0); impressed currents are the values of the current sources.
// Controlled sources are not impressed: their power is model-supplied.
func (n *Network) Invariants(sol *Solution, f *OrthogonalFrame, omega float64, dc bool) (Invariants, error) {
	b := len(n.branches)
	e := NewVec[Primitive](b)
	I := NewVec[Primitive](b)
	for k, br := range n.branches {
		src := br.ac
		if dc {
			src = complex(br.value, 0)
		}
		switch br.kind {
		case bV:
			e.Set(k, -src)
		case bI:
			I.Set(k, src)
		}
	}
	var iv Invariants
	iv.InputPrimitive = Dot(e, sol.I)
	iv.OutputPrimitive = Dot(sol.E, I)
	ic, _, err := f.Coordinates(sol.I)
	if err != nil {
		return iv, err
	}
	_, vo, err := f.Coordinates(sol.E)
	if err != nil {
		return iv, err
	}
	ep, err := TMulVec(f.Closed, e)
	if err != nil {
		return iv, err
	}
	Ip, err := TMulVec(f.Open, I)
	if err != nil {
		return iv, err
	}
	iv.InputFrame = Dot(ep, ic)
	iv.OutputFrame = Dot(vo, Ip)
	return iv, nil
}

// orthogonalResidual is max |xo| / max |x| of a current vector's open-path
// coordinates (0 when it satisfies KCL), used by the tests.
func orthogonalResidual(f *OrthogonalFrame, x Vec[Primitive]) (float64, error) {
	_, xo, err := f.Coordinates(x)
	if err != nil {
		return 0, err
	}
	m, s := 0.0, maxAbs(x.v)
	for _, v := range xo.v {
		if a := cmplx.Abs(v); a > m {
			m = a
		}
	}
	if s == 0 {
		return m, nil
	}
	return m / s, nil
}

func maxAbs(x []complex128) float64 {
	m := 0.0
	for _, v := range x {
		if a := cmplx.Abs(v); a > m {
			m = a
		}
	}
	return m
}
