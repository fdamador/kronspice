// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import "fmt"

// solveMesh is the mesh form. Unknowns: the mesh currents i′ (frame Mesh)
// and the junction-pair voltages across the branches whose relation has
// no impedance form (current sources, controlled current sources,
// capacitors open at DC, control coils; frame Junction — along these axes
// E is unknown and i known, Kron 1939 p. 411) — and the capacitors, the
// dual of the node form's inductors: 1/(jωC) diverges as ω → 0 and a
// capacitor voltage recovered from a difference of mesh currents would
// carry an error of ε·|i|/|ωC|, so the voltage is kept as an axis (the
// capacitor's state variable; at DC the capacitor is exactly an open
// junction pair). Every other branch's
// voltage is
//
//	E = Z·i + Mⱼ·Eⱼ − e        with i = C·i′,
//
// e being Kron's impressed emf (E + e = Z·i, Kron 1939 p. 417), so
// Kirchhoff's voltage law Cₜ·E = 0 becomes
//
//	Z′·i′ + Cₜ·Mⱼ·Eⱼ = e′,   Z′ = Cₜ·Z·C,  e′ = Cₜ·e   (Kron 1939 Ch. IV
//	pp. 101–104, eqs 4.23–4.37; Ch. V pp. 119–121)
//
// and each junction axis adds its branch's relation as an equation.
func (n *Network) solveMesh(rel []relation) (*Solution, error) {
	b := len(n.branches)
	C := n.g.loops()
	_, m := C.Dims()
	jn := make([]int, b) // branch → junction axis, or −1
	var jnBranches []int
	for k := range jn {
		jn[k] = -1
		if ce, _ := rel[k].own(k); ce == 0 || n.branches[k].kind == bC {
			jn[k] = len(jnBranches)
			jnBranches = append(jnBranches, k)
		}
	}
	nj := len(jnBranches)
	Z := NewMat[Primitive, Primitive](b, b)
	Mj := NewMat[Primitive, Junction](b, nj)
	e := NewVec[Primitive](b) // Kron's impressed emf: E + e = Z·i
	for k := 0; k < b; k++ {
		if jn[k] >= 0 {
			Mj.Set(k, jn[k], 1)
			continue
		}
		ce, _ := rel[k].own(k)
		for _, t := range rel[k].i {
			Z.Add(k, t.k, -t.c/ce)
		}
		for _, t := range rel[k].e {
			if t.k == k {
				continue
			}
			if jn[t.k] < 0 {
				return nil, fmt.Errorf("%w: branch %d senses the voltage of a branch that is not a junction axis", ErrUnsupported, k)
			}
			Mj.Add(k, jn[t.k], -t.c/ce)
		}
		e.Set(k, -rel[k].s/ce)
	}
	Zp, err := Transform(C, Z)
	if err != nil {
		return nil, err
	}
	CTMj, err := TMul(C, Mj)
	if err != nil {
		return nil, err
	}
	ep, err := TMulVec(C, e)
	if err != nil {
		return nil, err
	}
	size := m + nj
	sys := zeros(size, size)
	rhs := make([]complex128, size)
	for r := 0; r < m; r++ {
		for c := 0; c < m; c++ {
			sys[r][c] = Zp.At(r, c)
		}
		for c := 0; c < nj; c++ {
			sys[r][m+c] = CTMj.At(r, c)
		}
		rhs[r] = ep.At(r)
	}
	for a, k := range jnBranches {
		row := sys[m+a]
		for _, t := range rel[k].i {
			for c := 0; c < m; c++ {
				row[c] += t.c * C.At(t.k, c)
			}
		}
		for _, t := range rel[k].e {
			if jn[t.k] < 0 {
				return nil, fmt.Errorf("%w: branch %d senses the voltage of a branch that is not a junction axis", ErrUnsupported, k)
			}
			row[m+jn[t.k]] += t.c
		}
		rhs[m+a] = rel[k].s
	}
	x, err := solveSystem(sys, rhs)
	if err != nil {
		return nil, err
	}
	ip := Vec[Mesh]{v: x[:m]}
	ej := Vec[Junction]{v: x[m:]}
	I, err := MulVec(C, ip)
	if err != nil {
		return nil, err
	}
	ZI, err := MulVec(Z, I)
	if err != nil {
		return nil, err
	}
	ME, err := MulVec(Mj, ej)
	if err != nil {
		return nil, err
	}
	E := NewVec[Primitive](b)
	for k := 0; k < b; k++ {
		E.v[k] = ZI.v[k] + ME.v[k] - e.v[k]
	}
	sol := &Solution{Form: MeshForm, E: E, I: I,
		PrimitivePower: Dot(e, I), FramePower: Dot(ep, ip), Tellegen: Dot(E, I)}
	return sol, nil
}

// potentials returns the node potentials of a solution (node 0 = datum) by
// integrating the branch voltages along the spanning tree from ground. A
// sub-network that does not contain the datum has no defined potential.
func (n *Network) potentials(E Vec[Primitive]) ([]complex128, error) {
	g := n.g
	v := make([]complex128, g.n)
	if g.subNetworks() > 1 {
		return nil, fmt.Errorf("%w: %d sub-networks are not connected to ground", ErrSingular, g.subNetworks()-1)
	}
	// Breadth-first order: parents before children.
	order := make([]int, 0, g.n)
	order = append(order, 0)
	children := make([][]int, g.n)
	for x := 1; x < g.n; x++ {
		children[g.parent[x]] = append(children[g.parent[x]], x)
	}
	for i := 0; i < len(order); i++ {
		order = append(order, children[order[i]]...)
	}
	for _, x := range order[1:] {
		k, p := g.parentEdge[x], g.parent[x]
		if g.from[k] == x { // E = V(x) − V(p)
			v[x] = v[p] + E.v[k]
		} else { // E = V(p) − V(x)
			v[x] = v[p] - E.v[k]
		}
	}
	return v, nil
}
