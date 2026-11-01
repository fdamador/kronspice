// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import "fmt"

// solveNode is the node-pair form. Unknowns: the node potentials v (frame
// Node) and the currents of the apparent coils (frame Apparent) — the
// branches whose relation has no admittance form (voltage sources,
// controlled voltage sources, inductors shorted at DC; Kron 1939
// pp. 407–408) — and the inductors, whose admittance 1/(jωL) diverges as
// ω → 0: an inductor's current recovered from the difference of two node
// potentials would carry an error of ε·|v|/|ωL|, so the current is kept as
// an axis instead (it is the inductor's state variable, and at DC the
// inductor is exactly an apparent coil); likewise radiation edges (tiny
// resistances) and every branch whose current another relation senses (an
// F or H source's voltage source, the far end of a chord). Every other
// branch's current is
//
//	i = Y·E + Bₐ·iₐ + J        with E = A·v,
//
// so Kirchhoff's current law Aᵀ·i = 0 becomes
//
//	Y′·v + Aᵀ·Bₐ·iₐ = −Aᵀ·J,   Y′ = Aᵀ·Y·A   (Kron 1939 Ch. XIV p. 360)
//
// and each apparent coil adds its own relation as an equation.
//
// The products are formed from the sparse structure (each row of A has at
// most two entries and each branch relation a few terms), so the cost is
// proportional to the number of relation terms, not b²·N; the result is
// the same Y′ as Transform(A, Y) (TestNodeAssemblyEqualsTransform).
func (n *Network) solveNode(rel []relation) (*Solution, error) {
	ns, err := n.assembleNode(rel, nil)
	if err != nil {
		return nil, err
	}
	x, err := solveSystem(ns.dense(), ns.rhs)
	if err != nil {
		return nil, err
	}
	return ns.finish(x), nil
}

// nodeSystem is the assembled node-pair system and what is needed to turn
// its solution back into branch quantities.
type nodeSystem struct {
	rows       [][]entry // sparse rows, sorted by column
	rhs        []complex128
	nv, na     int
	from, to   []int    // each branch's node columns (−1 for the datum)
	yrow, brow [][]term // Y and Bₐ rows of the admittance-form branches (term.k: branch or apparent axis)
	J          []complex128
	app        []int // branch → apparent axis, or −1
}

// assembleNode builds the node-pair system. A branch marked in torn is kept
// in admittance form even if it is an inductor (diakoptics moves it to the
// interconnection, where its admittance 1/jωL is all that is needed).
func (n *Network) assembleNode(rel []relation, torn []bool) (*nodeSystem, error) {
	b := len(n.branches)
	nv := n.g.n - 1
	ns := &nodeSystem{nv: nv, from: make([]int, b), to: make([]int, b), yrow: make([][]term, b),
		brow: make([][]term, b), J: make([]complex128, b), app: make([]int, b)}
	for k, br := range n.branches {
		ns.from[k], ns.to[k] = br.from-1, br.to-1
	}
	// A branch whose current another relation senses (the voltage source of
	// an F or H, the far end of a chord) keeps its current as an axis.
	sensed := make([]bool, b)
	for k := range rel {
		for _, t := range rel[k].i {
			if t.k != k {
				sensed[t.k] = true
			}
		}
	}
	var appBranches []int
	for k := range ns.app {
		ns.app[k] = -1
		if torn != nil && torn[k] {
			continue
		}
		if _, ci := rel[k].own(k); ci == 0 || n.branches[k].kind == bL || n.branches[k].kind == bRadiation || sensed[k] {
			ns.app[k] = len(appBranches)
			appBranches = append(appBranches, k)
		}
	}
	na := len(appBranches)
	ns.na = na
	for k := 0; k < b; k++ {
		if ns.app[k] >= 0 {
			ns.brow[k] = []term{{ns.app[k], 1}}
			continue
		}
		_, ci := rel[k].own(k)
		for _, t := range rel[k].e {
			ns.yrow[k] = append(ns.yrow[k], term{t.k, -t.c / ci})
		}
		for _, t := range rel[k].i {
			if t.k == k {
				continue
			}
			if ns.app[t.k] < 0 {
				return nil, fmt.Errorf("%w: branch %d senses the current of a branch that is not an apparent coil", ErrUnsupported, k)
			}
			ns.brow[k] = append(ns.brow[k], term{ns.app[t.k], -t.c / ci})
		}
		ns.J[k] = rel[k].s / ci
	}
	size := nv + na
	sys := make([][]entry, size)
	rhs := make([]complex128, size)
	// Stamp twice: first counting each row's entries, then into one backing
	// array, so assembly makes a fixed number of allocations.
	count := make([]int, size)
	counting := true
	add := func(r, c int, v complex128) {
		if counting {
			count[r]++
			return
		}
		sys[r] = append(sys[r], entry{c, v})
	}
	for pass := 0; pass < 2; pass++ {
		if pass == 1 {
			total := 0
			for _, c := range count {
				total += c
			}
			back := make([]entry, total)
			at := 0
			for r, c := range count {
				sys[r] = back[at : at : at+c]
				at += c
			}
			counting = false
			for i := range rhs {
				rhs[i] = 0
			}
		}
		if err := ns.stamp(n, rel, appBranches, add, rhs); err != nil {
			return nil, err
		}
	}
	for r := range sys {
		sys[r] = compact(sys[r])
	}
	ns.rows, ns.rhs = sys, rhs
	return ns, nil
}

// stamp writes the node-pair system through add (KCL rows, then one row
// per apparent axis) and its right-hand side into rhs.
func (ns *nodeSystem) stamp(n *Network, rel []relation, appBranches []int, add func(r, c int, v complex128), rhs []complex128) error {
	b, nv := len(n.branches), ns.nv
	// KCL rows: Aᵀ·(Y·A·v + Bₐ·iₐ + J) = 0, row by row of A.
	kcl := func(k, r int, sign complex128) {
		for _, t := range ns.yrow[k] {
			if c := ns.from[t.k]; c >= 0 {
				add(r, c, sign*t.c)
			}
			if c := ns.to[t.k]; c >= 0 {
				add(r, c, -sign*t.c)
			}
		}
		for _, t := range ns.brow[k] {
			add(r, nv+t.k, sign*t.c)
		}
		rhs[r] -= sign * ns.J[k]
	}
	for k := 0; k < b; k++ {
		if r := ns.from[k]; r >= 0 {
			kcl(k, r, 1)
		}
		if r := ns.to[k]; r >= 0 {
			kcl(k, r, -1)
		}
	}
	for a, k := range appBranches {
		r := nv + a
		for _, t := range rel[k].e {
			if c := ns.from[t.k]; c >= 0 {
				add(r, c, t.c)
			}
			if c := ns.to[t.k]; c >= 0 {
				add(r, c, -t.c)
			}
		}
		for _, t := range rel[k].i {
			if ns.app[t.k] < 0 {
				return fmt.Errorf("%w: branch %d senses the current of a branch that is not an apparent coil", ErrUnsupported, k)
			}
			add(r, nv+ns.app[t.k], t.c)
		}
		rhs[nv+a] = rel[k].s
	}
	return nil
}

// entry is one stored coefficient of a sparse row.
type entry struct {
	j int
	v complex128
}

// compact merges a row's coefficients by column, adding them in the order
// they were stamped (so each sum is bit-identical to accumulating into a
// dense row) and sorting by column.
func compact(row []entry) []entry {
	// Stable insertion sort: rows hold a handful of entries.
	for i := 1; i < len(row); i++ {
		e := row[i]
		k := i
		for k > 0 && row[k-1].j > e.j {
			row[k] = row[k-1]
			k--
		}
		row[k] = e
	}
	out := row[:0]
	for _, e := range row {
		if n := len(out); n > 0 && out[n-1].j == e.j {
			out[n-1].v += e.v
		} else {
			out = append(out, e)
		}
	}
	return out
}

// dense returns the system matrix as a dense array.
func (ns *nodeSystem) dense() [][]complex128 {
	a := zeros(len(ns.rows), len(ns.rows))
	for i, row := range ns.rows {
		for _, e := range row {
			a[i][e.j] = e.v
		}
	}
	return a
}

// finish turns the solution x = (v, iₐ) into branch voltages and currents,
// with the frame invariants.
func (ns *nodeSystem) finish(x []complex128) *Solution {
	b := len(ns.from)
	pot := func(c int) complex128 {
		if c < 0 {
			return 0
		}
		return x[c]
	}
	E := NewVec[Primitive](b)
	for k := 0; k < b; k++ {
		E.v[k] = pot(ns.from[k]) - pot(ns.to[k])
	}
	I := NewVec[Primitive](b)
	for k := 0; k < b; k++ {
		s := ns.J[k]
		for _, t := range ns.yrow[k] {
			s += t.c * E.v[t.k]
		}
		for _, t := range ns.brow[k] {
			s += t.c * x[ns.nv+t.k]
		}
		I.v[k] = s
	}
	// The impressed currents J, as Kron's e·i in this frame: −J·E against
	// its transform −(AᵀJ)·v.
	J := Vec[Primitive]{v: ns.J}
	ATJv := complex128(0)
	for k := 0; k < b; k++ {
		if c := ns.from[k]; c >= 0 {
			ATJv += ns.J[k] * x[c]
		}
		if c := ns.to[k]; c >= 0 {
			ATJv -= ns.J[k] * x[c]
		}
	}
	return &Solution{Form: NodeForm, E: E, I: I,
		PrimitivePower: -Dot(J, E), FramePower: -ATJv, Tellegen: Dot(E, I)}
}
