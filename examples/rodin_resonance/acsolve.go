// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
	"math/cmplx"
)

// Complex AC solution of the sectional network with its drives and the
// output groups, by modified nodal analysis (Ho, Ruehli & Brennan 1975):
// unknowns are the node voltages and the currents of every inductive
// branch (sections and group coils), so mutual inductances of either sign
// enter directly. kronspice's K element takes 0 < k ≤ 1 only, and its
// N element (any inductance matrix) is transient-only and holds at most
// 64 windings, so the full network (288 sections) is solved here and
// twinned against a kronspice transient of a reduced network
// (TestACTwin).

// csolve solves a·x = b in place by Gaussian elimination with partial
// pivoting.
func csolve(a [][]complex128, b []complex128) ([]complex128, error) {
	n := len(a)
	for c := 0; c < n; c++ {
		p, best := c, cmplx.Abs(a[c][c])
		for r := c + 1; r < n; r++ {
			if v := cmplx.Abs(a[r][c]); v > best {
				p, best = r, v
			}
		}
		if best == 0 {
			return nil, errors.New("singular AC matrix")
		}
		a[c], a[p] = a[p], a[c]
		b[c], b[p] = b[p], b[c]
		piv := a[c][c]
		rowc := a[c]
		parallelFor(n-c-1, func(k int) {
			r := c + 1 + k
			f := a[r][c] / piv
			if f == 0 {
				return
			}
			row := a[r]
			for j := c; j < n; j++ {
				row[j] -= f * rowc[j]
			}
			b[r] -= f * b[c]
		})
	}
	x := make([]complex128, n)
	for r := n - 1; r >= 0; r-- {
		s := b[r]
		for j := r + 1; j < n; j++ {
			s -= a[r][j] * x[j]
		}
		x[r] = s / a[r][r]
	}
	return x, nil
}

// Source is a Thevenin drive (amplitude, phase in degrees, resistance)
// between node Node and ground.
type Source struct {
	Node     int
	V        float64
	PhaseDeg float64
	R        float64
}

// Group is an output coil coupled to the sections: its self inductance
// and resistance, the mutual to every section (same order as Net.L) and to
// the other groups, and its series tuning capacitance (0 = none).
type Group struct {
	L, R, Ctune float64
	MSec        []float64
}

// Combine selects how the groups feed the load(s).
type Combine int

const (
	Separate Combine = 1 // each group into its own load (own bridge)
	Series   Combine = 2 // all groups in series into one load (one bridge)
	Parallel Combine = 3 // all groups in parallel across one load (one bridge)
)

func (c Combine) String() string {
	return map[Combine]string{Separate: "separate (one bridge per group)", Series: "series (one bridge)", Parallel: "parallel (one bridge)"}[c]
}

// ACResult is the AC solution.
type ACResult struct {
	V       []complex128 // coil node voltages
	I       []complex128 // section currents
	IG      []complex128 // group currents
	ISrc    []complex128 // source currents (out of the source into the node)
	PSrc    []float64    // average power delivered by each source [W]
	PLoad   float64      // average power into the load(s) [W]
	PCopper float64      // sections + groups wire heat [W]
	PSrcR   float64      // heat in the source resistances [W]
	VLoad   []complex128 // load voltages
}

// SolveAC solves the network at f with the sources, the groups (with
// their mutual matrix MG among themselves), the combination and the load
// resistance.
func (n *Net) SolveAC(f float64, src []Source, groups []Group, mg [][]float64, comb Combine, rload float64) (*ACResult, error) {
	w := 2 * math.Pi * f
	jw := complex(0, w)
	nb := len(n.L)
	ng := len(groups)
	// output nodes: separate: one per group; series: ng−1 internal + top; parallel: one
	var nout int
	switch comb {
	case Separate:
		nout = ng
	case Series:
		nout = ng
	case Parallel:
		nout = 1
	}
	if ng == 0 {
		nout = 0
	}
	nn := n.Nodes + nout
	size := nn + nb + ng
	a := make([][]complex128, size)
	for i := range a {
		a[i] = make([]complex128, size)
	}
	b := make([]complex128, size)
	// coil capacitances
	for i := 0; i < n.Nodes; i++ {
		for j := 0; j < n.Nodes; j++ {
			a[i][j] += jw * complex(n.CNode[i][j], 0)
		}
	}
	// sources (Norton form)
	for _, s := range src {
		g := complex(1/s.R, 0)
		a[s.Node][s.Node] += g
		b[s.Node] += complex(s.V/s.R, 0) * cmplx.Rect(1, s.PhaseDeg*math.Pi/180)
	}
	// sections: branch g current is unknown nn+g, flows from Branches[g][0] to [1]
	for g := 0; g < nb; g++ {
		row := nn + g
		fr, to := n.Branches[g][0], n.Branches[g][1]
		if fr >= 0 {
			a[fr][row] += 1
			a[row][fr] += 1
		}
		if to >= 0 {
			a[to][row] -= 1
			a[row][to] -= 1
		}
		// −(R + jωL)·I − jω Σ M·I_other = 0  (v_fr − v_to − Z·I = 0)
		a[row][row] -= complex(n.R[g], 0)
		for h := 0; h < nb; h++ {
			a[row][nn+h] -= jw * complex(n.L[g][h], 0)
		}
		for k := 0; k < ng; k++ {
			a[row][nn+nb+k] -= jw * complex(groups[k].MSec[g], 0)
		}
	}
	// groups: branch k between output nodes; impedance R + jωL + 1/(jωC)
	gnode := func(k int) (int, int) { // (from, to), −1 = ground, output node ids offset n.Nodes
		switch comb {
		case Separate:
			return n.Nodes + k, -1
		case Series: // chain: ground → o1 → … → o_ng (top); group k from o_k to o_{k+1}, o_0 = ground
			fr := -1
			if k > 0 {
				fr = n.Nodes + k - 1
			}
			return n.Nodes + k, fr
		default:
			return n.Nodes, -1
		}
	}
	for k := 0; k < ng; k++ {
		row := nn + nb + k
		fr, to := gnode(k)
		if fr >= 0 {
			a[fr][row] += 1
			a[row][fr] += 1
		}
		if to >= 0 {
			a[to][row] -= 1
			a[row][to] -= 1
		}
		z := complex(groups[k].R, 0)
		if groups[k].Ctune > 0 {
			z += 1 / (jw * complex(groups[k].Ctune, 0))
		}
		a[row][row] -= z
		for h := 0; h < ng; h++ {
			a[row][nn+nb+h] -= jw * complex(mg[k][h], 0)
		}
		for g := 0; g < nb; g++ {
			a[row][nn+g] -= jw * complex(groups[k].MSec[g], 0)
		}
	}
	// loads
	var loads [][2]int
	if ng > 0 {
		switch comb {
		case Separate:
			for k := 0; k < ng; k++ {
				loads = append(loads, [2]int{n.Nodes + k, -1})
			}
		case Series:
			loads = append(loads, [2]int{n.Nodes + ng - 1, -1})
		case Parallel:
			loads = append(loads, [2]int{n.Nodes, -1})
		}
	}
	gl := complex(1/rload, 0)
	for _, l := range loads {
		a[l[0]][l[0]] += gl
	}
	// a tiny conductance on the series chain's internal nodes keeps them referenced
	if comb == Series {
		for k := 0; k+1 < ng; k++ {
			a[n.Nodes+k][n.Nodes+k] += 1e-12
		}
	}
	x, err := csolve(a, b)
	if err != nil {
		return nil, err
	}
	r := &ACResult{V: x[:n.Nodes], I: x[nn : nn+nb], IG: x[nn+nb:]}
	for _, s := range src {
		vs := complex(s.V, 0) * cmplx.Rect(1, s.PhaseDeg*math.Pi/180)
		i := (vs - x[s.Node]) / complex(s.R, 0)
		r.ISrc = append(r.ISrc, i)
		r.PSrc = append(r.PSrc, 0.5*real(vs*cmplx.Conj(i)))
		r.PSrcR += 0.5 * s.R * real(i*cmplx.Conj(i))
	}
	for g := 0; g < nb; g++ {
		r.PCopper += 0.5 * n.R[g] * real(r.I[g]*cmplx.Conj(r.I[g]))
	}
	for k := 0; k < ng; k++ {
		r.PCopper += 0.5 * groups[k].R * real(r.IG[k]*cmplx.Conj(r.IG[k]))
	}
	for _, l := range loads {
		v := x[l[0]]
		r.VLoad = append(r.VLoad, v)
		r.PLoad += 0.5 * real(v*cmplx.Conj(v)) / rload
	}
	return r, nil
}
