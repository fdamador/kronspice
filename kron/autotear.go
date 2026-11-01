// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// AutoTear chooses a cut set that splits c into about `pieces` torn
// pieces of similar size: the automatic graph partitioner of CLAUDE.md
// §3a rule 5 (manual cut sets stay available through Tearing).
//
// Only two-terminal R, C and L branches can be torn (an inductor only for
// AC: at DC it is a short, dc = true). Every other element keeps its
// nodes in one piece: sources, controlled sources together with their
// controlling nodes or controlling voltage source, nonlinear devices, and
// both inductors of a mutual coupling. Ground is the common reference of
// all pieces (the single-common-bus arrangement of Happ Ch. XII §4), so
// branches to ground are never cut.
//
// Method: nodes that must stay together are merged; the merged graph is
// ordered by breadth-first search from a pseudo-peripheral node (two
// sweeps; G. H. Golub, C. F. Van Loan, Matrix Computations §11.1, the
// Cuthill–McKee start), and the order is cut into contiguous runs of
// about equal node count. Torn are the tearable branches that join two
// runs. The result is a proposal: Tear and Solve report a piece that
// floats or is singular, and AutoTearSolvable retries with fewer pieces.
func AutoTear(c *circuit.Circuit, pieces int, dc bool) (Tearing, error) {
	n := len(c.NodeNames)
	if pieces < 2 || n < 3 {
		return Tearing{}, fmt.Errorf("%w: AutoTear needs at least 2 pieces and 2 nodes", ErrUnsupported)
	}
	par := make([]int, n)
	for i := range par {
		par[i] = i
	}
	var find func(int) int
	find = func(a int) int {
		for par[a] != a {
			par[a] = par[par[a]]
			a = par[a]
		}
		return a
	}
	join := func(a, b int) {
		if a == 0 || b == 0 {
			return // ground joins nothing: it is every piece's reference
		}
		par[find(a)] = find(b)
	}
	index := map[string]int{}
	for i, e := range c.Elements {
		index[strings.ToLower(e.Name)] = i
	}
	// At DC a capacitor is an open circuit: it joins nothing and carries
	// no current, so it is neither cut nor used to keep nodes together.
	open := func(e circuit.Element) bool { return dc && e.Kind == circuit.Capacitor }
	tearable := func(e circuit.Element) bool {
		switch e.Kind {
		case circuit.Resistor, circuit.Capacitor:
			return true
		case circuit.Inductor:
			return !dc
		}
		return false
	}
	coupled := map[string]bool{}
	for _, e := range c.Elements {
		if e.Kind == circuit.Mutual {
			coupled[strings.ToLower(e.Couples[0])] = true
			coupled[strings.ToLower(e.Couples[1])] = true
		}
	}
	nodesOf := func(e circuit.Element) []int {
		out := []int{e.Nodes[0], e.Nodes[1]}
		out = append(out, e.Extra...)
		if e.Kind.IsControlled() {
			if e.Kind.VoltageControlled() {
				out = append(out, e.CtrlNodes[0], e.CtrlNodes[1])
			} else if k, ok := index[strings.ToLower(e.CtrlSource)]; ok {
				out = append(out, c.Elements[k].Nodes[0], c.Elements[k].Nodes[1])
			}
		}
		if e.Kind == circuit.Mutual {
			for _, l := range e.Couples {
				if k, ok := index[strings.ToLower(l)]; ok {
					out = append(out, c.Elements[k].Nodes[0], c.Elements[k].Nodes[1])
				}
			}
		}
		if e.BCtrl == 'V' {
			out = append(out, e.CtrlNodes[0], e.CtrlNodes[1])
		}
		return out
	}
	for _, e := range c.Elements {
		if open(e) || tearable(e) && !coupled[strings.ToLower(e.Name)] {
			continue
		}
		ns := nodesOf(e)
		var first int
		for _, v := range ns {
			if v == 0 {
				continue
			}
			if first == 0 {
				first = v
				continue
			}
			join(first, v)
		}
	}
	// Merged graph over the groups of non-ground nodes.
	adj := map[int]map[int]bool{}
	weight := map[int]int{}
	for v := 1; v < n; v++ {
		r := find(v)
		weight[r]++
		if adj[r] == nil {
			adj[r] = map[int]bool{}
		}
	}
	for _, e := range c.Elements {
		if open(e) || !tearable(e) || coupled[strings.ToLower(e.Name)] {
			continue
		}
		a, b := e.Nodes[0], e.Nodes[1]
		if a == 0 || b == 0 {
			continue
		}
		ra, rb := find(a), find(b)
		if ra != rb {
			adj[ra][rb] = true
			adj[rb][ra] = true
		}
	}
	var groups []int
	for r := range adj {
		groups = append(groups, r)
	}
	sort.Ints(groups) // determinism
	if len(groups) < pieces {
		return Tearing{}, fmt.Errorf("%w: the elements that cannot be torn (sources, controlled sources with their controlling nodes, coupled inductors, nonlinear devices) bind the %d nodes into %d group(s); %d pieces asked", ErrUnsupported, n-1, len(groups), pieces)
	}
	neighbours := func(r int) []int {
		var out []int
		for s := range adj[r] {
			out = append(out, s)
		}
		sort.Ints(out)
		return out
	}
	bfs := func(start int, seen map[int]bool) []int {
		order := []int{start}
		seen[start] = true
		for i := 0; i < len(order); i++ {
			for _, s := range neighbours(order[i]) {
				if !seen[s] {
					seen[s] = true
					order = append(order, s)
				}
			}
		}
		return order
	}
	// Order every connected component, each from a pseudo-peripheral
	// start (the last node of a first sweep).
	seen := map[int]bool{}
	var order []int
	for _, g := range groups {
		if seen[g] {
			continue
		}
		probe := bfs(g, map[int]bool{})
		comp := bfs(probe[len(probe)-1], seen)
		order = append(order, comp...)
	}
	total := 0
	for _, g := range order {
		total += weight[g]
	}
	piece := map[int]int{}
	acc, k := 0, 0
	for _, g := range order {
		if acc >= (k+1)*total/pieces && k < pieces-1 {
			k++
		}
		piece[g] = k
		acc += weight[g]
	}
	var torn []string
	for _, e := range c.Elements {
		if open(e) || !tearable(e) || coupled[strings.ToLower(e.Name)] {
			continue
		}
		a, b := e.Nodes[0], e.Nodes[1]
		if a == 0 || b == 0 {
			continue
		}
		if piece[find(a)] != piece[find(b)] {
			torn = append(torn, e.Name)
		}
	}
	if len(torn) == 0 {
		return Tearing{}, fmt.Errorf("%w: no R, C or L branch between non-ground nodes joins the pieces (at DC capacitors are open and inductors shorts)", ErrUnsupported)
	}
	return Tearing{Torn: torn}, nil
}

// AutoTearSolvable tries AutoTear with `pieces`, then fewer, until the
// torn network solves at omega (dc for the operating point). It returns
// the tearing, the solver's statistics and the torn solution.
func AutoTearSolvable(c *circuit.Circuit, pieces int, omega float64, dc bool) (Tearing, *Diakoptics, *Solution, TearStats, error) {
	var last error
	for p := pieces; p >= 2; p-- {
		t, err := AutoTear(c, p, dc)
		if err != nil {
			last = err
			continue
		}
		d, err := Tear(c, t, nil)
		if err != nil {
			last = err
			continue
		}
		sol, st, err := d.Solve(omega, dc)
		if err != nil {
			last = err
			continue
		}
		return t, d, sol, st, nil
	}
	return Tearing{}, nil, nil, TearStats{}, fmt.Errorf("kron: no automatic tearing of %d down to 2 pieces solves: %w", pieces, last)
}

// ElementBranch returns the primitive branch of the named circuit element
// (case-insensitive), so a solution's E and I can be read per element.
func (n *Network) ElementBranch(name string) (int, bool) {
	for i, e := range n.c.Elements {
		if strings.EqualFold(e.Name, name) {
			b := n.elemBranch[i]
			return b, b >= 0
		}
	}
	return 0, false
}
