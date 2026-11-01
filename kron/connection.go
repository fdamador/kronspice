// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

// graph is the interconnection of the primitive branches: branch k runs
// from node From[k] to node To[k] (its current flows From → To through
// it, and its voltage is E = V(From) − V(To)). Node 0 is the datum.
type graph struct {
	n        int // nodes, including the datum
	from, to []int

	// Spanning forest, found once: tree[k] marks tree branches; for every
	// node, its parent and the tree branch to it (−1 at a root); comp is
	// the sub-network each node belongs to; roots[s] is sub-network s's
	// root (the datum for the sub-network that contains it).
	tree       []bool
	parent     []int
	parentEdge []int
	depth      []int
	comp       []int
	roots      []int
}

func newGraph(n int, from, to []int) *graph {
	g := &graph{n: n, from: from, to: to}
	g.forest()
	return g
}

// forest builds a breadth-first spanning forest, starting from the datum,
// so every node connected to ground gets its potential by a tree path
// from ground. Branch order decides ties, so the forest (and with it the
// mesh frame) is deterministic.
func (g *graph) forest() {
	b := len(g.from)
	adj := make([][]int, g.n)
	for k := 0; k < b; k++ {
		adj[g.from[k]] = append(adj[g.from[k]], k)
		if g.to[k] != g.from[k] {
			adj[g.to[k]] = append(adj[g.to[k]], k)
		}
	}
	g.tree = make([]bool, b)
	g.parent = make([]int, g.n)
	g.parentEdge = make([]int, g.n)
	g.depth = make([]int, g.n)
	g.comp = make([]int, g.n)
	for i := range g.comp {
		g.comp[i] = -1
	}
	queue := make([]int, 0, g.n)
	for r := 0; r < g.n; r++ {
		if g.comp[r] >= 0 {
			continue
		}
		s := len(g.roots)
		g.roots = append(g.roots, r)
		g.comp[r], g.parent[r], g.parentEdge[r] = s, -1, -1
		queue = append(queue[:0], r)
		for len(queue) > 0 {
			u := queue[0]
			queue = queue[1:]
			for _, k := range adj[u] {
				w := g.from[k]
				if w == u {
					w = g.to[k]
				}
				if g.comp[w] >= 0 {
					continue
				}
				g.comp[w], g.parent[w], g.parentEdge[w], g.depth[w] = s, u, k, g.depth[u]+1
				g.tree[k] = true
				queue = append(queue, w)
			}
		}
	}
}

// subNetworks is S, the number of separate sub-networks (connected
// components). A node that no branch touches counts as one.
func (g *graph) subNetworks() int { return len(g.roots) }

// incidence is the node-pair connection tensor A: E = A·v, one row per
// branch and one column per node other than the datum, with A[k][from] =
// +1 and A[k][to] = −1 (Kron's Cₜ⁻¹, Kron 1939 Ch. XIV p. 360; Tensors
// for Circuits eq 6.8). Kirchhoff's current law is Aᵀ·i = 0.
func (g *graph) incidence() Mat[Primitive, Node] {
	a := NewMat[Primitive, Node](len(g.from), g.n-1)
	for k := range g.from {
		if f := g.from[k]; f > 0 {
			a.Add(k, f-1, 1)
		}
		if t := g.to[k]; t > 0 {
			a.Add(k, t-1, -1)
		}
	}
	return a
}

// loops is the mesh connection tensor C: i = C·i′ (Kron 1939 eq 4.23),
// one column per fundamental loop. Each link (non-tree branch) closes one
// loop through the tree; the loop's direction is the link's, so C has +1
// on the link and ±1 on the tree branches of the loop, and Kirchhoff's
// current law holds identically: Aᵀ·C = 0. There are m = b − N + S
// columns (Kron 1939 p. 75, eqs 3.1–3.2). The columns span the null space
// of the junction constraints Aᵀ·i = 0, which is how Tensors for Circuits
// Ch. 7 (pp. 46–50) builds C.
func (g *graph) loops() Mat[Primitive, Mesh] {
	b := len(g.from)
	var links []int
	for k := 0; k < b; k++ {
		if !g.tree[k] {
			links = append(links, k)
		}
	}
	c := NewMat[Primitive, Mesh](b, len(links))
	for col, k := range links {
		c.Set(k, col, 1)
		// The loop runs from → to along the link, then back from `to` to
		// `from` through the tree: up from `to` to the common ancestor,
		// then down to `from`.
		u, w := g.from[k], g.to[k]
		for g.depth[w] > g.depth[u] {
			w = g.climb(c, col, w, true)
		}
		for g.depth[u] > g.depth[w] {
			u = g.climb(c, col, u, false)
		}
		for u != w {
			w = g.climb(c, col, w, true)
			u = g.climb(c, col, u, false)
		}
	}
	return c
}

// climb moves one tree step from node x to its parent, adding the tree
// branch to loop column col: up (x → parent) when the loop traverses it in
// that direction, otherwise down (parent → x). It returns the parent.
func (g *graph) climb(c Mat[Primitive, Mesh], col, x int, up bool) int {
	k, p := g.parentEdge[x], g.parent[x]
	sign := 1.0
	if up != (g.from[k] == x) {
		sign = -1 // the branch points the other way
	}
	c.Add(k, col, complex(sign, 0))
	return p
}
