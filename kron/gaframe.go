// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/hyper/ga"
)

// Kron's network in geometric algebra (CLAUDE.md §3b item 1, M9e).
//
// The branch space B has one basis vector per primitive branch; the
// primitive impedance (or admittance) tensor of a network of separate
// resistive coils is diagonal, so it is a metric on B: e_k·e_k = z_k. A
// mesh is a vector of B — column k of the connection tensor C, its
// components the branch currents of a unit mesh current (Kron 1939 eq
// 4.23) — and a node-pair axis is a column of A. The m mesh vectors span
// the loop space and the cut-set vectors its complement; that orthogonal
// pair of subspaces is Kron's orthogonal network in this language (Kron
// 1939 Ch. XVI; F. H. Branin, "The algebraic-topological basis for network
// analogies and the vector calculus," 1966; J. P. Roth, PNAS 41, 1955).
//
// The outer product of the mesh vectors is an m-blade of B. Its square is
// the Gram determinant of the vectors under the metric,
//
//	⟨B·B̃⟩ = det[cᵢ·z·cⱼ] = det(Cᵀ·Z·C) = det Z′,
//
// the determinant of Kron's transformed mesh impedance (D. Hestenes, G.
// Sobczyk, Clifford Algebra to Geometric Calculus, 1984, §1-3, the
// "outermorphism" determinant; L. Dorst et al., Geometric Algebra for
// Computer Science, 2007, §2.4 and §4.3). The same holds for the node-pair
// vectors and det(Aᵀ·Y·A). A vanishing blade is a dependent set of axes —
// a loop that has no impedance, or a group of nodes that no admittance
// holds — and the factor where the blade vanishes names it: the blade is
// built one vector at a time, B_k = B_{k−1} ∧ c_k, and |B_k| = |B_{k−1}|·
// |rej(c_k, B_{k−1})|, the rejection of c_k from the subspace already
// spanned (Dorst et al. §3.6). The ratio |rej|/|c_k| ∈ [0, 1] is how
// independent axis k is.

// GAFrame is one of the two frames of a linear network at DC, in the GA
// view: the axes as vectors of the branch space with its diagonal metric.
type GAFrame struct {
	// Kind is "mesh" (metric z: R for resistors, 0 for the inductors and
	// voltage sources that short at DC; capacitors and current sources,
	// open at DC, are Kron's junction pairs and not part of it) or "node"
	// (metric y = 1/R on resistors; inductors and voltage sources merge
	// their nodes; capacitors and current sources carry no admittance).
	Kind string
	// Branches names the branch axes (circuit element names).
	Branches []string
	// Metric is z_k (mesh) or y_k (node) of each branch axis.
	Metric []float64
	// Axes are the mesh or node-pair vectors, components on Branches.
	Axes [][]float64
	// AxisNames describes each axis: the branches of a loop, or the nodes
	// of a node-pair axis.
	AxisNames []string
}

// GAFrames returns the mesh and node-pair frames of the network at DC. It
// covers networks of resistors, inductors, capacitors and independent
// sources; controlled sources, transducers, mutual inductance and
// radiation edges are not a diagonal metric and give ErrUnsupported.
func (n *Network) GAFrames() (mesh, node *GAFrame, err error) {
	if len(n.chords) > 0 || len(n.radiation) > 0 {
		return nil, nil, fmt.Errorf("%w: GA frames need a diagonal primitive metric (no chords or radiation edges)", ErrUnsupported)
	}
	name := func(b branch) string { return n.c.Elements[b.elem].Name }
	var mb, nb []int // branches of each frame
	for k, b := range n.branches {
		switch b.kind {
		case bR:
			mb, nb = append(mb, k), append(nb, k)
		case bL, bV:
			mb = append(mb, k)
		case bC, bI:
		default:
			return nil, nil, fmt.Errorf("%w: GA frames: %s is not a passive element or independent source", ErrUnsupported, name(b))
		}
	}
	// Mesh frame: the loops of the graph of the branches that pass DC.
	mesh = &GAFrame{Kind: "mesh"}
	from, to := make([]int, len(mb)), make([]int, len(mb))
	for i, k := range mb {
		b := n.branches[k]
		from[i], to[i] = b.from, b.to
		mesh.Branches = append(mesh.Branches, name(b))
		if b.kind == bR {
			mesh.Metric = append(mesh.Metric, b.value)
		} else {
			mesh.Metric = append(mesh.Metric, 0)
		}
	}
	g := newGraph(n.nodes, from, to)
	c := g.loops()
	rows, cols := c.Dims()
	for j := 0; j < cols; j++ {
		v := make([]float64, rows)
		var parts []string
		for i := range v {
			v[i] = real(c.At(i, j))
			if v[i] != 0 {
				parts = append(parts, mesh.Branches[i])
			}
		}
		mesh.Axes = append(mesh.Axes, v)
		mesh.AxisNames = append(mesh.AxisNames, "loop "+strings.Join(parts, ", "))
	}
	// Node frame: merge the nodes joined by DC shorts, then the incidence
	// of the resistors on the merged nodes.
	parent := make([]int, n.nodes)
	for i := range parent {
		parent[i] = i
	}
	find := func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for _, b := range n.branches {
		if b.kind == bL || b.kind == bV {
			if a, c := find(b.from), find(b.to); a != c {
				if c == find(0) {
					a, c = c, a
				}
				parent[c] = a // keep the datum's class as the root
			}
		}
	}
	ground := find(0)
	index := map[int]int{}
	var classes [][]string
	for v := 0; v < n.nodes; v++ {
		r := find(v)
		if r == ground {
			continue
		}
		if _, ok := index[r]; !ok {
			index[r] = len(classes)
			classes = append(classes, nil)
		}
		if v < len(n.c.NodeNames) {
			classes[index[r]] = append(classes[index[r]], n.c.NodeNames[v])
		}
	}
	node = &GAFrame{Kind: "node"}
	for _, k := range nb {
		b := n.branches[k]
		node.Branches = append(node.Branches, name(b))
		node.Metric = append(node.Metric, 1/b.value)
	}
	for j := range classes {
		v := make([]float64, len(nb))
		for i, k := range nb {
			b := n.branches[k]
			if f := find(b.from); f != ground && index[f] == j {
				v[i]++
			}
			if t := find(b.to); t != ground && index[t] == j {
				v[i]--
			}
		}
		node.Axes = append(node.Axes, v)
		node.AxisNames = append(node.AxisNames, "node "+strings.Join(classes[j], "="))
	}
	return mesh, node, nil
}

// Gram returns the transformed tensor of the frame, [cᵢ·w·cⱼ]: Kron's
// Z′ = Cᵀ·Z·C (mesh) or Y′ = Aᵀ·Y·A (node).
func (f *GAFrame) Gram() [][]float64 {
	m := len(f.Axes)
	out := make([][]float64, m)
	for i := range out {
		out[i] = make([]float64, m)
		for j := range out[i] {
			out[i][j] = f.dot(f.Axes[i], f.Axes[j])
		}
	}
	return out
}

func (f *GAFrame) dot(u, v []float64) float64 {
	s := 0.0
	for k, w := range f.Metric {
		if u[k] != 0 && v[k] != 0 {
			s += u[k] * w * v[k]
		}
	}
	return s
}

// maxExplicitBranches bounds Blade: an explicit m-blade of a b-dimensional
// space can have C(b, m) terms.
const maxExplicitBranches = 24

// Blade returns the explicit outer product c₁∧c₂∧…∧cₘ of the axes in the
// branch algebra, and the algebra. The metric is carried by scaling each
// branch axis by √|w|: positive weights give the e_k² = +1 part of the
// signature, negative weights −1 and zero weights null vectors, so
// ⟨B·B̃⟩ = det(Gram) holds for any real diagonal metric. Limited to 24
// branches (the blade's size grows as C(b, m)); use Norm2 for larger frames.
func (f *GAFrame) Blade() (ga.Multivector, *ga.Algebra, error) {
	b := len(f.Branches)
	if b > maxExplicitBranches {
		return ga.Multivector{}, nil, fmt.Errorf("%w: explicit blade limited to %d branches, frame has %d", ErrUnsupported, maxExplicitBranches, b)
	}
	if b == 0 {
		return ga.Scalar(1), nil, nil
	}
	// Order the basis: positive, negative, then null weights.
	var pos, neg, null []int
	for k, w := range f.Metric {
		switch {
		case w > 0:
			pos = append(pos, k)
		case w < 0:
			neg = append(neg, k)
		default:
			null = append(null, k)
		}
	}
	alg, err := ga.New(len(pos), len(neg), len(null))
	if err != nil {
		return ga.Multivector{}, nil, err
	}
	order := append(append(append([]int(nil), pos...), neg...), null...)
	out := ga.Scalar(1)
	for _, ax := range f.Axes {
		c := make(map[ga.Blade]float64)
		for slot, k := range order {
			if ax[k] != 0 {
				s := math.Sqrt(math.Abs(f.Metric[k]))
				if s == 0 {
					s = 1
				}
				c[ga.Blade(1)<<slot] = ax[k] * s
			}
		}
		v, err := alg.FromBlades(c)
		if err != nil {
			return ga.Multivector{}, nil, err
		}
		out = alg.Outer(out, v)
	}
	return out, alg, nil
}

// Independence returns, for each axis in order, |rej(c_k, B_{k−1})|/|c_k|:
// the part of axis k not spanned by the earlier axes, as a fraction of its
// length, in the frame's metric (twice-iterated Gram–Schmidt rejection).
// It needs a non-negative metric (passive networks); for a negative weight
// it returns ok = false. The product of the squared rejections is
// ⟨B·B̃⟩ = det(Gram); a zero entry is an axis that depends on the earlier
// ones. Cost O(m²·b), about that of forming the Gram matrix.
func (f *GAFrame) Independence() (ratio []float64, norm2 float64, ok bool) {
	for _, w := range f.Metric {
		if w < 0 || math.IsNaN(w) {
			return nil, 0, false
		}
	}
	norm2 = 1
	var basis [][]float64 // orthonormalized rejections
	for _, ax := range f.Axes {
		r := append([]float64(nil), ax...)
		l0 := math.Sqrt(f.dot(ax, ax))
		for pass := 0; pass < 2; pass++ {
			for _, q := range basis {
				p := f.dot(r, q)
				for k := range r {
					r[k] -= p * q[k]
				}
			}
		}
		l := math.Sqrt(math.Max(f.dot(r, r), 0))
		rho := 0.0
		if l0 > 0 {
			rho = l / l0
		}
		// Treat a rejection at round-off level as zero, so an exactly
		// dependent axis reports 0 rather than 1e−16.
		if rho < 1e-13 {
			rho, l = 0, 0
		}
		ratio = append(ratio, rho)
		norm2 *= l * l
		if l > 0 {
			q := make([]float64, len(r))
			for k := range r {
				q[k] = r[k] / l
			}
			basis = append(basis, q)
		}
	}
	return ratio, norm2, true
}

// GramDet is det(Gram) by Gaussian elimination with partial pivoting: the
// conventional value the blade is checked against.
func GramDet(g [][]float64) float64 {
	m := len(g)
	a := make([][]float64, m)
	for i := range g {
		a[i] = append([]float64(nil), g[i]...)
	}
	d := 1.0
	for c := 0; c < m; c++ {
		p := c
		for r := c + 1; r < m; r++ {
			if math.Abs(a[r][c]) > math.Abs(a[p][c]) {
				p = r
			}
		}
		if a[p][c] == 0 {
			return 0
		}
		if p != c {
			a[c], a[p] = a[p], a[c]
			d = -d
		}
		d *= a[c][c]
		for r := c + 1; r < m; r++ {
			f := a[r][c] / a[c][c]
			for k := c; k < m; k++ {
				a[r][k] -= f * a[c][k]
			}
		}
	}
	return d
}

// PrecheckOptions selects what Precheck rejects.
type PrecheckOptions struct {
	// DC also rejects loops that have zero impedance only at DC (closed by
	// inductors or machine windings together with ideal voltage outputs):
	// set it when every analysis needs an operating point. Loops of ideal
	// voltage outputs alone are rejected always.
	DC bool
	// Conditioning also computes, for a linear circuit, how independent
	// the axes of its mesh and node frames are (PrecheckReport). It costs
	// about as much as a dense factorization of the frame's transformed
	// tensor, so it is off by default; the structural check alone is a
	// union–find pass over the elements.
	Conditioning bool
}

// PrecheckReport describes the frames of a circuit that passed.
type PrecheckReport struct {
	// Linear is true when conditioning was requested and the circuit had
	// GA frames (resistors, L, C, independent sources only); the fields
	// below are then set.
	Linear bool
	// MeshMin and NodeMin are the smallest axis independence ratios of
	// the mesh and node frames (1 = mutually orthogonal axes, → 0 = nearly
	// dependent: an ill-conditioned network); MeshAxis and NodeAxis name
	// those axes.
	MeshMin, NodeMin   float64
	MeshAxis, NodeAxis string
}

// Precheck rejects, before any matrix is assembled, a circuit whose mesh
// frame has an axis of zero impedance: a loop made only of ideal voltage
// outputs (V, E, H, voltage-output B elements) — and, with opt.DC, of
// inductors and windings too. In GA terms the mesh blade has a null factor
// (|rej| = 0 for every metric), so Kron's Z′ and the MNA matrix are
// singular whatever the element values. The error wraps ErrSingular and
// names the loop's elements. Floating nodes (a node-pair axis with no
// admittance) are circuit.Validate's job and are not repeated here.
//
// Precheck never rejects a circuit the simulator could solve: it reports
// only structural zeros, never ill-conditioning; the conditioning of a
// linear circuit is returned in the report for information.
func Precheck(c *circuit.Circuit, opt PrecheckOptions) (*PrecheckReport, error) {
	type edge struct {
		name     string
		from, to int
		dcOnly   bool
	}
	var zero []edge
	for _, e := range c.Elements {
		switch e.Kind {
		case circuit.VoltageSource, circuit.VCVS, circuit.CCVS:
			zero = append(zero, edge{e.Name, e.Nodes[0], e.Nodes[1], false})
		case circuit.Behavioral:
			if e.BOut == 'V' && !c.SelfControlled(e) {
				zero = append(zero, edge{e.Name, e.Nodes[0], e.Nodes[1], false})
			}
		case circuit.Inductor, circuit.NonlinearInductor:
			zero = append(zero, edge{e.Name, e.Nodes[0], e.Nodes[1], true})
		case circuit.Machine:
			w, _ := circuit.MachineWindings(e)
			for i, p := range w {
				zero = append(zero, edge{fmt.Sprintf("%s winding %d", e.Name, i+1), p[0], p[1], true})
			}
		}
	}
	// Voltage outputs first, so a loop of them is found even when DC-only
	// edges are also present.
	sort.SliceStable(zero, func(i, j int) bool { return !zero[i].dcOnly && zero[j].dcOnly })
	n := c.NumNodes()
	adj := make([][]int, n) // forest edges by node
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	find := func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for k, e := range zero {
		if e.dcOnly && !opt.DC {
			continue
		}
		if e.from < 0 || e.to < 0 || e.from >= n || e.to >= n {
			continue
		}
		if a, b := find(e.from), find(e.to); a != b {
			parent[a] = b
			adj[e.from] = append(adj[e.from], k)
			adj[e.to] = append(adj[e.to], k)
			continue
		}
		// The edge closes a loop through the forest: name it.
		path := forestPath(adj, e.from, e.to, func(k int) (int, int) { return zero[k].from, zero[k].to })
		names := []string{e.name}
		dc := e.dcOnly
		for _, p := range path {
			names = append(names, zero[p].name)
			dc = dc || zero[p].dcOnly
		}
		what := "ideal voltage outputs"
		if dc {
			what = "ideal voltage outputs and inductors (zero impedance at DC)"
		}
		return nil, fmt.Errorf("%w: loop of %s: %s — the mesh blade has a null factor, so Z′ and the MNA matrix are singular for any element values",
			ErrSingular, what, strings.Join(names, ", "))
	}
	rep := &PrecheckReport{MeshMin: 1, NodeMin: 1}
	if !opt.Conditioning {
		return rep, nil
	}
	if nw, err := Build(c); err == nil {
		if mesh, node, err := nw.GAFrames(); err == nil {
			rep.Linear = true
			for _, fr := range []*GAFrame{mesh, node} {
				r, _, ok := fr.Independence()
				if !ok {
					rep.Linear = false
					break
				}
				for i, v := range r {
					if fr.Kind == "mesh" && v < rep.MeshMin {
						rep.MeshMin, rep.MeshAxis = v, fr.AxisNames[i]
					}
					if fr.Kind == "node" && v < rep.NodeMin {
						rep.NodeMin, rep.NodeAxis = v, fr.AxisNames[i]
					}
				}
			}
		}
	}
	return rep, nil
}

// forestPath returns the forest edges on the path from u to v (both in the
// same tree), by breadth-first search.
func forestPath(adj [][]int, u, v int, ends func(k int) (int, int)) []int {
	prevEdge := make(map[int]int)
	prevNode := make(map[int]int)
	seen := map[int]bool{u: true}
	queue := []int{u}
	for len(queue) > 0 && !seen[v] {
		x := queue[0]
		queue = queue[1:]
		for _, k := range adj[x] {
			a, b := ends(k)
			y := a
			if y == x {
				y = b
			}
			if seen[y] {
				continue
			}
			seen[y], prevEdge[y], prevNode[y] = true, k, x
			queue = append(queue, y)
		}
	}
	var path []int
	for x := v; x != u; x = prevNode[x] {
		k, ok := prevEdge[x]
		if !ok {
			return nil
		}
		path = append(path, k)
	}
	return path
}
