// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"math/bits"
	"sort"
)

// Multi-dimensional networks: Kron's polyhedron of q-networks in any
// number of dimensions k.
//
// Sources. G. Kron, "Multidimensional curve-fitting with self-organizing
// automata," J. Math. Anal. Appl. 5, 46–69, 1962 ("K62"), pp. 53–55: a
// k-dimensional polyhedron whose p-networks are copies of his 1944 Maxwell
// network, even-dimensional ones magnetic and odd-dimensional ones
// dielectric, joined to a dual polyhedron with the roles exchanged, one
// complete wave spanning two neighbouring dimensions. The rule for q ≥ 2
// is Kron's 1959 space-filter diagram as reproduced by N. V.
// Balasubramanian, PhD thesis, Liverpool 1968, Fig. 9 p. 152 ("B68"): an
// orthogonal network per dimension, ε on odd and μ on even dimensions,
// the incidence matrices M^q_{q+1} as div and curl⁻¹ between them. The
// equations are H. Nicholson, ACSE Report 26, Sheffield 1974, pp. 2–6
// ("N74"): the orthogonal network E + e = Z(i + I) per q-network (eq. 3),
// its mesh and node solutions i = L(e − ZI), E = M(I − Ye) with
// L = C_c(C_cᵀZC_c)⁻¹C_cᵀ, M = A°((A°)ᵀYA°)⁻¹(A°)ᵀ = Z − ZLZ (p. 4), and
// the incidence relations M^{q−1}_q·M^q_{q+1} = 0, (C_c^{(q)})ᵀA°_{(q)} = 0
// (p. 6). Transcriptions in doc/kron-figure-catalogue.md §6–§7.
//
// Construction. A Complex is a chain complex: cells of every grade
// q = 0..k and the boundary ∂_q from q-cells to (q−1)-cells, with
// ∂_{q−1}∂_q = 0. Kron's M^q_{q+1} is the matrix of ∂_{q+1} (rows q-cells,
// columns (q+1)-cells), his "div" its transpose and "curl⁻¹" its
// pseudo-inverse (N74 p. 6). A Polyhedron puts a q-network on every
// grade: a positive weight W_q per cell, the cell's admittance (energy
// ½W_q u²), which for a cubical complex in a uniform medium is
// ε·h or h/μ with h = (dual-cell measure)/(cell measure) the diagonal
// Hodge star (E. Tonti, PIER 32, 1–44, 2001), ε on odd q and μ on even q
// for the primal polyhedron and the reverse for the dual (K62 p. 54). At
// k = 3 the grade-1 weight is Kron's 1944 dielectric capacitance
// εΔyΔz/Δx and the grade-2 weight the inverse of his magnetic inductance
// μΔxΔy/Δz (K62 Fig. 1; YeeGrid.Constitutive).
//
// The waves. With δ_q = ∂_{q+1}ᵀ (the coboundary) and its adjoint in the
// weights δ*_q = W_q⁻¹δ_qᵀW_{q+1}, the field u = (u_0 … u_k) obeys
//
//	du/dt = (δ − δ*)·u − W⁻¹(G·u + J)
//
// whose operator is skew in the energy ½Σ u_qᵀW_qu_q, so a lossless
// polyhedron keeps its energy and each neighbouring pair (q, q+1) carries
// one wave at 1/√(εμ): "one complete electromagnetic wave spans two
// neighbouring dimensions" (K62 p. 54). On grades 1–2 of a 3-cube this is
// Maxwell's equations in Yee's form (u_1 the edge voltages E·Δl, u_2 the
// face fluxes −B·ΔA); on all four grades it adds the two scalars that the
// standard theory removes, the electric one on the nodes and the magnetic
// one on the cubes (em/scalar, H47, H2), which stay zero while charge is
// conserved. Grades above 3 have no physical counterpart and are the
// estimator Kron used them as, never a field (CLAUDE.md H20).
//
// The second-order operator on a grade, Δ_q = δ*_qδ_q + δ_{q−1}δ*_{q−1}
// (the Hodge Laplacian in the weights), gives the polyhedron's natural
// frequencies, ω² its eigenvalues; its zero modes are the static
// (longitudinal) fields, as many as the free cells of the grade below
// plus the Betti number.

// ---------------------------------------------------------------------
// Complex

// Incident is one entry of a boundary: the face Cell, with orientation
// Sign = ±1.
type Incident struct{ Cell, Sign int }

// Complex is a finite chain complex of dimension Dim: Cells[q] cells of
// grade q and, for q ≥ 1, the boundary of every q-cell as a list of
// (q−1)-cells with signs.
type Complex struct {
	Dim   int
	Cells []int
	bnd   [][][]Incident // bnd[q][c]: boundary of q-cell c, q ≥ 1
}

// NewComplex makes a complex from cell counts and boundaries
// (boundaries[q][c] for q = 1..Dim) and checks it.
func NewComplex(cells []int, boundaries [][][]Incident) (*Complex, error) {
	k := len(cells) - 1
	if k < 0 || len(boundaries) != len(cells) {
		return nil, fmt.Errorf("%w: %d grades, %d boundary lists", ErrPolyhedron, len(cells), len(boundaries))
	}
	c := &Complex{Dim: k, Cells: append([]int(nil), cells...), bnd: boundaries}
	if err := c.Check(); err != nil {
		return nil, err
	}
	return c, nil
}

// Faces returns the boundary of q-cell c.
func (c *Complex) Faces(q, cell int) []Incident {
	if q < 1 || q > c.Dim {
		return nil
	}
	return c.bnd[q][cell]
}

// Check verifies the incidences and ∂_{q−1}∂_q = 0 exactly, which is Kron's
// M^{q−1}_q·M^q_{q+1} = 0 (N74 p. 6) and, at q = 2, Aᵀ·C = 0.
func (c *Complex) Check() error {
	if len(c.bnd) != c.Dim+1 {
		return fmt.Errorf("%w: %d boundary grades for dimension %d", ErrPolyhedron, len(c.bnd), c.Dim)
	}
	for q := 1; q <= c.Dim; q++ {
		if len(c.bnd[q]) != c.Cells[q] {
			return fmt.Errorf("%w: grade %d has %d cells and %d boundaries", ErrPolyhedron, q, c.Cells[q], len(c.bnd[q]))
		}
		for cell, fs := range c.bnd[q] {
			for _, f := range fs {
				if f.Cell < 0 || f.Cell >= c.Cells[q-1] || (f.Sign != 1 && f.Sign != -1) {
					return fmt.Errorf("%w: %d-cell %d has face %v", ErrPolyhedron, q, cell, f)
				}
			}
		}
	}
	for q := 2; q <= c.Dim; q++ {
		acc := map[int]int{}
		for cell, fs := range c.bnd[q] {
			clear(acc)
			for _, f := range fs {
				for _, g := range c.bnd[q-1][f.Cell] {
					acc[g.Cell] += f.Sign * g.Sign
				}
			}
			for g, v := range acc {
				if v != 0 {
					return fmt.Errorf("%w: ∂∂ ≠ 0 on %d-cell %d at %d-cell %d (%d)", ErrPolyhedron, q, cell, q-2, g, v)
				}
			}
		}
	}
	return nil
}

// Boundary applies ∂_q to a q-chain: a (q−1)-chain.
func (c *Complex) Boundary(q int, x []float64) []float64 {
	if q < 1 || q > c.Dim {
		return nil
	}
	out := make([]float64, c.Cells[q-1])
	for cell, fs := range c.bnd[q] {
		for _, f := range fs {
			out[f.Cell] += float64(f.Sign) * x[cell]
		}
	}
	return out
}

// Coboundary applies δ_q = ∂_{q+1}ᵀ to a q-cochain: a (q+1)-cochain (the
// generalised gradient, curl, … of B68 Fig. 8; Kron's e′ = Cᵀ·e at q = 1).
func (c *Complex) Coboundary(q int, x []float64) []float64 {
	if q < 0 || q >= c.Dim {
		return nil
	}
	out := make([]float64, c.Cells[q+1])
	for cell, fs := range c.bnd[q+1] {
		s := 0.0
		for _, f := range fs {
			s += float64(f.Sign) * x[f.Cell]
		}
		out[cell] = s
	}
	return out
}

// Incidence returns Kron's M^q_{q+1} as a dense matrix: rows q-cells,
// columns (q+1)-cells, the matrix of ∂_{q+1} (N74 Fig. 2, B68 Figs. 8–9).
func (c *Complex) Incidence(q int) ([][]float64, error) {
	if q < 0 || q >= c.Dim {
		return nil, fmt.Errorf("%w: no M^%d_%d in dimension %d", ErrPolyhedron, q, q+1, c.Dim)
	}
	m := make([][]float64, c.Cells[q])
	for i := range m {
		m[i] = make([]float64, c.Cells[q+1])
	}
	for cell, fs := range c.bnd[q+1] {
		for _, f := range fs {
			m[f.Cell][cell] += float64(f.Sign)
		}
	}
	return m, nil
}

// Euler is the Euler characteristic Σ(−1)^q·Cells[q].
func (c *Complex) Euler() int {
	e := 0
	for q, n := range c.Cells {
		if q%2 == 0 {
			e += n
		} else {
			e -= n
		}
	}
	return e
}

// Betti returns the Betti numbers b_q = Cells[q] − rank ∂_q − rank ∂_{q+1}
// (over the reals, by elimination); Σ(−1)^q·b_q = Euler(). Kron's
// p = N − S and m = b − N + S (K39 p. 75) are the q = 1 case.
func (c *Complex) Betti() []int {
	rank := make([]int, c.Dim+2)
	for q := 1; q <= c.Dim; q++ {
		m, _ := c.Incidence(q - 1)
		rank[q] = matRank(m)
	}
	b := make([]int, c.Dim+1)
	for q := range b {
		b[q] = c.Cells[q] - rank[q] - rank[q+1]
	}
	return b
}

// Simplicial builds the complex generated by the given simplices (each a
// list of vertex indices < vertices): every face of every simplex is a
// cell, a cell's orientation is its sorted vertex order and
// ∂[v₀…v_q] = Σ(−1)^i[v₀…v̂_i…v_q]. Cells of each grade are numbered in
// lexicographic order of their vertex lists; Vertices returns them.
func Simplicial(vertices int, simplices [][]int) (*SimplicialComplex, error) {
	if vertices <= 0 {
		return nil, fmt.Errorf("%w: %d vertices", ErrPolyhedron, vertices)
	}
	dim := 0
	set := map[string][]int{}
	var add func(s []int)
	add = func(s []int) {
		key := fmt.Sprint(s)
		if _, ok := set[key]; ok {
			return
		}
		set[key] = s
		if len(s) > 1 {
			for i := range s {
				f := append(append([]int(nil), s[:i]...), s[i+1:]...)
				add(f)
			}
		}
	}
	for _, s := range simplices {
		t := append([]int(nil), s...)
		sort.Ints(t)
		for i, v := range t {
			if v < 0 || v >= vertices || (i > 0 && v == t[i-1]) {
				return nil, fmt.Errorf("%w: simplex %v", ErrPolyhedron, s)
			}
		}
		if len(t) == 0 {
			continue
		}
		dim = max(dim, len(t)-1)
		add(t)
	}
	for v := 0; v < vertices; v++ {
		add([]int{v})
	}
	verts := make([][][]int, dim+1)
	for _, s := range set {
		verts[len(s)-1] = append(verts[len(s)-1], s)
	}
	index := make([]map[string]int, dim+1)
	for q := range verts {
		sort.Slice(verts[q], func(a, b int) bool { return lexLess(verts[q][a], verts[q][b]) })
		index[q] = map[string]int{}
		for i, s := range verts[q] {
			index[q][fmt.Sprint(s)] = i
		}
	}
	cells := make([]int, dim+1)
	bnd := make([][][]Incident, dim+1)
	for q := range verts {
		cells[q] = len(verts[q])
		if q == 0 {
			continue
		}
		bnd[q] = make([][]Incident, cells[q])
		for c, s := range verts[q] {
			for i := range s {
				f := append(append([]int(nil), s[:i]...), s[i+1:]...)
				sign := 1
				if i%2 == 1 {
					sign = -1
				}
				bnd[q][c] = append(bnd[q][c], Incident{index[q-1][fmt.Sprint(f)], sign})
			}
		}
	}
	cx, err := NewComplex(cells, bnd)
	if err != nil {
		return nil, err
	}
	return &SimplicialComplex{Complex: cx, verts: verts}, nil
}

// SimplicialComplex is a Complex built by Simplicial.
type SimplicialComplex struct {
	*Complex
	verts [][][]int
}

// Vertices returns the sorted vertex list of q-cell c.
func (s *SimplicialComplex) Vertices(q, c int) []int { return s.verts[q][c] }

func lexLess(a, b []int) bool {
	for i := range a {
		if i >= len(b) {
			return false
		}
		if a[i] != b[i] {
			return a[i] < b[i]
		}
	}
	return len(a) < len(b)
}

// ---------------------------------------------------------------------
// Cubical complex

// Cubical is the cubical complex of a rectangular grid of N[0]×…×N[k−1]
// vertices with spacing D in k dimensions: a q-cell is a vertex p and a
// set S of q axes with p_a < N_a − 1 for a ∈ S, spanning the box from p
// to p + Σ_{a∈S} e_a, oriented by the increasing order of its axes. Its
// boundary is ∂[p; S] = Σ_j (−1)^j ([p + e_{s_j}; S∖s_j] − [p; S∖s_j]),
// s_0 < s_1 < … the axes of S. At k = 3 the 1- and 2-cells are the edges
// and faces of YeeGrid, the faces with normal y carrying the opposite
// orientation (z∧x there, x∧z here).
type Cubical struct {
	*Complex
	N []int
	D []float64

	stride []int
	base   [][]int  // base[q][c]: vertex of q-cell c
	mask   [][]uint // mask[q][c]: axis set of q-cell c
	index  [][]int  // index[mask][vertex]: cell number in its grade, −1 if none
}

// NewCubical builds the grid (k = len(n) ≥ 1, every N ≥ 2, D > 0, k ≤ 16).
func NewCubical(n []int, d []float64) (*Cubical, error) {
	k := len(n)
	if k < 1 || k > 16 || len(d) != k {
		return nil, fmt.Errorf("%w: cubical grid %v spacing %v", ErrPolyhedron, n, d)
	}
	nv := 1
	stride := make([]int, k)
	for a := 0; a < k; a++ {
		if n[a] < 2 || !(d[a] > 0) {
			return nil, fmt.Errorf("%w: cubical grid %v spacing %v", ErrPolyhedron, n, d)
		}
		stride[a] = nv
		nv *= n[a]
	}
	g := &Cubical{N: append([]int(nil), n...), D: append([]float64(nil), d...), stride: stride,
		base: make([][]int, k+1), mask: make([][]uint, k+1), index: make([][]int, 1<<k)}
	coord := make([]int, k)
	for m := uint(0); m < 1<<k; m++ {
		q := bits.OnesCount(m)
		g.index[m] = make([]int, nv)
		for v := 0; v < nv; v++ {
			g.index[m][v] = -1
			g.coords(v, coord)
			ok := true
			for a := 0; a < k; a++ {
				if m&(1<<a) != 0 && coord[a] == n[a]-1 {
					ok = false
					break
				}
			}
			if ok {
				g.index[m][v] = len(g.base[q])
				g.base[q] = append(g.base[q], v)
				g.mask[q] = append(g.mask[q], m)
			}
		}
	}
	cells := make([]int, k+1)
	bnd := make([][][]Incident, k+1)
	for q := 0; q <= k; q++ {
		cells[q] = len(g.base[q])
		if q == 0 {
			continue
		}
		bnd[q] = make([][]Incident, cells[q])
		for c := range g.base[q] {
			v, m := g.base[q][c], g.mask[q][c]
			j := 0
			for a := 0; a < k; a++ {
				if m&(1<<a) == 0 {
					continue
				}
				fm := m &^ (1 << a)
				sign := 1
				if j%2 == 1 {
					sign = -1
				}
				bnd[q][c] = append(bnd[q][c],
					Incident{g.index[fm][v+stride[a]], sign},
					Incident{g.index[fm][v], -sign})
				j++
			}
		}
	}
	cx, err := NewComplex(cells, bnd)
	if err != nil {
		return nil, err
	}
	g.Complex = cx
	return g, nil
}

func (g *Cubical) coords(v int, out []int) {
	for a := range g.N {
		out[a] = v % g.N[a]
		v /= g.N[a]
	}
}

// Vertex is the vertex number of grid point p.
func (g *Cubical) Vertex(p []int) int {
	v := 0
	for a := range g.N {
		v += p[a] * g.stride[a]
	}
	return v
}

// Cell returns the number of the cell with lowest corner p and axis set
// axes (any order), its grade, and −1 if there is none.
func (g *Cubical) Cell(p []int, axes ...int) (cell, q int) {
	var m uint
	for _, a := range axes {
		if a < 0 || a >= len(g.N) {
			return -1, len(axes)
		}
		m |= 1 << a
	}
	return g.index[m][g.Vertex(p)], bits.OnesCount(m)
}

// Corner returns the lowest corner and the axis set (increasing) of
// q-cell c.
func (g *Cubical) Corner(q, c int) (p []int, axes []int) {
	p = make([]int, len(g.N))
	g.coords(g.base[q][c], p)
	for a := range g.N {
		if g.mask[q][c]&(1<<a) != 0 {
			axes = append(axes, a)
		}
	}
	return p, axes
}

// Centre is the midpoint of q-cell c [m], origin at vertex 0.
func (g *Cubical) Centre(q, c int) []float64 {
	p, _ := g.Corner(q, c)
	x := make([]float64, len(p))
	for a := range p {
		x[a] = float64(p[a]) * g.D[a]
		if g.mask[q][c]&(1<<a) != 0 {
			x[a] += g.D[a] / 2
		}
	}
	return x
}

// Hodge is the diagonal Hodge star of q-cell c: the measure of its dual
// (k−q)-cell over its own, Π_{a∉S}Δ_a / Π_{a∈S}Δ_a (full cells, as in
// YeeGrid.Constitutive; the half cells at a wall carry no free field).
func (g *Cubical) Hodge(q, c int) float64 {
	h := 1.0
	for a := range g.N {
		if g.mask[q][c]&(1<<a) != 0 {
			h /= g.D[a]
		} else {
			h *= g.D[a]
		}
	}
	return h
}

// OnWall reports whether q-cell c lies in the grid's boundary: some axis
// a not in the cell has p_a = 0 or N_a − 1.
func (g *Cubical) OnWall(q, c int) bool {
	v, m := g.base[q][c], g.mask[q][c]
	for a := range g.N {
		if m&(1<<a) != 0 {
			continue
		}
		pa := (v / g.stride[a]) % g.N[a]
		if pa == 0 || pa == g.N[a]-1 {
			return true
		}
	}
	return false
}

// ---------------------------------------------------------------------
// Polyhedron

// Polyhedron is Kron's polyhedron on a complex: a q-network on every grade
// with weight W[q][c] > 0 (the cell's admittance, energy ½W·u²), optional
// conductance G[q][c] ≥ 0, and cells held at zero (Fixed, a conducting
// wall). Only the grades Lo..Hi take part; the rest are absent (held at
// zero), so {1, 2} on a 3-cube is Maxwell's network and {0, …, k} Kron's
// full polyhedron.
type Polyhedron struct {
	C      *Complex
	W      [][]float64
	G      [][]float64
	Fixed  [][]bool
	Lo, Hi int
}

// KronPolyhedron returns the polyhedron of a cubical grid in a uniform
// medium with Kron's material rule: W = ε·h on odd grades and h/μ on even
// grades (primal), or the reverse (dual, K62 p. 54; B68 Fig. 9). Every
// grade is active and nothing is fixed.
func KronPolyhedron(g *Cubical, eps, mu float64, dual bool) (*Polyhedron, error) {
	if !(eps > 0) || !(mu > 0) {
		return nil, fmt.Errorf("%w: ε %g μ %g", ErrPolyhedron, eps, mu)
	}
	p := &Polyhedron{C: g.Complex, W: make([][]float64, g.Dim+1), G: make([][]float64, g.Dim+1),
		Fixed: make([][]bool, g.Dim+1), Lo: 0, Hi: g.Dim}
	for q := 0; q <= g.Dim; q++ {
		p.W[q] = make([]float64, g.Cells[q])
		p.G[q] = make([]float64, g.Cells[q])
		p.Fixed[q] = make([]bool, g.Cells[q])
		dielectric := q%2 == 1
		if dual {
			dielectric = !dielectric
		}
		for c := range p.W[q] {
			if dielectric {
				p.W[q][c] = eps * g.Hodge(q, c)
			} else {
				p.W[q][c] = g.Hodge(q, c) / mu
			}
		}
	}
	return p, nil
}

// Wall fixes every cell of the given grades that lies in the grid's
// boundary (a perfect electric wall on grades 0 and 1 of the primal
// 3-polyhedron: no tangential E, no potential).
func (p *Polyhedron) Wall(g *Cubical, grades ...int) {
	for _, q := range grades {
		for c := range p.Fixed[q] {
			if g.OnWall(q, c) {
				p.Fixed[q][c] = true
			}
		}
	}
}

// Check verifies the data.
func (p *Polyhedron) Check() error {
	k := p.C.Dim
	if p.Lo < 0 || p.Hi > k || p.Lo > p.Hi || len(p.W) != k+1 || len(p.Fixed) != k+1 {
		return fmt.Errorf("%w: grades %d..%d of a %d-complex", ErrPolyhedron, p.Lo, p.Hi, k)
	}
	for q := p.Lo; q <= p.Hi; q++ {
		if len(p.W[q]) != p.C.Cells[q] || len(p.Fixed[q]) != p.C.Cells[q] ||
			(p.G != nil && p.G[q] != nil && len(p.G[q]) != p.C.Cells[q]) {
			return fmt.Errorf("%w: grade %d data lengths", ErrPolyhedron, q)
		}
		for c, w := range p.W[q] {
			if !p.Fixed[q][c] && !(w > 0) {
				return fmt.Errorf("%w: %d-cell %d has W = %g", ErrPolyhedron, q, c, w)
			}
			if g := p.g(q, c); g < 0 || math.IsNaN(g) {
				return fmt.Errorf("%w: %d-cell %d has G = %g", ErrPolyhedron, q, c, g)
			}
		}
	}
	return nil
}

func (p *Polyhedron) g(q, c int) float64 {
	if p.G == nil || p.G[q] == nil {
		return 0
	}
	return p.G[q][c]
}

func (p *Polyhedron) active(q int) bool { return q >= p.Lo && q <= p.Hi }

func (p *Polyhedron) free(q, c int) bool { return p.active(q) && !p.Fixed[q][c] }

// Couple returns W_q·(Ku)_q = W_q·δ_{q−1}u_{q−1} − δ_qᵀ(W_{q+1}u_{q+1}) for
// grade q, zero on fixed cells and outside the active grades. Summed over
// grades, u_qᵀ·Couple(q) vanishes identically: the coupling between
// q-networks does no work, Kron's power invariant (K39 p. 409).
func (p *Polyhedron) Couple(u [][]float64, q int) []float64 {
	out := make([]float64, p.C.Cells[q])
	if !p.active(q) {
		return out
	}
	if q-1 >= p.Lo {
		d := p.C.Coboundary(q-1, u[q-1])
		for c := range out {
			out[c] = p.W[q][c] * d[c]
		}
	}
	if q+1 <= p.Hi {
		wu := make([]float64, p.C.Cells[q+1])
		for c := range wu {
			if p.free(q+1, c) {
				wu[c] = p.W[q+1][c] * u[q+1][c]
			}
		}
		d := p.C.Boundary(q+1, wu)
		for c := range out {
			out[c] -= d[c]
		}
	}
	for c := range out {
		if p.Fixed[q][c] {
			out[c] = 0
		}
	}
	return out
}

// Energy is ½Σ_q u_qᵀW_qu_q over the free cells of the active grades.
func (p *Polyhedron) Energy(u [][]float64) float64 {
	e := 0.0
	for q := p.Lo; q <= p.Hi; q++ {
		for c, w := range p.W[q] {
			if !p.Fixed[q][c] {
				e += 0.5 * w * u[q][c] * u[q][c]
			}
		}
	}
	return e
}

// Zero returns a zero field.
func (p *Polyhedron) Zero() [][]float64 {
	u := make([][]float64, p.C.Dim+1)
	for q := range u {
		u[q] = make([]float64, p.C.Cells[q])
	}
	return u
}

// Laplacian returns the second-order operator of grade q in the weights,
// Δ_q = δ*_qδ_q + δ_{q−1}δ*_{q−1} on its free cells, as the symmetric
// matrix W^{½}Δ_qW^{−½} (terms of inactive neighbouring grades dropped),
// with the free cells it is written on. Its eigenvalues are ω² of the
// polyhedron's natural modes on that grade.
func (p *Polyhedron) Laplacian(q int) ([][]float64, []int, error) {
	if err := p.Check(); err != nil {
		return nil, nil, err
	}
	if !p.active(q) {
		return nil, nil, fmt.Errorf("%w: grade %d is not active", ErrPolyhedron, q)
	}
	var cells []int
	pos := make([]int, p.C.Cells[q])
	for c := range pos {
		pos[c] = -1
		if !p.Fixed[q][c] {
			pos[c] = len(cells)
			cells = append(cells, c)
		}
	}
	n := len(cells)
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
	}
	// δ*_qδ_q: Σ over free (q+1)-cells τ of W_τ·(∂τ)_i(∂τ)_j/√(W_iW_j).
	if q+1 <= p.Hi {
		for t, fs := range p.C.bnd[q+1] {
			if p.Fixed[q+1][t] {
				continue
			}
			for _, f := range fs {
				i := pos[f.Cell]
				if i < 0 {
					continue
				}
				for _, h := range fs {
					j := pos[h.Cell]
					if j < 0 {
						continue
					}
					a[i][j] += p.W[q+1][t] * float64(f.Sign*h.Sign) /
						math.Sqrt(p.W[q][f.Cell]*p.W[q][h.Cell])
				}
			}
		}
	}
	// δ_{q−1}δ*_{q−1}: Σ over free (q−1)-cells σ of √(W_iW_j)·(δσ)_i(δσ)_j/W_σ.
	if q-1 >= p.Lo {
		co := make([][]Incident, p.C.Cells[q-1])
		for c, fs := range p.C.bnd[q] {
			if pos[c] < 0 {
				continue
			}
			for _, f := range fs {
				co[f.Cell] = append(co[f.Cell], Incident{c, f.Sign})
			}
		}
		for s, cs := range co {
			if p.Fixed[q-1][s] {
				continue
			}
			for _, f := range cs {
				for _, h := range cs {
					a[pos[f.Cell]][pos[h.Cell]] += math.Sqrt(p.W[q][f.Cell]*p.W[q][h.Cell]) *
						float64(f.Sign*h.Sign) / p.W[q-1][s]
				}
			}
		}
	}
	return a, cells, nil
}

// Modes returns the natural angular frequencies ω = √λ of grade q in
// increasing order (λ the eigenvalues of Laplacian(q), negative round-off
// set to 0).
func (p *Polyhedron) Modes(q int) ([]float64, error) {
	a, _, err := p.Laplacian(q)
	if err != nil {
		return nil, err
	}
	lam := symEigenvalues(a)
	w := make([]float64, len(lam))
	for i, l := range lam {
		w[i] = math.Sqrt(math.Max(l, 0))
	}
	sort.Float64s(w)
	return w, nil
}

// ---------------------------------------------------------------------
// Field (time domain)

// PolyField steps the polyhedron's field by the leapfrog: grades of parity
// IntParity live at the whole steps t = n·Dt, the others at the half
// steps, and each is updated from its neighbours, which have the other
// parity. On grades 1–2 of a 3-cube with IntParity 1 it is Yee's scheme
// (U1Field, em/fdtd.go). Losses are averaged over the update (as in
// em/fdtd.go) and sources J[q] (a function of grade and time, W·du/dt
// gains −J) enter at the update's midpoint, so the ledger
//
//	Energy(n) − Energy(0) = Work − Heat
//
// holds to round-off, with Energy the leapfrog's exact discrete energy
// (the half-step grades as the product of their values either side of the
// whole step, as U1Field.Energy).
type PolyField struct {
	P         *Polyhedron
	Dt        float64
	IntParity int
	// Source returns the impressed current on grade q at time t (nil: none).
	Source func(q int, t float64) []float64
	U      [][]float64 // whole-step grades at Time(); half-step grades at Time() + Dt/2
	Steps  int
	Work   float64 // energy delivered by the sources [J]
	Heat   float64 // energy dissipated in G [J]

	prev  [][]float64 // half-step grades at Time() − Dt/2
	tPrev [][]float64 // G·ū + J of the half-step update that made U from prev
}

// NewPolyField makes a field at rest.
func NewPolyField(p *Polyhedron, dt float64, intParity int) (*PolyField, error) {
	if err := p.Check(); err != nil {
		return nil, err
	}
	if !(dt > 0) || (intParity != 0 && intParity != 1) {
		return nil, fmt.Errorf("%w: dt %g parity %d", ErrPolyhedron, dt, intParity)
	}
	f := &PolyField{P: p, Dt: dt, IntParity: intParity, U: p.Zero(), prev: p.Zero(), tPrev: p.Zero()}
	return f, nil
}

// Time is the time of the whole-step grades.
func (f *PolyField) Time() float64 { return float64(f.Steps) * f.Dt }

func (f *PolyField) whole(q int) bool { return q%2 == f.IntParity }

// Set sets the field at Time(): the whole-step grades to u, and the
// half-step grades at Time() ± Dt/2 consistent with u lossless and
// source-free (half-step grades taken as u at Time() + Dt/2).
func (f *PolyField) Set(u [][]float64) {
	p := f.P
	for q := range f.U {
		copy(f.U[q], u[q])
		for c := range f.U[q] {
			if !p.free(q, c) {
				f.U[q][c] = 0
			}
		}
		clear(f.tPrev[q])
	}
	for q := p.Lo; q <= p.Hi; q++ {
		if f.whole(q) {
			continue
		}
		r := p.Couple(f.U, q)
		for c := range r {
			f.prev[q][c] = f.U[q][c]
			if p.free(q, c) {
				f.prev[q][c] -= f.Dt * r[c] / p.W[q][c]
			}
		}
	}
}

// Step advances one Dt: the whole-step grades from the half-step ones,
// then the half-step grades from the new whole-step ones.
func (f *PolyField) Step() {
	p, dt := f.P, f.Dt
	tMid := f.Time() + dt/2
	tNext := f.Time() + dt
	// Whole-step grades: W(u⁺ − u⁻)/dt = R − G(u⁺ + u⁻)/2 − J(t + dt/2).
	rs := make([][]float64, len(f.U))
	for q := p.Lo; q <= p.Hi; q++ {
		if f.whole(q) {
			rs[q] = p.Couple(f.U, q)
		}
	}
	for q := p.Lo; q <= p.Hi; q++ {
		if !f.whole(q) {
			continue
		}
		var j []float64
		if f.Source != nil {
			j = f.Source(q, tMid)
		}
		for c := range f.U[q] {
			if !p.free(q, c) {
				continue
			}
			w, g := p.W[q][c], p.g(q, c)
			jc := 0.0
			if j != nil {
				jc = j[c]
			}
			old := f.U[q][c]
			nu := ((w/dt-g/2)*old + rs[q][c] - jc) / (w/dt + g/2)
			bar := (old + nu) / 2
			f.Heat += dt * g * bar * bar
			f.Work -= dt * bar * jc
			f.U[q][c] = nu
		}
	}
	// Half-step grades, at time t + dt: T = G·ū + J(t + dt).
	for q := p.Lo; q <= p.Hi; q++ {
		if f.whole(q) {
			rs[q] = nil
			continue
		}
		rs[q] = p.Couple(f.U, q)
	}
	for q := p.Lo; q <= p.Hi; q++ {
		if f.whole(q) {
			continue
		}
		var j []float64
		if f.Source != nil {
			j = f.Source(q, tNext)
		}
		for c := range f.U[q] {
			if !p.free(q, c) {
				continue
			}
			w, g := p.W[q][c], p.g(q, c)
			jc := 0.0
			if j != nil {
				jc = j[c]
			}
			mid := f.U[q][c]
			nu := ((w/dt-g/2)*mid + rs[q][c] - jc) / (w/dt + g/2)
			tn := g*(mid+nu)/2 + jc
			// Ledger of this step's half-step part: −dt/2·u^{n+½}·(T_n + T_{n+1}),
			// split into heat (G terms) and work (J terms).
			gPrev := f.tPrev[q][c] - 0 // T_n holds G·ū + J of the previous update
			_ = gPrev
			f.ledgerHalf(q, c, mid, nu, jc)
			f.prev[q][c] = mid
			f.U[q][c] = nu
			f.tPrev[q][c] = tn
		}
	}
	f.Steps++
}

// ledgerHalf books −dt/2·u^{n+½}·(T_n + T_{n+1}) for one half-step cell,
// with T = G·ū + J split into heat and work. T_n's two parts are rebuilt
// from the stored previous values.
func (f *PolyField) ledgerHalf(q, c int, mid, next, jNext float64) {
	p, dt := f.P, f.Dt
	g := p.g(q, c)
	gPrevBar := g * (f.prev[q][c] + mid) / 2 // G·ū_n (prev is u^{n−½} until replaced)
	jPrev := f.tPrev[q][c] - gPrevBar        // J_n
	gNextBar := g * (mid + next) / 2
	f.Heat += dt / 2 * mid * (gPrevBar + gNextBar)
	f.Work -= dt / 2 * mid * (jPrev + jNext)
}

// Energy is the leapfrog's discrete energy at Time(): ½uᵀWu on the
// whole-step grades and ½u⁻ᵀWu⁺ on the half-step grades (u∓ at
// Time() ∓ Dt/2).
func (f *PolyField) Energy() float64 {
	p := f.P
	e := 0.0
	for q := p.Lo; q <= p.Hi; q++ {
		for c, w := range p.W[q] {
			if !p.free(q, c) {
				continue
			}
			if f.whole(q) {
				e += 0.5 * w * f.U[q][c] * f.U[q][c]
			} else {
				e += 0.5 * w * f.prev[q][c] * f.U[q][c]
			}
		}
	}
	return e
}

// CourantStep is the leapfrog's stability limit for a cubical grid in a
// uniform medium, Dt < 1/(v·√Σ_a Δ_a⁻²) with v = 1/√(εμ): the largest
// natural frequency of every grade is 2v·√ΣΔ_a⁻² (that of the Yee cell
// at k = 3).
func CourantStep(g *Cubical, eps, mu float64) float64 {
	s := 0.0
	for _, d := range g.D {
		s += 1 / (d * d)
	}
	return math.Sqrt(eps*mu) / math.Sqrt(s)
}

// ---------------------------------------------------------------------
// q-networks: Nicholson's orthogonal solution

// QNetwork is the q-network of a polyhedron as one of Kron's orthogonal
// networks (N74 eqs. 1–3): branches are the free q-cells, junctions the
// free (q−1)-cells, primitive impedance Z = 1/W. Kirchhoff's current law
// is ∂_q i = 0 at the free junctions; the closed paths C_c (columns, a
// basis of the cycles) and the open paths A° (columns, a basis of the
// coboundaries, the row space of ∂_q) together span the branch space, and
// C_cᵀA° = 0 (N74 p. 6).
type QNetwork struct {
	Q        int
	Branches []int       // q-cell of each branch
	Z        []float64   // primitive impedance per branch
	Cc       [][]float64 // b × m, closed paths
	Ao       [][]float64 // b × p, open paths
}

// QNetwork returns the q-network of grade q.
func (p *Polyhedron) QNetwork(q int) (*QNetwork, error) {
	if err := p.Check(); err != nil {
		return nil, err
	}
	if q < 0 || q > p.C.Dim {
		return nil, fmt.Errorf("%w: no grade %d", ErrPolyhedron, q)
	}
	n := &QNetwork{Q: q}
	pos := make([]int, p.C.Cells[q])
	for c := range pos {
		pos[c] = -1
		if !p.Fixed[q][c] {
			pos[c] = len(n.Branches)
			n.Branches = append(n.Branches, c)
			n.Z = append(n.Z, 1/p.W[q][c])
		}
	}
	b := len(n.Branches)
	// ∂_q restricted: rows free (q−1)-cells, columns branches.
	var d [][]float64
	if q >= 1 {
		row := make([]int, p.C.Cells[q-1])
		for s := range row {
			row[s] = -1
			if !p.Fixed[q-1][s] {
				row[s] = len(d)
				d = append(d, make([]float64, b))
			}
		}
		for c, fs := range p.C.bnd[q] {
			if pos[c] < 0 {
				continue
			}
			for _, f := range fs {
				if r := row[f.Cell]; r >= 0 {
					d[r][pos[c]] += float64(f.Sign)
				}
			}
		}
	}
	n.Cc = nullSpace(d, b)
	n.Ao = rowSpace(d, b)
	if len(n.Cc[0])+len(n.Ao[0]) != b && b > 0 {
		return nil, fmt.Errorf("%w: q-network %d: %d closed + %d open paths for %d branches",
			ErrPolyhedron, q, len(n.Cc[0]), len(n.Ao[0]), b)
	}
	return n, nil
}

// Counts returns the numbers of branches b, closed paths m and open paths
// p (m + p = b; K39 p. 75 at q = 1).
func (n *QNetwork) Counts() (b, m, p int) {
	b = len(n.Branches)
	if b == 0 {
		return 0, 0, 0
	}
	return b, len(n.Cc[0]), len(n.Ao[0])
}

// MeshOperator is N74's L = C_c(C_cᵀZC_c)⁻¹C_cᵀ (p. 4), b × b.
func (n *QNetwork) MeshOperator() ([][]float64, error) {
	return projector(n.Cc, n.Z)
}

// NodeOperator is N74's M = A°((A°)ᵀYA°)⁻¹(A°)ᵀ with Y = Z⁻¹ (p. 4).
func (n *QNetwork) NodeOperator() ([][]float64, error) {
	y := make([]float64, len(n.Z))
	for i, z := range n.Z {
		y[i] = 1 / z
	}
	return projector(n.Ao, y)
}

// Solve returns the orthogonal network's response to impressed branch
// voltages e and impressed currents I: the mesh current i = L(e − ZI) and
// the junction-pair voltage E = M(I − Ye) (N74 p. 4), which satisfy
// E + e = Z(i + I), ∂i = 0, C_cᵀE = 0 and Eᵀi = 0.
func (n *QNetwork) Solve(e, I []float64) (i, E []float64, err error) {
	b := len(n.Branches)
	if len(e) != b || len(I) != b {
		return nil, nil, fmt.Errorf("%w: %d/%d sources for %d branches", ErrPolyhedron, len(e), len(I), b)
	}
	l, err := n.MeshOperator()
	if err != nil {
		return nil, nil, err
	}
	m, err := n.NodeOperator()
	if err != nil {
		return nil, nil, err
	}
	u, w := make([]float64, b), make([]float64, b)
	for k := range u {
		u[k] = e[k] - n.Z[k]*I[k]
		w[k] = I[k] - e[k]/n.Z[k]
	}
	return matVec(l, u), matVec(m, w), nil
}

// projector returns B(BᵀDB)⁻¹Bᵀ for diagonal D.
func projector(b [][]float64, d []float64) ([][]float64, error) {
	n := len(b)
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, n)
	}
	if n == 0 || len(b[0]) == 0 {
		return out, nil
	}
	m := len(b[0])
	g := make([][]float64, m)
	for r := range g {
		g[r] = make([]float64, m)
		for c := range g[r] {
			s := 0.0
			for k := 0; k < n; k++ {
				s += b[k][r] * d[k] * b[k][c]
			}
			g[r][c] = s
		}
	}
	// X = (BᵀDB)⁻¹Bᵀ, column by column of Bᵀ.
	bt := make([][]float64, m)
	for r := range bt {
		bt[r] = make([]float64, n)
		for k := 0; k < n; k++ {
			bt[r][k] = b[k][r]
		}
	}
	x, err := solveReal(g, bt)
	if err != nil {
		return nil, err
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			s := 0.0
			for r := 0; r < m; r++ {
				s += b[i][r] * x[r][j]
			}
			out[i][j] = s
		}
	}
	return out, nil
}

// ---------------------------------------------------------------------
// AC solve, direct and torn into its q-networks

// ACSolve returns the polyhedron's sinusoidal steady state at ω for the
// impressed currents J (per grade, phasors): (jωW + G − WK)u = −J on the
// free cells of the active grades, solved as one system.
func (p *Polyhedron) ACSolve(omega float64, j [][]complex128) ([][]complex128, error) {
	idx, n, err := p.acIndex()
	if err != nil {
		return nil, err
	}
	a := zeros(n, n)
	rhs := make([]complex128, n)
	for q := p.Lo; q <= p.Hi; q++ {
		for c, i := range idx[q] {
			if i < 0 {
				continue
			}
			a[i][i] = complex(p.g(q, c), omega*p.W[q][c])
			if j != nil && j[q] != nil {
				rhs[i] = -j[q][c]
			}
		}
	}
	p.acCouple(idx, func(r, c int, v float64) { a[r][c] += complex(v, 0) })
	x, err := solveSystem(a, rhs)
	if err != nil {
		return nil, err
	}
	return p.acScatter(idx, x), nil
}

// TornStats reports a torn AC solve: the pieces (one per q-network of the
// torn parity, each cell its own block because W and G are diagonal) and
// the order of the interconnection system (the free cells of the other
// parity's networks).
type TornStats struct {
	Pieces            int
	PieceCells        int
	Interconnection   int
	TornPower, Power  complex128
}

// ACTorn solves the same system by tearing the polyhedron into its
// q-networks (K62 pp. 54–55: "the polyhedron is torn into isolated
// q-networks"): the networks of parity tear are solved on their own,
// u_a = D_a⁻¹(s_a − X·u_b), and joined through the networks of the other
// parity by the interconnection system (D_b − Y·D_a⁻¹·X)u_b = s_b −
// Y·D_a⁻¹s_a (Kron's tearing as a change of frame, Q. Appl. Math. 13,
// 147, 1955; the block form of Aitchison 1983 with the coupling entries
// the incidence matrices). Its matrix is the second-order operator of the
// interconnecting grades, jωW + G + W·Δ·W/(jωW + G) in effect: Nicholson's
// "propagation from the 0-points up through the dimensions" (N74 p. 6).
// The complex power Σ ū·J of the torn solution, summed piece by piece and
// over the interconnection, equals that of the untorn one (Roth 1959 §5).
func (p *Polyhedron) ACTorn(omega float64, j [][]complex128, tear int) ([][]complex128, TornStats, error) {
	var st TornStats
	if tear != 0 && tear != 1 {
		return nil, st, fmt.Errorf("%w: tear parity %d", ErrPolyhedron, tear)
	}
	idx, n, err := p.acIndex()
	if err != nil {
		return nil, st, err
	}
	d := make([]complex128, n)
	s := make([]complex128, n)
	torn := make([]bool, n)
	for q := p.Lo; q <= p.Hi; q++ {
		if q%2 == tear {
			st.Pieces++
		}
		for c, i := range idx[q] {
			if i < 0 {
				continue
			}
			d[i] = complex(p.g(q, c), omega*p.W[q][c])
			if j != nil && j[q] != nil {
				s[i] = -j[q][c]
			}
			torn[i] = q%2 == tear
		}
	}
	// Coupling entries: X (torn rows ← interconnection columns) and Y.
	type ent struct {
		r, c int
		v float64
	}
	var xs, ys []ent
	p.acCouple(idx, func(r, c int, v float64) {
		if torn[r] {
			xs = append(xs, ent{r, c, v})
		} else {
			ys = append(ys, ent{r, c, v})
		}
	})
	// Number the interconnection unknowns.
	bi := make([]int, n)
	nb := 0
	for i := range bi {
		bi[i] = -1
		if torn[i] {
			st.PieceCells++
		} else {
			bi[i] = nb
			nb++
		}
	}
	st.Interconnection = nb
	if nb == 0 || st.PieceCells == 0 {
		return nil, st, fmt.Errorf("%w: nothing to tear (%d torn, %d joining cells)", ErrPolyhedron, st.PieceCells, nb)
	}
	// Pieces: D_a⁻¹ s_a and D_a⁻¹ X (sparse, per torn row).
	for i := range d {
		if torn[i] && d[i] == 0 {
			return nil, st, fmt.Errorf("%w: torn piece cell %d is singular at ω = %g", ErrSingular, i, omega)
		}
	}
	dax := make([][]ent, n) // per torn row: entries of D_a⁻¹X
	for _, e := range xs {
		dax[e.r] = append(dax[e.r], ent{e.r, e.c, e.v})
	}
	a := zeros(nb, nb)
	rhs := make([]complex128, nb)
	for i := range d {
		if !torn[i] {
			a[bi[i]][bi[i]] += d[i]
			rhs[bi[i]] += s[i]
		}
	}
	for _, e := range ys { // row e.r (interconnection) ← column e.c (torn)
		if !torn[e.c] {
			a[bi[e.r]][bi[e.c]] += complex(e.v, 0)
			continue
		}
		yv := complex(e.v, 0)
		rhs[bi[e.r]] -= yv * s[e.c] / d[e.c]
		for _, x := range dax[e.c] {
			a[bi[e.r]][bi[x.c]] -= yv * complex(x.v, 0) / d[e.c]
		}
	}
	ub, err := solveSystem(a, rhs)
	if err != nil {
		return nil, st, err
	}
	x := make([]complex128, n)
	for i := range x {
		if !torn[i] {
			x[i] = ub[bi[i]]
		}
	}
	for i := range x {
		if torn[i] {
			v := s[i]
			for _, e := range dax[i] {
				v -= complex(e.v, 0) * x[e.c]
			}
			x[i] = v / d[i]
		}
	}
	for i := range x {
		pw := -complexConj(x[i]) * s[i] // Σ ū·J with J = −s
		if torn[i] {
			st.TornPower += pw
		}
		st.Power += pw
	}
	return p.acScatter(idx, x), st, nil
}

func complexConj(z complex128) complex128 { return complex(real(z), -imag(z)) }

func (p *Polyhedron) acIndex() ([][]int, int, error) {
	if err := p.Check(); err != nil {
		return nil, 0, err
	}
	idx := make([][]int, p.C.Dim+1)
	n := 0
	for q := range idx {
		idx[q] = make([]int, p.C.Cells[q])
		for c := range idx[q] {
			idx[q][c] = -1
			if p.free(q, c) {
				idx[q][c] = n
				n++
			}
		}
	}
	if n == 0 {
		return nil, 0, fmt.Errorf("%w: no free cells", ErrPolyhedron)
	}
	return idx, n, nil
}

// acCouple emits the entries of −WK: row (q, c) gets −W_q·δ_{q−1} from
// grade q−1 and +δ_qᵀW_{q+1} from grade q+1.
func (p *Polyhedron) acCouple(idx [][]int, add func(r, c int, v float64)) {
	for q := p.Lo + 1; q <= p.Hi; q++ {
		for t, fs := range p.C.bnd[q] {
			it := idx[q][t]
			if it < 0 {
				continue
			}
			for _, f := range fs {
				is := idx[q-1][f.Cell]
				if is < 0 {
					continue
				}
				sg := float64(f.Sign)
				add(it, is, -p.W[q][t]*sg) // row q: −W_q(δ u_{q−1})
				add(is, it, p.W[q][t]*sg)  // row q−1: +(∂ W_q u_q)
			}
		}
	}
}

func (p *Polyhedron) acScatter(idx [][]int, x []complex128) [][]complex128 {
	u := make([][]complex128, p.C.Dim+1)
	for q := range u {
		u[q] = make([]complex128, p.C.Cells[q])
		for c, i := range idx[q] {
			if i >= 0 {
				u[q][c] = x[i]
			}
		}
	}
	return u
}

// ---------------------------------------------------------------------
// Dense real helpers (small systems only)

const elimTol = 1e-10

// rref returns the reduced row echelon form of a copy of a (r × n) and its
// pivot columns.
func rref(a [][]float64, n int) ([][]float64, []int) {
	m := make([][]float64, len(a))
	for i := range a {
		m[i] = append([]float64(nil), a[i]...)
	}
	var piv []int
	row := 0
	for col := 0; col < n && row < len(m); col++ {
		best, bv := -1, elimTol
		for r := row; r < len(m); r++ {
			if v := math.Abs(m[r][col]); v > bv {
				best, bv = r, v
			}
		}
		if best < 0 {
			continue
		}
		m[row], m[best] = m[best], m[row]
		pv := m[row][col]
		for c := range m[row] {
			m[row][c] /= pv
		}
		for r := range m {
			if r != row && m[r][col] != 0 {
				f := m[r][col]
				for c := range m[r] {
					m[r][c] -= f * m[row][c]
				}
			}
		}
		piv = append(piv, col)
		row++
	}
	return m[:row], piv
}

func matRank(a [][]float64) int {
	if len(a) == 0 {
		return 0
	}
	_, piv := rref(a, len(a[0]))
	return len(piv)
}

// nullSpace returns a basis of {x : a·x = 0} as the columns of an n × m
// matrix (a has n columns; nil a means every vector).
func nullSpace(a [][]float64, n int) [][]float64 {
	r, piv := rref(a, n)
	isPiv := make([]bool, n)
	for _, c := range piv {
		isPiv[c] = true
	}
	var freeCols []int
	for c := 0; c < n; c++ {
		if !isPiv[c] {
			freeCols = append(freeCols, c)
		}
	}
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, len(freeCols))
	}
	for k, fc := range freeCols {
		out[fc][k] = 1
		for i, pc := range piv {
			out[pc][k] = -r[i][fc]
		}
	}
	return out
}

// rowSpace returns a basis of the row space of a as the columns of an
// n × rank matrix (independent rows of a, transposed).
func rowSpace(a [][]float64, n int) [][]float64 {
	r, _ := rref(a, n)
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, len(r))
		for k := range r {
			out[i][k] = r[k][i]
		}
	}
	return out
}

// solveReal solves g·X = B (g square) by Gaussian elimination with partial
// pivoting; B and X are given row-wise (g's order × columns).
func solveReal(g [][]float64, b [][]float64) ([][]float64, error) {
	n := len(g)
	a := make([][]float64, n)
	for i := range g {
		a[i] = append(append([]float64(nil), g[i]...), b[i]...)
	}
	w := len(a[0])
	scale := 0.0
	for i := range g {
		for _, v := range g[i] {
			scale = math.Max(scale, math.Abs(v))
		}
	}
	for col := 0; col < n; col++ {
		best, bv := col, math.Abs(a[col][col])
		for r := col + 1; r < n; r++ {
			if v := math.Abs(a[r][col]); v > bv {
				best, bv = r, v
			}
		}
		if bv <= 1e-14*scale {
			return nil, fmt.Errorf("%w: pivot %d", ErrSingular, col)
		}
		a[col], a[best] = a[best], a[col]
		for r := col + 1; r < n; r++ {
			f := a[r][col] / a[col][col]
			if f == 0 {
				continue
			}
			for c := col; c < w; c++ {
				a[r][c] -= f * a[col][c]
			}
		}
	}
	x := make([][]float64, n)
	for i := range x {
		x[i] = make([]float64, w-n)
	}
	for c := 0; c < w-n; c++ {
		for r := n - 1; r >= 0; r-- {
			s := a[r][n+c]
			for k := r + 1; k < n; k++ {
				s -= a[r][k] * x[k][c]
			}
			x[r][c] = s / a[r][r]
		}
	}
	return x, nil
}

func matVec(a [][]float64, x []float64) []float64 {
	out := make([]float64, len(a))
	for i, row := range a {
		s := 0.0
		for j, v := range row {
			s += v * x[j]
		}
		out[i] = s
	}
	return out
}

// symEigenvalues returns the eigenvalues of a symmetric matrix by cyclic
// Jacobi rotations (Golub & Van Loan, Matrix Computations, §8.5).
func symEigenvalues(a0 [][]float64) []float64 {
	n := len(a0)
	a := make([][]float64, n)
	for i := range a0 {
		a[i] = append([]float64(nil), a0[i]...)
	}
	for sweep := 0; sweep < 100; sweep++ {
		off, diag := 0.0, 0.0
		for i := 0; i < n; i++ {
			diag += a[i][i] * a[i][i]
			for j := i + 1; j < n; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off <= 1e-30*(diag+1e-300) {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if a[p][q] == 0 {
					continue
				}
				theta := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := 1 / (math.Abs(theta) + math.Sqrt(theta*theta+1))
				if theta < 0 {
					t = -t
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p] = c*akp - s*akq
					a[k][q] = s*akp + c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k] = c*apk - s*aqk
					a[q][k] = s*apk + c*aqk
				}
			}
		}
	}
	out := make([]float64, n)
	for i := range out {
		out[i] = a[i][i]
	}
	return out
}
