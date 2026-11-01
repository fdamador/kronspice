// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
)

// Quasi-static extraction of lumped parameters from geometry (CLAUDE.md
// §3c "What em/ does" item 1, milestone M9i).
//
// Cross-sections (2-D): the electrostatic potential of a uniform line
// obeys div(ε grad φ) = 0 between conductors — the scalar part of
// Maxwell's ∇𝔇 with no free charge (Treatise art. 83 in SI, doc.go). It is
// solved by finite volumes on a uniform grid inside a grounded rectangular
// box, conductors held at fixed potentials. Conductor boundaries cut grid
// edges at their true position (Shortley–Weller: the flux across an edge
// cut at fraction θ is ε·Δφ/θ), so curved conductors converge at second
// order instead of as a staircase. The charge per unit length on a
// conductor is the discrete flux of D out of it (Gauss), which makes the
// capacitance matrix exactly symmetric. For lines whose field is TEM in a
// homogeneous medium (or quasi-TEM), the inductance matrix per unit length
// follows from the vacuum capacitance matrix C₀: L = μ₀ε₀·C₀⁻¹ (C. R. Paul,
// Analysis of Multiconductor Transmission Lines, 2nd ed., Wiley 2008,
// §3.2.3), so no separate magnetostatic solve is needed for non-magnetic
// lines.
//
// Coils (3-D, filamentary): mutual inductance by Neumann's formula
// M = μ₀/4π ∮∮ dl₁·dl₂/|r₁ − r₂| (Maxwell, Treatise Vol. II art. 524),
// exact for coaxial circles in elliptic integrals (art. 701), and a
// solenoid's current sheet by integrating that M over its length.

// ErrSolve is a field solve that failed to converge.
var ErrSolve = errors.New("em: field solve did not converge")

// Conductor2D is a conductor of the cross-section: Inside reports whether
// a point lies in it. A zero-thickness strip (Strip) instead occupies the
// grid row nearest its height between X0 and X1, its ends located exactly.
type Conductor2D struct {
	Name   string
	Inside func(x, y float64) bool

	IsStrip   bool
	X0, X1, Y float64
}

// Section is a 2-D problem: a grounded box [0, W] × [0, H] with grid
// spacing h, conductors, and a relative permittivity map (nil = vacuum).
type Section struct {
	W, H, Step float64
	Conductors []Conductor2D
	EpsR       func(x, y float64) float64
}

// Disk is a round conductor of radius r centred at (cx, cy).
func Disk(name string, cx, cy, r float64) Conductor2D {
	return Conductor2D{Name: name, Inside: func(x, y float64) bool { return (x-cx)*(x-cx)+(y-cy)*(y-cy) <= r*r }}
}

// Rect is a rectangular conductor [x0, x1] × [y0, y1].
func Rect(name string, x0, y0, x1, y1 float64) Conductor2D {
	return Conductor2D{Name: name, Inside: func(x, y float64) bool { return x >= x0 && x <= x1 && y >= y0 && y <= y1 }}
}

// Strip is a zero-thickness strip from x0 to x1 at height y (on the
// nearest grid row; place y on a row for exact results).
func Strip(name string, x0, x1, y float64) Conductor2D {
	return Conductor2D{Name: name, IsStrip: true, X0: x0, X1: x1, Y: y,
		Inside: func(xx, yy float64) bool { return false }}
}

// Outside is the complement of a disk: the outer conductor of a coax of
// inner radius r (everything outside it is metal).
func Outside(name string, cx, cy, r float64) Conductor2D {
	return Conductor2D{Name: name, Inside: func(x, y float64) bool { return (x-cx)*(x-cx)+(y-cy)*(y-cy) >= r*r }}
}

// SectionResult holds the per-unit-length matrices of a cross-section.
type SectionResult struct {
	Names []string
	C     [][]float64 // capacitance matrix [F/m] (Maxwell's coefficients)
	C0    [][]float64 // with every dielectric replaced by vacuum
	L     [][]float64 // inductance matrix μ₀ε₀·C₀⁻¹ [H/m]
}

// Solve extracts the capacitance and inductance matrices. Conductors
// named with the prefix "gnd" are held at 0 V with the box.
func (s Section) Solve() (*SectionResult, error) {
	g, err := s.grid()
	if err != nil {
		return nil, err
	}
	var sig []int
	res := &SectionResult{}
	for k, c := range s.Conductors {
		if len(c.Name) >= 3 && c.Name[:3] == "gnd" {
			continue
		}
		sig = append(sig, k)
		res.Names = append(res.Names, c.Name)
	}
	solveAll := func(vac bool) ([][]float64, error) {
		n := len(sig)
		m := make([][]float64, n)
		for i := range m {
			m[i] = make([]float64, n)
		}
		for col, k := range sig {
			v := make([]float64, len(s.Conductors))
			v[k] = 1
			q, err := g.charges(v, vac)
			if err != nil {
				return nil, err
			}
			for row, j := range sig {
				m[row][col] = q[j]
			}
		}
		// Symmetrize round-off (the discrete operator is symmetric).
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				a := 0.5 * (m[i][j] + m[j][i])
				m[i][j], m[j][i] = a, a
			}
		}
		return m, nil
	}
	if res.C, err = solveAll(false); err != nil {
		return nil, err
	}
	if s.EpsR == nil {
		res.C0 = res.C
	} else if res.C0, err = solveAll(true); err != nil {
		return nil, err
	}
	inv, err := invert(res.C0)
	if err != nil {
		return nil, err
	}
	res.L = make([][]float64, len(inv))
	for i := range inv {
		res.L[i] = make([]float64, len(inv))
		for j := range inv[i] {
			res.L[i][j] = Mu0 * Eps0 * inv[i][j]
		}
	}
	return res, nil
}

// grid2 is the discretized section.
type grid2 struct {
	s      Section
	nx, ny int
	h      float64
	owner  []int     // conductor index of each node, −1 free, −2 box
	eps    []float64 // relative permittivity of each cell (nx−1)·(ny−1)
	// For every edge from a free node to a fixed node: the fraction θ of
	// the edge at which the conductor boundary lies.
}

func (s Section) grid() (*grid2, error) {
	if !(s.Step > 0) || s.W <= 4*s.Step || s.H <= 4*s.Step {
		return nil, fmt.Errorf("%w: section %g × %g with step %g", ErrGrid, s.W, s.H, s.Step)
	}
	nx, ny := int(math.Round(s.W/s.Step))+1, int(math.Round(s.H/s.Step))+1
	if nx*ny > 4_000_000 {
		return nil, fmt.Errorf("%w: %d × %d grid", ErrGrid, nx, ny)
	}
	g := &grid2{s: s, nx: nx, ny: ny, h: s.Step, owner: make([]int, nx*ny), eps: make([]float64, (nx-1)*(ny-1))}
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			n := i + nx*j
			g.owner[n] = -1
			if i == 0 || j == 0 || i == nx-1 || j == ny-1 {
				g.owner[n] = -2
				continue
			}
			x, y := float64(i)*g.h, float64(j)*g.h
			for k, c := range s.Conductors {
				in := c.Inside(x, y)
				if c.IsStrip {
					in = j == int(math.Round(c.Y/g.h)) && x >= c.X0-1e-12*g.h && x <= c.X1+1e-12*g.h
				}
				if in {
					g.owner[n] = k
					break
				}
			}
		}
	}
	for j := 0; j < ny-1; j++ {
		for i := 0; i < nx-1; i++ {
			e := 1.0
			if s.EpsR != nil {
				e = s.EpsR((float64(i)+0.5)*g.h, (float64(j)+0.5)*g.h)
			}
			g.eps[i+(nx-1)*j] = e
		}
	}
	return g, nil
}

// edgeEps is ε_r on the dual face of the edge from node a to its neighbour
// along axis ax (0 = x, 1 = y): the mean of the two cells the face crosses.
func (g *grid2) edgeEps(i, j, ax int) float64 {
	cell := func(ci, cj int) (float64, bool) {
		if ci < 0 || cj < 0 || ci >= g.nx-1 || cj >= g.ny-1 {
			return 0, false
		}
		return g.eps[ci+(g.nx-1)*cj], true
	}
	var a, b float64
	var oa, ob bool
	if ax == 0 { // edge (i,j)–(i+1,j): cells (i, j−1) and (i, j)
		a, oa = cell(i, j-1)
		b, ob = cell(i, j)
	} else { // edge (i,j)–(i,j+1): cells (i−1, j) and (i, j)
		a, oa = cell(i-1, j)
		b, ob = cell(i, j)
	}
	switch {
	case oa && ob:
		return 0.5 * (a + b)
	case oa:
		return a
	case ob:
		return b
	}
	return 1
}

// fraction is where, along the edge from free point (x0,y0) to fixed point
// (x1,y1), the conductor boundary lies (θ ∈ (0, 1]), by bisection on the
// conductor's Inside test; the box wall and strips are taken at the node.
func (g *grid2) fraction(owner int, x0, y0, x1, y1 float64) float64 {
	if owner < 0 {
		return 1
	}
	c := g.s.Conductors[owner]
	if c.IsStrip {
		// Along the strip's row the boundary is its end; across it, the
		// strip itself.
		if y0 != y1 {
			return 1
		}
		end := c.X0
		if x1 > x0 {
			end = c.X0
		} else {
			end = c.X1
		}
		return math.Min(1, math.Max(math.Abs(end-x0)/math.Abs(x1-x0), 1e-3))
	}
	lo, hi := 0.0, 1.0
	for it := 0; it < 40; it++ {
		m := 0.5 * (lo + hi)
		if c.Inside(x0+m*(x1-x0), y0+m*(y1-y0)) {
			hi = m
		} else {
			lo = m
		}
	}
	return math.Max(hi, 1e-3)
}

// charges solves for the potential with conductor k at v[k] and returns
// the charge per unit length on every conductor [C/m].
func (g *grid2) charges(v []float64, vacuum bool) ([]float64, error) {
	nx, ny, h := g.nx, g.ny, g.h
	// Index the free nodes.
	idx := make([]int, nx*ny)
	nf := 0
	for n, o := range g.owner {
		if o == -1 {
			idx[n] = nf
			nf++
		} else {
			idx[n] = -1
		}
	}
	type nb struct {
		j int     // free index, or −1
		w float64 // coupling ε/θ
		v float64 // fixed potential when j = −1
		o int     // owner of the fixed node
	}
	nbs := make([][4]nb, nf)
	diag := make([]float64, nf)
	rhs := make([]float64, nf)
	pot := func(o int) float64 {
		if o < 0 {
			return 0
		}
		return v[o]
	}
	for n, o := range g.owner {
		if o != -1 {
			continue
		}
		i, j := n%nx, n/nx
		f := idx[n]
		x0, y0 := float64(i)*h, float64(j)*h
		dirs := [4][3]int{{1, 0, 0}, {-1, 0, 0}, {0, 1, 1}, {0, -1, 1}}
		for d, dd := range dirs {
			ii, jj := i+dd[0], j+dd[1]
			m := ii + nx*jj
			// ε on the edge's dual face.
			var e float64
			if dd[2] == 0 {
				e = g.edgeEps(min(i, ii), j, 0)
			} else {
				e = g.edgeEps(i, min(j, jj), 1)
			}
			if vacuum {
				e = 1
			}
			if g.owner[m] == -1 {
				nbs[f][d] = nb{j: idx[m], w: e}
				diag[f] += e
				continue
			}
			th := g.fraction(g.owner[m], x0, y0, float64(ii)*h, float64(jj)*h)
			w := e / th
			nbs[f][d] = nb{j: -1, w: w, v: pot(g.owner[m]), o: g.owner[m]}
			diag[f] += w
			rhs[f] += w * pot(g.owner[m])
		}
	}
	// Conjugate gradient on A·φ = rhs (A symmetric positive definite),
	// preconditioned by symmetric successive over-relaxation (SSOR, ω = 1.6;
	// Y. Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., SIAM
	// 2003, §10.2).
	apply := func(x, out []float64) {
		for f := range out {
			s := diag[f] * x[f]
			for _, b := range nbs[f] {
				if b.j >= 0 {
					s -= b.w * x[b.j]
				}
			}
			out[f] = s
		}
	}
	const om = 1.6
	precond := func(r, z []float64) {
		// Forward: (D + ωL)·y = r.
		for f := 0; f < nf; f++ {
			s := r[f]
			for _, b := range nbs[f] {
				if b.j >= 0 && b.j < f {
					s += om * b.w * z[b.j]
				}
			}
			z[f] = s / diag[f]
		}
		// u = D·y, then backward: (D + ωU)·z = u.
		for f := nf - 1; f >= 0; f-- {
			s := diag[f] * z[f]
			for _, b := range nbs[f] {
				if b.j > f {
					s += om * b.w * z[b.j]
				}
			}
			z[f] = s / diag[f]
		}
		for f := range z {
			z[f] *= om * (2 - om)
		}
	}
	x := make([]float64, nf)
	r := append([]float64(nil), rhs...)
	z := make([]float64, nf)
	precond(r, z)
	p := append([]float64(nil), z...)
	ap := make([]float64, nf)
	rz := dotv(r, z)
	bnorm := math.Sqrt(dotv(rhs, rhs))
	if bnorm == 0 {
		return make([]float64, len(g.s.Conductors)), nil
	}
	conv := false
	for it := 0; it < 20*nf+100; it++ {
		apply(p, ap)
		alpha := rz / dotv(p, ap)
		for f := range x {
			x[f] += alpha * p[f]
			r[f] -= alpha * ap[f]
		}
		if math.Sqrt(dotv(r, r)) <= 1e-12*bnorm {
			conv = true
			break
		}
		precond(r, z)
		rz2 := dotv(r, z)
		beta := rz2 / rz
		rz = rz2
		for f := range p {
			p[f] = z[f] + beta*p[f]
		}
	}
	if !conv {
		return nil, ErrSolve
	}
	// Charge on each conductor: flux of D out of it through the edges to
	// free nodes, ε₀·Σ w·(V_c − φ_f) (per unit length; h cancels in 2-D).
	q := make([]float64, len(g.s.Conductors))
	for f := range nbs {
		for _, b := range nbs[f] {
			if b.j == -1 && b.o >= 0 {
				q[b.o] += Eps0 * b.w * (b.v - x[f])
			}
		}
	}
	// Conductors touching each other directly (conductor–conductor edges)
	// carry flux too: add it.
	for n, o := range g.owner {
		if o < 0 {
			continue
		}
		i, j := n%nx, n/nx
		for _, dd := range [2][2]int{{1, 0}, {0, 1}} {
			ii, jj := i+dd[0], j+dd[1]
			if ii >= nx || jj >= ny {
				continue
			}
			o2 := g.owner[ii+nx*jj]
			if o2 == -1 || o2 == o {
				continue
			}
			var e float64
			if dd[0] == 1 {
				e = g.edgeEps(i, j, 0)
			} else {
				e = g.edgeEps(i, j, 1)
			}
			if vacuum {
				e = 1
			}
			dq := Eps0 * e * (pot(o) - pot(o2))
			q[o] += dq
			if o2 >= 0 {
				q[o2] -= dq
			}
		}
	}
	return q, nil
}

func dotv(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}

// invert inverts a small dense matrix (Gauss–Jordan, partial pivoting).
func invert(a [][]float64) ([][]float64, error) {
	n := len(a)
	m := make([][]float64, n)
	for i := range a {
		m[i] = make([]float64, 2*n)
		copy(m[i], a[i])
		m[i][n+i] = 1
	}
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		if m[p][c] == 0 {
			return nil, fmt.Errorf("em: singular capacitance matrix")
		}
		m[c], m[p] = m[p], m[c]
		piv := m[c][c]
		for k := range m[c] {
			m[c][k] /= piv
		}
		for r := 0; r < n; r++ {
			if r == c || m[r][c] == 0 {
				continue
			}
			f := m[r][c]
			for k := range m[r] {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	out := make([][]float64, n)
	for i := range out {
		out[i] = m[i][n:]
	}
	return out, nil
}
