// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
)

// Gauge theory on Kron's network (CLAUDE.md §1, register entry H22).
//
// A lattice gauge theory lives on a graph: a group element on every link
// (the parallel transporter; K. G. Wilson, Phys. Rev. D 10, 2445, 1974), a
// gauge transformation at every node, and gauge-invariant quantities on
// closed loops (Wilson loops). Kron's network already has the three
// objects: the branch, the junction and the mesh. This file writes the
// correspondence out for the abelian group U(1), where it is Maxwell's
// electrodynamics on Kron's 1944 field network (Proc. IRE 32, 289, 1944;
// E. Tonti, PIER 32, 1, 2001) and, in the lumped limit, Kirchhoff's laws;
// and for SU(3), the gauge group of the quark model, in the Hamiltonian
// lattice form of J. Kogut and L. Susskind (Phys. Rev. D 11, 395, 1975), on
// the same graph, so that the two can be set side by side (H22 step 2).
//
// # Cell complex
//
// Nodes (0-cells), links (1-cells, oriented from → to) and plaquettes
// (2-cells, oriented closed paths of links). Two operators,
//
//	Grad χ : nodes → links,      (Grad χ)_k = χ_from − χ_to    (Kron's E = A·v)
//	Curl a : links → plaquettes, (Curl a)_p = Σ_{k∈∂p} ±a_k   (Kron's e′ = Cᵀ·e)
//
// with Curl∘Grad = 0: the identity Aᵀ·C = 0 of connection.go, d∘d = 0 of
// the exterior calculus (Roth, PNAS 41, 518, 1955; Branin 1966). Their
// transposes are Div = Gradᵀ (Kron's Aᵀ·i, Kirchhoff's current law) and
// Curlᵀ (Kron's i = C·i′, branch currents from mesh currents).
//
// # U(1)
//
// Link variable a_k = ∫_k A·dl [V·s] in the temporal gauge A₀ = 0; the
// phase of the parallel transporter is θ_k = q·a_k/ħ (Peierls). Conjugate
// variable Q_k [C], the electric flux through the dual face of link k,
// which is the charge on the link's capacitance. Hamiltonian
//
//	H = Σ_k Q_k²/(2C_k) + Σ_p (Curl a)_p²/(2L_p)
//
// whose two terms are the ledger's ½Cv² (v = Q/C) and ½L i² (i_p = Φ_p/L_p,
// the mesh current; Φ_p = (Curl a)_p the flux through plaquette p).
// Hamilton's equations ȧ = ∂H/∂(−Q) and Q̇ = ∂H/∂a are Faraday's and
// Ampère–Maxwell's laws on the lattice — Yee's scheme on a rectangular
// grid (K. S. Yee, IEEE Trans. Antennas Propag. 14, 302, 1966) — and the
// generator of gauge transformations,
//
//	G_n = (Div Q)_n − ρ_n,
//
// is Gauss's law. Ġ = Div Q̇ = Div Curlᵀ(Φ/L) = 0 because Curl∘Grad = 0;
// with Q̇ = i that statement is Kirchhoff's current law. A pure gauge
// a = Grad χ has zero holonomy: it links no flux and induces no EMF,
// which is Kirchhoff's voltage law. (Flux–charge pair as in U. Vool and
// M. Devoret, Int. J. Circuit Theory Appl. 45, 897, 2017, whose branch
// flux is −a_k.)
//
// # SU(3)
//
// Link variable U_k ∈ SU(3); gauge transformation U_k → g_from·U_k·g_to†;
// Wilson loop W_p = P∏_{∂p} U (path ordered, U† when a link is traversed
// against its orientation); Gauss generator
//
//	G_n = Σ_{from=n} E_k − Σ_{to=n} U_k†·E_k·U_k − ρ_n
//
// with Lie-algebra-valued electric fields E_k (Kogut & Susskind §II: the
// field at the far end of a link is the near-end field parallel-
// transported, with the sign of the orientation). G_n transforms
// covariantly, G_n → g_n·G_n·g_n†, not invariantly: the "currents" meeting
// at a node add as matrices after transport, and the result depends on
// the transport. Two loops at one node compose by matrix multiplication,
// W(p₁p₂) ≠ W(p₂p₁) in general; in U(1) they compose by adding fluxes,
// which is the additivity every Kirchhoff sum rests on.

// ErrLattice reports an inconsistent cell complex or field.
var ErrLattice = errors.New("kron: invalid lattice")

// Step is one link of an oriented path: Dir = +1 along from → to, −1
// against it.
type Step struct{ Link, Dir int }

// Lattice is a cell complex: nodes, oriented links and oriented
// plaquettes (closed paths of links). It is the structure a gauge theory
// needs and nothing else; the constitutive data (capacitance per link,
// inductance per plaquette) belong to a U1Field on it.
type Lattice struct {
	Nodes      int
	From, To   []int
	Plaquettes [][]Step
	// Free marks links whose variables evolve. A link that is not free
	// holds a = 0 and Q = 0 (a conducting wall: no tangential E, no
	// flux through a face in the wall). nil means every link is free.
	Free []bool
	// Names labels the links (circuit element names), optional.
	Names []string
}

// Links is the number of links b.
func (l *Lattice) Links() int { return len(l.From) }

// free reports whether link k evolves.
func (l *Lattice) free(k int) bool { return l.Free == nil || l.Free[k] }

// Check verifies the complex: links join existing nodes, plaquettes are
// closed paths of consecutive links, and Curl∘Grad = 0 on a random
// potential (the Kirchhoff identity Aᵀ·C = 0).
func (l *Lattice) Check() error {
	if l.Nodes <= 0 || len(l.From) != len(l.To) {
		return fmt.Errorf("%w: %d nodes, %d/%d link ends", ErrLattice, l.Nodes, len(l.From), len(l.To))
	}
	for k := range l.From {
		if l.From[k] < 0 || l.From[k] >= l.Nodes || l.To[k] < 0 || l.To[k] >= l.Nodes {
			return fmt.Errorf("%w: link %d joins nodes %d–%d", ErrLattice, k, l.From[k], l.To[k])
		}
	}
	for p, path := range l.Plaquettes {
		if len(path) == 0 {
			return fmt.Errorf("%w: plaquette %d is empty", ErrLattice, p)
		}
		for i, s := range path {
			if s.Link < 0 || s.Link >= len(l.From) || (s.Dir != 1 && s.Dir != -1) {
				return fmt.Errorf("%w: plaquette %d step %d: link %d dir %d", ErrLattice, p, i, s.Link, s.Dir)
			}
		}
		start, at := l.pathStart(path[0]), l.pathStart(path[0])
		for i, s := range path {
			if l.pathStart(s) != at {
				return fmt.Errorf("%w: plaquette %d is not a path at step %d", ErrLattice, p, i)
			}
			at = l.pathEnd(s)
		}
		if at != start {
			return fmt.Errorf("%w: plaquette %d does not close", ErrLattice, p)
		}
	}
	r := rand.New(rand.NewPCG(7, 11))
	chi := make([]float64, l.Nodes)
	scale := 0.0
	for n := range chi {
		chi[n] = r.NormFloat64()
		scale = math.Max(scale, math.Abs(chi[n]))
	}
	for p, v := range l.Curl(l.Grad(chi)) {
		if math.Abs(v) > 1e-12*scale*float64(len(l.Plaquettes[p])) {
			return fmt.Errorf("%w: Curl∘Grad ≠ 0 on plaquette %d (%g)", ErrLattice, p, v)
		}
	}
	return nil
}

func (l *Lattice) pathStart(s Step) int {
	if s.Dir > 0 {
		return l.From[s.Link]
	}
	return l.To[s.Link]
}

func (l *Lattice) pathEnd(s Step) int {
	if s.Dir > 0 {
		return l.To[s.Link]
	}
	return l.From[s.Link]
}

// Grad is the coboundary of a node function: (Grad χ)_k = χ_from − χ_to,
// Kron's E = A·v (connection.go incidence; Kron 1939 Ch. XIV p. 360).
func (l *Lattice) Grad(chi []float64) []float64 {
	out := make([]float64, len(l.From))
	for k := range out {
		out[k] = chi[l.From[k]] - chi[l.To[k]]
	}
	return out
}

// Div is the transpose of Grad: (Div q)_n = Σ_{from=n} q_k − Σ_{to=n} q_k,
// the net flux leaving node n. Kron's Aᵀ·i; Kirchhoff's current law is
// Div i = 0 (continuity is Div i + ρ̇ = 0).
func (l *Lattice) Div(q []float64) []float64 {
	out := make([]float64, l.Nodes)
	for k := range l.From {
		out[l.From[k]] += q[k]
		out[l.To[k]] -= q[k]
	}
	return out
}

// Curl is the coboundary of a link function: (Curl a)_p = Σ_{∂p} ±a_k,
// the holonomy of a around plaquette p. On Kron's network it is e′ = Cᵀ·e,
// the loop sum of branch quantities (Kron 1939 eq 4.26).
func (l *Lattice) Curl(a []float64) []float64 {
	out := make([]float64, len(l.Plaquettes))
	for p, path := range l.Plaquettes {
		s := 0.0
		for _, st := range path {
			s += float64(st.Dir) * a[st.Link]
		}
		out[p] = s
	}
	return out
}

// CurlT is the transpose of Curl: the link quantity made by plaquette
// quantities, (Curlᵀ i′)_k = Σ_{p∋k} ±i′_p. On Kron's network it is
// i = C·i′, the branch currents of the mesh currents (Kron 1939 eq 4.23).
func (l *Lattice) CurlT(ip []float64) []float64 {
	out := make([]float64, len(l.From))
	for p, path := range l.Plaquettes {
		for _, st := range path {
			out[st.Link] += float64(st.Dir) * ip[p]
		}
	}
	return out
}

// Gauge applies the U(1) gauge transformation generated by χ: a → a +
// Grad χ (A → A + ∇χ along every link). Holonomies are unchanged because
// Curl∘Grad = 0.
func (l *Lattice) Gauge(a, chi []float64) []float64 {
	g := l.Grad(chi)
	out := make([]float64, len(a))
	for k := range out {
		out[k] = a[k] + g[k]
	}
	return out
}

// WilsonLoop is the U(1) Wilson loop of plaquette p for a charge q [C]:
// exp(i·q·Φ_p/ħ) with Φ_p = (Curl a)_p the flux the loop links. For a
// Cooper pair (q = 2e) the phase is 2π·Φ_p/Φ₀ with Φ₀ = h/2e.
func (l *Lattice) WilsonLoop(a []float64, p int, q float64) complex128 {
	return cmplx.Exp(complex(0, q*l.Curl(a)[p]/HBar))
}

// Rebase rotates a closed path so that it starts and ends at node n; ok
// is false when the path does not pass through n.
func (l *Lattice) Rebase(path []Step, n int) (out []Step, ok bool) {
	for i, s := range path {
		if l.pathStart(s) == n {
			out = append(out, path[i:]...)
			out = append(out, path[:i]...)
			return out, true
		}
	}
	return nil, false
}

// Physical constants (CODATA 2018, exact since the 2019 SI).
const (
	// ElementaryCharge e [C].
	ElementaryCharge = 1.602176634e-19
	// Planck h [J·s].
	Planck = 6.62607015e-34
	// HBar ħ = h/2π [J·s].
	HBar = Planck / (2 * math.Pi)
	// FluxQuantum Φ₀ = h/2e [Wb], the flux that turns a Cooper pair's
	// Wilson loop through 2π.
	FluxQuantum = Planck / (2 * ElementaryCharge)
	// JosephsonConstant K_J = 2e/h [Hz/V]: the volt from a frequency.
	JosephsonConstant = 2 * ElementaryCharge / Planck
	// VonKlitzingConstant R_K = h/e² [Ω]: the ohm from the quantum Hall
	// plateaus.
	VonKlitzingConstant = Planck / (ElementaryCharge * ElementaryCharge)
)

// LatticeFromNetwork is the cell complex of a Kron network: its junctions,
// its branches as links, and its fundamental meshes as plaquettes, each
// an ordered closed path (the link, then back through the spanning tree),
// so that the signed incidence of plaquette p is column p of the
// connection tensor C (connection.go loops).
func LatticeFromNetwork(n *Network) (*Lattice, error) {
	g := n.g
	l := &Lattice{Nodes: g.n, From: append([]int(nil), g.from...), To: append([]int(nil), g.to...)}
	l.Names = make([]string, len(n.branches))
	for k, b := range n.branches {
		switch {
		case b.elem >= 0:
			l.Names[k] = n.c.Elements[b.elem].Name
		case b.kind == bControl:
			l.Names[k] = fmt.Sprintf("ctrl%d", k)
		default:
			l.Names[k] = fmt.Sprintf("branch%d", k)
		}
	}
	for k := range g.from {
		if g.tree[k] {
			continue
		}
		l.Plaquettes = append(l.Plaquettes, g.loopPath(k))
	}
	if err := l.Check(); err != nil {
		return nil, err
	}
	return l, nil
}

// loopPath is the fundamental loop closed by link k as an ordered path:
// along k from its From node to its To node, then back through the tree
// (up from To to the common ancestor, down to From). Same loop and
// direction as column k of loops().
func (g *graph) loopPath(k int) []Step {
	path := []Step{{k, 1}}
	u, w := g.from[k], g.to[k]
	var up, down []Step // w → ancestor, u → ancestor
	climb := func(x int) (Step, int) {
		pe, p := g.parentEdge[x], g.parent[x]
		s := Step{pe, 1}
		if g.from[pe] != x {
			s.Dir = -1
		}
		return s, p
	}
	for g.depth[w] > g.depth[u] {
		var s Step
		s, w = climb(w)
		up = append(up, s)
	}
	for g.depth[u] > g.depth[w] {
		var s Step
		s, u = climb(u)
		down = append(down, s)
	}
	for u != w {
		var s Step
		s, w = climb(w)
		up = append(up, s)
		s, u = climb(u)
		down = append(down, s)
	}
	path = append(path, up...)
	for i := len(down) - 1; i >= 0; i-- {
		path = append(path, Step{down[i].Link, -down[i].Dir})
	}
	return path
}

// U1Field is U(1) electrodynamics on a lattice: a link variable a and a
// link charge Q on every link, with the constitutive data C per link and
// L per plaquette, stepped by the leapfrog that Hamilton's equations
// give (Yee's scheme on a grid). a lives at t + Dt/2 and Q at t.
type U1Field struct {
	Lat *Lattice
	C   []float64 // capacitance per link [F] (electric energy Q²/2C)
	L   []float64 // inductance per plaquette [H] (magnetic energy Φ²/2L)
	A   []float64 // link variables a_k [V·s] at Time() + Dt/2
	Q   []float64 // link charges [C] at Time()
	Dt  float64
	// Steps completed.
	Steps int

	aPrev []float64 // a at Time() − Dt/2, for the exact discrete energy
}

// NewU1Field makes a field on l with the given constitutive data and time
// step, every variable zero. Links that are not free need no C; plaquettes
// need L > 0.
func NewU1Field(l *Lattice, c, lp []float64, dt float64) (*U1Field, error) {
	if err := l.Check(); err != nil {
		return nil, err
	}
	if len(c) != l.Links() || len(lp) != len(l.Plaquettes) || !(dt > 0) {
		return nil, fmt.Errorf("%w: %d C for %d links, %d L for %d plaquettes, dt %g",
			ErrLattice, len(c), l.Links(), len(lp), len(l.Plaquettes), dt)
	}
	for k, v := range c {
		if l.free(k) && !(v > 0) {
			return nil, fmt.Errorf("%w: link %d has C = %g", ErrLattice, k, v)
		}
	}
	for p, v := range lp {
		if !(v > 0) {
			return nil, fmt.Errorf("%w: plaquette %d has L = %g", ErrLattice, p, v)
		}
	}
	b := l.Links()
	return &U1Field{Lat: l, C: c, L: lp, A: make([]float64, b), Q: make([]float64, b),
		aPrev: make([]float64, b), Dt: dt}, nil
}

// SetVoltage sets the link charges from branch voltages at Time(), Q = C·v,
// on the free links, and makes a(−Dt/2) consistent with them.
func (f *U1Field) SetVoltage(v []float64) {
	for k := range f.Q {
		if f.Lat.free(k) {
			f.Q[k] = f.C[k] * v[k]
		}
	}
	f.begin()
}

// begin sets a(−Dt/2) consistent with a(+Dt/2) and Q(0): a^{n+½} =
// a^{n−½} − Dt·v^n.
func (f *U1Field) begin() {
	for k := range f.A {
		f.aPrev[k] = f.A[k]
		if f.Lat.free(k) {
			f.aPrev[k] += f.Dt * f.Q[k] / f.C[k]
		}
	}
}

// Time is the time of Q.
func (f *U1Field) Time() float64 { return float64(f.Steps) * f.Dt }

// Voltage is v = Q/C on every link (0 on fixed links).
func (f *U1Field) Voltage() []float64 {
	v := make([]float64, len(f.Q))
	for k := range v {
		if f.Lat.free(k) {
			v[k] = f.Q[k] / f.C[k]
		}
	}
	return v
}

// Flux is Φ_p = (Curl a)_p, the flux through every plaquette at
// Time() + Dt/2.
func (f *U1Field) Flux() []float64 { return f.Lat.Curl(f.A) }

// MeshCurrent is i_p = Φ_p/L_p, the mesh (plaquette) current.
func (f *U1Field) MeshCurrent() []float64 {
	ip := f.Flux()
	for p := range ip {
		ip[p] /= f.L[p]
	}
	return ip
}

// Step advances one Dt: Ampère–Maxwell (Q by the currents the plaquettes
// drive through each link, Q̇ = Curlᵀ(Φ/L)), then Faraday (ȧ = −v).
func (f *U1Field) Step() {
	ik := f.Lat.CurlT(f.MeshCurrent())
	for k := range f.Q {
		if f.Lat.free(k) {
			f.Q[k] += f.Dt * ik[k]
		}
	}
	for k := range f.A {
		f.aPrev[k] = f.A[k]
		if f.Lat.free(k) {
			f.A[k] -= f.Dt * f.Q[k] / f.C[k]
		}
	}
	f.Steps++
}

// Energy returns the Hamiltonian's two terms [J]: electric Σ Q²/2C
// (= Σ ½Cv²) at Time(), and magnetic Σ Φ⁻·Φ⁺/2L with the fluxes half a
// step either side of Time() (= Σ ½L i² in the leapfrog's exact discrete
// form, the same product em/fdtd.go's ledger uses; it is conserved to
// round-off, whereas Σ Φ⁺²/2L oscillates at O(Dt²)).
func (f *U1Field) Energy() (electric, magnetic float64) {
	for k, q := range f.Q {
		if f.Lat.free(k) {
			electric += q * q / (2 * f.C[k])
		}
	}
	prev, next := f.Lat.Curl(f.aPrev), f.Lat.Curl(f.A)
	for p := range prev {
		magnetic += prev[p] * next[p] / (2 * f.L[p])
	}
	return electric, magnetic
}

// Gauss is the gauge generator G_n = (Div Q)_n − ρ_n at every node, zero
// where Gauss's law holds with the declared charge ρ (nil: no charge).
// On a conducting wall its value is the induced surface charge.
func (f *U1Field) Gauss(rho []float64) []float64 {
	g := f.Lat.Div(f.Q)
	if rho != nil {
		for n := range g {
			g[n] -= rho[n]
		}
	}
	return g
}

// YeeGrid is a rectangular lattice of N[0]×N[1]×N[2] nodes with spacing D,
// closed by conducting walls: links are the grid edges (tangential edges
// on the walls are fixed), plaquettes the grid faces, oriented by the
// right-hand rule about their normal. It is Kron's 1944 field network and
// the Yee cell of em/fdtd.go, so a U1Field on it is the FDTD solve.
type YeeGrid struct {
	*Lattice
	N    [3]int
	D    [3]float64
	link [3][]int
	face [3][]int
}

// NewYeeGrid builds the grid (every N ≥ 3).
func NewYeeGrid(n [3]int, d [3]float64) (*YeeGrid, error) {
	for a := 0; a < 3; a++ {
		if n[a] < 3 || !(d[a] > 0) {
			return nil, fmt.Errorf("%w: Yee grid %v spacing %v", ErrLattice, n, d)
		}
	}
	g := &YeeGrid{Lattice: &Lattice{Nodes: n[0] * n[1] * n[2]}, N: n, D: d}
	idx := func(i, j, k int) int { return i + n[0]*(j+n[1]*k) }
	for c := 0; c < 3; c++ {
		g.link[c] = make([]int, g.Nodes)
		g.face[c] = make([]int, g.Nodes)
		for i := range g.link[c] {
			g.link[c][i], g.face[c][i] = -1, -1
		}
	}
	for c := 0; c < 3; c++ {
		for k := 0; k < n[2]; k++ {
			for j := 0; j < n[1]; j++ {
				for i := 0; i < n[0]; i++ {
					p := [3]int{i, j, k}
					if p[c] == n[c]-1 {
						continue
					}
					q := p
					q[c]++
					g.link[c][idx(i, j, k)] = len(g.From)
					g.From = append(g.From, idx(i, j, k))
					g.To = append(g.To, idx(q[0], q[1], q[2]))
					free := true
					for a := 0; a < 3; a++ {
						if a != c && (p[a] == 0 || p[a] == n[a]-1) {
							free = false
						}
					}
					g.Free = append(g.Free, free)
				}
			}
		}
	}
	// Faces: normal c, spanned by axes a = c+1 and b = c+2 (cyclic), path
	// +a at the corner, +b at corner+a, −a at corner+b, −b at the corner.
	for c := 0; c < 3; c++ {
		a, b := (c+1)%3, (c+2)%3
		for k := 0; k < n[2]; k++ {
			for j := 0; j < n[1]; j++ {
				for i := 0; i < n[0]; i++ {
					p := [3]int{i, j, k}
					if p[a] == n[a]-1 || p[b] == n[b]-1 {
						continue
					}
					pa, pb := p, p
					pa[a]++
					pb[b]++
					g.face[c][idx(i, j, k)] = len(g.Plaquettes)
					g.Plaquettes = append(g.Plaquettes, []Step{
						{g.link[a][idx(p[0], p[1], p[2])], 1},
						{g.link[b][idx(pa[0], pa[1], pa[2])], 1},
						{g.link[a][idx(pb[0], pb[1], pb[2])], -1},
						{g.link[b][idx(p[0], p[1], p[2])], -1},
					})
				}
			}
		}
	}
	if err := g.Check(); err != nil {
		return nil, err
	}
	return g, nil
}

// Link is the link index of the edge along axis c starting at node
// (i, j, k), −1 if there is none.
func (g *YeeGrid) Link(c, i, j, k int) int { return g.link[c][i+g.N[0]*(j+g.N[1]*k)] }

// Face is the plaquette index of the face with normal c whose lowest
// corner is node (i, j, k), −1 if there is none.
func (g *YeeGrid) Face(c, i, j, k int) int { return g.face[c][i+g.N[0]*(j+g.N[1]*k)] }

// Axis is the axis of link k.
func (g *YeeGrid) Axis(k int) int {
	n := g.From[k]
	for c := 0; c < 3; c++ {
		if g.link[c][n] == k {
			return c
		}
	}
	return -1
}

// Position is the midpoint of link k [m] (origin at node 0).
func (g *YeeGrid) Position(k int) [3]float64 {
	n := g.From[k]
	i, j, l := n%g.N[0], (n/g.N[0])%g.N[1], n/(g.N[0]*g.N[1])
	pos := [3]float64{float64(i) * g.D[0], float64(j) * g.D[1], float64(l) * g.D[2]}
	if c := g.Axis(k); c >= 0 {
		pos[c] += g.D[c] / 2
	}
	return pos
}

// Constitutive returns the capacitance of every edge, C = ε·(dual face)/
// (edge length), and the inductance of every face, L = μ·(face)/(dual edge
// length), for a uniform medium: Kron's 1944 network constants, which
// make Q²/2C = ½εE²·dV and Φ²/2L = B²/2μ·dV cell by cell.
func (g *YeeGrid) Constitutive(eps, mu float64) (c, lp []float64) {
	c = make([]float64, g.Links())
	for k := range c {
		ax := g.Axis(k)
		a, b := (ax+1)%3, (ax+2)%3
		c[k] = eps * g.D[a] * g.D[b] / g.D[ax]
	}
	lp = make([]float64, len(g.Plaquettes))
	for n := range g.face[0] {
		for ax := 0; ax < 3; ax++ {
			if p := g.face[ax][n]; p >= 0 {
				a, b := (ax+1)%3, (ax+2)%3
				lp[p] = mu * g.D[a] * g.D[b] / g.D[ax]
			}
		}
	}
	return c, lp
}

// ---------------------------------------------------------------------
// SU(3)

// SU3 is a 3×3 complex matrix; the group elements are the unitary ones
// of determinant 1, the Lie algebra the traceless Hermitian ones (so that
// exp(i·h) is in the group).
type SU3 [3][3]complex128

// Identity3 is the unit matrix.
func Identity3() SU3 { return SU3{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} }

// Mul is the matrix product a·b.
func (a SU3) Mul(b SU3) SU3 {
	var c SU3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j]
		}
	}
	return c
}

// Dagger is the conjugate transpose.
func (a SU3) Dagger() SU3 {
	var c SU3
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			c[i][j] = cmplx.Conj(a[j][i])
		}
	}
	return c
}

// Add is a + b.
func (a SU3) Add(b SU3) SU3 {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			a[i][j] += b[i][j]
		}
	}
	return a
}

// Sub is a − b.
func (a SU3) Sub(b SU3) SU3 { return a.Add(b.Scale(-1)) }

// Scale is s·a.
func (a SU3) Scale(s complex128) SU3 {
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			a[i][j] *= s
		}
	}
	return a
}

// Trace is Σ a_ii.
func (a SU3) Trace() complex128 { return a[0][0] + a[1][1] + a[2][2] }

// Det is the determinant.
func (a SU3) Det() complex128 {
	return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) -
		a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) +
		a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
}

// Norm is the Frobenius norm.
func (a SU3) Norm() float64 {
	s := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			s += real(a[i][j])*real(a[i][j]) + imag(a[i][j])*imag(a[i][j])
		}
	}
	return math.Sqrt(s)
}

// Commutator is a·b − b·a.
func (a SU3) Commutator(b SU3) SU3 { return a.Mul(b).Sub(b.Mul(a)) }

// UnitaryDefect is ‖a†a − 1‖ + |det a − 1|, zero for an SU(3) element.
func (a SU3) UnitaryDefect() float64 {
	return a.Dagger().Mul(a).Sub(Identity3()).Norm() + cmplx.Abs(a.Det()-1)
}

// GellMann are the eight generators λ₁…λ₈ of su(3) (M. Gell-Mann, Phys.
// Rev. 125, 1067, 1962), traceless Hermitian with Tr(λ_a λ_b) = 2δ_ab;
// λ₃ and λ₈ span the Cartan subalgebra (the diagonal, commuting part).
var GellMann = [8]SU3{
	{{0, 1, 0}, {1, 0, 0}, {0, 0, 0}},
	{{0, -1i, 0}, {1i, 0, 0}, {0, 0, 0}},
	{{1, 0, 0}, {0, -1, 0}, {0, 0, 0}},
	{{0, 0, 1}, {0, 0, 0}, {1, 0, 0}},
	{{0, 0, -1i}, {0, 0, 0}, {1i, 0, 0}},
	{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}},
	{{0, 0, 0}, {0, 0, -1i}, {0, 1i, 0}},
	{{complex(1/math.Sqrt(3), 0), 0, 0}, {0, complex(1/math.Sqrt(3), 0), 0}, {0, 0, complex(-2/math.Sqrt(3), 0)}},
}

// Algebra is Σ c_a·λ_a/2, a traceless Hermitian matrix with the usual
// normalization T_a = λ_a/2.
func Algebra(c [8]float64) SU3 {
	var h SU3
	for a := range c {
		h = h.Add(GellMann[a].Scale(complex(c[a]/2, 0)))
	}
	return h
}

// Expm is the matrix exponential exp(m) by scaling and squaring with a
// Taylor series (Moler & Van Loan, SIAM Rev. 45, 3, 2003, method 3);
// for a traceless Hermitian h, Expm(i·h) is in SU(3).
func Expm(m SU3) SU3 {
	s := 0
	if nm := m.Norm(); nm > 0.5 {
		s = int(math.Ceil(math.Log2(nm/0.5))) + 1
	}
	m = m.Scale(complex(math.Ldexp(1, -s), 0))
	sum, term := Identity3(), Identity3()
	for n := 1; n <= 40; n++ {
		term = term.Mul(m).Scale(complex(1/float64(n), 0))
		sum = sum.Add(term)
		if term.Norm() < 1e-18 {
			break
		}
	}
	for ; s > 0; s-- {
		sum = sum.Mul(sum)
	}
	return sum
}

// Exp is the group element exp(i·h) of an algebra element h.
func Exp(h SU3) SU3 { return Expm(h.Scale(1i)) }

// RandomAlgebra is a random algebra element with Gaussian coefficients
// of the given scale on all eight generators; cartan restricts it to λ₃
// and λ₈, the commuting U(1)×U(1) part, which is the abelian case.
func RandomAlgebra(r *rand.Rand, scale float64, cartan bool) SU3 {
	var c [8]float64
	for a := range c {
		if !cartan || a == 2 || a == 7 {
			c[a] = scale * r.NormFloat64()
		}
	}
	return Algebra(c)
}

// RandomSU3 is exp(i·h) for a random algebra element (see RandomAlgebra).
func RandomSU3(r *rand.Rand, scale float64, cartan bool) SU3 {
	return Exp(RandomAlgebra(r, scale, cartan))
}

// SU3Field is an SU(3) link variable on every link of a lattice.
type SU3Field struct {
	Lat *Lattice
	U   []SU3
}

// NewSU3Field makes the trivial field (every link the identity).
func NewSU3Field(l *Lattice) *SU3Field {
	f := &SU3Field{Lat: l, U: make([]SU3, l.Links())}
	for k := range f.U {
		f.U[k] = Identity3()
	}
	return f
}

// Randomize sets every link to a random group element (RandomSU3).
func (f *SU3Field) Randomize(r *rand.Rand, scale float64, cartan bool) {
	for k := range f.U {
		f.U[k] = RandomSU3(r, scale, cartan)
	}
}

// Gauge applies the gauge transformation g (one group element per node):
// U_k → g_from·U_k·g_to† (Wilson 1974; Kogut & Susskind 1975 §II).
func (f *SU3Field) Gauge(g []SU3) {
	for k := range f.U {
		f.U[k] = g[f.Lat.From[k]].Mul(f.U[k]).Mul(g[f.Lat.To[k]].Dagger())
	}
}

// Transport is the path-ordered product along a path: the parallel
// transporter from the path's first node to its last, with U† for a link
// traversed against its orientation.
func (f *SU3Field) Transport(path []Step) SU3 {
	w := Identity3()
	for _, s := range path {
		if s.Dir > 0 {
			w = w.Mul(f.U[s.Link])
		} else {
			w = w.Mul(f.U[s.Link].Dagger())
		}
	}
	return w
}

// WilsonLoop is the holonomy of plaquette p, based at the plaquette's
// first node. Its trace is gauge invariant; the matrix itself transforms
// by conjugation with g at the base node.
func (f *SU3Field) WilsonLoop(p int) SU3 { return f.Transport(f.Lat.Plaquettes[p]) }

// Gauss is the non-abelian Gauss generator at every node,
//
//	G_n = Σ_{from=n} E_k − Σ_{to=n} U_k†·E_k·U_k − ρ_n,
//
// for left electric fields E_k (algebra elements at each link's From
// end) and charges ρ_n (nil: none). It transforms covariantly, G_n →
// g_n·G_n·g_n†; for diagonal E, U and ρ it reduces to the U(1) generator
// Div Q − ρ on each diagonal entry.
func (f *SU3Field) Gauss(e []SU3, rho []SU3) []SU3 {
	g := make([]SU3, f.Lat.Nodes)
	for k := range f.U {
		g[f.Lat.From[k]] = g[f.Lat.From[k]].Add(e[k])
		g[f.Lat.To[k]] = g[f.Lat.To[k]].Sub(f.U[k].Dagger().Mul(e[k]).Mul(f.U[k]))
	}
	if rho != nil {
		for n := range g {
			g[n] = g[n].Sub(rho[n])
		}
	}
	return g
}

// Compose returns the holonomies of plaquettes p and q based at their
// common node n and the composition defect ‖W_p·W_q − W_q·W_p‖/3 (the
// Frobenius norm of the commutator over that of the identity). The
// defect is zero for every abelian field, so a circuit quantity that
// composes as these loops do would have to show a non-zero defect;
// Kirchhoff's sums show none. ok is false if the plaquettes do not share
// node n.
func (f *SU3Field) Compose(p, q, n int) (wp, wq SU3, defect float64, ok bool) {
	pp, ok1 := f.Lat.Rebase(f.Lat.Plaquettes[p], n)
	pq, ok2 := f.Lat.Rebase(f.Lat.Plaquettes[q], n)
	if !ok1 || !ok2 {
		return SU3{}, SU3{}, 0, false
	}
	wp, wq = f.Transport(pp), f.Transport(pq)
	return wp, wq, wp.Commutator(wq).Norm() / 3, true
}
