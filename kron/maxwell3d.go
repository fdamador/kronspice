// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// Kron's three-dimensional equivalent circuit of Maxwell's equations,
// with its ideal transformers, as a kronspice netlist (stage B of
// doc/kron-figure-catalogue.md §10; H20 step 1, CLAUDE.md §1).
//
// Source: G. Kron, J. Math. Anal. Appl. 5, 46–69, 1962, Fig. 1 (p. 49),
// the 1944 Proc. IRE network with every value written out; its working
// is N. V. Balasubramanian, Liverpool thesis 1968, pp. 153–158 and
// Figs. 10–11 (read from the page images).
//
// The drawing, as built here:
//
//   - At every lattice node a "sphere": one attachment point for each
//     lattice edge that meets the node (up to six), and between every two
//     attachment points on perpendicular axes an arc carrying a quarter
//     inductor (up to twelve arcs, one for each face that has a corner at
//     the node).
//   - Along every lattice edge a dielectric branch, capacitor ∥
//     conductance, from the attachment point on one sphere to the
//     attachment point on the next: C = ε·(dual face)/(edge),
//     G = σ·(dual face)/(edge). Its voltage is E·Δ, its current
//     (J + Ḋ)·(dual face).
//   - Every face is a "large mesh" (Balasubramanian's ABCDEFGH): its four
//     dielectric branches and the four arcs, one in each corner sphere,
//     joining them. The four quarter inductors each take ¼ of the face's
//     flux rate, L/4 with L = μ·(face)/(dual edge).
//   - Ideal 1:1 transformers (L = M = ∞) in the arcs, poled so that round
//     the large mesh each primary voltage cancels its secondary. They are
//     invisible to the large mesh and offer infinite impedance to any
//     current in the "small" meshes (loops of arcs within the spheres), so
//     every arc of a face carries the same current H·Δ (Balasubramanian
//     p. 157).
//
// Then Kirchhoff's voltage law round a large mesh is Faraday's law,
// Σ E·Δ = −L·dI/dt = −dΦ/dt, and Kirchhoff's current law at an attachment
// point is Ampère's: the dielectric branch's current is the sum of the
// four face currents round its edge, *(J + Ḋ) = dH. These are exactly the
// space-discrete equations of Yee's scheme (em/fdtd.go) and of the U(1)
// lattice field of gauge.go on the same grid.
//
// How the ties are drawn. Balasubramanian (p. 153) says each large mesh
// holds two ideal transformers, primaries and secondaries "placed
// diagonally opposite": arcs at corners 0 and 2 tied, 1 and 3 tied
// (KronTies). Within one face that forces only two pairs of equal
// currents, so ChainTies, which adds the third tie (corners 0–1, 1–2,
// 2–3), is built beside it. TestMaxwellNetworkIsYee finds both exact in
// the limit L = M = ∞: in Kron's arrangement the lattice of spheres
// fixes the remaining arc currents, as his drawing implies.
//
// Ideal transformer. Kron writes it as two coils with L = M = ∞ (Fig. 1,
// "L = M = ∞"), and so it is built: two windings of inductance L_w,
// coupled with k = 1 (kronspice K element), the secondary poled against
// the circulation. Then, with I_c and I_c′ the two arcs' currents along
// the circulation, both winding voltages are jωL_w·(I_c − I_c′): they
// cancel round the large mesh at any L_w, and the arcs' current
// difference (the magnetizing current) is that voltage over jωL_w, which
// vanishes as L_w → ∞. L_w = Ideal × (L/4); the departure from Maxwell's
// equations falls as 1/Ideal (TestMaxwellNetworkIsYee). The windings
// are coils, so the spheres stay galvanically whole, as in the drawing;
// a realization by gyrators would cut every link into an island whose
// potential no equation fixes.
//
// Conducting walls: edges lying in a wall carry no tangential E; their
// dielectric branch is replaced by a 0 V source (the wall short), and
// faces lying entirely in a wall (normal B = 0 there) are left out.

// TieRule selects how the ideal transformers tie a face's four arcs.
type TieRule int

const (
	// KronTies (the default): two transformers per face, corners 0–2 and
	// 1–3, the arrangement Balasubramanian describes for Kron's drawing
	// (p. 153).
	KronTies TieRule = iota
	// ChainTies: three transformers per face, corners 0–1, 1–2, 2–3, so
	// each face ties all four of its quarter inductors itself.
	ChainTies
)

// MaxwellNetwork is Kron's 1944/1962 three-dimensional field network on
// a YeeGrid (conducting walls on the box faces).
type MaxwellNetwork struct {
	Grid  *YeeGrid
	Eps   float64 // permittivity [F/m]
	Mu    float64 // permeability [H/m]
	Sigma float64 // conductivity [S/m]
	Ties  TieRule
	// Ideal is the winding inductance of the ideal transformers in units
	// of the quarter inductance (default 1e7): Kron's L = M = ∞ as a limit.
	// The departure from Yee falls as 1/Ideal, but the attachment points of
	// different links are tied to one another only through the windings'
	// magnetizing inductance, so the network's matrix has a near-null space
	// of strength ~1/Ideal: in AC the admittance 1/(jωL_w), in a transient
	// step the conductance h/(2L_w). Transient runs need a smaller Ideal
	// (1e3 for steps of a fifth of a cell's transit time) for the LU to
	// keep those pivots.
	Ideal float64
	// Leakage, when > 0, couples the windings with k = 1 − Leakage instead
	// of 1, and takes each winding's leakage inductance Leakage·L_w off the
	// quarter inductor of its arc, so the network is unchanged (the
	// T-equivalent of a coupled pair is the leakage in series with each
	// winding, an ideal 1:1 transformer and the magnetizing inductance
	// k·L_w across it). With k = 1 the two branch rows of a pair are equal
	// in their current columns and the time-step matrix, whose entries are
	// 2L_w/h, can cancel one of them to round-off; Leakage of 1e-5 keeps the
	// pair invertible. It must leave every quarter inductor positive.
	Leakage float64
	// The DC operating point kronspice computes before an AC or transient
	// run needs paths Kron's lossless network does not have: a small
	// Sigma joins the spheres at DC (MaxwellYee includes it exactly);
	// SeriesR, when > 0, is a resistor in series with every quarter
	// inductor (a face loop of inductors and transformers is otherwise a
	// superconducting ring, its DC current undetermined; MaxwellYee
	// includes it exactly); Leak, when > 0, is a resistor to ground from
	// the first free link's start point, the network's reference (no current flows in
	// it: nothing else touches ground).
	Leak    float64
	SeriesR float64
}

// Conductor marks a perfect conductor inside the grid: every link whose
// two end nodes satisfy in(i, j, k) is held at E·Δ = 0, as the wall links
// are (the inner conductor of a coaxial line). Faces lying wholly in the
// conductor then drop out of a MaxwellNetwork.
func (g *YeeGrid) Conductor(in func(i, j, k int) bool) {
	at := func(n int) bool { return in(n%g.N[0], (n/g.N[0])%g.N[1], n/(g.N[0]*g.N[1])) }
	for k := range g.Free {
		if at(g.From[k]) && at(g.To[k]) {
			g.Free[k] = false
		}
	}
}

// Point is the node name of the attachment point of link k on the sphere
// at its start (end = 0) or its end (end = 1).
func (m *MaxwellNetwork) Point(k, end int) string { return fmt.Sprintf("p%d_%d", k, end) }

// Branch is the element name of link k's dielectric branch (its
// capacitor; the conductance is "G" + the same suffix, the wall short
// "VW" + it).
func (m *MaxwellNetwork) Branch(k int) string { return fmt.Sprintf("C%d", k) }

// Arc is the element name of the quarter inductor at corner c (0–3) of
// face p; its current, from the attachment point of the face's step c to
// that of step c+1, is the face's mesh current H·Δ in the face's own
// circulation.
func (m *MaxwellNetwork) Arc(p, c int) string { return fmt.Sprintf("L%d_%d", p, c) }

// Faces reports whether face p is in the network (not lying in a wall).
func (m *MaxwellNetwork) Faces(p int) bool {
	for _, s := range m.Grid.Plaquettes[p] {
		if m.Grid.free(s.Link) {
			return true
		}
	}
	return false
}

// Constants returns the capacitance of every link and the full inductance
// of every face (each quarter inductor is a quarter of it).
func (m *MaxwellNetwork) Constants() (c, l []float64) { return m.Grid.Constitutive(m.Eps, m.Mu) }

// Build writes the network into a builder.
func (m *MaxwellNetwork) Build(b *circuit.Builder) error {
	g := m.Grid
	if g == nil || !(m.Eps > 0 && m.Mu > 0) || m.Sigma < 0 || math.IsNaN(m.Eps+m.Mu+m.Sigma) {
		return fmt.Errorf("%w: Maxwell network ε %g μ %g σ %g", ErrPolyhedron, m.Eps, m.Mu, m.Sigma)
	}
	if m.Ties != ChainTies && m.Ties != KronTies {
		return fmt.Errorf("%w: tie rule %d", ErrPolyhedron, m.Ties)
	}
	ideal := m.Ideal
	if ideal == 0 {
		ideal = 1e7
	}
	if !(m.Leakage >= 0 && m.Leakage < 1) {
		return fmt.Errorf("%w: transformer leakage %g", ErrPolyhedron, m.Leakage)
	}
	if !(ideal > 0) {
		return fmt.Errorf("%w: ideal-transformer ratio %g", ErrPolyhedron, ideal)
	}
	capv, ind := m.Constants()
	use := func(n string) string { return n }
	// Dielectric branches and wall shorts (a wall link on no face of the
	// network, along an edge of the box, joins nothing and is left out).
	onFace := make([]bool, g.Links())
	for p, steps := range g.Plaquettes {
		if m.Faces(p) {
			for _, s := range steps {
				onFace[s.Link] = true
			}
		}
	}
	for k := 0; k < g.Links(); k++ {
		a, z := use(m.Point(k, 0)), use(m.Point(k, 1))
		if !g.free(k) {
			if !onFace[k] {
				continue
			}
			b.VoltageSource(fmt.Sprintf("VW%d", k), a, z, 0)
			continue
		}
		b.Capacitor(m.Branch(k), a, z, capv[k])
		if m.Sigma > 0 {
			b.Resistor(fmt.Sprintf("G%d", k), a, z, m.Eps/(m.Sigma*capv[k]))
		}
	}
	// Large meshes.
	for p, steps := range g.Plaquettes {
		if !m.Faces(p) {
			continue
		}
		// Corner c joins step c to step c+1 at the node where step c
		// ends; the arc runs from link(step c)'s attachment point there to
		// link(step c+1)'s.
		var from, to [4]string
		for c := 0; c < 4; c++ {
			s, t := steps[c], steps[(c+1)%4]
			from[c] = m.Point(s.Link, endIndex(s, true))
			to[c] = m.Point(t.Link, endIndex(t, false))
		}
		// Windings per corner: +1 along the circulation (primary), −1
		// against it (secondary); each entry names its transformer.
		type winding struct {
			tr   int
			sign int
		}
		var w [4][]winding
		switch m.Ties {
		case ChainTies:
			for t := 0; t < 3; t++ {
				w[t] = append(w[t], winding{t, 1})
				w[t+1] = append(w[t+1], winding{t, -1})
			}
		case KronTies:
			w[0] = append(w[0], winding{0, 1})
			w[2] = append(w[2], winding{0, -1})
			w[1] = append(w[1], winding{1, 1})
			w[3] = append(w[3], winding{1, -1})
		}
		ports := map[int][2][2]string{} // transformer → [primary, secondary] (p, n)
		for c := 0; c < 4; c++ {
			x := use(fmt.Sprintf("f%d_%d_l", p, c))
			lq := ind[p]/4 - float64(len(w[c]))*m.Leakage*ideal*ind[p]/4
			if !(lq > 0) {
				return fmt.Errorf("%w: leakage %g × Ideal %g leaves no quarter inductance", ErrPolyhedron, m.Leakage, ideal)
			}
			b.Inductor(m.Arc(p, c), from[c], x, lq)
			if m.SeriesR > 0 {
				y := use(fmt.Sprintf("f%d_%d_r", p, c))
				b.Resistor(fmt.Sprintf("RL%d_%d", p, c), x, y, m.SeriesR)
				x = y
			}
			for i, wd := range w[c] {
				next := to[c]
				if i < len(w[c])-1 {
					next = use(fmt.Sprintf("f%d_%d_w%d", p, c, i))
				}
				pr := ports[wd.tr]
				if wd.sign > 0 {
					pr[0] = [2]string{x, next} // along the circulation: p = entry
				} else {
					pr[1] = [2]string{next, x} // against it: p = exit
				}
				ports[wd.tr] = pr
				x = next
			}
			if len(w[c]) == 0 {
				// no winding at this corner: the arc ends at the point
				b.VoltageSource(fmt.Sprintf("VA%d_%d", p, c), x, to[c], 0)
			}
		}
		lw := ideal * ind[p] / 4
		for t := 0; t < len(ports); t++ {
			pr := ports[t]
			wa, wb := fmt.Sprintf("WA%d_%d", p, t), fmt.Sprintf("WB%d_%d", p, t)
			b.Inductor(wa, pr[0][0], pr[0][1], lw)
			b.Inductor(wb, pr[1][0], pr[1][1], lw)
			b.Mutual(fmt.Sprintf("KT%d_%d", p, t), wa, wb, 1-m.Leakage)
		}
	}
	if m.Leak > 0 {
		for k := 0; k < g.Links(); k++ {
			if g.free(k) {
				b.Resistor("RK", m.Point(k, 0), "0", m.Leak)
				break
			}
		}
	}
	return nil
}

// endIndex is the end of a step's link at the step's start (first) or
// end: a step with Dir = +1 runs from the link's start (0) to its end
// (1).
func endIndex(s Step, last bool) int {
	if (s.Dir > 0) == last {
		return 1
	}
	return 0
}

// MaxwellYee solves the space-discrete Yee equations that the network
// stands for, at angular frequency ω, with impressed edge currents j
// (one per link, along the link; zero on wall links):
//
//	(G_e + jωC_e)·V_e = Σ_p K_pe·I_p − j_e     (Ampère, per free link)
//	(jωL_p + 4R)·I_p = −Σ_e K_pe·V_e            (Faraday, per face)
//
// K the face–link incidence (±1 along the face's circulation), R the
// SeriesR of each quarter inductor. It returns the link voltages and
// face currents: the twin against which the circuit is checked.
func (m *MaxwellNetwork) MaxwellYee(omega float64, j []float64) (v, cur []complex128, err error) {
	g := m.Grid
	nl, np := g.Links(), len(g.Plaquettes)
	if len(j) != nl || !(omega > 0) {
		return nil, nil, fmt.Errorf("%w: Yee twin ω %g, %d sources for %d links", ErrPolyhedron, omega, len(j), nl)
	}
	capv, ind := m.Constants()
	idx := make([]int, nl)
	n := 0
	for k := range idx {
		idx[k] = -1
		if g.free(k) {
			idx[k] = n
			n++
		}
	}
	zf := make([]complex128, np)
	for p := range zf {
		zf[p] = complex(4*m.SeriesR, omega*ind[p])
	}
	a := make([][]complex128, n)
	for i := range a {
		a[i] = make([]complex128, n+1)
	}
	for k := 0; k < nl; k++ {
		if i := idx[k]; i >= 0 {
			a[i][i] += complex(m.Sigma*capv[k]/m.Eps, omega*capv[k])
			a[i][n] = complex(-j[k], 0)
		}
	}
	// Kᵀ·diag(1/Z)·K V moved to the left: (Y + KᵀZ⁻¹K)V = −j.
	for p, steps := range g.Plaquettes {
		if !m.Faces(p) {
			continue
		}
		for _, s := range steps {
			i := idx[s.Link]
			if i < 0 {
				continue
			}
			for _, t := range steps {
				if jj := idx[t.Link]; jj >= 0 {
					a[i][jj] += complex(float64(s.Dir*t.Dir), 0) / zf[p]
				}
			}
		}
	}
	x, err := solveComplex(a)
	if err != nil {
		return nil, nil, err
	}
	v = make([]complex128, nl)
	for k := range v {
		if i := idx[k]; i >= 0 {
			v[k] = x[i]
		}
	}
	cur = make([]complex128, np)
	for p, steps := range g.Plaquettes {
		if !m.Faces(p) {
			continue
		}
		var s complex128
		for _, st := range steps {
			s += complex(float64(st.Dir), 0) * v[st.Link]
		}
		cur[p] = -s / zf[p]
	}
	return v, cur, nil
}

// solveComplex solves the augmented system a (n × n+1) by Gaussian
// elimination with partial pivoting.
func solveComplex(a [][]complex128) ([]complex128, error) {
	n := len(a)
	for c := 0; c < n; c++ {
		piv, best := c, 0.0
		for r := c; r < n; r++ {
			if v := cmplxAbs(a[r][c]); v > best {
				piv, best = r, v
			}
		}
		if best == 0 {
			return nil, fmt.Errorf("%w: singular Yee system", ErrPolyhedron)
		}
		a[c], a[piv] = a[piv], a[c]
		for r := c + 1; r < n; r++ {
			f := a[r][c] / a[c][c]
			if f == 0 {
				continue
			}
			for k := c; k <= n; k++ {
				a[r][k] -= f * a[c][k]
			}
		}
	}
	x := make([]complex128, n)
	for r := n - 1; r >= 0; r-- {
		s := a[r][n]
		for k := r + 1; k < n; k++ {
			s -= a[r][k] * x[k]
		}
		x[r] = s / a[r][r]
	}
	return x, nil
}

func cmplxAbs(z complex128) float64 { return math.Hypot(real(z), imag(z)) }
