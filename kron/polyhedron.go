// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// Kron's equivalent circuit of Maxwell's field equations as a kronspice
// netlist (H20 step 1, CLAUDE.md §1; doc/kron-figure-catalogue.md §§1, 6).
//
// The cell is Fig. 1 of G. Kron, "Multidimensional curve-fitting with
// self-organizing automata," J. Math. Anal. Appl. 5, 46–69, 1962, p. 49
// (the 1944 Proc. IRE network with every element's impedance written
// out): along each axis a dielectric branch (capacitor) in parallel with
// a conductive branch (resistor), magnetic branches (inductors) round the
// faces, and ideal transformers at the nodes. In SI, for the branch along
// x of a cell Δx × Δy × Δz with scale factors h₁, h₂, h₃ (N. V.
// Balasubramanian, Liverpool thesis 1968, eqs. 5-3, 5-4 and p. 156):
//
//	C₁ = (h₂h₃/h₁)(ΔyΔz/Δx)·ε     G₁ = (h₂h₃/h₁)(ΔyΔz/Δx)·σ
//	L₂ = ¼ (h₃h₁/h₂)(ΔzΔx/Δy)·μ   (one quarter of the face flux on each
//	                               of the four inductors round a face)
//
// The gauge lattice of gauge.go (YeeGrid, U1Field) is the same cell in
// field variables; this file writes it as circuit elements so that the
// circuit engine, Kron's frames and the energy ledger solve it.
//
// Two-dimensional reduction (Balasubramanian Fig. 13, p. 165): when
// nothing varies along y, Maxwell's equations split into two independent
// sets (his eqs. 6-2) and the TE set is a plain lattice with no ideal
// transformers: at every node a shunt capacitor and conductance carrying
// Ḋ²ΔzΔx + J²ΔzΔx with the node voltage E₂Δy, and series inductors
// carrying H₃Δz along x (voltage Ḃ³ΔxΔy) and H₁Δx along z (voltage
// Ḃ¹ΔyΔz). Conducting walls are nodes held at zero (no tangential E).
// Kron's 1945 cavity (J. Appl. Phys. 16, 172, Fig. 14) is the dual
// assignment on the same lattice: inductor currents E_x, E_z and
// capacitor voltages H_y, with the inductors lying in a wall removed.

// ErrPolyhedron reports an invalid field-network specification.
var ErrPolyhedron = errors.New("kron: invalid field network")

// FieldCell is one cell of the lattice: its sides, scale factors (1 for
// Cartesian coordinates) and material.
type FieldCell struct {
	Dx, Dy, Dz float64
	H1, H2, H3 float64 // scale factors; zero means 1
	Eps, Mu    float64 // permittivity, permeability [F/m, H/m]
	Sigma      float64 // conductivity [S/m]
}

func (c FieldCell) h() (float64, float64, float64) {
	h := func(v float64) float64 {
		if v == 0 {
			return 1
		}
		return v
	}
	return h(c.H1), h(c.H2), h(c.H3)
}

func (c FieldCell) check() error {
	if !(c.Dx > 0 && c.Dy > 0 && c.Dz > 0) || !(c.Eps > 0 && c.Mu > 0) || c.Sigma < 0 || math.IsNaN(c.Dx+c.Dy+c.Dz+c.Eps+c.Mu+c.Sigma) {
		return fmt.Errorf("%w: cell %+v", ErrPolyhedron, c)
	}
	return nil
}

// Capacitance is the dielectric branch along axis a (0, 1, 2 for x, y,
// z), Kron 1962 Fig. 1 in SI.
func (c FieldCell) Capacitance(a int) float64 { return c.axial(a) * c.Eps }

// Conductance is the conductive branch along axis a.
func (c FieldCell) Conductance(a int) float64 { return c.axial(a) * c.Sigma }

// axial is (h_b h_c / h_a)(Δb Δc / Δa) for axis a.
func (c FieldCell) axial(a int) float64 {
	h1, h2, h3 := c.h()
	hh := [3]float64{h1, h2, h3}
	d := [3]float64{c.Dx, c.Dy, c.Dz}
	b, cc := (a+1)%3, (a+2)%3
	return hh[b] * hh[cc] / hh[a] * d[b] * d[cc] / d[a]
}

// FaceInductance is the magnetic branch of a face normal to axis a
// carrying the whole face flux (voltage Ḃ_a·Δb·Δc for the current H_a·Δa):
// (h_b h_c / h_a)(Δb Δc / Δa)·μ. Kron's three-dimensional cell puts a
// quarter of it on each of the four inductors round the face.
func (c FieldCell) FaceInductance(a int) float64 { return c.axial(a) * c.Mu }

// QuarterInductance is the inductor of Kron 1962 Fig. 1, one quarter of
// the face inductance.
func (c FieldCell) QuarterInductance(a int) float64 { return 0.25 * c.FaceInductance(a) }

// TELattice is Balasubramanian's Fig. 13: a rectangular TE field network
// of nx × nz cells in the x–z plane, Δy deep, with the field E_y as the
// node voltage E_y·Δy. Nodes are named "e_i_j" for i = 0..nx, j = 0..nz;
// the walls (i = 0, nx or j = 0, nz) are ground when Walls is true, so
// only interior nodes exist in the circuit. Drive and termination are
// added by the caller on the returned builder.
type TELattice struct {
	Nx, Nz int
	Cell   FieldCell
	Walls  bool
	// Lossless drops the conductances even when Sigma > 0 (for Kron's
	// lossless lattices). Shunt adds a resistor to ground at every node
	// and SeriesR a resistor in series with every inductor (0 means none):
	// numerical leaks that give a lossless lattice a DC operating point
	// for the MNA analyses (a lattice of pure inductors is a set of
	// zero-resistance loops, singular at DC); Kron's own AC solve
	// (Network.Solve) needs neither.
	Lossless bool
	Shunt    float64
	SeriesR  float64
	// HalfEndRow halves the shunt admittance of the last row (j = Nz),
	// its capacitors and conductances and the transverse inductors along
	// it: a node holds half its shunt from each neighbouring section, so a
	// lattice cut there and terminated keeps only its own half.
	HalfEndRow bool
}

// inductor adds L between a and b, through SeriesR when it is set.
func (l TELattice) inductor(b *circuit.Builder, name, from, to string, henries float64) {
	if l.SeriesR > 0 {
		mid := "m" + name
		b.Inductor(name, from, mid, henries)
		b.Resistor("R"+name, mid, to, l.SeriesR)
		return
	}
	b.Inductor(name, from, to, henries)
}

// Node is the circuit node name of lattice point (i, j), "0" on a wall.
func (l TELattice) Node(i, j int) string {
	if l.Walls && (i == 0 || i == l.Nx || j == 0 || j == l.Nz) {
		return "0"
	}
	return fmt.Sprintf("e_%d_%d", i, j)
}

// Build writes the lattice into a builder: shunt C (and G) at every
// non-wall node, series inductors along x (H_z·Δz, L = μΔxΔy/Δz) and
// along z (H_x·Δx, L = μΔyΔz/Δx) between neighbouring nodes.
func (l TELattice) Build(b *circuit.Builder) error {
	if l.Nx < 1 || l.Nz < 1 {
		return fmt.Errorf("%w: %d × %d cells", ErrPolyhedron, l.Nx, l.Nz)
	}
	if err := l.Cell.check(); err != nil {
		return err
	}
	c := l.Cell
	cap := c.Capacitance(1) // shunt along y: ε ΔzΔx/Δy
	g := c.Conductance(1)
	lx := c.FaceInductance(2) // link along x carries H_z·Δz: μ ΔxΔy/Δz
	lz := c.FaceInductance(0) // link along z carries H_x·Δx: μ ΔyΔz/Δx
	for i := 0; i <= l.Nx; i++ {
		for j := 0; j <= l.Nz; j++ {
			n := l.Node(i, j)
			if n == "0" {
				continue
			}
			share := 1.0
			if l.HalfEndRow && j == l.Nz {
				share = 0.5
			}
			b.Capacitor(fmt.Sprintf("C_%d_%d", i, j), n, "0", share*cap)
			if g > 0 && !l.Lossless {
				b.Resistor(fmt.Sprintf("G_%d_%d", i, j), n, "0", 1/(share*g))
			}
			if l.Shunt > 0 {
				b.Resistor(fmt.Sprintf("RS_%d_%d", i, j), n, "0", l.Shunt)
			}
		}
	}
	for i := 0; i <= l.Nx; i++ {
		for j := 0; j <= l.Nz; j++ {
			n := l.Node(i, j)
			if i < l.Nx {
				if m := l.Node(i+1, j); !(n == "0" && m == "0") {
					lxv := lx
					if l.HalfEndRow && j == l.Nz {
						lxv = 2 * lx // half the admittance: the row's share of the shunt
					}
					l.inductor(b, fmt.Sprintf("LX_%d_%d", i, j), n, m, lxv)
				}
			}
			if j < l.Nz {
				if m := l.Node(i, j+1); !(n == "0" && m == "0") {
					l.inductor(b, fmt.Sprintf("LZ_%d_%d", i, j), n, m, lz)
				}
			}
		}
	}
	return nil
}

// Dispersion returns the lattice's own resonance of the (p, q) mode of a
// closed box of nx × nz cells against the continuum value: on this
// lattice ω² = (4/(μεΔ²))·[sin²(pπΔx/2a)·(Δ/Δx)² + sin²(qπΔz/2b)·(Δ/Δz)²]
// (the spatial part of the Yee dispersion relation, which TestFDTDCavity
// and TestU1YeeCavity use), and in the continuum ω = c·π·√((p/a)² + (q/b)²).
func (l TELattice) Dispersion(p, q int) (lattice, continuum float64) {
	c := l.Cell
	a, bb := float64(l.Nx)*c.Dx, float64(l.Nz)*c.Dz
	v := 1 / math.Sqrt(c.Mu*c.Eps)
	sx := 2 / c.Dx * math.Sin(float64(p)*math.Pi*c.Dx/(2*a))
	sz := 2 / c.Dz * math.Sin(float64(q)*math.Pi*c.Dz/(2*bb))
	return v * math.Hypot(sx, sz), v * math.Pi * math.Hypot(float64(p)/a, float64(q)/bb)
}

// Kron1945Cavity is Fig. 14 of G. Kron, "Numerical solution of ordinary
// and partial differential equations by means of equivalent circuits,"
// J. Appl. Phys. 16, 172–186, 1945: "a square wave guide or cavity
// resonator, by a conducting surface that extends to infinity in the y
// direction, be divided into sixteen squares. The net of coils in each
// square consists of two inductances (z = 0.8 or y = 1.25) and a ground
// capacitor (z = 1.3 or y = 0.77). The currents in the inductance
// represent E_x and E_z, while the voltages across the capacitors
// represent H_y." The drawing (read at 260 dpi) has five rows and five
// columns of junctions, a capacitor to ground at every junction and a
// coil on every one of the forty links, horizontal coils carrying E_z
// (arrow to the left) and vertical coils E_x (arrow down); the conducting
// boundary on three sides is where the lattice stops. At the open (top)
// end "a sinusoidal current distribution in space, E_x, is impressed":
// vertical arrows into the five top junctions labelled 50 and 35.35
// pointing down on the left, 35.35 and 50 pointing up on the right and
// nothing at the centre, so the currents I₀·cos(πi/4) = 50, 35.35, 0,
// −35.35, −50 enter the top-row nodes. Solved at the analyzer frequency
// taken as ω = 1 rad/s, so L = 0.8 H and C = 1/1.3 F give Kron's
// reactances. Nodes are "h_i_j", i = 0..4 across, j = 0..4 down. leak, when
// > 0, is a resistor from every node to ground, a numerical path for the
// DC operating point (1e9 Ω against the capacitors' 0.77 S).
func Kron1945Cavity(leak float64) (*circuit.Circuit, error) {
	const n = 5
	b := circuit.NewBuilder("Kron 1945 Fig. 14: square cavity, Hy at the junctions")
	node := func(i, j int) string { return fmt.Sprintf("h_%d_%d", i, j) }
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			b.Capacitor(fmt.Sprintf("C_%d_%d", i, j), node(i, j), "0", 1/1.3)
			if leak > 0 {
				b.Resistor(fmt.Sprintf("RL_%d_%d", i, j), node(i, j), "0", leak)
			}
			if i < n-1 {
				b.Inductor(fmt.Sprintf("LZ_%d_%d", i, j), node(i, j), node(i+1, j), 0.8)
			}
			if j < n-1 {
				b.Inductor(fmt.Sprintf("LX_%d_%d", i, j), node(i, j), node(i, j+1), 0.8)
			}
		}
	}
	for i := 0; i < n; i++ {
		amp := 50 * math.Cos(math.Pi*float64(i)/4)
		if math.Abs(amp) < 1e-12 {
			continue
		}
		b.Add(circuit.Element{Name: fmt.Sprintf("IX_%d", i), Kind: circuit.CurrentSource, ACMag: amp}, "0", node(i, 0))
	}
	return b.Build()
}

// Kron1945Potentials are the junction potentials H_y of Fig. 14 as Kron
// read them on the network analyzer, rows top to bottom, columns left to
// right; Fig. 15 gives their unbalanced currents, the largest 7.93 A
// against a largest coil current of 102.8 A ("8 percent inaccuracy").
var Kron1945Potentials = [5][5]float64{
	{84.8, 53.68, 0.24, -53.28, -84},
	{25.66, 16.34, 0.24, -15.46, -25},
	{-44.72, -27.86, 0.44, 28.46, 46.2},
	{-100.2, -60.7, 0.48, 63.4, 101.8},
	{-130.2, -80.2, 0.7, 83.7, 133.7},
}

// Kron1945Residuals are the unbalanced currents i_u of Fig. 15 [A] at the
// same junctions.
var Kron1945Residuals = [5][5]float64{
	{2.4, 2.15, -0.02, -0.75, -2.2},
	{2.12, -4.36, 0.63, 3.41, 7.93},
	{5.26, 6.97, -0.21, 2.15, -6.1},
	{4.1, -2.87, 2.69, -0.1, 0.9},
	{-0.2, 1.17, 2.97, -2.2, 0.4},
}

// Waveguide is Balasubramanian's rectangular TE₁₀ waveguide computation
// (thesis §6-2–6-3, Figs. 15–16): a guide of width A and height Dy, the
// width divided into N cells of Δx = A/N and the axis into Sections of
// Δz, TE lattice with the side walls grounded, driven at z = 0 by current
// sources with the sin(πx/a) profile and terminated at the far end by a
// resistance R per node. The matching impedance is the R for which the
// voltage along the axis changes in phase only, and the theoretical value
// per section is η/(βc/ω)·Δy/Δx with β = √((ω/c)² − (π/a)²) (his eqs. 6-3
// to 6-5; at 2 GHz in the 11 cm guide βc/ω = 0.731955 and 2.3395×10⁴·Δy Ω
// per section at Δx = 2.2 cm).
type Waveguide struct {
	A, Dy, Dz float64
	N         int     // cells across the width
	Sections  int     // cells along the axis
	Freq      float64 // [Hz]
	Eps, Mu   float64
}

// Beta returns βc/ω for the TE₁₀ mode.
func (w Waveguide) Beta() float64 {
	k := 2 * math.Pi * w.Freq * math.Sqrt(w.Mu*w.Eps)
	kc := math.Pi / w.A
	if k <= kc {
		return 0
	}
	return math.Sqrt(k*k-kc*kc) / k
}

// Theory is the theoretical matching impedance per section (eq. 6-5).
func (w Waveguide) Theory() float64 {
	return math.Sqrt(w.Mu/w.Eps) / w.Beta() * w.Dy / (w.A / float64(w.N))
}

// LatticeImpedance is the lattice's own matching impedance per node for
// the TE₁₀ profile, from the dispersion of the ladder each column reduces
// to: series Z = jωL_z, shunt Y = jωC + 2(1 − cos(π/N))/(jωL_x), and the
// π-section iterative impedance Z_π = √(Z/Y)/√(1 + ZY/4) (the impedance
// seen at a node whose shunt is shared with the continuation). Real in
// the pass band.
func (w Waveguide) LatticeImpedance() float64 {
	om := 2 * math.Pi * w.Freq
	cell := FieldCell{Dx: w.A / float64(w.N), Dy: w.Dy, Dz: w.Dz, Eps: w.Eps, Mu: w.Mu}
	z := om * cell.FaceInductance(0)
	y := om*cell.Capacitance(1) - 2*(1-math.Cos(math.Pi/float64(w.N)))/(om*cell.FaceInductance(2))
	zy := -z * y
	return math.Sqrt(z/y) / math.Sqrt(1+zy/4)
}

// Circuit builds the driven lattice terminated by r at every interior
// node of the last row. The last row's capacitors are halved: a node
// holds half its shunt from each neighbouring section, so the end node
// keeps only the half that belongs to the lattice and the termination
// stands for the other half plus the continuation. (With the full
// capacitor no resistance can match the lattice, and the residual
// standing wave is what sets the "matching" found by minimum ripple.)
func (w Waveguide) Circuit(r float64) (*circuit.Circuit, *TELattice, error) {
	if w.N < 2 || w.Sections < 2 || w.Beta() == 0 {
		return nil, nil, fmt.Errorf("%w: waveguide %+v", ErrPolyhedron, w)
	}
	cell := FieldCell{Dx: w.A / float64(w.N), Dy: w.Dy, Dz: w.Dz, Eps: w.Eps, Mu: w.Mu}
	l := &TELattice{Nx: w.N, Nz: w.Sections, Cell: cell, Walls: false, Lossless: true, HalfEndRow: true}
	b := circuit.NewBuilder("TE10 waveguide on Kron's lattice")
	if err := l.Build(b); err != nil {
		return nil, nil, err
	}
	for i := 0; i <= w.N; i++ {
		for j := 0; j <= w.Sections; j++ {
			n := l.Node(i, j)
			if i == 0 || i == w.N {
				// side walls: E_y = 0
				b.VoltageSource(fmt.Sprintf("VW_%d_%d", i, j), n, "0", 0)
				continue
			}
			if j == 0 {
				amp := math.Sin(math.Pi * float64(i) / float64(w.N))
				b.Add(circuit.Element{Name: fmt.Sprintf("IS_%d", i), Kind: circuit.CurrentSource, ACMag: amp}, "0", n)
			}
			if j == w.Sections {
				b.Resistor(fmt.Sprintf("RT_%d", i), n, "0", r)
			}
		}
	}
	c, err := b.Build()
	return c, l, err
}
