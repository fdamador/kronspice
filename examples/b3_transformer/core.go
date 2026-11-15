// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"fmt"
	"math"
)

const mu0 = 4e-7 * math.Pi

// Material is a declared core material (linear, small signal; no
// hysteresis, no eddy currents).
type Material struct {
	Name   string
	Mur    float64 // relative permeability
	Rho    float64 // resistivity [Ω·m]
	Eddy   float64 // half-thickness of the conducting part (wire radius, or half a tape) [m]; 0 = the wire radius
	Bsat   float64 // saturation flux density [T]
	Source string
}

// Materials are the three declared cores. Galfenol's μr is the value
// declared in examples/magnetostrictive_transformer (Etrema: 75–100);
// silicon steel and nanocrystalline tape are typical datasheet values.
var Materials = []Material{
	{"Galfenol wire", 87, 6e-7, 0, 1.5, "Etrema Galfenol data: μr 75–100, Bs 1.5 T (as in magnetostrictive_transformer); ρ ≈ 0.6 µΩ·m declared"},
	{"silicon-steel wire", 4000, 4.7e-7, 0, 1.6, "non-oriented 3 % Si electrical steel: μr ≈ 2000–8000, ρ 0.47 µΩ·m, Bs ≈ 1.6 T (declared)"},
	{"nanocrystalline tape", 1e5, 1.15e-6, 9e-6, 1.2, "Fe-based nanocrystalline tape (Finemet type), 18 µm: μr ≈ 5e4–1.5e5 at low frequency, ρ 1.15 µΩ·m, Bs 1.2 T (declared)"},
}

// SkinDepth is √(ρ/(π f μ)) in the core material at f.
func (m Material) SkinDepth(f float64) float64 {
	return math.Sqrt(m.Rho / (math.Pi * f * mu0 * m.Mur))
}

// EddyOK reports whether eddy currents may be left out at f: the skin
// depth is at least the conducting half-thickness.
func (s Spring) EddyOK(f float64) bool {
	h := s.Mat.Eddy
	if h == 0 {
		h = s.A
	}
	return s.Mat.SkinDepth(f) >= h
}

// Spring is the N-turn helical core closed by a straight return leg.
type Spring struct {
	N     int     // core turns
	R     float64 // helix radius (to the wire centre) [m]
	A     float64 // wire radius [m]
	P     float64 // pitch [m]
	Reach float64 // radial reach of the return leg beyond the helix [m]
	Mat   Material
	S     int // segments per core turn
	// RetArea is the return leg's cross-section as a multiple of the
	// wire's (a thick return keeps neighbouring turns at nearly the same
	// magnetic potential, which is what removes the gap leakage).
	RetArea float64
}

// DefaultSpring is the design: 10 turns, Ø20 mm, Ø4 mm core wire,
// 12 mm pitch (8 mm gaps), return leg 16× the wire's cross-section
// (Ø16 mm). The thick core, wide gaps and thick return are what keep the
// flux off the gaps (section 1 of the run).
func DefaultSpring(m Material) Spring {
	return Spring{N: 10, R: 10e-3, A: 2e-3, P: 12e-3, Reach: 10e-3, Mat: m, S: 24, RetArea: 16}
}

// TurnLength is the length of core wire in one turn.
func (s Spring) TurnLength() float64 { return math.Hypot(2*math.Pi*s.R, s.P) }

// Area is the core wire's cross-section.
func (s Spring) Area() float64 { return math.Pi * s.A * s.A }

// SegReluctance is one segment's reluctance [A/Wb].
func (s Spring) SegReluctance() float64 {
	return s.TurnLength() / float64(s.S) / (mu0 * s.Mat.Mur * s.Area())
}

// ReturnLength is the core length of the return leg: up the spring's
// height and out and back across the reach.
func (s Spring) ReturnLength() float64 { return float64(s.N)*s.P + 2*s.Reach }

// ReturnReluctance is the return leg's reluctance.
func (s Spring) ReturnReluctance() float64 {
	ra := s.RetArea
	if ra <= 0 {
		ra = 1
	}
	return s.ReturnLength() / (mu0 * s.Mat.Mur * s.Area() * ra)
}

// LeakPermeance is the air permeance between the same azimuth on two
// neighbouring turns, per segment: two parallel cylinders of radius a at
// centre distance p carry μ0·π/acosh(p/2a) per unit length (the
// two-wire capacitance formula with ε0 → μ0; the core surface is
// nearly a magnetic equipotential when μr ≫ 1).
func (s Spring) LeakPermeance() float64 {
	return mu0 * math.Pi / math.Acosh(s.P/(2*s.A)) * (2 * math.Pi * s.R / float64(s.S))
}

// HelixNumber is leakage permeance per turn × reluctance per turn,
// (2πR)²/(μr·a²·acosh(p/2a)). Flux follows the wire round the helix only
// when it is ≪ 1; above 1 it short-cuts axially between neighbouring
// turns.
func (s Spring) HelixNumber() float64 {
	return s.LeakPermeance() * float64(s.S) * s.SegReluctance() * float64(s.S)
}

// Winding is a coil wrapped round the spring's wall at one azimuth, each of
// its Turns enclosing the core passes First … First+Passes−1 (or, for the
// axis pickup, a straight wire up the bore closed outside, Turns = 1,
// enclosing every pass at its azimuth).
type Winding struct {
	Name    string
	Azimuth float64 // [rad]; the helix starts at 90° (where the return leg joins)
	First   int
	Passes  int
	Turns   int
}

// segments returns the core segment crossing the winding's azimuth in each
// enclosed pass.
func (s Spring) segments(w Winding) ([]int, error) {
	if w.Passes < 1 || w.First < 0 || w.First+w.Passes > s.N {
		return nil, fmt.Errorf("winding %s: passes %d…%d outside the %d-turn core", w.Name, w.First, w.First+w.Passes-1, s.N)
	}
	off := math.Mod(w.Azimuth-math.Pi/2, 2*math.Pi)
	if off < 0 {
		off += 2 * math.Pi
	}
	o := int(off / (2 * math.Pi) * float64(s.S))
	seg := make([]int, w.Passes)
	for i := range seg {
		seg[i] = (w.First+i)*s.S + o
	}
	return seg, nil
}

// Field is the magnetic network's solution for one set of MMFs.
type Field struct {
	Seg    []float64 // flux in each core segment along the wire [Wb]
	Return float64   // flux down the return leg [Wb]
	Leak   []float64 // flux from node j to node j+S through the air [Wb]
}

// Solve finds the fluxes for the MMF mmf[j] (A-turns) in core segment j,
// acting along the wire. Nodes are the magnetic potentials at segment
// boundaries 0…M (M = N·S); node 0 is the reference; segments join j to
// j+1, leakage permeances j to j+S, the return leg M to 0. The nodal
// matrix is a symmetric positive-definite band of half-width S, solved by
// band Cholesky.
func (s Spring) Solve(mmf []float64) (Field, error) {
	m := s.N * s.S
	if len(mmf) != m {
		return Field{}, errors.New("mmf length")
	}
	gs := 1 / s.SegReluctance()
	gl := s.LeakPermeance()
	gr := 1 / s.ReturnReluctance()
	bw := s.S
	n := m // unknowns u_1…u_M, index i = node−1
	band := make([][]float64, n)
	for i := range band {
		band[i] = make([]float64, bw+1) // band[i][d] = A[i][i+d]
	}
	rhs := make([]float64, n)
	add := func(a, b int, g float64) { // conductance between nodes a, b (0 = reference)
		if a > 0 {
			band[a-1][0] += g
		}
		if b > 0 {
			band[b-1][0] += g
		}
		if a > 0 && b > 0 {
			i, j := a-1, b-1
			if i > j {
				i, j = j, i
			}
			band[i][j-i] -= g
		}
	}
	for j := 0; j < m; j++ {
		add(j, j+1, gs)
		// MMF F in the segment drives flux F·gs from j to j+1:
		// a source injecting F·gs into node j+1 and drawing it from j.
		if j+1 > 0 {
			rhs[j] += mmf[j] * gs
		}
		if j > 0 {
			rhs[j-1] -= mmf[j] * gs
		}
	}
	for j := 0; j+s.S <= m; j++ {
		add(j, j+s.S, gl)
	}
	add(m, 0, gr)
	// Band Cholesky: A = Lᵀ-stored upper factor U, A = UᵀU.
	for i := 0; i < n; i++ {
		d := band[i][0]
		if !(d > 0) {
			return Field{}, fmt.Errorf("magnetic network not positive definite at node %d", i+1)
		}
		d = math.Sqrt(d)
		band[i][0] = d
		for k := 1; k <= bw && i+k < n; k++ {
			band[i][k] /= d
		}
		for k := 1; k <= bw && i+k < n; k++ {
			uik := band[i][k]
			if uik == 0 {
				continue
			}
			for l := k; l <= bw && i+l < n; l++ {
				band[i+k][l-k] -= uik * band[i][l]
			}
		}
	}
	y := make([]float64, n) // Uᵀy = rhs
	for i := 0; i < n; i++ {
		v := rhs[i]
		for k := 1; k <= bw && i-k >= 0; k++ {
			v -= band[i-k][k] * y[i-k]
		}
		y[i] = v / band[i][0]
	}
	u := make([]float64, n+1) // node potentials, u[0] = 0
	for i := n - 1; i >= 0; i-- {
		v := y[i]
		for k := 1; k <= bw && i+k < n; k++ {
			v -= band[i][k] * u[i+k+1]
		}
		u[i+1] = v / band[i][0]
	}
	f := Field{Seg: make([]float64, m), Leak: make([]float64, m-s.S+1)}
	for j := 0; j < m; j++ {
		f.Seg[j] = (u[j] - u[j+1] + mmf[j]) * gs
	}
	for j := range f.Leak {
		f.Leak[j] = (u[j] - u[j+s.S]) * gl
	}
	f.Return = (u[m] - u[0]) * gr
	return f, nil
}

// unitField is the field of 1 A in winding w.
func (s Spring) unitField(w Winding) (Field, error) {
	seg, err := s.segments(w)
	if err != nil {
		return Field{}, err
	}
	mmf := make([]float64, s.N*s.S)
	for _, j := range seg {
		mmf[j] += float64(w.Turns)
	}
	return s.Solve(mmf)
}

// Linkage is winding v's flux linkage in field f: Turns × Σ (flux in each
// enclosed pass).
func (s Spring) Linkage(v Winding, f Field) (float64, error) {
	seg, err := s.segments(v)
	if err != nil {
		return 0, err
	}
	l := 0.0
	for _, j := range seg {
		l += f.Seg[j]
	}
	return float64(v.Turns) * l, nil
}

// Inductance is the inductance matrix of the windings from the magnetic
// network (core and turn-to-turn leakage only).
func (s Spring) Inductance(ws []Winding) ([][]float64, []Field, error) {
	L := make([][]float64, len(ws))
	fs := make([]Field, len(ws))
	for j, w := range ws {
		f, err := s.unitField(w)
		if err != nil {
			return nil, nil, err
		}
		fs[j] = f
		L[j] = make([]float64, len(ws))
	}
	for i, v := range ws {
		for j := range ws {
			l, err := s.Linkage(v, fs[j])
			if err != nil {
				return nil, nil, err
			}
			L[i][j] = l
		}
	}
	return L, fs, nil
}

// PassFlux is the flux in each core pass at azimuth 0 (one value per turn).
func (s Spring) PassFlux(f Field) []float64 {
	seg, _ := s.segments(Winding{Azimuth: 0, First: 0, Passes: s.N, Turns: 1})
	out := make([]float64, len(seg))
	for i, j := range seg {
		out[i] = f.Seg[j]
	}
	return out
}

// Ratio is the open-circuit voltage ratio v_out/v_in with `in` driven:
// L_out,in / L_in,in (dΦ/dt is common to both, so the ratio of linkages).
func (s Spring) Ratio(in, out Winding) (float64, error) {
	f, err := s.unitField(in)
	if err != nil {
		return 0, err
	}
	li, err := s.Linkage(in, f)
	if err != nil {
		return 0, err
	}
	lo, err := s.Linkage(out, f)
	if err != nil {
		return 0, err
	}
	return lo / li, nil
}

// IdealRatio is (N_out·n_out)/(N_in·n_in): every pass carrying the same flux.
func IdealRatio(in, out Winding) float64 {
	return float64(out.Turns*out.Passes) / float64(in.Turns*in.Passes)
}

// Leak is the leakage summary of a spring with the design windings.
type Leak struct {
	Spread float64 // (max−min)/max of pass flux, 10-pass primary driven
	K      float64 // coupling of the 10-pass and 2-pass windings
	Swap   float64 // v(10-pass)/v(2-pass) with the 2-pass winding driven (ideal 5)
}

// Leakage measures both gap leaks for the design windings.
func (s Spring) Leakage() (Leak, error) {
	pri, sec := designWindings(s.N)
	L, fs, err := s.Inductance([]Winding{pri, sec})
	if err != nil {
		return Leak{}, err
	}
	pf := s.PassFlux(fs[0])
	mn, mx := pf[0], pf[0]
	for _, v := range pf {
		mn, mx = math.Min(mn, v), math.Max(mx, v)
	}
	return Leak{Spread: (mx - mn) / mx, K: L[0][1] / math.Sqrt(L[0][0]*L[1][1]), Swap: L[0][1] / L[1][1]}, nil
}
