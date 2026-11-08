// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 8 in the field solver (em/fdtd.go): the triplet m = +1, 0, −1 on
// the dielectric torus of Parts 5b and 7b (εr = 100, a = 10 mm,
// r = 2.5 mm, lossless). m = 0 is the E-ring (TE01δ); m = ±1 are the
// turning patterns of Part 7b, (cos φ ± i·sin φ)/√2 of the lowest odd
// lattice mode, the sin φ copy being the cos φ field turned 90° on the
// grid. Each state is a filtered phasor (Ê, B̂) at its own frequency
// (spinfield.go: FindMode), normalized to a time-averaged energy of 1 J,
//
//	Ū = ¼·Σ(ε|Ê|² + |B̂|²/μ0)·dV.
//
// From the phasors alone (no further time stepping):
//
//   - the overlaps ⟨a|b⟩ = ¼·Σ(ε Ê_a*·Ê_b + B̂_a*·B̂_b/μ0)·dV (8a);
//   - over one period, at phases θ: U_E, U_B, J_z (Abraham ε0∫ρ(E×B)_φ
//     and Minkowski ∫ρ(D×B)_φ, spinfield.go's integrals), the magnetic
//     dipole m(t) = Re(m̂·e^{iθ}) with m̂ = ½·Σ r × Ĵ·dV and
//     Ĵ = iω(ε − ε0)Ê the polarization current (the only current; the
//     torus has no conduction), and the electric dipole p̂ = Σ(ε − ε0)Ê·dV
//     of the neutral polarized body (8d);
//   - the field's moment of inertia per unit energy, ⟨ρ²⟩/c² with
//     ⟨ρ²⟩ = Σ ρ²·ū·dV / Σ ū·dV (8e).
//
// A single-frequency field has m(t) = Re(m̂·e^{iωt}): its time average is
// zero for every state, standing or turning. What differs is how the
// dipoles move: the E-ring (m = 0) has a magnetic dipole along the axis
// that passes through zero twice a period and no electric dipole; the
// turning m = ±1 patterns have no magnetic dipole but an in-plane
// electric dipole of constant size that turns with the pattern.

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// TState is one triplet state in the field solver.
type TState struct {
	Name     string
	M        int // +1, 0, −1 (labelled by the sign of the computed J_z)
	F        float64
	Ms       []Mode
	A        []complex128 // coefficients, normalized to Ū = 1 J
	MHat     [3]complex128
	PHat     [3]complex128
	Rho2     float64 // ⟨ρ²⟩ [m²]
	Theta    []float64
	UE, UB   []float64 // [J]
	JzA, JzM []float64 // [J·s]
	Mt, Pt   [][3]float64
	JzAMean  float64
	JzMMean  float64
	UBMin    float64 // smallest U_B/Ū over the period
	UBMax    float64
	MRange   [2]float64 // min, max over a period of |m(t)| / max|m|
	PRange   [2]float64 // the same for |p(t)|
}

// TripletFieldResult is Part 8 on one grid.
type TripletFieldResult struct {
	Cells   float64
	Dx      float64
	FE, F1  float64 // E-ring and m = 1 frequencies [Hz]
	States  [3]TState
	Gram    [3][3]complex128
	GramEig [3]float64
	Seconds float64
	grid    spinGrid
}

// tripletGridArrays holds ε per E edge for the grid.
func (g spinGrid) epsEdges(f *em.FDTD) [3][]float64 {
	var eps [3][]float64
	for c := 0; c < 3; c++ {
		eps[c] = make([]float64, len(f.E[c]))
		for k := 0; k < g.n; k++ {
			for j := 0; j < g.n; j++ {
				for i := 0; i < g.n; i++ {
					p := [3]float64{float64(i), float64(j), float64(k)}
					p[c] += 0.5
					er, _ := g.mat(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx)
					eps[c][g.idx(i, j, k)] = Eps0 * er
				}
			}
		}
	}
	return eps
}

// phasor returns the state's Ê, B̂.
func (s TState) phasor() (E, B [3][]complex128) {
	for c := 0; c < 3; c++ {
		E[c] = make([]complex128, len(s.Ms[0].E[c]))
		B[c] = make([]complex128, len(s.Ms[0].B[c]))
		for k, m := range s.Ms {
			for n := range E[c] {
				E[c][n] += s.A[k] * m.E[c][n]
			}
			for n := range B[c] {
				B[c][n] += s.A[k] * m.B[c][n]
			}
		}
	}
	return
}

// inner is ¼·Σ(ε Êa*·Êb + B̂a*·B̂b/μ0)·dV.
func (g spinGrid) inner(eps [3][]float64, ea, ba, eb, bb [3][]complex128) complex128 {
	var s complex128
	for c := 0; c < 3; c++ {
		for n := range ea[c] {
			s += complex(eps[c][n], 0) * cmplx.Conj(ea[c][n]) * eb[c][n]
		}
		for n := range ba[c] {
			s += cmplx.Conj(ba[c][n]) * bb[c][n] / complex(Mu0, 0)
		}
	}
	return s * complex(0.25*g.dx*g.dx*g.dx, 0)
}

// RunTripletField finds the E-ring and m = 1 modes on sf's grid (as Part
// 7b) and evaluates the triplet. nTheta phases per period.
func RunTripletField(sf SpinField, nTheta int, log func(string, ...any)) (TripletFieldResult, error) {
	sf.Log = log
	cells := sf.Ring.Cells
	g := sf.grid()
	res := TripletFieldResult{Cells: cells, Dx: g.dx, grid: g}
	mE, err := sf.FindMode("E-ring (TE01δ, m = 0)", "E", 1.0e9, 4e9)
	if err != nil {
		return res, err
	}
	m1, err := sf.FindMode("m = 1 (cos φ)", "m1", 1.0e9, 6e9)
	if err != nil {
		return res, err
	}
	m1s := g.Rot90(m1)
	res.FE, res.F1 = mE.F, m1.F
	f, err := g.fdtd()
	if err != nil {
		return res, err
	}
	eps := g.epsEdges(f)
	s2 := complex(math.Sqrt2/2, 0)
	i := complex(0, 1)
	sts := []TState{
		{Name: "m = 0 (E-ring, TE01δ)", F: mE.F, Ms: []Mode{mE}, A: []complex128{1}},
		{Name: "turning (cos φ + i sin φ)", F: m1.F, Ms: []Mode{m1, m1s}, A: []complex128{s2, i * s2}},
		{Name: "turning (cos φ − i sin φ)", F: m1.F, Ms: []Mode{m1, m1s}, A: []complex128{s2, -i * s2}},
	}
	dv := g.dx * g.dx * g.dx
	zero := [3][]float64{make([]float64, len(f.E[0])), make([]float64, len(f.E[1])), make([]float64, len(f.E[2]))}
	zeroB := [3][]float64{make([]float64, len(f.B[0])), make([]float64, len(f.B[1])), make([]float64, len(f.B[2]))}
	var phE, phB [3][3][]complex128
	for q := range sts {
		s := &sts[q]
		E, B := s.phasor()
		nrm := real(g.inner(eps, E, B, E, B))
		for k := range s.A {
			s.A[k] /= complex(math.Sqrt(nrm), 0)
		}
		E, B = s.phasor()
		phE[q], phB[q] = E, B
		w := 2 * math.Pi * s.F
		// dipoles and ⟨ρ²⟩
		var r2num, r2den float64
		for c := 0; c < 3; c++ {
			for k := 0; k < g.n; k++ {
				for j := 0; j < g.n; j++ {
					for ii := 0; ii < g.n; ii++ {
						n := g.idx(ii, j, k)
						p := [3]float64{float64(ii), float64(j), float64(k)}
						p[c] += 0.5
						x, y, z := -g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx
						ue := 0.25 * eps[c][n] * sq(cmplx.Abs(E[c][n]))
						r2num += (x*x + y*y) * ue
						r2den += ue
						if de := eps[c][n] - Eps0; de != 0 {
							pc := complex(de*dv, 0) * E[c][n]
							s.PHat[c] += pc
							// m̂ += ½·r × Ĵ dV, Ĵ = iω·P̂
							jc := complex(0, w) * pc
							r := [3]float64{x, y, z}
							a1, a2 := (c+1)%3, (c+2)%3
							s.MHat[a1] += 0.5 * complex(r[a2], 0) * jc
							s.MHat[a2] -= 0.5 * complex(r[a1], 0) * jc
						}
						pb := [3]float64{float64(ii) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
						pb[c] -= 0.5
						if n < len(B[c]) {
							xb, yb := -g.half+pb[0]*g.dx, -g.half+pb[1]*g.dx
							ub := 0.25 * sq(cmplx.Abs(B[c][n])) / Mu0
							r2num += (xb*xb + yb*yb) * ub
							r2den += ub
						}
					}
				}
			}
		}
		s.Rho2 = r2num / r2den
		// one period
		D := [3][]float64{make([]float64, len(f.E[0])), make([]float64, len(f.E[1])), make([]float64, len(f.E[2]))}
		s.UBMin, s.UBMax = math.Inf(1), math.Inf(-1)
		for k := 0; k < nTheta; k++ {
			th := 2 * math.Pi * float64(k) / float64(nTheta)
			ph := cmplx.Exp(complex(0, th))
			var ue, ub float64
			for c := 0; c < 3; c++ {
				for n := range f.E[c] {
					f.E[c][n] = real(E[c][n] * ph)
					D[c][n] = eps[c][n] * f.E[c][n]
					ue += 0.5 * eps[c][n] * f.E[c][n] * f.E[c][n] * dv
				}
				for n := range f.B[c] {
					f.B[c][n] = real(B[c][n] * ph)
					ub += 0.5 * f.B[c][n] * f.B[c][n] / Mu0 * dv
				}
			}
			jz, jzm, _, _, _ := g.integrals(f.E, f.B, D, zero, zeroB)
			s.Theta, s.UE, s.UB = append(s.Theta, th), append(s.UE, ue), append(s.UB, ub)
			s.JzA, s.JzM = append(s.JzA, jz), append(s.JzM, jzm)
			var mt, pt [3]float64
			for c := 0; c < 3; c++ {
				mt[c], pt[c] = real(s.MHat[c]*ph), real(s.PHat[c]*ph)
			}
			s.Mt, s.Pt = append(s.Mt, mt), append(s.Pt, pt)
			s.JzAMean += jz / float64(nTheta)
			s.JzMMean += jzm / float64(nTheta)
			s.UBMin, s.UBMax = math.Min(s.UBMin, ub), math.Max(s.UBMax, ub)
		}
		rng := func(v [][3]float64) [2]float64 {
			var mx float64
			for _, x := range v {
				mx = math.Max(mx, math.Sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))
			}
			r := [2]float64{math.Inf(1), 1}
			for _, x := range v {
				r[0] = math.Min(r[0], math.Sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/mx)
			}
			return r
		}
		s.MRange, s.PRange = rng(s.Mt), rng(s.Pt)
		if log != nil {
			log("%s: J_z Minkowski mean %.4g J·s, U_B/Ū [%.3f, %.3f]", s.Name, s.JzMMean, s.UBMin, s.UBMax)
		}
	}
	// label m by the sign of J_z (Minkowski) for the turning states
	sts[0].M = 0
	for q := 1; q < 3; q++ {
		if sts[q].JzMMean > 0 {
			sts[q].M = 1
		} else {
			sts[q].M = -1
		}
	}
	// order +1, 0, −1
	var plus, minus int
	if sts[1].M == 1 {
		plus, minus = 1, 2
	} else {
		plus, minus = 2, 1
	}
	order := [3]int{plus, 0, minus}
	for a, q := range order {
		res.States[a] = sts[q]
	}
	for a := 0; a < 3; a++ {
		for b := 0; b < 3; b++ {
			res.Gram[a][b] = g.inner(eps, phE[order[a]], phB[order[a]], phE[order[b]], phB[order[b]])
		}
	}
	res.GramEig = hermEig3(res.Gram)
	return res, nil
}

// MomentPerQuantum is the amplitude of the magnetic dipole when the state
// holds one quantum ħω: |m̂|·√(ħω/1 J) (the phasors are at Ū = 1 J), in
// A·m², split into the axial and in-plane parts.
func (s TState) MomentPerQuantum() (axial, perp float64) {
	k := math.Sqrt(Hbar * 2 * math.Pi * s.F)
	return cmplx.Abs(s.MHat[2]) * k, math.Hypot(cmplx.Abs(s.MHat[0]), cmplx.Abs(s.MHat[1])) * k
}

// DipolePerQuantum is the electric dipole amplitude per quantum [C·m].
func (s TState) DipolePerQuantum() float64 {
	k := math.Sqrt(Hbar * 2 * math.Pi * s.F)
	return math.Sqrt(sq(cmplx.Abs(s.PHat[0]))+sq(cmplx.Abs(s.PHat[1]))+sq(cmplx.Abs(s.PHat[2]))) * k
}

// Sample returns E and B of a state at phase θ at the given points
// (trilinear, em.FDTD.FieldAt).
func (r TripletFieldResult) Sample(q int, theta float64, pts [][3]float64) (Es, Bs [][3]float64, err error) {
	g := r.grid
	f, err := g.fdtd()
	if err != nil {
		return nil, nil, err
	}
	E, B := r.States[q].phasor()
	ph := cmplx.Exp(complex(0, theta))
	for c := 0; c < 3; c++ {
		for n := range f.E[c] {
			f.E[c][n] = real(E[c][n] * ph)
		}
		for n := range f.B[c] {
			f.B[c][n] = real(B[c][n] * ph)
		}
	}
	for _, p := range pts {
		e, b := f.FieldAt(p[0], p[1], p[2])
		Es, Bs = append(Es, e), append(Bs, b)
	}
	return
}
