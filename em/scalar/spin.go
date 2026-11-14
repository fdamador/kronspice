// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"math"

	"github.com/fdamador/kronspice/em"
)

// SpinConnection is ECE2's spin connection in vacuum, declared through its
// wave numbers κ₀(x), κ(x) [1/m] (Eckardt, *ECE Unified Field Theory* Vol. 1,
// 6.164–6.169; see ece.go). It acts on the field as the geometric current
//
//	J_ω = (2/μ₀)(κ₀E/c + κ×B)
//
// in Ampère's law: a conductivity σ = 2κ₀/(μ₀c), which drains the field for
// κ₀ > 0 and feeds it for κ₀ < 0, and a Hall-like term (2/μ₀)κ×B. Its work
// −∫J_ω·E is booked in Ledger.SpinWork, the hypothetical column: when the
// field gains energy the spin connection is the source booked for it, as a
// negative resistance is model-supplied in CLAUDE.md §3c. Hypothetical
// mode (H48); κ₀ = κ = 0 is Maxwell exactly.
type SpinConnection struct {
	Kappa0 func(x, y, z float64) float64    // κ₀ [1/m]; nil = 0
	Kappa  func(x, y, z float64) [3]float64 // κ [1/m]; nil = 0
	// Envelope scales both in time (a spin connection that is switched or
	// modulated); nil = 1.
	Envelope func(t float64) float64
}

type spinState struct {
	k0 [3][]float64    // κ₀ at each E edge, per component
	kv [3][][3]float64 // κ at each E edge
	g  func(t float64) float64
}

// setupSpin samples κ₀ and κ at the E edges.
func (r *Region) setupSpin(s *SpinConnection) {
	n := r.N[0] * r.N[1] * r.N[2]
	st := &spinState{g: s.Envelope}
	for c := 0; c < 3; c++ {
		st.k0[c] = make([]float64, n)
		st.kv[c] = make([][3]float64, n)
	}
	r.each(func(i, j, k, m int) {
		for c := 0; c < 3; c++ {
			if !r.eFree(c, i, j, k) {
				continue
			}
			off := [3]float64{}
			off[c] = 0.5
			x, y, z := r.pos(off, i, j, k)
			if s.Kappa0 != nil {
				st.k0[c][m] = s.Kappa0(x, y, z)
			}
			if s.Kappa != nil {
				st.kv[c][m] = s.Kappa(x, y, z)
			}
		}
	})
	r.spin = st
}

// bAtEdge is B^{n+½} interpolated to E edge (c, m): components c lie on the
// same line; the other two are averaged over the two faces that straddle
// the edge.
func (r *Region) bAtEdge(c, m int) [3]float64 {
	step := [3]int{1, r.N[0], r.N[0] * r.N[1]}
	var b [3]float64
	for d := 0; d < 3; d++ {
		if d == c {
			continue
		}
		e := 3 - c - d
		b[d] = 0.5 * (r.B[d][m] + r.B[d][m-step[e]])
	}
	return b
}

// applySpin replaces E* (the Ampère update without the spin current) by
// E^{n+1} with J_ω^{n+½} = (2/μ₀)(κ₀Ē/c + κ×B^{n+½}), Ē = (E^n + E^{n+1})/2,
// solved exactly per edge (semi-implicit, so κ₀ < 0 stays stable for
// dt·|κ₀|c < 1). B^{n+½} is r.B as it stands before Faraday. Returns the
// current on each edge for the work booking.
func (r *Region) applySpin(eOld [3][]float64, tHalf float64) (jw [3][]float64) {
	st := r.spin
	g := 1.0
	if st.g != nil {
		g = st.g(tHalf)
	}
	cl := em.SpeedOfLight
	cb := r.Dt / em.Eps0
	for c := 0; c < 3; c++ {
		jw[c] = make([]float64, len(r.E[c]))
	}
	r.each(func(i, j, k, m int) {
		for c := 0; c < 3; c++ {
			if !r.eFree(c, i, j, k) {
				continue
			}
			k0 := g * st.k0[c][m]
			kv := st.kv[c][m]
			var hall float64
			if kv != ([3]float64{}) {
				b := r.bAtEdge(c, m)
				d, e := (c+1)%3, (c+2)%3
				hall = g * 2 * (kv[d]*b[e] - kv[e]*b[d]) / em.Mu0
			}
			if k0 == 0 && hall == 0 {
				continue
			}
			// σ/ε₀·dt/2 = (2κ₀/(μ₀c))·dt/(2ε₀) = κ₀ c dt.
			a := k0 * cl * r.Dt
			en := (r.E[c][m] - a*eOld[c][m] - cb*hall) / (1 + a)
			r.E[c][m] = en
			jw[c][m] = 2*k0/(em.Mu0*cl)*0.5*(eOld[c][m]+en) + hall
		}
	})
	return jw
}

// CoulombDefect checks ECE2's Coulomb law in vacuum, ∇·E = −2κ·E (6.168),
// on the field the evolution produced: max over interior nodes of
// |∇·E + 2κ·E| relative to the larger of max|∇·E| and max|2κ·E|. Ampère's
// law with J_ω carries charge ε₀∇·E by continuity; the Coulomb law fixes
// it again from κ. The two agree only for special fields, and this shows
// how far apart they are. Each κ_cE_c at a node is the mean of its two
// c-edges.
func (r *Region) CoulombDefect() float64 {
	if r.spin == nil {
		return 0
	}
	step := [3]int{1, r.N[0], r.N[0] * r.N[1]}
	var worst, sd, sk float64
	r.each(func(i, j, k, m int) {
		if !r.interior(i, j, k) {
			return
		}
		var ke float64
		for c := 0; c < 3; c++ {
			ec := 0.5 * (r.E[c][m] + r.E[c][m-step[c]])
			kc := 0.5 * (r.spin.kv[c][m][c] + r.spin.kv[c][m-step[c]][c])
			ke += kc * ec
		}
		div := r.divE(m)
		worst = math.Max(worst, math.Abs(div+2*ke))
		sd = math.Max(sd, math.Abs(div))
		sk = math.Max(sk, math.Abs(2*ke))
	})
	s := math.Max(sd, sk)
	if s == 0 {
		return 0
	}
	return worst / s
}
