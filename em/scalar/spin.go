// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"fmt"
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
//
// KappaL0, KappaL add the homogeneous current of UFT 454 eqs. 24–27 (UFT 456
// eqs. 28–29, Eckardt & Goel 2024) to Faraday's law, ∂B/∂t + ∇×E = V with
//
//	V = 2(cκ_Λ0 B − κ_Λ×E)        [V/m²]
//
// The κ_Λ0 part is a magnetic conductivity σ* = −η₀ with UFT 456's
// η₀ = 2cκ_Λ0μ₀ (eq. 39): κ_Λ0 > 0 is σ* < 0, a gain medium that feeds the
// field η₀H² per unit volume (UFT 456 eq. 57), κ_Λ0 < 0 a magnetic loss (the
// PML's σ*, Taflove & Hagness ch. 7). Its work is booked in
// Ledger.SpinWorkHom, exact on the grid. Kappa0 = KappaL0 is UFT 456's
// "equilibrium" medium (eq. 109): E loses and H gains at the same rate, and
// the dispersion ω² = c²k² − (2cκ₀)² grows every mode with k < 2κ₀
// (ECE2Dispersion). Heaviside's distortionless medium is KappaL0 = −Kappa0.
type SpinConnection struct {
	Kappa0 func(x, y, z float64) float64    // κ₀ [1/m]; nil = 0
	Kappa  func(x, y, z float64) [3]float64 // κ [1/m]; nil = 0
	// Homogeneous (Faraday) part, UFT 456 eq. 29.
	KappaL0 func(x, y, z float64) float64    // κ_Λ0 [1/m]; nil = 0
	KappaL  func(x, y, z float64) [3]float64 // κ_Λ [1/m]; nil = 0
	// Envelope scales all four in time (a spin connection that is switched
	// or modulated, |Envelope| ≤ 1); nil = 1.
	Envelope func(t float64) float64
}

type spinState struct {
	k0 [3][]float64    // κ₀ at each E edge, per component
	kv [3][][3]float64 // κ at each E edge
	g  func(t float64) float64

	hom   bool            // homogeneous part present
	kL0   [3][]float64    // κ_Λ0 at each B face
	kL    [3][][3]float64 // κ_Λ at each B face
	vPrev [3][]float64    // V^n on each face, for the work booking
}

// setupSpin samples κ₀ and κ at the E edges and κ_Λ0, κ_Λ at the B faces.
// It refuses a gain κ_Λ0 > 0 (or a κ₀ < 0) whose semi-implicit update would
// change sign within one step: c|κ|dt must stay below 1.
func (r *Region) setupSpin(s *SpinConnection) error {
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
	cl := em.SpeedOfLight
	for c := 0; c < 3; c++ {
		for _, k0 := range st.k0[c] {
			if math.Abs(k0)*cl*r.Dt >= 1 {
				return fmt.Errorf("%w: c|κ₀|dt = %.3g ≥ 1", ErrConfig, math.Abs(k0)*cl*r.Dt)
			}
		}
	}
	if s.KappaL0 != nil || s.KappaL != nil {
		st.hom = true
		for c := 0; c < 3; c++ {
			n3 := len(r.B[c])
			st.kL0[c] = make([]float64, n3)
			st.kL[c] = make([][3]float64, n3)
			st.vPrev[c] = make([]float64, n3)
			imax, jmax, kmax := r.faceRange(c)
			for k := 0; k < kmax; k++ {
				for j := 0; j < jmax; j++ {
					for i := 0; i < imax; i++ {
						off := [3]float64{0.5, 0.5, 0.5}
						off[c] = 0
						x, y, z := r.pos(off, i, j, k)
						m := r.index(i, j, k)
						if s.KappaL0 != nil {
							st.kL0[c][m] = s.KappaL0(x, y, z)
							if math.Abs(st.kL0[c][m])*cl*r.Dt >= 1 {
								return fmt.Errorf("%w: c|κ_Λ0|dt = %.3g ≥ 1", ErrConfig, math.Abs(st.kL0[c][m])*cl*r.Dt)
							}
						}
						if s.KappaL != nil {
							st.kL[c][m] = s.KappaL(x, y, z)
						}
					}
				}
			}
		}
	}
	r.spin = st
	return nil
}

// eAtFace is E^{n+1} interpolated to B face (c, m): the two components
// tangential to the face, each the mean of the two edges that straddle it.
func (r *Region) eAtFace(c, i, j, k, m int) [3]float64 {
	step := [3]int{1, r.N[0], r.N[0] * r.N[1]}
	ijk := [3]int{i, j, k}
	var e [3]float64
	for d := 0; d < 3; d++ {
		if d == c {
			continue
		}
		o := 3 - c - d // the other tangential axis: E_d is averaged along it
		lo := r.E[d][m]
		if ijk[o]+1 < r.N[o] {
			e[d] = 0.5 * (lo + r.E[d][m+step[o]])
		} else {
			e[d] = 0.5 * lo
		}
	}
	return e
}

// homFaraday is the Faraday update of face (c, m) with UFT 456's
// V = 2(cκ_Λ0 B̄ − κ_Λ×E^{n+1}), B̄ = (B^{n+½} + B^{n+3/2})/2, solved exactly
// (semi-implicit): rest is the rest of the update, −dt(curl E + …). It
// returns B^{n+3/2} and the V it applied.
func (r *Region) homFaraday(c, i, j, k, m int, old, rest, g float64) (nw, v float64) {
	st := r.spin
	a := g * st.kL0[c][m] * em.SpeedOfLight * r.Dt
	var lor float64
	if kl := st.kL[c][m]; kl != ([3]float64{}) {
		e := r.eAtFace(c, i, j, k, m)
		d, f := (c+1)%3, (c+2)%3
		lor = g * 2 * (kl[d]*e[f] - kl[f]*e[d])
	}
	if a == 0 && lor == 0 {
		return old + rest, 0
	}
	nw = ((1+a)*old + rest - r.Dt*lor) / (1 - a)
	v = (nw - old - rest) / r.Dt
	return nw, v
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
