// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Step 1b: the ether under tension: n twisted at a uniform rate k.
//
// With φ = a·e^{ikx} (n turning in a plane) the ether thins to the local
// minimum a0(k) of V + κ2k²a² (Q: a0² = 1 − 2ξ²k²) and carries the twist
// stress j = 2κ2·a0²·k (the "supercurrent" of a superfluid, the torque of a
// twisted rod). The uniform state is stable against long-wave amplitude–
// phase perturbations only while (Bloch analysis of the second variation,
// q → 0; Eckhaus 1965; Langer & Ambegaokar, Phys. Rev. 164, 498, 1967)
//
//	½ V″(a0) > 3 κ2 k²     ⇔     dj/dk > 0,
//
// so the largest stress the ether can carry is j_c = max_k j(k): its
// tensile strength. Q: k_c = 1/(√6 ξ), a0² = 2/3, j_c = (4/3√6)√(κ2λ).
// Past k_c the ether tears: the amplitude falls to zero at a point (a phase
// slip, the 1-D cavity) and one turn of twist is released. Below k_c but
// with the cavity phase favoured (S, κ2k² > ε) it can tear only by
// nucleation, over a barrier (radial.go): the binodal and the spinodal of a
// liquid.
//
// An S² ether (n ∈ S², as Faddeev's) does not hold a uniform twist at all:
// tilting n out of the plane lowers κ2a²k²cos²χ, so the twist unwinds over
// the pole with a ≈ 1, no cavity. Tension in an S² ether needs topology
// that forbids the escape (antipodal ends, reversal.go; a Hopf knot).

import (
	"math"
	"math/rand/v2"
)

// TwistPoint is the uniform twisted state at rate k.
type TwistPoint struct {
	K, A0, J float64
	Stable   bool // Eckhaus criterion
	Exists   bool
}

// UniformTwist evaluates the twisted state of e at rate k.
func UniformTwist(e EOS, k float64) TwistPoint {
	t := TwistPoint{K: k}
	a0, ok := e.Twist(k).Ether()
	if !ok {
		return t
	}
	t.Exists, t.A0, t.J = true, a0, 2*kappa2*a0*a0*k
	// V″ of the potential without the twist term (the twist is in the
	// gradient here): D2 of the twisted EOS minus 2κ2k².
	v2 := e.Twist(k).D2(a0) - 2*kappa2*k*k
	t.Stable = 0.5*v2 > 3*kappa2*k*k
	return t
}

// Spinodal returns k_c by bisection on the Eckhaus criterion, and k_j, the
// maximum of j(k) on a fine grid (the two must agree).
func Spinodal(e EOS) (kc, kj, jc float64) {
	lo, hi := 0.0, 10/e.Xi()
	for it := 0; it < 100; it++ {
		m := (lo + hi) / 2
		if p := UniformTwist(e, m); p.Exists && p.Stable {
			lo = m
		} else {
			hi = m
		}
	}
	kc = (lo + hi) / 2
	const n = 200000
	for i := 1; i < n; i++ {
		k := 2 * kc * float64(i) / n
		if p := UniformTwist(e, k); p.Exists && p.J > jc {
			jc, kj = p.J, k
		}
	}
	return kc, kj, jc
}

// Binodal (S only, ε > 0): the twist at which the cavity phase (a = 0,
// energy ε) costs the same as the twisted ether.
func Binodal(e EOS) float64 {
	f := func(k float64) float64 {
		a0, ok := e.Twist(k).Ether()
		if !ok {
			return -1
		}
		return e.Eps - e.Twist(k).V(a0)
	}
	lo, hi := 0.0, 10/e.Xi()
	for it := 0; it < 100; it++ {
		m := (lo + hi) / 2
		if f(m) > 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return (lo + hi) / 2
}

// RingRun is a relaxation (gradient-flow) run of a twisted ring.
type RingRun struct {
	W0, W1        int     // winding at start and end
	MinA          float64 // least |φ| seen
	SlipT         float64 // first time |φ| < 0.05 (−1: never)
	E0, E1, Diss  float64 // energy start, end, dissipated
	LedgerRel     float64 // |E0 − E1 − Diss| / E0
	MinAOutOfLine float64 // S² run: least |φ| while the twist unwound
}

// Ring relaxes a ring of m nodes, spacing h, wound w times in the xy
// plane, with seeded noise of amplitude noise; comps = 2 (U(1)) or 3 (S²).
// Flow: ∂φ/∂t = −(1/W)∂E/∂φ (continuum ∂φ/∂t = 2κ2φ″ − 2V_u φ), explicit
// Euler with dt = 0.2h². Ledger: dE = −Σ W|∂φ/∂t|² dt per step (first
// order in dt; the residual is reported).
func Ring(e EOS, m int, h float64, w, comps int, noise, T float64, seed uint64) RingRun {
	rng := rand.New(rand.NewPCG(seed, 1))
	L := float64(m) * h
	k := 2 * math.Pi * float64(w) / L
	a0, _ := e.Twist(k).Ether()
	p := make([][3]float64, m)
	for i := range p {
		x := float64(i) * h
		p[i] = [3]float64{a0 * math.Cos(k*x), a0 * math.Sin(k*x), 0}
		for c := 0; c < comps; c++ {
			p[i][c] += noise * (2*rng.Float64() - 1)
		}
	}
	energy := func() float64 {
		var s float64
		for i := range p {
			j := (i + 1) % m
			var d2, u float64
			for c := 0; c < 3; c++ {
				d := p[j][c] - p[i][c]
				d2 += d * d
				u += p[i][c] * p[i][c]
			}
			s += kappa2*d2/h + h*e.V(math.Sqrt(u))
		}
		return s
	}
	winding := func() int {
		var s float64
		for i := range p {
			j := (i + 1) % m
			s += math.Atan2(p[i][0]*p[j][1]-p[i][1]*p[j][0], p[i][0]*p[j][0]+p[i][1]*p[j][1])
		}
		return int(math.Round(s / (2 * math.Pi)))
	}
	r := RingRun{W0: winding(), MinA: math.Inf(1), SlipT: -1, E0: energy()}
	dt := 0.2 * h * h
	vel := make([][3]float64, m)
	steps := int(T / dt)
	for s := 0; s < steps; s++ {
		var diss float64
		for i := range p {
			im, ip := (i+m-1)%m, (i+1)%m
			var u float64
			for c := 0; c < 3; c++ {
				u += p[i][c] * p[i][c]
			}
			_, vu, _ := e.Vu(u)
			for c := 0; c < comps; c++ {
				g := 2*h*vu*p[i][c] + 2*kappa2*(2*p[i][c]-p[im][c]-p[ip][c])/h
				vel[i][c] = -g / h
				diss += h * vel[i][c] * vel[i][c] * dt
			}
		}
		for i := range p {
			for c := 0; c < comps; c++ {
				p[i][c] += dt * vel[i][c]
			}
			a := math.Sqrt(p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2])
			if a < r.MinA {
				r.MinA = a
				if a < 0.05 && r.SlipT < 0 {
					r.SlipT = float64(s) * dt
				}
			}
		}
		r.Diss += diss
	}
	r.W1, r.E1 = winding(), energy()
	r.LedgerRel = math.Abs(r.E0-r.E1-r.Diss) / r.E0
	return r
}
