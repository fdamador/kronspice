// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Step 1c: a reversal squeezed into a slab (see wall.go for the physics).
// φ(0) = (1, 0), φ(L) = (−1, 0). Two branches:
//
//	straight: φ = (f(x), 0), passing through zero: a cavity sheet;
//	turned:   φ = a(x)(cos θ, sin θ), θ: 0 → π, a > 0: the ether bends.
//
// L_c: the straight branch loses stability (out-of-line) as L grows.
// L_t: the turned branch loses stability or ceases to exist as L shrinks;
// below L_t the ether forced to reverse must tear. L_t = L_c means a
// continuous (pitchfork) onset; L_t < L_c means hysteresis: a sheet once torn
// stays torn until the slab is opened past L_c.

import "math"

// Branch is one stationary reversal.
type Branch struct {
	Energy, MinA float64
	Neg          int  // negative Hessian eigenvalues (0 = stable)
	Turned       bool // max |φ_y| > 1e-6
	P            []v2
}

func slabLadder(L, h float64) (Ladder, int) {
	n := int(math.Round(L/h)) + 1
	if n < 7 {
		n = 7
	}
	return Planar(n, L/float64(n-1)), n
}

func branchOf(l Ladder, e EOS, p []v2) (Branch, error) {
	_, neg, err := l.Solve2(e, p, 1e-10, 400)
	if err != nil {
		return Branch{}, err
	}
	var my float64
	for _, q := range p {
		my = math.Max(my, math.Abs(q[1]))
	}
	return Branch{Energy: l.Energy2(e, p), MinA: minAmp(p), Neg: neg, Turned: my > 1e-6, P: p}, nil
}

// Straight solves the through-zero reversal at length L.
func Straight(e EOS, L, h float64) (Branch, error) {
	l, n := slabLadder(L, h)
	p := make([]v2, n)
	for i := range p {
		x := float64(i)*l.H - L/2
		p[i] = v2{-math.Tanh(x / (2 * e.Xi())), 0}
	}
	p[0], p[n-1] = v2{1, 0}, v2{-1, 0}
	return branchOf(l, e, p)
}

// Turned solves the bending reversal at L from a resampled previous
// solution (nil: a half circle at a = 1).
func Turned(e EOS, L, h float64, prev []v2) (Branch, error) {
	l, n := slabLadder(L, h)
	p := make([]v2, n)
	for i := range p {
		s := float64(i) / float64(n-1)
		if prev == nil {
			p[i] = v2{math.Cos(math.Pi * s), math.Sin(math.Pi * s)}
			continue
		}
		x := s * float64(len(prev)-1)
		j := int(x)
		if j >= len(prev)-1 {
			j = len(prev) - 2
		}
		t := x - float64(j)
		p[i] = v2{prev[j][0]*(1-t) + prev[j+1][0]*t, prev[j][1]*(1-t) + prev[j+1][1]*t}
	}
	p[0], p[n-1] = v2{1, 0}, v2{-1, 0}
	return branchOf(l, e, p)
}

// ReversalResult holds the two critical lengths.
type ReversalResult struct {
	Lc, Lt float64 // in ℓ
	// Energies at L = Lt (turned, straight) for the record.
	ETurned, EStraight float64
}

// CriticalLengths finds L_c by bisection on the straight branch's
// stability, and L_t by continuation of the turned branch down from 12ξ.
func CriticalLengths(e EOS, h float64) (ReversalResult, error) {
	var r ReversalResult
	xi := e.Xi()
	lo, hi := 0.2*xi, 30*xi
	for it := 0; it < 40; it++ {
		m := (lo + hi) / 2
		b, err := Straight(e, m, h)
		if err != nil {
			return r, err
		}
		if b.Neg == 0 {
			lo = m
		} else {
			hi = m
		}
	}
	r.Lc = (lo + hi) / 2
	good := func(b Branch) bool { return b.Turned && b.Neg == 0 && b.MinA > 1e-3 }
	L := 12 * xi
	b, err := Turned(e, L, h, nil)
	if err != nil || !good(b) {
		return r, errNoConverge
	}
	step := 0.5 * xi
	last := b
	for step > 1e-5*xi && L > 0.05*xi {
		nb, err := Turned(e, L-step, h, last.P)
		if err == nil && good(nb) {
			L -= step
			last = nb
		} else {
			step /= 2
		}
	}
	r.Lt = L
	r.ETurned = last.Energy
	if s, err := Straight(e, L, h); err == nil {
		r.EStraight = s.Energy
	}
	return r, nil
}
