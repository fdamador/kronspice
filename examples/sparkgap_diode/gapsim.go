// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// A spark gap is a switch with memory: it strikes at a threshold and stays
// on while the arc carries current. A static current-controlled curve
// (kronspice's S-type B element) has no such memory, and the snap from
// 2.4 kV to 50 V at breakdown stalls the transient, so the gap rectifiers
// are stepped here with their own small nodal solver (backward Euler on
// the load capacitor), checked against kronspice in TestGapSolverTwin.

type gapState struct {
	on  bool
	sgn float64 // +1 forward (point positive), −1 reverse
	t   float64 // time it struck
	ext float64 // time the last arc went out (−1: never)
}

// RunGaps simulates a rectifier whose arms are PolarGaps (one per arm).
// Nodes: ac1, p, n (n = ground for half-wave).
func (r Rectifier) RunGaps(g PolarGap) RectResult {
	type arm struct{ a, k int } // node indices: 0 ground, 1 ac1, 2 p, 3 n
	var arms []arm
	if r.Bridge {
		arms = []arm{{1, 2}, {0, 2}, {3, 1}, {3, 0}}
	} else {
		arms = []arm{{1, 2}}
	}
	out := 3
	if !r.Bridge {
		out = 0
	}
	per := 1 / r.F
	dt := per / 20000
	steps := r.Cycles * 20000
	st := make([]gapState, len(arms))
	for i := range st {
		st[i].ext = -1
	}
	gl := 1 / g.Rleak()
	v := make([]float64, 4) // node voltages, v[0] = 0
	vcPrev := 0.0
	var res RectResult
	t0 := float64(r.Skip) * per
	lo, hi := math.Inf(1), math.Inf(-1)
	var T, sv, spl, sps float64
	for s := 1; s <= steps; s++ {
		t := float64(s) * dt
		vs := r.Vp * math.Sin(2*math.Pi*r.F*t)
		// assemble G·x = I for x = (v1, v2, v3)
		var G [3][3]float64
		var I [3]float64
		stamp := func(a, b int, gg, cur float64) { // conductance gg from a to b, plus current cur injected a→b
			if a > 0 {
				G[a-1][a-1] += gg
				I[a-1] -= cur
			}
			if b > 0 {
				G[b-1][b-1] += gg
				I[b-1] += cur
			}
			if a > 0 && b > 0 {
				G[a-1][b-1] -= gg
				G[b-1][a-1] -= gg
			}
		}
		// source through Rs into ac1
		G[0][0] += 1 / r.Rs
		I[0] += vs / r.Rs
		gc := r.Cload / dt
		stamp(2, out, gc+1/r.Rload, -gc*vcPrev)
		if r.Bridge {
			stamp(3, 0, 1e-9, 0)
		}
		for i, a := range arms {
			if st[i].on {
				// i = (v − sgn·Varc)/Rarc
				stamp(a.a, a.k, 1/g.Rarc, -st[i].sgn*g.Varc/g.Rarc)
			} else {
				stamp(a.a, a.k, gl, 0)
			}
		}
		if !r.Bridge {
			G[2][2] = 1 // n is unused: pin it to 0
			I[2] = 0
		}
		x := solve3(G, I)
		v[1], v[2], v[3] = x[0], x[1], x[2]
		if !r.Bridge {
			v[3] = 0
		}
		// switching
		for i, a := range arms {
			vg := v[a.a] - v[a.k]
			if st[i].on {
				ig := (vg - st[i].sgn*g.Varc) / g.Rarc
				if t-st[i].t >= g.Hold && ig*st[i].sgn < g.Iext {
					st[i].on = false
					st[i].ext = t
				}
			} else {
				k := 1.0
				if g.Recovery > 0 && st[i].ext >= 0 {
					k = -math.Expm1(-(t - st[i].ext) / g.Recovery)
				}
				fwd, rev := g.restrike(k, g.Vpos()), g.restrike(k, g.Vneg())
				if vg >= fwd {
					st[i] = gapState{true, 1, t, st[i].ext}
					res.Forward++
				} else if vg <= -rev {
					st[i] = gapState{true, -1, t, st[i].ext}
					res.Backward++
				}
			}
			res.ReversePeak = math.Max(res.ReversePeak, -vg)
		}
		vo := v[2] - v[out]
		is := (vs - v[1]) / r.Rs
		res.EnergyIn += vs * is * dt
		vcPrev = vo
		res.Time = append(res.Time, t)
		res.Vout = append(res.Vout, vo)
		if t < t0 {
			continue
		}
		lo, hi = math.Min(lo, vo), math.Max(hi, vo)
		T += dt
		sv += vo * dt
		spl += vo * vo / r.Rload * dt
		sps += vs * is * dt
	}
	res.Vdc = sv / T
	res.Ripple = hi - lo
	res.Idc = res.Vdc / r.Rload
	res.Pload, res.Psource = spl/T, sps/T
	res.Efficiency = res.Pload / res.Psource
	res.Stored = 0.5 * r.Cload * vcPrev * vcPrev
	res.Heat = res.EnergyIn - res.Stored // backward Euler: whatever is not stored has been dissipated
	return res
}

// solve3 solves a 3×3 system by Gaussian elimination with partial pivoting.
func solve3(A [3][3]float64, b [3]float64) [3]float64 {
	for c := 0; c < 3; c++ {
		p := c
		for r := c + 1; r < 3; r++ {
			if math.Abs(A[r][c]) > math.Abs(A[p][c]) {
				p = r
			}
		}
		A[c], A[p] = A[p], A[c]
		b[c], b[p] = b[p], b[c]
		for r := c + 1; r < 3; r++ {
			f := A[r][c] / A[c][c]
			for k := c; k < 3; k++ {
				A[r][k] -= f * A[c][k]
			}
			b[r] -= f * b[c]
		}
	}
	var x [3]float64
	for r := 2; r >= 0; r-- {
		s := b[r]
		for k := r + 1; k < 3; k++ {
			s -= A[r][k] * x[k]
		}
		x[r] = s / A[r][r]
	}
	return x
}
