// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 9c field twin: two (or N) dielectric tori of Part 5b (εr = 100,
// a = 10 mm, r = 2.5 mm, lossless) stacked on a common axis, centres D
// apart, in em/fdtd.go. Each run starts every torus with the azimuthal
// displacement bump of Part 5b, D_φ = s_i·D0·cos²(πs/2r) (s_i the ring's
// weight: +1, +1 in phase; +1, −1 antiphase; the band-edge weights for a
// chain), B = 0, and rings down into the absorbing layer. The mirror
// symmetry z → −z keeps an in-phase start in the even modes and an
// antiphase start in the odd ones, and the azimuthal symmetry keeps both
// in the m = 0 (TE0) family; the slowest of the excited modes is what
// is left after the start-up, and its Q is read from the decay of the
// field energy inside a surface 2a outside the rings, its frequency from
// the zero crossings of E_φ at (a, 0, z_0). The field ledger (field +
// absorbed = start) is checked.
//
// Approximations: staircased tori on r/Cells cells, absorber 3a from the
// outermost ring (the subradiant pair's near field reaches it weakly),
// no dielectric loss (compare with the circuit's radiation-only Q).

import (
	"math"
	"math/cmplx"
	"time"

	"github.com/fdamador/kronspice/em"
)

// PairField is one field run.
type PairField struct {
	A, R, EpsR float64
	Cells      float64   // cells per tube radius
	D          float64   // centre spacing [m]
	W          []float64 // weights of the rings (their number is the chain length)
	HalfA      float64   // vacuum margin beyond the outermost ring, in a (absorber outside it)
	SurfA      float64   // energy surface margin, in a
	Absorb     int
	Periods    float64 // run length, lumped periods
	FitFrom    float64 // fit window start, lumped periods
	FLo, FHi   float64 // band searched for the mode [Hz]
	Log        func(format string, args ...any)
}

// DefaultPairField is the r/2 grid.
func DefaultPairField(d float64, w []float64) PairField {
	return PairField{A: 10e-3, R: 2.5e-3, EpsR: 100, Cells: 2, D: d, W: w, HalfA: 3, SurfA: 2, Absorb: 12, Periods: 40, FitFrom: 4, FLo: 1.0e9, FHi: 4.5e9}
}

// PairFieldResult is what one run measured.
type PairFieldResult struct {
	Cells, D float64
	Nodes    int
	F, Q     float64 // mode frequency and Q (demodulated probe)
	QEnergy  float64 // from the decay of all the field energy inside the surface
	AT, Amp  []float64
	Residual float64 // ledger residual / start energy
	T, U     []float64
	Seconds  float64
	Weights  []float64
}

func (p PairField) zc(i int) float64 { return (float64(i) - float64(len(p.W)-1)/2) * p.D }

func (p PairField) inRing(x, y, z float64) int {
	rho := math.Hypot(x, y)
	for i := range p.W {
		if math.Hypot(rho-p.A, z-p.zc(i)) < p.R {
			return i
		}
	}
	return -1
}

// bump is −∫_{a−r}^{ρ} cos²(πs/2r) dρ' at height z (s from the tube centre).
func (p PairField) bump(rho, z float64) float64 {
	lo := math.Max(0, p.A-p.R)
	if rho <= lo || math.Abs(z) >= p.R {
		return 0
	}
	hi := math.Min(rho, p.A+p.R)
	const m = 64
	h := (hi - lo) / m
	var s float64
	for q := 0; q < m; q++ {
		r := lo + (float64(q)+0.5)*h
		sd := math.Hypot(r-p.A, z)
		if sd < p.R {
			c := math.Cos(math.Pi * sd / (2 * p.R))
			s += c * c * h
		}
	}
	return -s
}

// Run performs the ring-down.
func (p PairField) Run() (PairFieldResult, error) {
	t0 := time.Now()
	dx := p.R / p.Cells
	mxy := int(math.Ceil((1 + p.HalfA) * p.A / dx))
	zmax := p.zc(len(p.W)-1) + p.HalfA*p.A + p.R
	mz := int(math.Ceil(zmax / dx))
	nxy, nz := 2*mxy+2*p.Absorb+1, 2*mz+2*p.Absorb+1
	hxy, hz := float64(nxy-1)/2*dx, float64(nz-1)/2*dx
	mat := func(x, y, z float64) (float64, float64) {
		if p.inRing(x, y, z) >= 0 {
			return p.EpsR, 0
		}
		return 1, 0
	}
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{nxy, nxy, nz}, D: [3]float64{dx, dx, dx},
		Origin: [3]float64{-hxy, -hxy, -hz}, Absorber: p.Absorb, Ledger: true, Material: mat})
	if err != nil {
		return PairFieldResult{}, err
	}
	f.SetEFromPotential(func(x, y, z float64) [3]float64 {
		rho := math.Hypot(x, y)
		var c float64
		for i, w := range p.W {
			c += w * p.bump(rho, z-p.zc(i))
		}
		return [3]float64{0, 0, c}
	})
	// D was set: divide by εr on every E edge
	for c := 0; c < 3; c++ {
		for k := 0; k < nz; k++ {
			for j := 0; j < nxy; j++ {
				for i := 0; i < nxy; i++ {
					q := [3]float64{float64(i), float64(j), float64(k)}
					q[c] += 0.5
					er, _ := mat(-hxy+q[0]*dx, -hxy+q[1]*dx, -hz+q[2]*dx)
					f.E[c][i+nxy*(j+nxy*k)] /= er
				}
			}
		}
	}
	bxy := (1+p.SurfA)*p.A + 0.25*dx
	bz := p.zc(len(p.W)-1) + p.SurfA*p.A + 0.25*dx
	inside := func(x, y, z float64) bool { return math.Abs(x) < bxy && math.Abs(y) < bxy && math.Abs(z) < bz }
	W0 := f.RegionEnergy(func(x, y, z float64) bool { return true })
	e := DefaultERing()
	Tl := 2 * math.Pi / e.Omega0()
	every := int(math.Max(1, math.Round(Tl/16/f.Dt)))
	res := PairFieldResult{Cells: p.Cells, D: p.D, Nodes: nxy * nxy * nz, Weights: p.W}
	z0 := p.zc(0)
	var pt, pv []float64 // probe E_φ at (a, 0, z0) every step
	for f.Time() < p.Periods*Tl {
		f.Step()
		E, _ := f.FieldAt(p.A, 0, z0)
		pt, pv = append(pt, f.Time()), append(pv, E[1])
		if f.Steps%every == 0 {
			res.T = append(res.T, f.Time()/Tl)
			res.U = append(res.U, f.RegionEnergy(inside))
		}
	}
	res.Residual = f.Ledger().Residual() / W0
	// Energy decay of everything inside the surface (all excited modes).
	var xs, ys []float64
	for k := range res.T {
		if res.T[k] >= p.FitFrom && res.U[k] >= 1e-7*res.U[0] {
			xs, ys = append(xs, res.T[k]*Tl), append(ys, math.Log(res.U[k]))
		}
	}
	if len(xs) > 3 {
		res.QEnergy = -e.Omega0() / fitSlope(xs, ys)
	}
	// The mode: the strongest spectral peak of the probe between FLo and
	// FHi (after 2 periods, Hann window), then its complex amplitude by
	// demodulation over a sliding 3-period Hann window, Q = ω/(2σ) from
	// the decay of |amplitude| (other modes are filtered out).
	var ts, ps []float64
	for k := range pt {
		if pt[k] > 2*Tl {
			ts, ps = append(ts, pt[k]), append(ps, pv[k])
		}
	}
	spec := func(fr float64, lo, hi int) complex128 {
		var s complex128
		for k := lo; k < hi; k++ {
			w := 0.5 - 0.5*math.Cos(2*math.Pi*float64(k-lo)/float64(hi-lo-1))
			s += complex(w*ps[k], 0) * cmplx.Exp(complex(0, -2*math.Pi*fr*ts[k]))
		}
		return s
	}
	best, bf := -1.0, p.FLo
	for fr := p.FLo; fr <= p.FHi; fr += (p.FHi - p.FLo) / 400 {
		if v := cmplx.Abs(spec(fr, 0, len(ts))); v > best {
			best, bf = v, fr
		}
	}
	lo, hi := bf-(p.FHi-p.FLo)/400, bf+(p.FHi-p.FLo)/400
	gr := (math.Sqrt(5) - 1) / 2
	for it := 0; it < 50; it++ {
		m1, m2 := hi-gr*(hi-lo), lo+gr*(hi-lo)
		if cmplx.Abs(spec(m1, 0, len(ts))) > cmplx.Abs(spec(m2, 0, len(ts))) {
			hi = m2
		} else {
			lo = m1
		}
	}
	res.F = (lo + hi) / 2
	dtp := ts[1] - ts[0]
	win := int(math.Round(3 / res.F / dtp))
	slide := int(math.Max(1, math.Round(0.25/res.F/dtp)))
	var ax, ay []float64
	var amax float64
	for s0 := 0; s0+win < len(ts); s0 += slide {
		a := cmplx.Abs(spec(res.F, s0, s0+win))
		tm := ts[s0+win/2]
		res.AT = append(res.AT, tm/Tl)
		res.Amp = append(res.Amp, a)
		amax = math.Max(amax, a)
	}
	for k := range res.AT {
		if res.AT[k] >= p.FitFrom && res.Amp[k] > 1e-5*amax {
			ax, ay = append(ax, res.AT[k]*Tl), append(ay, math.Log(res.Amp[k]))
		}
	}
	if len(ax) > 3 {
		res.Q = -2 * math.Pi * res.F / (2 * fitSlope(ax, ay))
	}
	res.Seconds = time.Since(t0).Seconds()
	if p.Log != nil {
		p.Log("field: %d ring(s), d = %.3g a, weights %v, r/%g: f = %.5g GHz, Q = %.4g, residual %.1e (%d nodes, %.0f s)",
			len(p.W), p.D/p.A, p.W, p.Cells, res.F/1e9, res.Q, res.Residual, res.Nodes, res.Seconds)
	}
	return res, nil
}
