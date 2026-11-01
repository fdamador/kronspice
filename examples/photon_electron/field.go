// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/em"
)

// Part 3: the paper's central assumption under linear Maxwell equations.
//
// Williamson & van der Mark place one wavelength of light on a closed path
// that runs twice round a torus and assume it stays there. Here that light
// is set up in the field solver (em/fdtd.go) and left alone in empty space:
//
//   - path: the (2, 1) torus knot, twice round the ring (radius R) while
//     once round the tube (radius R/2); R = λ/4π as in the paper;
//   - field: a travelling wave cos(ks − ωt) along the path, k = 2π/L for
//     path length L, in a tube of Gaussian width w round it, with E and B
//     across the path so the energy flows along it (E×B along the
//     tangent);
//   - E and B are built as discrete curls of sampled potentials, so
//     div E = div B = 0 exactly: the start is a valid charge-free field.
//
// The solver then measures how much energy is still within 2R of the
// centre as time goes on, and how much left through the absorbing walls.
// Linear Maxwell equations have no term that holds light in a loop, so the
// standard prediction is that it leaves within about one period; the
// model needs it to stay for ever. Lengths scale out (linear equations):
// results are in units of R, λ and the period T = λ/c.

// KnotField sets up the light on the double loop.
type KnotField struct {
	Lambda  float64 // wavelength [m]
	CellsR  float64 // grid cells per ring radius R
	WidthR  float64 // Gaussian width of the tube of light, in R
	HalfR   float64 // half-width of the region inside the absorber, in R
	Periods float64 // how long to follow it, in T
}

// DefaultKnot: λ = 1 m, 12 cells per R, tube width R/3, region ±6R
// (about λ/2, so the absorber sits outside the near field), 4 periods.
func DefaultKnot() KnotField {
	return KnotField{Lambda: 1, CellsR: 12, WidthR: 1.0 / 3, HalfR: 6, Periods: 4}
}

// KnotSample is the energy at one time.
type KnotSample struct {
	T        float64 // time in periods
	Inside   float64 // fraction of the starting energy within 2R of the centre
	InRegion float64 // fraction anywhere in the region
	Absorbed float64 // fraction that left through the walls
}

// KnotRun is the outcome.
type KnotRun struct {
	R, Dx, PathLength, W0 float64
	Samples               []KnotSample
	Residual              float64 // ledger residual / W0
	DivB                  float64 // max |div B| / (max|B|/dx), inside the absorber
}

func sq(x float64) float64 { return x * x }

// knot is the (2, 1) torus knot at parameter u ∈ [0, 4π): φ = u (twice
// round the ring), θ = u/2 (once round the tube). Returns point and unit
// tangent.
func (k KnotField) knot(u float64) (p, t [3]float64) {
	R := k.Lambda / (4 * math.Pi)
	a := R / 2
	th := u / 2
	rr := R + a*math.Cos(th)
	p = [3]float64{rr * math.Cos(u), rr * math.Sin(u), a * math.Sin(th)}
	d := [3]float64{
		-a/2*math.Sin(th)*math.Cos(u) - rr*math.Sin(u),
		-a/2*math.Sin(th)*math.Sin(u) + rr*math.Cos(u),
		a / 2 * math.Cos(th),
	}
	n := math.Sqrt(sq(d[0]) + sq(d[1]) + sq(d[2]))
	t = [3]float64{d[0] / n, d[1] / n, d[2] / n}
	return
}

// Run builds the region, launches the light and follows it.
func (k KnotField) Run() (KnotRun, error) {
	R := k.Lambda / (4 * math.Pi)
	dx := R / k.CellsR
	abs := 12
	n := 2*int(math.Ceil(k.HalfR*k.CellsR)) + 2*abs + 1
	half := float64(n-1) / 2 * dx
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx},
		Origin: [3]float64{-half, -half, -half}, Absorber: abs, Ledger: true})
	if err != nil {
		return KnotRun{}, err
	}
	// Path samples with arc length s and the frame (t, n1, n2).
	const M = 800
	type pt struct {
		p, n1, n2 [3]float64
		s         float64
	}
	pts := make([]pt, M)
	s := 0.0
	prev, _ := k.knot(0)
	for i := 0; i < M; i++ {
		p, t := k.knot(4 * math.Pi * float64(i) / M)
		s += math.Sqrt(sq(p[0]-prev[0]) + sq(p[1]-prev[1]) + sq(p[2]-prev[2]))
		prev = p
		n2 := [3]float64{-t[2] * t[0], -t[2] * t[1], 1 - t[2]*t[2]} // ẑ off the tangent
		nn := math.Sqrt(sq(n2[0]) + sq(n2[1]) + sq(n2[2]))
		n2 = [3]float64{n2[0] / nn, n2[1] / nn, n2[2] / nn}
		n1 := [3]float64{n2[1]*t[2] - n2[2]*t[1], n2[2]*t[0] - n2[0]*t[2], n2[0]*t[1] - n2[1]*t[0]} // n2 × t
		pts[i] = pt{p: p, n1: n1, n2: n2, s: s}
	}
	first, _ := k.knot(0)
	L := s + math.Sqrt(sq(first[0]-prev[0])+sq(first[1]-prev[1])+sq(first[2]-prev[2]))
	ds := L / M
	kw := 2 * math.Pi / L
	w := k.WidthR * R
	c := em.SpeedOfLight
	omega := c * kw
	// Potentials as sums of Gaussian blobs along the path:
	// C = Σ g·cos(ks)·n1/k gives E = curl C ≈ −g sin(ks)·n2,
	// A = Σ g·cos(ks)·n2/(kc) gives B = curl A ≈ g sin(ks)·n1/c,
	// so E×B points along the tangent.
	pot := func(x, y, z float64, useN1 bool, phase, scale float64) [3]float64 {
		var out [3]float64
		if math.Hypot(math.Hypot(x, y)-R, z) > R/2+4*w {
			return out
		}
		for _, q := range pts {
			r2 := sq(x-q.p[0]) + sq(y-q.p[1]) + sq(z-q.p[2])
			if r2 > 16*w*w {
				continue
			}
			g := math.Exp(-r2/(2*w*w)) * math.Cos(kw*q.s-phase) * ds * scale
			v := q.n2
			if useN1 {
				v = q.n1
			}
			for i := 0; i < 3; i++ {
				out[i] += g * v[i]
			}
		}
		return out
	}
	f.SetEFromPotential(func(x, y, z float64) [3]float64 { return pot(x, y, z, true, 0, 1/kw) })
	f.SetBFromPotential(func(x, y, z float64) [3]float64 { return pot(x, y, z, false, omega*f.Dt/2, 1/(kw*c)) })
	inside := func(x, y, z float64) bool { return x*x+y*y+z*z < 4*R*R }
	all := func(x, y, z float64) bool { return true }
	run := KnotRun{R: R, Dx: dx, PathLength: L}
	run.W0 = f.RegionEnergy(all)
	T := k.Lambda / c
	for _, frac := range []float64{0, 0.125, 0.25, 0.5, 1, 1.5, 2, 3, 4} {
		if frac > k.Periods {
			break
		}
		f.Run(frac * T)
		run.Samples = append(run.Samples, KnotSample{T: frac, Inside: f.RegionEnergy(inside) / run.W0,
			InRegion: f.RegionEnergy(all) / run.W0, Absorbed: f.Ledger().Absorbed / run.W0})
	}
	run.Residual = f.Ledger().Residual() / run.W0
	wd, sc := f.DivB()
	run.DivB = wd / sc
	return run, nil
}
