// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"sync"
)

// The Mandelbulb (White–Nylander power-n map), the user's challenger
// cavity (H46, proposed). v ↦ vⁿ + c in spherical form, starting at v = 0.
//
// Two angle conventions are in use (Hvidtfeldt 2011) and they give
// different solids:
//   - Polar (Wikipedia's): θ = arccos(z/r),
//     vⁿ = rⁿ(sin nθ cos nφ, sin nθ sin nφ, cos nθ).
//   - Latitude: θ = arcsin(z/r),
//     vⁿ = rⁿ(cos nθ cos nφ, cos nθ sin nφ, sin nθ).
//     Its equatorial slice z = 0 is exactly the degree-n Multibrot set
//     zⁿ + c (tested).
// Both have (n − 1)-fold symmetry about z: turning c by 2π/(n − 1) turns
// every iterate by the same angle, because n·α ≡ α (tested).

// Bulb is one Mandelbulb.
type Bulb struct {
	N        float64 // power
	Latitude bool    // latitude convention (default: polar)
	MaxIt    int     // iteration depth: how much fine detail exists
}

// Escape returns the iteration at which |v| > 2, or MaxIt if it never does.
func (b Bulb) Escape(cx, cy, cz float64) int {
	x, y, z := 0.0, 0.0, 0.0
	for i := 0; i < b.MaxIt; i++ {
		r := math.Sqrt(x*x + y*y + z*z)
		if r > 2 {
			return i
		}
		if r == 0 {
			x, y, z = cx, cy, cz
			continue
		}
		rn := math.Pow(r, b.N)
		ph := math.Atan2(y, x) * b.N
		if b.Latitude {
			th := math.Asin(z/r) * b.N
			x = rn*math.Cos(th)*math.Cos(ph) + cx
			y = rn*math.Cos(th)*math.Sin(ph) + cy
			z = rn*math.Sin(th) + cz
		} else {
			th := math.Acos(z/r) * b.N
			x = rn*math.Sin(th)*math.Cos(ph) + cx
			y = rn*math.Sin(th)*math.Sin(ph) + cy
			z = rn*math.Cos(th) + cz
		}
	}
	return b.MaxIt
}

// In reports membership.
func (b Bulb) In(x, y, z float64) bool { return b.Escape(x, y, z) >= b.MaxIt }

// MultibrotIn is the degree-n Multibrot set zⁿ + c in the plane, for the
// equatorial-slice test.
func MultibrotIn(n float64, cx, cy float64, maxIt int) bool {
	x, y := 0.0, 0.0
	for i := 0; i < maxIt; i++ {
		r := math.Hypot(x, y)
		if r > 2 {
			return false
		}
		if r == 0 {
			x, y = cx, cy
			continue
		}
		rn := math.Pow(r, n)
		a := math.Atan2(y, x) * n
		x, y = rn*math.Cos(a)+cx, rn*math.Sin(a)+cy
	}
	return true
}

// Radius returns the outer surface radius along the ray at latitude lat and
// azimuth phi: the largest r ≤ 1.5 with the point inside, found by marching
// in and bisecting (0 if the ray never meets the set).
func (b Bulb) Radius(lat, phi float64) float64 {
	dx, dy, dz := math.Cos(lat)*math.Cos(phi), math.Cos(lat)*math.Sin(phi), math.Sin(lat)
	const step = 0.004
	out := 1.5
	r := out
	for ; r > 0; r -= step {
		if b.In(r*dx, r*dy, r*dz) {
			break
		}
		out = r
	}
	if r <= 0 {
		return 0
	}
	lo, hi := r, out
	for k := 0; k < 30; k++ {
		m := 0.5 * (lo + hi)
		if b.In(m*dx, m*dy, m*dz) {
			lo = m
		} else {
			hi = m
		}
	}
	return lo
}

// Profile samples the outer radius round a latitude ring.
func (b Bulb) Profile(lat float64, samples int) []float64 {
	p := make([]float64, samples)
	for i := range p {
		p[i] = b.Radius(lat, 2*math.Pi*float64(i)/float64(samples))
	}
	return p
}

// Valley is one cavity cut from a profile, in units of the bulb.
type Valley struct {
	Width, Depth, Area float64 // arc at the outer radius, deepest point, cross-section area
}

// Valleys splits a closed profile, scaled so its largest radius is 1, into
// cavities: the stretches lying deeper than res below the outer circle.
// Each stretch is one cavity, whatever fine structure is inside it.
func Valleys(p []float64, res float64) []Valley {
	n := len(p)
	rmax := 0.0
	for _, r := range p {
		rmax = math.Max(rmax, r)
	}
	if rmax == 0 {
		return nil
	}
	d := make([]float64, n)
	start := -1
	for i, r := range p {
		d[i] = 1 - r/rmax
		if d[i] < res && start < 0 {
			start = i // a point on the outer circle to begin from
		}
	}
	if start < 0 {
		return nil
	}
	dphi := 2 * math.Pi / float64(n)
	var out []Valley
	var cur *Valley
	for k := 1; k <= n; k++ {
		i := (start + k) % n
		if d[i] >= res {
			if cur == nil {
				cur = &Valley{}
			}
			cur.Width += dphi
			cur.Depth = math.Max(cur.Depth, d[i])
			// annular sector between the outer circle and the profile
			ri := 1 - d[i]
			cur.Area += 0.5 * (1 - ri*ri) * dphi
		} else if cur != nil {
			out = append(out, *cur)
			cur = nil
		}
	}
	if cur != nil {
		out = append(out, *cur)
	}
	return out
}

// Lobes counts the stretches of a profile within res of the outer circle
// (the lobes' tips) on the ring.
func Lobes(p []float64, res float64) int {
	return len(Valleys(p, res))
}

var bulbCache sync.Map // key bulbKey → [][]float64 (profiles)

type bulbKey struct {
	n      float64
	lat    bool
	rows   int
	latM   float64
	maxIt  int
	sample int
}

// BulbProfiles returns the outer-radius profiles of rows latitude rings
// between ±latMax, cached (they do not depend on the scale).
func BulbProfiles(b Bulb, rows int, latMax float64) [][]float64 {
	const samples = 720
	k := bulbKey{b.N, b.Latitude, rows, latMax, b.MaxIt, samples}
	if v, ok := bulbCache.Load(k); ok {
		return v.([][]float64)
	}
	out := make([][]float64, rows)
	for i := range out {
		lat := 0.0
		if rows > 1 {
			lat = -latMax + 2*latMax*float64(i)/float64(rows-1)
		}
		out[i] = b.Profile(lat, samples)
	}
	bulbCache.Store(k, out)
	return out
}

// BulbRows returns the valleys of each ring at depth resolution res (a
// fraction of the bulb's outer radius).
func BulbRows(b Bulb, rows int, latMax, res float64) [][]Valley {
	prof := BulbProfiles(b, rows, latMax)
	out := make([][]Valley, rows)
	for i, p := range prof {
		out[i] = Valleys(p, res)
	}
	return out
}

// BulbMachine puts a Mandelbulb band on each drum: rows latitude rings of
// the bulb, scaled so the bulb's outer radius is the band radius r. Each
// valley becomes a cavity: mouth W = arc width × r, depth × r, volume =
// cross-section × r² × axial extent (taken as min(W, row length)),
// projected area = depth × axial extent. Detail smaller than res (metres)
// does not exist in the part.
func BulbMachine(b Bulb, r float64, rows int, latMax, res, ph, gap, rpm, pOut float64) Machine {
	rowLen := 2 * latMax * r / float64(rows)
	if rows == 1 {
		rowLen = 0.2 * r
	}
	m := Machine{Name: "Mandelbulb", R: r, Rows: rows, RowLen: rowLen, Gap: gap, RPM: rpm, POut: pOut,
		Hub: 0.010, Drums: 2}
	vr := BulbRows(b, rows, latMax, res/r)
	m.Rotor = make([][]Cavity, rows)
	for i, row := range vr {
		for _, v := range row {
			w := v.Width * r
			t := math.Min(w, rowLen)
			m.Rotor[i] = append(m.Rotor[i], Cavity{W: w, V: v.Area * r * r * t, A: v.Depth * r * t, Count: 1})
		}
	}
	if ph > 0 {
		// housing dimples Ø 10 mm hemispheres at pitch ph·W
		w := 0.01
		m.Housing = []Cavity{Hemicap(w, w/2, int(2*math.Pi*r/(ph*w)))}
	}
	return m
}
