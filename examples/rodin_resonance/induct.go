// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/fdamador/kronspice/em"
)

// Partial-inductance (segment) method: a wire is a chain of straight
// segments; L_ij = μ0/4π ∫∫ dl_i·dl_j/|r_i − r_j| (Neumann) for i ≠ j and
// the straight-wire self term μ0ℓ/2π·(ln(2ℓ/r) − 3/4) (low frequency,
// internal inductance included; F. W. Grover, Inductance Calculations,
// 1946, ch. 6; A. E. Ruehli, IBM J. Res. Dev. 16, 470, 1972). Near pairs
// are split recursively until each piece is short against its distance;
// far pairs use Gauss–Legendre quadrature (8 or 2 points per segment).

var gx8 = [8]float64{-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363}
var gw8 = [8]float64{0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
	0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763}
var gx2 = [2]float64{-0.5773502691896258, 0.5773502691896258}

func sub(a, b P3) P3         { return P3{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func add(a, b P3) P3         { return P3{a[0] + b[0], a[1] + b[1], a[2] + b[2]} }
func scl(a P3, s float64) P3 { return P3{a[0] * s, a[1] * s, a[2] * s} }
func dot(a, b P3) float64    { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
func norm(a P3) float64      { return math.Sqrt(dot(a, a)) }
func cross(a, b P3) P3 {
	return P3{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}
func mid(a, b P3) P3 { return scl(add(a, b), 0.5) }
func unit(a P3) P3   { return scl(a, 1/norm(a)) }

// segMutual is the mutual partial inductance of segments a→b and c→d
// divided by μ0/4π.
func segMutual(a, b, c, d P3, depth int) float64 {
	d1, d2 := sub(b, a), sub(d, c)
	dd := dot(d1, d2)
	if dd == 0 {
		return 0
	}
	l1, l2 := norm(d1), norm(d2)
	lm := math.Max(l1, l2)
	sep := norm(sub(mid(a, b), mid(c, d))) - 0.5*(l1+l2)
	if sep < lm && depth < 14 {
		if l1 >= l2 {
			m := mid(a, b)
			return segMutual(a, m, c, d, depth+1) + segMutual(m, b, c, d, depth+1)
		}
		m := mid(c, d)
		return segMutual(a, b, c, m, depth+1) + segMutual(a, b, m, d, depth+1)
	}
	s := 0.0
	if sep > 4*lm {
		for p := range gx2 {
			x := add(a, scl(d1, 0.5*(1+gx2[p])))
			for q := range gx2 {
				y := add(c, scl(d2, 0.5*(1+gx2[q])))
				s += 0.25 / norm(sub(x, y))
			}
		}
		return s * dd
	}
	for p := range gx8 {
		x := add(a, scl(d1, 0.5*(1+gx8[p])))
		for q := range gx8 {
			y := add(c, scl(d2, 0.5*(1+gx8[q])))
			s += 0.25 * gw8[p] * gw8[q] / norm(sub(x, y))
		}
	}
	return s * dd
}

// selfSeg is a straight wire's partial self inductance / (μ0/4π).
func selfSeg(l, r float64) float64 { return 2 * l * (math.Log(2*l/r) - 0.75) }

// Segs is a list of segments with their section labels.
type Segs struct {
	A, B []P3 // endpoints
	Sec  []int
	R    []float64 // wire radius per segment
}

// AddPolyline appends the segments of an open polyline, labelled by
// sections numbered from base (nsec sections of equal segment count).
func (s *Segs) AddPolyline(p []P3, base, nsec int, r float64) {
	n := len(p) - 1
	for i := 0; i < n; i++ {
		s.A = append(s.A, p[i])
		s.B = append(s.B, p[i+1])
		s.Sec = append(s.Sec, base+i*nsec/n)
		s.R = append(s.R, r)
	}
}

// SectionMatrix returns the section inductance matrix [H] (nsec × nsec),
// summed over segment pairs in parallel.
func (s *Segs) SectionMatrix(nsec int) [][]float64 {
	n := len(s.A)
	workers := runtime.NumCPU()
	parts := make([][][]float64, workers)
	var wg sync.WaitGroup
	next := make(chan int, 64)
	go func() {
		for i := 0; i < n; i++ {
			next <- i
		}
		close(next)
	}()
	for w := 0; w < workers; w++ {
		m := newMat(nsec)
		parts[w] = m
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range next {
				si := s.Sec[i]
				m[si][si] += selfSeg(norm(sub(s.B[i], s.A[i])), s.R[i])
				for j := i + 1; j < n; j++ {
					v := segMutual(s.A[i], s.B[i], s.A[j], s.B[j], 0)
					sj := s.Sec[j]
					m[si][sj] += v
					m[sj][si] += v
				}
			}
		}()
	}
	wg.Wait()
	out := newMat(nsec)
	for _, m := range parts {
		for i := range m {
			for j := range m[i] {
				out[i][j] += m[i][j] * em.Mu0 / (4 * math.Pi)
			}
		}
	}
	return out
}

// MutualTo is the mutual inductance [H] between every section of s and a
// closed loop (polygon) carrying the same current in all its turns.
func (s *Segs) MutualTo(loop []P3, turns float64, nsec int) []float64 {
	out := make([]float64, nsec)
	nl := len(loop)
	for i := range s.A {
		v := 0.0
		for k := 0; k < nl; k++ {
			v += segMutual(s.A[i], s.B[i], loop[k], loop[(k+1)%nl], 0)
		}
		out[s.Sec[i]] += v * turns * em.Mu0 / (4 * math.Pi)
	}
	return out
}

func newMat(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	return m
}

// FieldOf is the field [T] at p of a polyline carrying current i.
func FieldOf(pl []P3, i float64, p P3) P3 {
	var b P3
	for k := 0; k+1 < len(pl); k++ {
		b = add(b, em.SegmentB(pl[k], pl[k+1], p, i))
	}
	return b
}

// LoopInductance is a closed loop's self inductance [H] by the segment
// method (wire radius r).
func LoopInductance(loop []P3, r float64) float64 {
	var s Segs
	p := append(append([]P3(nil), loop...), loop[0])
	s.AddPolyline(p, 0, 1, r)
	return s.SectionMatrix(1)[0][0]
}

// CircleRing is the low-frequency self inductance of a thin round ring
// (Kirchhoff; Grover 1946, ch. 13): μ0R(ln(8R/r) − 7/4).
func CircleRing(R, r float64) float64 { return em.Mu0 * R * (math.Log(8*R/r) - 1.75) }

func numWorkers() int { return runtime.NumCPU() }
