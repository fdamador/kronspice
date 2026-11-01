// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// Magnetostatic fields and forces of filamentary currents (Biot–Savart;
// J. D. Jackson, Classical Electrodynamics, 3rd ed., §5.2–5.3), used by
// the momentum ledger, the homopolar and tether examples and the
// rotating-field study (M9n).

// SegmentB is the field at p of a straight segment a→b carrying current i
// (closed form: μ₀i/(4πd)·(cos θ₁ − cos θ₂) around the segment). Points on
// the segment's line give 0.
func SegmentB(a, b, p [3]float64, i float64) [3]float64 {
	l := sub3(b, a)
	ll := math.Sqrt(dot(l, l))
	if ll == 0 {
		return [3]float64{}
	}
	u := scale3(l, 1/ll)
	ap, bp := sub3(p, a), sub3(p, b)
	perp := sub3(ap, scale3(u, dot(ap, u)))
	d := math.Sqrt(dot(perp, perp))
	if d < 1e-12*ll {
		return [3]float64{}
	}
	c1 := dot(ap, u) / math.Sqrt(dot(ap, ap))
	c2 := dot(bp, u) / math.Sqrt(dot(bp, bp))
	return scale3(cross(u, scale3(perp, 1/d)), Mu0*i/(4*math.Pi*d)*(c1-c2))
}

// Loop is a closed polygon of points carrying current I (the last point
// joins the first).
type Loop struct {
	P [][3]float64
	I float64
}

// B is the loop's field at p.
func (l Loop) B(p [3]float64) [3]float64 {
	var b [3]float64
	n := len(l.P)
	for k := 0; k < n; k++ {
		b = add3(b, SegmentB(l.P[k], l.P[(k+1)%n], p, l.I))
	}
	return b
}

// ForceIn is the force on the loop in the external field ext,
// ∮ I dl × B, by 8-point Gauss–Legendre quadrature on each segment.
func (l Loop) ForceIn(ext func(p [3]float64) [3]float64) [3]float64 {
	var f [3]float64
	n := len(l.P)
	for k := 0; k < n; k++ {
		a, b := l.P[k], l.P[(k+1)%n]
		dl := sub3(b, a)
		for g := range gl8x {
			p := add3(a, scale3(dl, 0.5*(1+gl8x[g])))
			f = add3(f, scale3(cross(dl, ext(p)), l.I*0.5*gl8w[g]))
		}
	}
	return f
}

// Circle returns an n-gon approximating a circle of radius r about centre
// c in the plane normal to axis (unit), traversed right-handed.
func Circle(c, axis [3]float64, r float64, n int) [][3]float64 {
	a := [3]float64{1, 0, 0}
	if math.Abs(axis[0]) > 0.9 {
		a = [3]float64{0, 1, 0}
	}
	u := cross(axis, a)
	u = scale3(u, 1/math.Sqrt(dot(u, u)))
	w := cross(axis, u)
	pts := make([][3]float64, n)
	for k := range pts {
		s, co := math.Sincos(2 * math.Pi * float64(k) / float64(n))
		pts[k] = add3(c, add3(scale3(u, r*co), scale3(w, r*s)))
	}
	return pts
}

var gl8x = [8]float64{-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363}
var gl8w = [8]float64{0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
	0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763}

func sub3(a, b [3]float64) [3]float64 { return [3]float64{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func add3(a, b [3]float64) [3]float64 { return [3]float64{a[0] + b[0], a[1] + b[1], a[2] + b[2]} }
func scale3(a [3]float64, s float64) [3]float64 {
	return [3]float64{a[0] * s, a[1] * s, a[2] * s}
}
