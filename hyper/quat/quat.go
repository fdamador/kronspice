// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package quat is Hamilton's quaternions: q = w + x·i + y·j + z·k with
// i² = j² = k² = ijk = −1 (W. R. Hamilton, Elements of Quaternions, 1866,
// Book II, Ch. I, arts. 181–183). A pure quaternion x·i + y·j + z·k is a
// 3-vector; the product of two pure quaternions is
//
//	u·v = −u·v (dot) + u×v,
//
// so its scalar part is minus the dot product and its vector part the cross
// product — the sign convention Maxwell uses (Treatise, Vol. II, art. 11
// and art. 25). A unit quaternion r rotates a vector v by r·v·r* (Hamilton,
// Elements, art. 179): r = cos(θ/2) + sin(θ/2)·n turns v by θ about the unit
// axis n, right-handed.
//
// The package depends on the standard library only (the project keeps no
// third-party dependencies; gonum's num/quat has the same algebra).
package quat

import "math"

// Q is a quaternion W + X·i + Y·j + Z·k.
type Q struct{ W, X, Y, Z float64 }

// Vec is a pure quaternion (scalar part 0) from a 3-vector.
func Vec(x, y, z float64) Q { return Q{0, x, y, z} }

// Scalar is the quaternion with only a scalar part.
func Scalar(w float64) Q { return Q{W: w} }

// Add returns p + q.
func (p Q) Add(q Q) Q { return Q{p.W + q.W, p.X + q.X, p.Y + q.Y, p.Z + q.Z} }

// Sub returns p − q.
func (p Q) Sub(q Q) Q { return Q{p.W - q.W, p.X - q.X, p.Y - q.Y, p.Z - q.Z} }

// Scale returns s·p.
func (p Q) Scale(s float64) Q { return Q{s * p.W, s * p.X, s * p.Y, s * p.Z} }

// Mul is Hamilton's product p·q (non-commutative).
func (p Q) Mul(q Q) Q {
	return Q{
		p.W*q.W - p.X*q.X - p.Y*q.Y - p.Z*q.Z,
		p.W*q.X + p.X*q.W + p.Y*q.Z - p.Z*q.Y,
		p.W*q.Y - p.X*q.Z + p.Y*q.W + p.Z*q.X,
		p.W*q.Z + p.X*q.Y - p.Y*q.X + p.Z*q.W,
	}
}

// Conj is the conjugate p* = W − X·i − Y·j − Z·k (Hamilton's K.q).
func (p Q) Conj() Q { return Q{p.W, -p.X, -p.Y, -p.Z} }

// Norm2 is p·p* = W² + X² + Y² + Z² (the square of Hamilton's tensor T.q).
func (p Q) Norm2() float64 { return p.W*p.W + p.X*p.X + p.Y*p.Y + p.Z*p.Z }

// Abs is the tensor |p|, computed without overflow for large components.
func (p Q) Abs() float64 { return math.Hypot(math.Hypot(p.W, p.X), math.Hypot(p.Y, p.Z)) }

// Inv returns p⁻¹ = p*/|p|², and false for p = 0.
func (p Q) Inv() (Q, bool) {
	n := p.Norm2()
	if n == 0 || math.IsNaN(n) {
		return Q{}, false
	}
	return p.Conj().Scale(1 / n), true
}

// Unit returns p/|p| (Hamilton's versor U.q), and false for p = 0.
func (p Q) Unit() (Q, bool) {
	a := p.Abs()
	if a == 0 || math.IsNaN(a) {
		return Q{}, false
	}
	return p.Scale(1 / a), true
}

// Real is the scalar part (Hamilton's S.q).
func (p Q) Real() float64 { return p.W }

// Imag is the vector part (Hamilton's V.q) as x, y, z.
func (p Q) Imag() (x, y, z float64) { return p.X, p.Y, p.Z }

// Dot is the Euclidean inner product of the four components,
// W·W′ + X·X′ + Y·Y′ + Z·Z′ (= S.(p·q*)).
func (p Q) Dot(q Q) float64 { return p.W*q.W + p.X*q.X + p.Y*q.Y + p.Z*q.Z }

// Exp returns e^p = e^W·(cos|v| + sin|v|·v/|v|), v the vector part.
func (p Q) Exp() Q {
	vn := math.Sqrt(p.X*p.X + p.Y*p.Y + p.Z*p.Z)
	ew := math.Exp(p.W)
	if vn == 0 {
		return Q{W: ew}
	}
	s := ew * math.Sin(vn) / vn
	return Q{ew * math.Cos(vn), s * p.X, s * p.Y, s * p.Z}
}

// Log returns the principal logarithm ln|p| + (v/|v|)·atan2(|v|, W), and
// false for p = 0. For a negative real p the axis is undefined; i is used.
func (p Q) Log() (Q, bool) {
	a := p.Abs()
	if a == 0 || math.IsNaN(a) {
		return Q{}, false
	}
	vn := math.Sqrt(p.X*p.X + p.Y*p.Y + p.Z*p.Z)
	if vn == 0 {
		if p.W > 0 {
			return Q{W: math.Log(a)}, true
		}
		return Q{math.Log(a), math.Pi, 0, 0}, true
	}
	s := math.Atan2(vn, p.W) / vn
	return Q{math.Log(a), s * p.X, s * p.Y, s * p.Z}, true
}

// Rotor returns the unit quaternion turning vectors by angle θ (radians,
// right-handed) about the axis (x, y, z), and false for a zero axis.
func Rotor(x, y, z, theta float64) (Q, bool) {
	n, ok := Vec(x, y, z).Unit()
	if !ok {
		return Q{}, false
	}
	s, c := math.Sincos(theta / 2)
	return Q{c, s * n.X, s * n.Y, s * n.Z}, true
}

// Rotate returns r·v·r* for the vector (x, y, z); r should be a unit
// quaternion (for any non-zero r this is a rotation scaled by |r|²).
func (r Q) Rotate(x, y, z float64) (float64, float64, float64) {
	v := r.Mul(Vec(x, y, z)).Mul(r.Conj())
	return v.X, v.Y, v.Z
}

// Matrix returns the 3×3 rotation matrix M (row-major) with M·v = r·v·r*
// for a unit quaternion r.
func (r Q) Matrix() [3][3]float64 {
	w, x, y, z := r.W, r.X, r.Y, r.Z
	return [3][3]float64{
		{w*w + x*x - y*y - z*z, 2 * (x*y - w*z), 2 * (x*z + w*y)},
		{2 * (x*y + w*z), w*w - x*x + y*y - z*z, 2 * (y*z - w*x)},
		{2 * (x*z - w*y), 2 * (y*z + w*x), w*w - x*x - y*y + z*z},
	}
}

// FromMatrix returns a unit quaternion r (with W ≥ 0) whose rotation matrix
// is m, a proper orthogonal matrix, by Shepperd's method (S. W. Shepperd,
// "Quaternion from rotation matrix," J. Guidance and Control 1(3), 1978).
func FromMatrix(m [3][3]float64) Q {
	tr := m[0][0] + m[1][1] + m[2][2]
	var q Q
	switch {
	case tr >= m[0][0] && tr >= m[1][1] && tr >= m[2][2]:
		s := 2 * math.Sqrt(1+tr)
		q = Q{s / 4, (m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s}
	case m[0][0] >= m[1][1] && m[0][0] >= m[2][2]:
		s := 2 * math.Sqrt(1+m[0][0]-m[1][1]-m[2][2])
		q = Q{(m[2][1] - m[1][2]) / s, s / 4, (m[0][1] + m[1][0]) / s, (m[0][2] + m[2][0]) / s}
	case m[1][1] >= m[2][2]:
		s := 2 * math.Sqrt(1+m[1][1]-m[0][0]-m[2][2])
		q = Q{(m[0][2] - m[2][0]) / s, (m[0][1] + m[1][0]) / s, s / 4, (m[1][2] + m[2][1]) / s}
	default:
		s := 2 * math.Sqrt(1+m[2][2]-m[0][0]-m[1][1])
		q = Q{(m[1][0] - m[0][1]) / s, (m[0][2] + m[2][0]) / s, (m[1][2] + m[2][1]) / s, s / 4}
	}
	if q.W < 0 {
		q = q.Scale(-1)
	}
	return q
}

// Slerp interpolates between unit quaternions p and q along the shorter
// great arc, t ∈ [0, 1] (K. Shoemake, "Animating rotation with quaternion
// curves," SIGGRAPH 1985).
func Slerp(p, q Q, t float64) Q {
	d := p.Dot(q)
	if d < 0 {
		q, d = q.Scale(-1), -d
	}
	if d > 1-1e-12 {
		r, _ := p.Scale(1 - t).Add(q.Scale(t)).Unit()
		return r
	}
	th := math.Acos(d)
	s := math.Sin(th)
	return p.Scale(math.Sin((1-t)*th) / s).Add(q.Scale(math.Sin(t*th) / s))
}
