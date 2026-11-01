// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"

	"github.com/fdamador/kronspice/hyper/ga"
	"github.com/fdamador/kronspice/hyper/quat"
)

// Frame changes of the field (CLAUDE.md §3c rule 3: parity and moving
// observer). Both are computed twice — by the textbook vector rules and by
// geometric algebra — so a frame bookkeeping error cannot pass as a
// difference between the quaternion/GA form and Maxwell's equations.

// TransformVector applies an orthogonal 3×3 matrix (proper or improper) to
// E and B by the vector-calculus rules: E is a polar vector, E → M·E; B is
// an axial vector (a cross product), B → det(M)·M·B.
func TransformVector(m [3][3]float64, e, b [3]float64) ([3]float64, [3]float64) {
	d := det3(m)
	me, mb := mulVec(m, e), mulVec(m, b)
	for k := range mb {
		mb[k] *= d
	}
	return me, mb
}

// TransformGA applies the same orthogonal matrix in G3, where B is carried
// as the bivector I·B: M is factored into reflections (Householder), each
// unit normal n acting as v → −n·v·n on vectors and X → n·X·n on bivectors
// (Dorst et al. 2007, §7.1; Doran & Lasenby §2.7), and B is read back as
// −I·(I·B)′. The determinant never appears: it is the parity of the number
// of reflections, carried by the grade.
func TransformGA(m [3][3]float64, e, b [3]float64) ([3]float64, [3]float64) {
	E := ga.Vector3(e[0], e[1], e[2])
	I := ga.G3{ga.G3I: 1}
	Bb := I.Mul(ga.Vector3(b[0], b[1], b[2])) // I·B, a bivector
	for _, n := range reflections(m) {
		nv := ga.Vector3(n[0], n[1], n[2])
		E = nv.Mul(E).Mul(nv).Scale(-1)
		Bb = nv.Mul(Bb).Mul(nv)
	}
	ex, ey, ez := E.Vec()
	bx, by, bz := I.Mul(Bb).Scale(-1).Vec()
	return [3]float64{ex, ey, ez}, [3]float64{bx, by, bz}
}

// reflections factors an orthogonal M into unit normals n₁ … n_m with
// M = R(n₁)·R(n₂)·…·R(n_m), R(n) = 1 − 2nnᵀ, returned in the order they act
// on a vector (n_m first).
func reflections(m [3][3]float64) [][3]float64 {
	// Householder: H_k … H_1·M = D (diagonal ±1, M being orthogonal), so
	// M = H_1·…·H_k·D, and D is a product of axis reflections.
	a := m
	var hs [][3]float64
	for col := 0; col < 2; col++ {
		var x [3]float64
		for r := col; r < 3; r++ {
			x[r] = a[r][col]
		}
		alpha := -math.Copysign(math.Sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]), x[col])
		v := x
		v[col] -= alpha
		nv := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
		if nv < 1e-14 {
			continue
		}
		for r := range v {
			v[r] /= nv
		}
		// a ← (1 − 2vvᵀ)·a
		for c := 0; c < 3; c++ {
			dot := v[0]*a[0][c] + v[1]*a[1][c] + v[2]*a[2][c]
			for r := 0; r < 3; r++ {
				a[r][c] -= 2 * v[r] * dot
			}
		}
		hs = append(hs, v)
	}
	// M = H_1·H_2·D: acting on a vector, D first, then H_2, then H_1.
	var out [][3]float64
	for k := 0; k < 3; k++ {
		if a[k][k] < 0 {
			var n [3]float64
			n[k] = 1
			out = append(out, n)
		}
	}
	for k := len(hs) - 1; k >= 0; k-- {
		out = append(out, hs[k])
	}
	return out
}

func det3(m [3][3]float64) float64 {
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
}

func mulVec(m [3][3]float64, v [3]float64) [3]float64 {
	return [3]float64{
		m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2],
		m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2],
		m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2],
	}
}

// BoostVector returns the fields seen by an observer moving with velocity
// v (|v| < c), by the textbook transformation (J. D. Jackson, Classical
// Electrodynamics, 3rd ed., eq. 11.149):
//
//	E′ = γ(E + v×B) − γ²/(γ+1)·v(v·E)/c²
//	B′ = γ(B − v×E/c²) − γ²/(γ+1)·v(v·B)/c²
func BoostVector(e, b, v [3]float64) ([3]float64, [3]float64) {
	c2 := SpeedOfLight * SpeedOfLight
	v2 := v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
	g := 1 / math.Sqrt(1-v2/c2)
	vxb, vxe := cross(v, b), cross(v, e)
	ve, vb := dot(v, e), dot(v, b)
	var e2, b2 [3]float64
	for k := 0; k < 3; k++ {
		e2[k] = g*(e[k]+vxb[k]) - g*g/(g+1)*v[k]*ve/c2
		b2[k] = g*(b[k]-vxe[k]/c2) - g*g/(g+1)*v[k]*vb/c2
	}
	return e2, b2
}

// Boost returns the fields seen by an observer moving with velocity v by
// the spacetime algebra: F′ = R̃·F·R with the boost rotor
// R = exp(α·v̂/2), v̂ = Σ v̂ₖσₖ the unit relative vector, tanh α = |v|/c
// (Doran & Lasenby §5.4: a boost is a rotor whose plane is timelike;
// the passive change of observer is the reverse of the active boost).
func (s *STA) Boost(e, b, v [3]float64) ([3]float64, [3]float64) {
	a := s.Alg
	speed := math.Sqrt(dot(v, v))
	if speed == 0 {
		return e, b
	}
	alpha := math.Atanh(speed / SpeedOfLight)
	gen := ga.Multivector{}
	for k := 0; k < 3; k++ {
		gen = ga.Add(gen, s.Sigma[k].Scale(v[k]/speed*alpha/2))
	}
	R := a.Exp(gen)
	F := s.Field(vecQ(e), vecQ(b))
	F2 := a.Mul(a.Mul(R.Reverse(), F), R)
	var e2, b2 [3]float64
	for k := 0; k < 3; k++ {
		e2[k] = a.ScalarProduct(F2, s.Sigma[k]) / a.ScalarProduct(s.Sigma[k], s.Sigma[k])
		b2[k] = a.ScalarProduct(F2, s.ISig[k]) / a.ScalarProduct(s.ISig[k], s.ISig[k]) / SpeedOfLight
	}
	return e2, b2
}

func cross(a, b [3]float64) [3]float64 {
	return [3]float64{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

func dot(a, b [3]float64) float64 { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }

func vecQ(v [3]float64) quat.Q { return quat.Vec(v[0], v[1], v[2]) }
