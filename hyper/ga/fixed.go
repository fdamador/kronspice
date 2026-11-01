// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package ga

import (
	"math"

	"github.com/fdamador/kronspice/hyper/quat"
)

// G2 is a multivector of the Euclidean plane Cl(2,0,0), indexed by blade
// mask: [1, e₁, e₂, e₁₂]. Its even part S + B·e₁₂ is the complex number
// S + B·i (e₁₂² = −1). Values, allocation-free.
type G2 [4]float64

// G3 is a multivector of Euclidean 3-space Cl(3,0,0), indexed by blade
// mask: [1, e₁, e₂, e₁₂, e₃, e₁₃, e₂₃, e₁₂₃]. Values, allocation-free.
type G3 [8]float64

// Blade indices of G3, by mask.
const (
	G3S   = 0
	G3E1  = 1
	G3E2  = 2
	G3E12 = 3
	G3E3  = 4
	G3E13 = 5
	G3E23 = 6
	G3I   = 7
)

var sign2 [4][4]float64
var sign3 [8][8]float64

func init() {
	e2, _ := New(2, 0, 0)
	for i := range sign2 {
		for j := range sign2[i] {
			sign2[i][j], _ = e2.bladeMul(Blade(i), Blade(j))
		}
	}
	e3, _ := New(3, 0, 0)
	for i := range sign3 {
		for j := range sign3[i] {
			sign3[i][j], _ = e3.bladeMul(Blade(i), Blade(j))
		}
	}
}

// Mul is the geometric product.
func (x G2) Mul(y G2) G2 {
	var z G2
	for i, a := range x {
		if a == 0 {
			continue
		}
		for j, b := range y {
			z[i^j] += sign2[i][j] * a * b
		}
	}
	return z
}

// Reverse is x̃ (e₁₂ → −e₁₂).
func (x G2) Reverse() G2 { return G2{x[0], x[1], x[2], -x[3]} }

// Complex is the even part as a complex number S + B·i.
func (x G2) Complex() complex128 { return complex(x[0], x[3]) }

// G2FromComplex is the even multivector re + im·e₁₂.
func G2FromComplex(z complex128) G2 { return G2{real(z), 0, 0, imag(z)} }

// Multivector converts to the sparse form in Cl(2,0,0).
func (x G2) Multivector() Multivector { return fromArray(x[:]) }

// Mul is the geometric product.
func (x G3) Mul(y G3) G3 {
	var z G3
	for i, a := range x {
		if a == 0 {
			continue
		}
		for j, b := range y {
			z[i^j] += sign3[i][j] * a * b
		}
	}
	return z
}

// Add returns x + y.
func (x G3) Add(y G3) G3 {
	for i := range x {
		x[i] += y[i]
	}
	return x
}

// Scale returns s·x.
func (x G3) Scale(s float64) G3 {
	for i := range x {
		x[i] *= s
	}
	return x
}

// Reverse is x̃: bivectors and the trivector change sign.
func (x G3) Reverse() G3 {
	return G3{x[0], x[1], x[2], -x[3], x[4], -x[5], -x[6], -x[7]}
}

// Grade returns ⟨x⟩ₖ.
func (x G3) Grade(k int) G3 {
	var z G3
	for i, v := range x {
		if Blade(i).Grade() == k {
			z[i] = v
		}
	}
	return z
}

// Outer is x∧y.
func (x G3) Outer(y G3) G3 {
	var z G3
	for i, a := range x {
		for j, b := range y {
			if i&j == 0 {
				z[i^j] += sign3[i][j] * a * b
			}
		}
	}
	return z
}

// Vector3 is the vector x·e₁ + y·e₂ + z·e₃.
func Vector3(x, y, z float64) G3 { return G3{G3E1: x, G3E2: y, G3E3: z} }

// Vec returns the vector part as x, y, z.
func (x G3) Vec() (float64, float64, float64) { return x[G3E1], x[G3E2], x[G3E3] }

// Rotor3 returns the rotor R = cos(θ/2) − sin(θ/2)·I·n that turns vectors by
// θ about the axis n = (x, y, z) through R·v·R̃, and false for a zero axis.
// I·n = n₁e₂₃ + n₂e₃₁ + n₃e₁₂ = n₁e₂₃ − n₂e₁₃ + n₃e₁₂.
func Rotor3(x, y, z, theta float64) (G3, bool) {
	n := math.Sqrt(x*x + y*y + z*z)
	if n == 0 || math.IsNaN(n) {
		return G3{}, false
	}
	s, c := math.Sincos(theta / 2)
	s /= n
	return G3{G3S: c, G3E23: -s * x, G3E13: s * y, G3E12: -s * z}, true
}

// Sandwich returns r·x·r̃.
func (r G3) Sandwich(x G3) G3 { return r.Mul(x).Mul(r.Reverse()) }

// Quat returns the even part as a quaternion under i ↔ −e₂₃, j ↔ e₁₃,
// k ↔ −e₁₂ (package documentation). The odd part is ignored.
func (x G3) Quat() quat.Q {
	return quat.Q{W: x[G3S], X: -x[G3E23], Y: x[G3E13], Z: -x[G3E12]}
}

// G3FromQuat is the even multivector of a quaternion (the inverse of Quat).
func G3FromQuat(q quat.Q) G3 {
	return G3{G3S: q.W, G3E23: -q.X, G3E13: q.Y, G3E12: -q.Z}
}

// Multivector converts to the sparse form in Cl(3,0,0).
func (x G3) Multivector() Multivector { return fromArray(x[:]) }

// G3From converts a sparse Cl(3,0,0) multivector; false if it has blades
// outside the algebra.
func G3From(m Multivector) (G3, bool) {
	var x G3
	for _, t := range m.t {
		if t.b > 7 {
			return G3{}, false
		}
		x[t.b] = t.c
	}
	return x, true
}

func fromArray(c []float64) Multivector {
	var m Multivector
	for i, v := range c {
		if v != 0 {
			m.t = append(m.t, term{Blade(i), v})
		}
	}
	return m
}
