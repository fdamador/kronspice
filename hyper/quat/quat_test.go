// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package quat

import (
	"math"
	"math/rand"
	"testing"
)

func randQ(r *rand.Rand) Q {
	return Q{r.NormFloat64(), r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
}

func near(a, b Q, tol float64) bool {
	return a.Sub(b).Abs() <= tol*math.Max(1, math.Max(a.Abs(), b.Abs()))
}

// Hamilton's defining relations (Elements, art. 183): i² = j² = k² = ijk = −1,
// ij = k = −ji, jk = i = −kj, ki = j = −ik.
func TestHamiltonRelations(t *testing.T) {
	i, j, k, one := Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, 1), Scalar(1)
	m1 := one.Scale(-1)
	for _, c := range []struct {
		name      string
		got, want Q
	}{
		{"i²", i.Mul(i), m1}, {"j²", j.Mul(j), m1}, {"k²", k.Mul(k), m1},
		{"ijk", i.Mul(j).Mul(k), m1},
		{"ij", i.Mul(j), k}, {"ji", j.Mul(i), k.Scale(-1)},
		{"jk", j.Mul(k), i}, {"kj", k.Mul(j), i.Scale(-1)},
		{"ki", k.Mul(i), j}, {"ik", i.Mul(k), j.Scale(-1)},
	} {
		if c.got != c.want {
			t.Errorf("%s = %v, want %v", c.name, c.got, c.want)
		}
	}
}

// Algebra laws on random quaternions: associativity, distributivity,
// |pq| = |p||q|, (pq)* = q*p*, p·p⁻¹ = 1, and for pure quaternions
// uv = −u·v + u×v.
func TestQuaternionLaws(t *testing.T) {
	r := rand.New(rand.NewSource(1))
	for n := 0; n < 2000; n++ {
		p, q, s := randQ(r), randQ(r), randQ(r)
		if !near(p.Mul(q).Mul(s), p.Mul(q.Mul(s)), 1e-14) {
			t.Fatalf("associativity fails for %v %v %v", p, q, s)
		}
		if !near(p.Mul(q.Add(s)), p.Mul(q).Add(p.Mul(s)), 1e-14) {
			t.Fatalf("distributivity fails")
		}
		if d := math.Abs(p.Mul(q).Abs() - p.Abs()*q.Abs()); d > 1e-14*p.Abs()*q.Abs() {
			t.Fatalf("|pq| ≠ |p||q| by %g", d)
		}
		if !near(p.Mul(q).Conj(), q.Conj().Mul(p.Conj()), 1e-15) {
			t.Fatalf("(pq)* ≠ q*p*")
		}
		pi, _ := p.Inv()
		if !near(p.Mul(pi), Scalar(1), 1e-15) || !near(pi.Mul(p), Scalar(1), 1e-15) {
			t.Fatalf("p·p⁻¹ ≠ 1")
		}
		u, v := Vec(p.X, p.Y, p.Z), Vec(q.X, q.Y, q.Z)
		dot := u.X*v.X + u.Y*v.Y + u.Z*v.Z
		cross := Vec(u.Y*v.Z-u.Z*v.Y, u.Z*v.X-u.X*v.Z, u.X*v.Y-u.Y*v.X)
		if !near(u.Mul(v), Scalar(-dot).Add(cross), 1e-15) {
			t.Fatalf("uv ≠ −u·v + u×v")
		}
	}
	if _, ok := (Q{}).Inv(); ok {
		t.Error("0 has an inverse")
	}
}

// Rotations: a rotor turns i into j by π/2 about k; r·v·r* agrees with the
// rotation matrix, preserves length, and the matrix round-trips through
// FromMatrix; exp and log are inverse; Slerp ends at its endpoints and keeps
// unit length.
func TestRotors(t *testing.T) {
	r90, _ := Rotor(0, 0, 1, math.Pi/2)
	x, y, z := r90.Rotate(1, 0, 0)
	if math.Abs(x) > 1e-15 || math.Abs(y-1) > 1e-15 || math.Abs(z) > 1e-15 {
		t.Errorf("rotating i by π/2 about k gives (%g, %g, %g), want j", x, y, z)
	}
	rng := rand.New(rand.NewSource(2))
	for n := 0; n < 1000; n++ {
		q, _ := randQ(rng).Unit()
		v := [3]float64{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		a, b, c := q.Rotate(v[0], v[1], v[2])
		m := q.Matrix()
		for i, w := range [3]float64{a, b, c} {
			mv := m[i][0]*v[0] + m[i][1]*v[1] + m[i][2]*v[2]
			if math.Abs(mv-w) > 1e-14*(1+math.Abs(w)) {
				t.Fatalf("matrix and r·v·r* differ: %g vs %g", mv, w)
			}
		}
		if d := math.Abs(math.Sqrt(a*a+b*b+c*c) - math.Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])); d > 1e-14*(1+math.Abs(a)) {
			t.Fatalf("rotation changed the length by %g", d)
		}
		back := FromMatrix(m)
		if q.W < 0 {
			q = q.Scale(-1)
		}
		if !near(back, q, 1e-14) {
			t.Fatalf("FromMatrix(Matrix(q)) = %v, want %v", back, q)
		}
		p := randQ(rng)
		if l, ok := p.Log(); !ok || !near(l.Exp(), p, 1e-13) {
			t.Fatalf("exp(log p) ≠ p for %v", p)
		}
		q2, _ := randQ(rng).Unit()
		if s := Slerp(q, q2, 0); !near(s, q, 1e-14) {
			t.Fatalf("Slerp(·,·,0) = %v, want %v", s, q)
		}
		if s := Slerp(q, q2, 0.37); math.Abs(s.Abs()-1) > 1e-14 {
			t.Fatalf("Slerp left the unit sphere: %g", s.Abs())
		}
	}
}
