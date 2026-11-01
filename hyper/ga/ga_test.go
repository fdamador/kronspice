// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package ga

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/hyper/quat"
)

func mustNew(t *testing.T, p, q, r int) *Algebra {
	t.Helper()
	a, err := New(p, q, r)
	if err != nil {
		t.Fatal(err)
	}
	return a
}

// randMV is a random multivector with every blade present.
func randMV(r *rand.Rand, a *Algebra) Multivector {
	c := map[Blade]float64{}
	for b := Blade(0); b <= a.full; b++ {
		c[b] = r.NormFloat64()
	}
	m, _ := a.FromBlades(c)
	return m
}

func closeMV(x, y Multivector, tol float64) bool {
	return Sub(x, y).MaxAbs() <= tol*math.Max(1, math.Max(x.MaxAbs(), y.MaxAbs()))
}

// Algebra laws in several signatures, including degenerate ones:
// associativity, distributivity, the reverse of a product (xy)~ = ỹx̃,
// eₖ² = metric, and anticommutation of distinct basis vectors.
func TestAlgebraLaws(t *testing.T) {
	r := rand.New(rand.NewSource(3))
	for _, sig := range [][3]int{{2, 0, 0}, {3, 0, 0}, {0, 2, 0}, {1, 3, 0}, {3, 0, 1}, {2, 2, 1}, {4, 1, 0}} {
		a := mustNew(t, sig[0], sig[1], sig[2])
		for n := 0; n < 40; n++ {
			x, y, z := randMV(r, a), randMV(r, a), randMV(r, a)
			if !closeMV(a.Mul(a.Mul(x, y), z), a.Mul(x, a.Mul(y, z)), 1e-12) {
				t.Fatalf("Cl%v: associativity fails", sig)
			}
			if !closeMV(a.Mul(x, Add(y, z)), Add(a.Mul(x, y), a.Mul(x, z)), 1e-13) {
				t.Fatalf("Cl%v: distributivity fails", sig)
			}
			if !closeMV(a.Mul(x, y).Reverse(), a.Mul(y.Reverse(), x.Reverse()), 1e-13) {
				t.Fatalf("Cl%v: (xy)~ ≠ ỹx̃", sig)
			}
			if !closeMV(a.Mul(x, y).Conjugate(), a.Mul(y.Conjugate(), x.Conjugate()), 1e-13) {
				t.Fatalf("Cl%v: conjugate is not an anti-automorphism", sig)
			}
			if !closeMV(a.Mul(x, y).GradeInvolution(), a.Mul(x.GradeInvolution(), y.GradeInvolution()), 1e-13) {
				t.Fatalf("Cl%v: grade involution is not an automorphism", sig)
			}
		}
		for i := 1; i <= a.Dim(); i++ {
			ei, _ := a.Basis(i)
			want := a.sq[i-1]
			if got := a.Mul(ei, ei); !closeMV(got, Scalar(want), 0) {
				t.Fatalf("Cl%v: e%d² = %v, want %g", sig, i, got, want)
			}
			for j := i + 1; j <= a.Dim(); j++ {
				ej, _ := a.Basis(j)
				if !closeMV(a.Mul(ei, ej), a.Mul(ej, ei).Scale(-1), 0) {
					t.Fatalf("Cl%v: e%de%d ≠ −e%de%d", sig, i, j, j, i)
				}
			}
		}
	}
	if _, err := New(0, 0, 0); err == nil {
		t.Error("Cl(0,0,0) accepted")
	}
	if _, err := New(40, 30, 0); err == nil {
		t.Error("n = 70 accepted")
	}
}

// Outer-product rules: grade(A∧B) = grade A + grade B for blades (or 0),
// a∧a = 0 for vectors, and the wedge of n vectors in Cl(n) is det[v]·I —
// the identity kron/gaframe builds on.
func TestOuterProduct(t *testing.T) {
	r := rand.New(rand.NewSource(4))
	a := mustNew(t, 5, 0, 0)
	vec := func() Multivector {
		c := make([]float64, 5)
		for i := range c {
			c[i] = r.NormFloat64()
		}
		v, _ := a.Vector(c...)
		return v
	}
	for n := 0; n < 200; n++ {
		u, v, w := vec(), vec(), vec()
		if g := a.OuterAll(u, v).Grades(); len(g) != 1 || g[0] != 2 {
			t.Fatalf("u∧v has grades %v", g)
		}
		if g := a.OuterAll(u, v, w).Grades(); len(g) != 1 || g[0] != 3 {
			t.Fatalf("u∧v∧w has grades %v", g)
		}
		if x := a.Outer(u, u); x.MaxAbs() > 1e-15 {
			t.Fatalf("u∧u = %v", x)
		}
		// Antisymmetry and associativity of ∧.
		if !closeMV(a.Outer(u, v), a.Outer(v, u).Scale(-1), 1e-15) ||
			!closeMV(a.Outer(a.Outer(u, v), w), a.Outer(u, a.Outer(v, w)), 1e-14) {
			t.Fatal("∧ is not antisymmetric and associative")
		}
		// a⌋(b∧c) = (a·b)c − (a·c)b (Dorst et al. eq. 3.17).
		ab, ac := a.ScalarProduct(u, v), a.ScalarProduct(u, w)
		if !closeMV(a.LeftContract(u, a.Outer(v, w)), Sub(w.Scale(ab), v.Scale(ac)), 1e-13) {
			t.Fatal("a⌋(b∧c) ≠ (a·b)c − (a·c)b")
		}
		// For vectors uv = u⌋v + u∧v.
		if !closeMV(a.Mul(u, v), Add(a.LeftContract(u, v), a.Outer(u, v)), 1e-14) {
			t.Fatal("uv ≠ u·v + u∧v")
		}
		// Wedge of 5 vectors = det·I.
		var m [5][5]float64
		vs := make([]Multivector, 5)
		for i := range vs {
			vs[i] = vec()
			for j := 0; j < 5; j++ {
				m[j][i] = vs[i].Coef(Blade(1) << j)
			}
		}
		got := a.OuterAll(vs...)
		if d := det5(m); math.Abs(got.Coef(a.full)-d) > 1e-12*math.Max(1, math.Abs(d)) || got.Len() > 1 {
			t.Fatalf("v₁∧…∧v₅ = %v, det = %g", got, d)
		}
	}
}

func det5(m [5][5]float64) float64 {
	d := 1.0
	for c := 0; c < 5; c++ {
		p := c
		for r := c + 1; r < 5; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		if p != c {
			m[c], m[p] = m[p], m[c]
			d = -d
		}
		d *= m[c][c]
		for r := c + 1; r < 5; r++ {
			f := m[r][c] / m[c][c]
			for k := c; k < 5; k++ {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	return d
}

// Twin test against a matrix representation: Cl(3,0,0) ≅ 2×2 complex
// matrices with eₖ ↦ Pauli σₖ, and Cl(1,3,0) ≅ 4×4 complex matrices with
// eₖ ↦ Dirac γ⁰…γ³. Random products in the algebra map to matrix products.
func TestMatrixRepresentations(t *testing.T) {
	s1 := [][]complex128{{0, 1}, {1, 0}}
	s2 := [][]complex128{{0, -1i}, {1i, 0}}
	s3 := [][]complex128{{1, 0}, {0, -1}}
	checkRep(t, mustNew(t, 3, 0, 0), [][][]complex128{s1, s2, s3}, 5)
	// Dirac representation: γ⁰ = diag(I, −I), γᵏ = [[0, σₖ], [−σₖ, 0]].
	g := make([][][]complex128, 4)
	g[0] = [][]complex128{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}
	for k, s := range [][][]complex128{s1, s2, s3} {
		m := [][]complex128{make([]complex128, 4), make([]complex128, 4), make([]complex128, 4), make([]complex128, 4)}
		for i := 0; i < 2; i++ {
			for j := 0; j < 2; j++ {
				m[i][2+j] = s[i][j]
				m[2+i][j] = -s[i][j]
			}
		}
		g[k+1] = m
	}
	checkRep(t, mustNew(t, 1, 3, 0), g, 6)
}

func checkRep(t *testing.T, a *Algebra, gen [][][]complex128, seed int64) {
	t.Helper()
	n := len(gen[0])
	id := make([][]complex128, n)
	for i := range id {
		id[i] = make([]complex128, n)
		id[i][i] = 1
	}
	matMul := func(x, y [][]complex128) [][]complex128 {
		z := make([][]complex128, n)
		for i := range z {
			z[i] = make([]complex128, n)
			for k := 0; k < n; k++ {
				for j := 0; j < n; j++ {
					z[i][j] += x[i][k] * y[k][j]
				}
			}
		}
		return z
	}
	blade := func(b Blade) [][]complex128 {
		m := id
		for k := 0; k < a.Dim(); k++ {
			if b&(Blade(1)<<k) != 0 {
				m = matMul(m, gen[k])
			}
		}
		return m
	}
	rep := func(x Multivector) [][]complex128 {
		z := make([][]complex128, n)
		for i := range z {
			z[i] = make([]complex128, n)
		}
		for _, tt := range x.t {
			bm := blade(tt.b)
			for i := range z {
				for j := range z[i] {
					z[i][j] += complex(tt.c, 0) * bm[i][j]
				}
			}
		}
		return z
	}
	r := rand.New(rand.NewSource(seed))
	for k := 0; k < 50; k++ {
		x, y := randMV(r, a), randMV(r, a)
		got, want := rep(a.Mul(x, y)), matMul(rep(x), rep(y))
		for i := range got {
			for j := range got[i] {
				if cmplx.Abs(got[i][j]-want[i][j]) > 1e-12*(1+cmplx.Abs(want[i][j])) {
					t.Fatalf("matrix representation differs at (%d,%d): %v vs %v", i, j, got[i][j], want[i][j])
				}
			}
		}
	}
}

// Quaternions: Cl(0,2,0) ≅ ℍ (i = e₁, j = e₂, k = e₁e₂), and the even
// subalgebra of G3 ≅ ℍ under Quat / G3FromQuat: the conversion round-trips
// exactly and products agree to round-off; rotors from either side rotate
// vectors identically.
func TestQuaternionIsomorphisms(t *testing.T) {
	r := rand.New(rand.NewSource(7))
	h := mustNew(t, 0, 2, 0)
	toH := func(q quat.Q) Multivector {
		m, _ := h.FromBlades(map[Blade]float64{0: q.W, 1: q.X, 2: q.Y, 3: q.Z})
		return m
	}
	for n := 0; n < 1000; n++ {
		p := quat.Q{W: r.NormFloat64(), X: r.NormFloat64(), Y: r.NormFloat64(), Z: r.NormFloat64()}
		q := quat.Q{W: r.NormFloat64(), X: r.NormFloat64(), Y: r.NormFloat64(), Z: r.NormFloat64()}
		pq := p.Mul(q)
		if !closeMV(h.Mul(toH(p), toH(q)), toH(pq), 1e-15) {
			t.Fatalf("Cl(0,2) product ≠ quaternion product")
		}
		if back := G3FromQuat(p).Quat(); back != p {
			t.Fatalf("Quat(G3FromQuat(p)) = %v, want %v exactly", back, p)
		}
		got := G3FromQuat(p).Mul(G3FromQuat(q)).Quat()
		if got.Sub(pq).Abs() > 1e-15*(1+pq.Abs()) {
			t.Fatalf("G3 even product %v ≠ quaternion product %v", got, pq)
		}
		axis := [3]float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		th := r.Float64() * 2 * math.Pi
		rq, _ := quat.Rotor(axis[0], axis[1], axis[2], th)
		rg, _ := Rotor3(axis[0], axis[1], axis[2], th)
		if d := rg.Quat().Sub(rq).Abs(); d > 1e-15 {
			t.Fatalf("Rotor3 and quat.Rotor differ by %g", d)
		}
		v := [3]float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		x1, y1, z1 := rq.Rotate(v[0], v[1], v[2])
		x2, y2, z2 := rg.Sandwich(Vector3(v[0], v[1], v[2])).Vec()
		if math.Abs(x1-x2)+math.Abs(y1-y2)+math.Abs(z1-z2) > 1e-14*(1+math.Abs(x1)+math.Abs(y1)+math.Abs(z1)) {
			t.Fatalf("rotations differ: (%g,%g,%g) vs (%g,%g,%g)", x1, y1, z1, x2, y2, z2)
		}
		// R R̃ = 1.
		rr := rg.Mul(rg.Reverse())
		rr[0]--
		if maxAbs3(rr) > 1e-15 {
			t.Fatalf("R R̃ = %v", rr)
		}
	}
}

func maxAbs3(x G3) float64 {
	m := 0.0
	for _, v := range x {
		m = math.Max(m, math.Abs(v))
	}
	return m
}

// The fixed G2 and G3 agree with the sparse algebra; G2's even part is the
// complex numbers.
func TestFixedTypes(t *testing.T) {
	r := rand.New(rand.NewSource(8))
	a2, a3 := mustNew(t, 2, 0, 0), mustNew(t, 3, 0, 0)
	for n := 0; n < 500; n++ {
		var x, y G3
		for i := range x {
			x[i], y[i] = r.NormFloat64(), r.NormFloat64()
		}
		if !closeMV(x.Mul(y).Multivector(), a3.Mul(x.Multivector(), y.Multivector()), 1e-15) {
			t.Fatal("G3.Mul differs from the sparse product")
		}
		if !closeMV(x.Outer(y).Multivector(), a3.Outer(x.Multivector(), y.Multivector()), 1e-15) {
			t.Fatal("G3.Outer differs from the sparse product")
		}
		if back, ok := G3From(x.Multivector()); !ok || back != x {
			t.Fatal("G3 ↔ Multivector does not round-trip")
		}
		var u, v G2
		for i := range u {
			u[i], v[i] = r.NormFloat64(), r.NormFloat64()
		}
		if !closeMV(u.Mul(v).Multivector(), a2.Mul(u.Multivector(), v.Multivector()), 1e-15) {
			t.Fatal("G2.Mul differs from the sparse product")
		}
		z, w := complex(r.NormFloat64(), r.NormFloat64()), complex(r.NormFloat64(), r.NormFloat64())
		if d := cmplx.Abs(G2FromComplex(z).Mul(G2FromComplex(w)).Complex() - z*w); d > 1e-15*(1+cmplx.Abs(z*w)) {
			t.Fatalf("G2 even product ≠ complex product by %g", d)
		}
	}
	if testing.AllocsPerRun(100, func() {
		var x, y G3
		x[1], y[2] = 1, 2
		_ = x.Mul(y).Sandwich(y)
	}) != 0 {
		t.Error("G3 products allocate")
	}
}

// Inverse, dual and exponential.
func TestInverseDualExp(t *testing.T) {
	r := rand.New(rand.NewSource(9))
	for _, sig := range [][3]int{{3, 0, 0}, {1, 3, 0}, {2, 1, 0}, {4, 0, 0}} {
		a := mustNew(t, sig[0], sig[1], sig[2])
		ok := 0
		for n := 0; n < 30; n++ {
			x := randMV(r, a)
			inv, err := a.Inverse(x)
			if err != nil {
				continue // a random multivector can be singular
			}
			ok++
			if !closeMV(a.Mul(x, inv), Scalar(1), 1e-10) {
				t.Fatalf("Cl%v: x·x⁻¹ ≠ 1", sig)
			}
		}
		if ok < 25 {
			t.Errorf("Cl%v: only %d of 30 random multivectors inverted", sig, ok)
		}
		// A vector with v² ≠ 0 inverts to v/v².
		v, _ := a.Vector(1, 2, 0.5)
		inv, err := a.Inverse(v)
		v2 := a.Mul(v, v).Coef(0)
		if err != nil || !closeMV(inv, v.Scale(1/v2), 1e-15) {
			t.Fatalf("Cl%v: v⁻¹ = %v, want v/v²", sig, inv)
		}
	}
	// A null vector is not invertible; a degenerate algebra has no dual.
	m := mustNew(t, 1, 1, 0)
	null, _ := m.Vector(1, 1)
	if _, err := m.Inverse(null); err == nil {
		t.Error("null vector inverted")
	}
	pga := mustNew(t, 3, 0, 1)
	if _, err := pga.Dual(Scalar(1)); err == nil {
		t.Error("dual in a degenerate algebra")
	}
	// G3 dual: e₁* = e₁·I⁻¹ = −e₂₃ (I⁻¹ = −I).
	a3 := mustNew(t, 3, 0, 0)
	e1, _ := a3.Basis(1)
	d, _ := a3.Dual(e1)
	want, _ := a3.Blade(-1, 2, 3)
	if !closeMV(d, want, 0) {
		t.Errorf("e₁* = %v, want −e₂e₃", d)
	}
	// Exp: e^(e₁₂θ) = cos θ + e₁₂ sin θ; e^(e₁θ) = cosh θ + e₁ sinh θ; the
	// series path (a non-blade argument) agrees with a product of commuting
	// factors e^(A+B) = e^A e^B for A = e₁₂, B = scalar.
	e12, _ := a3.Blade(0.7, 1, 2)
	if got, wantE := a3.Exp(e12), Add(Scalar(math.Cos(0.7)), e12.Scale(math.Sin(0.7)/0.7)); !closeMV(got, wantE, 1e-15) {
		t.Errorf("exp(0.7e₁₂) = %v", got)
	}
	mixed := Add(e12, Scalar(0.3))
	if got, wantE := a3.Exp(mixed), a3.Mul(a3.Exp(Scalar(0.3)), a3.Exp(e12)); !closeMV(got, wantE, 1e-13) {
		t.Errorf("exp(0.3 + 0.7e₁₂) = %v, want %v", got, wantE)
	}
	// Rotor in a general algebra: rotate e₁ by π/2 in the e₁e₂ plane → e₂.
	b, _ := a3.Blade(1, 1, 2)
	rot, err := a3.Rotor(b, math.Pi/2)
	if err != nil {
		t.Fatal(err)
	}
	e2, _ := a3.Basis(2)
	if got := a3.Sandwich(rot, e1); !closeMV(got, e2, 1e-15) {
		t.Errorf("rotating e₁ by π/2 in e₁e₂ gives %v", got)
	}
	// A Lorentz boost in Cl(1,3): e^(e₁e₂·φ/2) sandwich preserves v² for
	// a timelike vector (the spacetime interval).
	st := mustNew(t, 1, 3, 0)
	boostPlane, _ := st.Blade(0.8, 1, 2)
	R := st.Exp(boostPlane)
	tv, _ := st.Vector(2, 0.3, 0.1, -0.4)
	tv2 := st.Mul(tv, tv).Coef(0)
	if got := st.Mul(st.Sandwich(R, tv), st.Sandwich(R, tv)).Coef(0); math.Abs(got-tv2) > 1e-13 {
		t.Errorf("boost changed the interval: %g → %g", tv2, got)
	}
}

// Determinism: the same inputs give bit-identical products.
func TestDeterministic(t *testing.T) {
	a := mustNew(t, 4, 2, 1)
	r := rand.New(rand.NewSource(10))
	x, y := randMV(r, a), randMV(r, a)
	first := a.Mul(x, y)
	for i := 0; i < 20; i++ {
		got := a.Mul(x, y)
		if got.Len() != first.Len() {
			t.Fatal("length differs")
		}
		for k := range got.t {
			if got.t[k] != first.t[k] {
				t.Fatal("products differ between runs")
			}
		}
	}
}

func BenchmarkG3Mul(b *testing.B) {
	x, y := G3{1, 2, 3, 4, 5, 6, 7, 8}, G3{8, 7, 6, 5, 4, 3, 2, 1}
	for i := 0; i < b.N; i++ {
		x = x.Mul(y).Scale(0.01)
	}
	_ = x
}

func BenchmarkSparseMulCl6(b *testing.B) {
	a, _ := New(6, 0, 0)
	r := rand.New(rand.NewSource(1))
	x, y := randMV(r, a), randMV(r, a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_ = a.Mul(x, y)
	}
}
