// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package ga

import (
	"errors"
	"fmt"
	"math"
	"math/bits"
	"sort"
	"strings"
)

// Errors returned by the algebra.
var (
	ErrSignature     = errors.New("ga: invalid signature")
	ErrBasis         = errors.New("ga: basis vector out of range")
	ErrNotInvertible = errors.New("ga: multivector is not invertible")
	ErrDegenerate    = errors.New("ga: pseudoscalar is not invertible (degenerate metric)")
)

// Blade is a basis blade as a bitmask (bit k−1 ↔ eₖ), see the package
// documentation.
type Blade uint64

// Grade is the number of basis vectors in the blade.
func (b Blade) Grade() int { return bits.OnesCount64(uint64(b)) }

// Algebra is Cl(p, q, r).
type Algebra struct {
	p, q, r int
	sq      []float64 // eₖ²
	full    Blade     // the pseudoscalar
}

// New returns Cl(p, q, r); n = p + q + r must be 1 … 64.
func New(p, q, r int) (*Algebra, error) {
	n := p + q + r
	if p < 0 || q < 0 || r < 0 || n < 1 || n > 64 {
		return nil, fmt.Errorf("%w: Cl(%d,%d,%d)", ErrSignature, p, q, r)
	}
	a := &Algebra{p: p, q: q, r: r, sq: make([]float64, n)}
	for k := range a.sq {
		switch {
		case k < p:
			a.sq[k] = 1
		case k < p+q:
			a.sq[k] = -1
		}
	}
	if n == 64 {
		a.full = ^Blade(0)
	} else {
		a.full = Blade(1)<<n - 1
	}
	return a, nil
}

// Dim is the number of basis vectors n.
func (a *Algebra) Dim() int { return len(a.sq) }

// Signature returns (p, q, r).
func (a *Algebra) Signature() (int, int, int) { return a.p, a.q, a.r }

// Pseudoscalar is the blade e₁e₂…eₙ.
func (a *Algebra) Pseudoscalar() Blade { return a.full }

// term is one coefficient on a canonical blade.
type term struct {
	b Blade
	c float64
}

// Multivector is a sparse element of an algebra: non-zero coefficients on
// canonical blades, sorted by blade. The zero value is 0. Multivectors are
// values: operations never modify their arguments.
type Multivector struct{ t []term }

// Scalar returns the multivector s.
func Scalar(s float64) Multivector {
	if s == 0 {
		return Multivector{}
	}
	return Multivector{[]term{{0, s}}}
}

// Basis returns the vector eₖ (k = 1 … n).
func (a *Algebra) Basis(k int) (Multivector, error) {
	if k < 1 || k > a.Dim() {
		return Multivector{}, fmt.Errorf("%w: e%d in Cl(%d,%d,%d)", ErrBasis, k, a.p, a.q, a.r)
	}
	return Multivector{[]term{{Blade(1) << (k - 1), 1}}}, nil
}

// Vector returns Σ cₖ·eₖ for c = (c₁, …, cₘ), m ≤ n.
func (a *Algebra) Vector(c ...float64) (Multivector, error) {
	if len(c) > a.Dim() {
		return Multivector{}, fmt.Errorf("%w: %d components in Cl(%d,%d,%d)", ErrBasis, len(c), a.p, a.q, a.r)
	}
	var m Multivector
	for k, v := range c {
		if v != 0 {
			m.t = append(m.t, term{Blade(1) << k, v})
		}
	}
	return m, nil
}

// Blade returns c·eᵢ eⱼ … (factors in the given order, any order, repeated
// factors contracted by the metric).
func (a *Algebra) Blade(c float64, factors ...int) (Multivector, error) {
	m := Scalar(c)
	for _, k := range factors {
		e, err := a.Basis(k)
		if err != nil {
			return Multivector{}, err
		}
		m = a.Mul(m, e)
	}
	return m, nil
}

// FromBlades builds a multivector from blade → coefficient pairs; blades
// outside the algebra are an error.
func (a *Algebra) FromBlades(c map[Blade]float64) (Multivector, error) {
	var m Multivector
	for b, v := range c {
		if b&^a.full != 0 {
			return Multivector{}, fmt.Errorf("%w: blade %b", ErrBasis, b)
		}
		if v != 0 {
			m.t = append(m.t, term{b, v})
		}
	}
	sort.Slice(m.t, func(i, j int) bool { return m.t[i].b < m.t[j].b })
	return m, nil
}

// Coef returns the coefficient on blade b.
func (m Multivector) Coef(b Blade) float64 {
	i := sort.Search(len(m.t), func(i int) bool { return m.t[i].b >= b })
	if i < len(m.t) && m.t[i].b == b {
		return m.t[i].c
	}
	return 0
}

// Blades returns the blades with non-zero coefficients, in increasing order.
func (m Multivector) Blades() []Blade {
	out := make([]Blade, len(m.t))
	for i, t := range m.t {
		out[i] = t.b
	}
	return out
}

// Len is the number of non-zero coefficients.
func (m Multivector) Len() int { return len(m.t) }

// IsZero reports whether every coefficient is zero.
func (m Multivector) IsZero() bool { return len(m.t) == 0 }

// String writes the multivector as "1.5 + 2·e1 − e1e3".
func (m Multivector) String() string {
	if len(m.t) == 0 {
		return "0"
	}
	var sb strings.Builder
	for i, t := range m.t {
		c := t.c
		if i > 0 {
			if c < 0 {
				sb.WriteString(" − ")
				c = -c
			} else {
				sb.WriteString(" + ")
			}
		}
		sb.WriteString(fmt.Sprintf("%g", c))
		if t.b != 0 {
			sb.WriteString("·")
			for k := 0; k < 64; k++ {
				if t.b&(Blade(1)<<k) != 0 {
					sb.WriteString(fmt.Sprintf("e%d", k+1))
				}
			}
		}
	}
	return sb.String()
}

// reorderSign is (−1)^(transpositions needed to sort eA·eB into canonical
// order): for each factor of A, the factors of B with a lower index
// (Dorst et al. 2007, §19.1.2).
func reorderSign(x, y Blade) float64 {
	x >>= 1
	n := 0
	for x != 0 {
		n += bits.OnesCount64(uint64(x & y))
		x >>= 1
	}
	if n&1 != 0 {
		return -1
	}
	return 1
}

// bladeMul returns the sign (including the metric) and blade of eA·eB.
func (a *Algebra) bladeMul(x, y Blade) (float64, Blade) {
	s := reorderSign(x, y)
	common := x & y
	for common != 0 {
		k := bits.TrailingZeros64(uint64(common))
		if k >= len(a.sq) {
			return math.NaN(), x ^ y // a blade from another algebra
		}
		s *= a.sq[k]
		if s == 0 {
			return 0, 0
		}
		common &= common - 1
	}
	return s, x ^ y
}

// accumulate collects products into a sorted multivector, summing in a
// deterministic order (inputs are sorted, so each blade's sum is formed in
// the same sequence every run).
type accumulator struct {
	idx map[Blade]int
	t   []term
}

func (ac *accumulator) add(b Blade, c float64) {
	if ac.idx == nil {
		ac.idx = make(map[Blade]int)
	}
	if i, ok := ac.idx[b]; ok {
		ac.t[i].c += c
		return
	}
	ac.idx[b] = len(ac.t)
	ac.t = append(ac.t, term{b, c})
}

func (ac *accumulator) result() Multivector {
	out := ac.t[:0]
	for _, t := range ac.t {
		if t.c != 0 {
			out = append(out, t)
		}
	}
	sort.Slice(out, func(i, j int) bool { return out[i].b < out[j].b })
	return Multivector{out}
}

// product is the common kernel: keep(x, y) selects which blade pairs
// contribute.
func (a *Algebra) product(x, y Multivector, keep func(bx, by Blade) bool) Multivector {
	var ac accumulator
	for _, s := range x.t {
		for _, t := range y.t {
			if !keep(s.b, t.b) {
				continue
			}
			sg, b := a.bladeMul(s.b, t.b)
			if sg != 0 {
				ac.add(b, sg*s.c*t.c)
			}
		}
	}
	return ac.result()
}

// Mul is the geometric product xy.
func (a *Algebra) Mul(x, y Multivector) Multivector {
	return a.product(x, y, func(Blade, Blade) bool { return true })
}

// Outer is the outer (wedge) product x∧y.
func (a *Algebra) Outer(x, y Multivector) Multivector {
	return a.product(x, y, func(bx, by Blade) bool { return bx&by == 0 })
}

// LeftContract is x⌋y: blades of x contained in blades of y.
func (a *Algebra) LeftContract(x, y Multivector) Multivector {
	return a.product(x, y, func(bx, by Blade) bool { return bx&by == bx })
}

// RightContract is x⌊y: blades of y contained in blades of x.
func (a *Algebra) RightContract(x, y Multivector) Multivector {
	return a.product(x, y, func(bx, by Blade) bool { return bx&by == by })
}

// ScalarProduct is ⟨xy⟩₀.
func (a *Algebra) ScalarProduct(x, y Multivector) float64 {
	return a.product(x, y, func(bx, by Blade) bool { return bx == by }).Coef(0)
}

// Add returns x + y.
func Add(x, y Multivector) Multivector { return combine(x, y, 1) }

// Sub returns x − y.
func Sub(x, y Multivector) Multivector { return combine(x, y, -1) }

func combine(x, y Multivector, s float64) Multivector {
	out := make([]term, 0, len(x.t)+len(y.t))
	i, j := 0, 0
	for i < len(x.t) || j < len(y.t) {
		switch {
		case j == len(y.t) || (i < len(x.t) && x.t[i].b < y.t[j].b):
			out = append(out, x.t[i])
			i++
		case i == len(x.t) || y.t[j].b < x.t[i].b:
			out = append(out, term{y.t[j].b, s * y.t[j].c})
			j++
		default:
			if c := x.t[i].c + s*y.t[j].c; c != 0 {
				out = append(out, term{x.t[i].b, c})
			}
			i++
			j++
		}
	}
	return Multivector{out}
}

// Scale returns s·x.
func (x Multivector) Scale(s float64) Multivector {
	if s == 0 {
		return Multivector{}
	}
	out := make([]term, len(x.t))
	for i, t := range x.t {
		out[i] = term{t.b, s * t.c}
	}
	return Multivector{out}
}

// mapGrades applies a sign that depends on the grade.
func (x Multivector) mapGrades(f func(k int) float64) Multivector {
	out := make([]term, len(x.t))
	for i, t := range x.t {
		out[i] = term{t.b, f(t.b.Grade()) * t.c}
	}
	return Multivector{out}
}

// Reverse is x̃: (−1)^(k(k−1)/2) on grade k.
func (x Multivector) Reverse() Multivector {
	return x.mapGrades(func(k int) float64 {
		if (k*(k-1)/2)&1 != 0 {
			return -1
		}
		return 1
	})
}

// GradeInvolution is x̂: (−1)^k on grade k.
func (x Multivector) GradeInvolution() Multivector {
	return x.mapGrades(func(k int) float64 {
		if k&1 != 0 {
			return -1
		}
		return 1
	})
}

// Conjugate is the Clifford conjugate x̄: (−1)^(k(k+1)/2) on grade k.
func (x Multivector) Conjugate() Multivector {
	return x.mapGrades(func(k int) float64 {
		if (k*(k+1)/2)&1 != 0 {
			return -1
		}
		return 1
	})
}

// Grade returns the grade-k part ⟨x⟩ₖ.
func (x Multivector) Grade(k int) Multivector {
	var out []term
	for _, t := range x.t {
		if t.b.Grade() == k {
			out = append(out, t)
		}
	}
	return Multivector{out}
}

// Even returns the even-grade part.
func (x Multivector) Even() Multivector {
	var out []term
	for _, t := range x.t {
		if t.b.Grade()&1 == 0 {
			out = append(out, t)
		}
	}
	return Multivector{out}
}

// Grades returns the sorted list of grades present.
func (x Multivector) Grades() []int {
	seen := map[int]bool{}
	var g []int
	for _, t := range x.t {
		if k := t.b.Grade(); !seen[k] {
			seen[k] = true
			g = append(g, k)
		}
	}
	sort.Ints(g)
	return g
}

// Norm2 is ⟨x·x̃⟩₀.
func (a *Algebra) Norm2(x Multivector) float64 { return a.ScalarProduct(x, x.Reverse()) }

// Norm is √|⟨x·x̃⟩₀|.
func (a *Algebra) Norm(x Multivector) float64 { return math.Sqrt(math.Abs(a.Norm2(x))) }

// MaxAbs is the largest coefficient magnitude (a norm that is positive for
// every non-zero multivector, whatever the signature).
func (x Multivector) MaxAbs() float64 {
	m := 0.0
	for _, t := range x.t {
		m = math.Max(m, math.Abs(t.c))
	}
	return m
}

// Inverse returns x⁻¹ with x·x⁻¹ = x⁻¹·x = 1. When x·x̃ is a non-zero
// scalar (every versor and blade) it is x̃/(x·x̃); otherwise, for algebras
// of up to 10 basis vectors, the linear system L_x·y = 1 of the left
// multiplication map is solved (dimension 2ⁿ). ErrNotInvertible otherwise.
func (a *Algebra) Inverse(x Multivector) (Multivector, error) {
	if x.IsZero() {
		return Multivector{}, ErrNotInvertible
	}
	xr := x.Reverse()
	p := a.Mul(x, xr)
	scale := x.MaxAbs() * x.MaxAbs()
	if s := p.Coef(0); s != 0 && Sub(p, Scalar(s)).MaxAbs() <= 1e-13*math.Max(scale, math.Abs(s)) {
		inv := xr.Scale(1 / s)
		if a.isInverse(x, inv) {
			return inv, nil
		}
	}
	n := a.Dim()
	if n > 10 {
		return Multivector{}, fmt.Errorf("%w: general inverse limited to 10 basis vectors", ErrNotInvertible)
	}
	dim := 1 << n
	m := make([][]float64, dim)
	for i := range m {
		m[i] = make([]float64, dim+1)
	}
	for col := 0; col < dim; col++ {
		for _, t := range x.t {
			s, b := a.bladeMul(t.b, Blade(col))
			m[b][col] += s * t.c
		}
	}
	m[0][dim] = 1
	y, ok := solveDense(m)
	if !ok {
		return Multivector{}, ErrNotInvertible
	}
	c := make(map[Blade]float64)
	for i, v := range y {
		if v != 0 {
			c[Blade(i)] = v
		}
	}
	inv, _ := a.FromBlades(c)
	if !a.isInverse(x, inv) {
		return Multivector{}, ErrNotInvertible
	}
	return inv, nil
}

func (a *Algebra) isInverse(x, inv Multivector) bool {
	tol := 1e-10
	return Sub(a.Mul(x, inv), Scalar(1)).MaxAbs() <= tol && Sub(a.Mul(inv, x), Scalar(1)).MaxAbs() <= tol
}

// solveDense solves an augmented system by Gaussian elimination with
// partial pivoting; false when singular to working precision.
func solveDense(m [][]float64) ([]float64, bool) {
	n := len(m)
	big := 0.0
	for _, r := range m {
		for _, v := range r[:n] {
			big = math.Max(big, math.Abs(v))
		}
	}
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		if math.Abs(m[p][c]) <= 1e-13*big {
			return nil, false
		}
		m[c], m[p] = m[p], m[c]
		for r := c + 1; r < n; r++ {
			f := m[r][c] / m[c][c]
			if f == 0 {
				continue
			}
			for k := c; k <= n; k++ {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	x := make([]float64, n)
	for r := n - 1; r >= 0; r-- {
		s := m[r][n]
		for k := r + 1; k < n; k++ {
			s -= m[r][k] * x[k]
		}
		x[r] = s / m[r][r]
	}
	return x, true
}

// Dual returns x·I⁻¹; ErrDegenerate when the metric has null vectors.
func (a *Algebra) Dual(x Multivector) (Multivector, error) {
	if a.r > 0 {
		return Multivector{}, ErrDegenerate
	}
	i := Multivector{[]term{{a.full, 1}}}
	ii, err := a.Inverse(i)
	if err != nil {
		return Multivector{}, err
	}
	return a.Mul(x, ii), nil
}

// Exp returns e^x. When x² is a scalar (every blade, and in particular a
// simple bivector) the closed form cos/cosh is used; otherwise the Taylor
// series with scaling and squaring.
func (a *Algebra) Exp(x Multivector) Multivector {
	if x.IsZero() {
		return Scalar(1)
	}
	x2 := a.Mul(x, x)
	scale := x.MaxAbs() * x.MaxAbs()
	if s := x2.Coef(0); Sub(x2, Scalar(s)).MaxAbs() <= 1e-14*math.Max(scale, 1e-300) {
		switch {
		case s < 0:
			th := math.Sqrt(-s)
			return Add(Scalar(math.Cos(th)), x.Scale(math.Sin(th)/th))
		case s > 0:
			th := math.Sqrt(s)
			return Add(Scalar(math.Cosh(th)), x.Scale(math.Sinh(th)/th))
		default:
			return Add(Scalar(1), x)
		}
	}
	// Scaling and squaring: e^x = (e^(x/2^k))^(2^k), with ‖x/2^k‖ ≤ 1/2.
	k := 0
	for m := x.MaxAbs() * float64(len(x.t)); m > 0.5 && k < 60; m /= 2 {
		k++
	}
	y := x.Scale(math.Ldexp(1, -k))
	sum, term := Scalar(1), Scalar(1)
	for n := 1; n <= 30; n++ {
		term = a.Mul(term, y).Scale(1 / float64(n))
		sum = Add(sum, term)
		if term.MaxAbs() <= 1e-17*sum.MaxAbs() {
			break
		}
	}
	for ; k > 0; k-- {
		sum = a.Mul(sum, sum)
	}
	return sum
}

// Sandwich returns r·x·r̃ (a rotor or versor applied to x).
func (a *Algebra) Sandwich(r, x Multivector) Multivector {
	return a.Mul(a.Mul(r, x), r.Reverse())
}

// Rotor returns e^(−B·θ/2) for the plane B (normalized to |B²| = 1) and
// angle θ; it rotates vectors in that plane by θ (the sense of B) through
// Sandwich. B must be a 2-blade squaring to a negative scalar.
func (a *Algebra) Rotor(b Multivector, theta float64) (Multivector, error) {
	b2 := a.Mul(b, b)
	s := b2.Coef(0)
	if !(s < 0) || Sub(b2, Scalar(s)).MaxAbs() > 1e-12*math.Abs(s) {
		return Multivector{}, fmt.Errorf("ga: rotor plane must be a blade with B² < 0, got B² = %v", b2)
	}
	return a.Exp(b.Scale(-theta / (2 * math.Sqrt(-s)))), nil
}

// OuterAll returns v₁∧v₂∧…∧vₖ (1 for no vectors).
func (a *Algebra) OuterAll(v ...Multivector) Multivector {
	out := Scalar(1)
	for _, x := range v {
		out = a.Outer(out, x)
	}
	return out
}
