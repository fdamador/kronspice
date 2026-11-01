// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// stamper is satisfied by every real backend.
type stamper interface {
	Add(i, j int, v float64)
	Clear()
	Solve(b, x []float64) error
}

// gridMNA stamps an MNA system of a w×h resistor grid (conductances from
// rng) with nv grounded voltage sources on the first nodes, each adding a
// branch row with a zero diagonal (so pivoting is required). It returns the
// size and the right-hand side.
func gridMNA(m stamper, w, h, nv int, rng *rand.Rand) []float64 {
	nn := w * h
	n := nn + nv
	m.Clear()
	g := func(a, b int, v float64) {
		m.Add(a, a, v)
		m.Add(b, b, v)
		m.Add(a, b, -v)
		m.Add(b, a, -v)
	}
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			i := y*w + x
			if x+1 < w {
				g(i, i+1, 1/(10+990*rng.Float64()))
			}
			if y+1 < h {
				g(i, i+w, 1/(10+990*rng.Float64()))
			}
			m.Add(i, i, 1e-9) // leak, so no node floats
		}
	}
	b := make([]float64, n)
	for k := 0; k < nv; k++ {
		node, br := (k*7919)%nn, nn+k
		m.Add(node, br, 1)
		m.Add(br, node, 1)
		b[br] = 1 + rng.Float64()
	}
	return b
}

func backwardError(d *Dense, x, b []float64) float64 {
	n := len(b)
	an, xn, bn, rn := 0.0, 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		r := b[i]
		for j := 0; j < n; j++ {
			r -= d.At(i, j) * x[j]
			an = math.Max(an, math.Abs(d.At(i, j)))
		}
		rn = math.Max(rn, math.Abs(r))
		xn = math.Max(xn, math.Abs(x[i]))
		bn = math.Max(bn, math.Abs(b[i]))
	}
	return rn / (an*xn*float64(n) + bn)
}

func relDiff(a, b []float64) float64 {
	worst, scale := 0.0, 0.0
	for i := range a {
		worst = math.Max(worst, math.Abs(a[i]-b[i]))
		scale = math.Max(scale, math.Abs(b[i]))
	}
	return worst / scale
}

// Sparse and dense backends agree on MNA-shaped systems of many sizes.
func TestSparseMatchesDense(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, c := range []struct{ w, h, nv int }{{1, 1, 1}, {2, 2, 1}, {5, 1, 2}, {6, 5, 3}, {12, 10, 6}, {25, 20, 10}} {
		n := c.w*c.h + c.nv
		d, _ := NewDense(n)
		s, _ := NewSparse(n)
		seed := rng.Int63()
		b := gridMNA(d, c.w, c.h, c.nv, rand.New(rand.NewSource(seed)))
		gridMNA(s, c.w, c.h, c.nv, rand.New(rand.NewSource(seed)))
		xd, xs := make([]float64, n), make([]float64, n)
		if err := d.Solve(b, xd); err != nil {
			t.Fatal(err)
		}
		if err := s.Solve(b, xs); err != nil {
			t.Fatal(err)
		}
		if r := relDiff(xs, xd); r > 1e-9 {
			t.Errorf("%dx%d+%d: sparse vs dense %g", c.w, c.h, c.nv, r)
		}
	}
}

// Random unsymmetric sparse matrices with a dominant diagonal on only some
// rows: checks the residual, not just agreement.
func TestSparseRandomResidual(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	for trial := 0; trial < 50; trial++ {
		n := 1 + rng.Intn(60)
		s, _ := NewSparse(n)
		d, _ := NewDense(n)
		for i := 0; i < n; i++ {
			// A random permutation matrix plus noise is nonsingular with
			// high probability and needs off-diagonal pivots.
			j := rng.Intn(n)
			for _, m := range []stamper{s, d} {
				m.Add(i, (i+j)%n, 1)
			}
		}
		for k := 0; k < 3*n; k++ {
			i, j, v := rng.Intn(n), rng.Intn(n), rng.NormFloat64()*0.1
			s.Add(i, j, v)
			d.Add(i, j, v)
		}
		b := make([]float64, n)
		for i := range b {
			b[i] = rng.NormFloat64()
		}
		xd, xs := make([]float64, n), make([]float64, n)
		errD, errS := d.Solve(b, xd), s.Solve(b, xs)
		if (errD == nil) != (errS == nil) {
			t.Fatalf("trial %d: dense %v, sparse %v", trial, errD, errS)
		}
		if errD != nil {
			continue
		}
		// Normwise backward error |b − A·x| / (|A|·|x| + |b|), the
		// quantity partial pivoting keeps near machine precision; compared
		// with the dense backend's on the same system.
		if es, ed := backwardError(d, xs, b), backwardError(d, xd, b); es > 1e-13 {
			t.Fatalf("trial %d (n=%d): backward error sparse %g, dense %g", trial, n, es, ed)
		}
	}
}

// A voltage source across a single node: zero diagonal on the branch row.
func TestSparseZeroDiagonal(t *testing.T) {
	s, _ := NewSparse(2)
	s.Add(0, 0, 1e-3)
	s.Add(0, 1, 1)
	s.Add(1, 0, 1)
	x := make([]float64, 2)
	if err := s.Solve([]float64{0, 5}, x); err != nil {
		t.Fatal(err)
	}
	if math.Abs(x[0]-5) > 1e-12 || math.Abs(x[1]+5e-3) > 1e-15 {
		t.Fatalf("x = %v", x)
	}
}

func TestSparseSingular(t *testing.T) {
	cases := map[string]func(*Sparse){
		"empty":       func(*Sparse) {},
		"zero column": func(s *Sparse) { s.Add(0, 0, 1); s.Add(1, 0, 1); s.Add(2, 2, 1) },
		"dependent rows": func(s *Sparse) {
			for j := 0; j < 3; j++ {
				s.Add(0, j, 1)
				s.Add(1, j, 2)
			}
			s.Add(2, 2, 1)
		},
	}
	for name, f := range cases {
		s, _ := NewSparse(3)
		f(s)
		if err := s.Solve(make([]float64, 3), make([]float64, 3)); !errors.Is(err, ErrSingularMatrix) {
			t.Errorf("%s: %v", name, err)
		}
	}
	s, _ := NewSparse(2)
	s.Add(0, 0, math.NaN())
	s.Add(1, 1, 1)
	if err := s.Solve(make([]float64, 2), make([]float64, 2)); !errors.Is(err, ErrNonFinite) {
		t.Errorf("NaN: %v", err)
	}
	if _, err := NewSparse(0); !errors.Is(err, ErrBadSize) {
		t.Errorf("size 0: %v", err)
	}
	if err := s.Solve(make([]float64, 1), make([]float64, 2)); !errors.Is(err, ErrBadSize) {
		t.Errorf("bad rhs: %v", err)
	}
}

// Ground indices are ignored, as in Dense.
func TestSparseIgnoresGround(t *testing.T) {
	s, _ := NewSparse(1)
	s.Add(-1, 0, 5)
	s.Add(0, -1, 5)
	s.Add(0, 0, 2)
	x := []float64{0}
	if err := s.Solve([]float64{4}, x); err != nil || x[0] != 2 {
		t.Fatalf("x %v err %v", x, err)
	}
}

// After the pattern is known, re-stamping and solving allocate nothing.
func TestSparseSolveDoesNotAllocate(t *testing.T) {
	s, _ := NewSparse(12*10 + 6)
	rng := rand.New(rand.NewSource(3))
	b := gridMNA(s, 12, 10, 6, rng)
	x := make([]float64, len(b))
	if err := s.Solve(b, x); err != nil {
		t.Fatal(err)
	}
	allocs := testing.AllocsPerRun(20, func() {
		gridMNA2(s, 12, 10, 6)
		if err := s.Solve(b, x); err != nil {
			t.Fatal(err)
		}
	})
	if allocs != 0 {
		t.Fatalf("%v allocations per re-stamp + solve", allocs)
	}
}

// gridMNA2 restamps the same pattern as gridMNA without allocating.
func gridMNA2(m stamper, w, h, nv int) {
	nn := w * h
	m.Clear()
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			i := y*w + x
			for _, j := range [2]int{i + 1, i + w} {
				if (j == i+1 && x+1 < w) || (j == i+w && y+1 < h) {
					m.Add(i, i, 1e-2)
					m.Add(j, j, 1e-2)
					m.Add(i, j, -1e-2)
					m.Add(j, i, -1e-2)
				}
			}
			m.Add(i, i, 1e-9)
		}
	}
	for k := 0; k < nv; k++ {
		node, br := (k*7919)%nn, nn+k
		m.Add(node, br, 1)
		m.Add(br, node, 1)
	}
}

// A new entry after the first solve (a device that stamps only in some
// regions) grows the pattern, and the next solve is still right.
func TestSparsePatternGrowth(t *testing.T) {
	s, _ := NewSparse(2)
	d, _ := NewDense(2)
	for _, m := range []stamper{s, d} {
		m.Add(0, 0, 2)
		m.Add(1, 1, 3)
	}
	x := make([]float64, 2)
	if err := s.Solve([]float64{2, 3}, x); err != nil || x[0] != 1 || x[1] != 1 {
		t.Fatalf("x %v err %v", x, err)
	}
	for _, m := range []stamper{s, d} {
		m.Add(0, 1, 1)
		m.Add(1, 0, -1)
	}
	xd := make([]float64, 2)
	_ = d.Solve([]float64{2, 3}, xd)
	if err := s.Solve([]float64{2, 3}, x); err != nil || relDiff(x, xd) > 1e-14 {
		t.Fatalf("x %v want %v err %v", x, xd, err)
	}
}

// Fill stays near-linear for a grid under minimum degree, far below n².
func TestSparseFillIsModest(t *testing.T) {
	s, _ := NewSparse(40*40 + 4)
	b := gridMNA(s, 40, 40, 4, rand.New(rand.NewSource(5)))
	if err := s.Solve(b, make([]float64, len(b))); err != nil {
		t.Fatal(err)
	}
	a, l, u := s.NNZ()
	n := len(b)
	if l+u > 30*a || l+u > n*n/20 {
		t.Fatalf("fill: A %d, L %d, U %d (n %d)", a, l, u, n)
	}
	t.Logf("n %d: nnz A %d, L %d, U %d", n, a, l, u)
}

func TestCSparseMatchesCDense(t *testing.T) {
	rng := rand.New(rand.NewSource(11))
	for trial := 0; trial < 30; trial++ {
		n := 1 + rng.Intn(40)
		s, _ := NewCSparse(n)
		d, _ := NewCDense(n)
		add := func(i, j int, v complex128) { s.Add(i, j, v); d.Add(i, j, v) }
		for i := 0; i < n; i++ {
			add(i, (i+rng.Intn(n))%n, complex(1+rng.Float64(), rng.Float64()))
		}
		for k := 0; k < 2*n; k++ {
			add(rng.Intn(n), rng.Intn(n), complex(rng.NormFloat64(), rng.NormFloat64())*0.1)
		}
		b := make([]complex128, n)
		for i := range b {
			b[i] = complex(rng.NormFloat64(), rng.NormFloat64())
		}
		xs, xd := make([]complex128, n), make([]complex128, n)
		errS, errD := s.Solve(b, xs), d.Solve(b, xd)
		if (errS == nil) != (errD == nil) {
			t.Fatalf("trial %d: sparse %v dense %v", trial, errS, errD)
		}
		if errS != nil {
			continue
		}
		worst, scale := 0.0, 0.0
		for i := range xs {
			worst = math.Max(worst, cmplx.Abs(xs[i]-xd[i]))
			scale = math.Max(scale, cmplx.Abs(xd[i]))
		}
		if worst > 1e-9*scale {
			t.Fatalf("trial %d: diff %g", trial, worst/scale)
		}
	}
}

func FuzzSparseMatchesDense(f *testing.F) {
	f.Add(int64(1), uint8(5))
	f.Add(int64(42), uint8(30))
	f.Fuzz(func(t *testing.T, seed int64, size uint8) {
		n := 1 + int(size)%40
		rng := rand.New(rand.NewSource(seed))
		s, _ := NewSparse(n)
		d, _ := NewDense(n)
		for k := 0; k < 3*n; k++ {
			i, j, v := rng.Intn(n), rng.Intn(n), rng.NormFloat64()
			if rng.Intn(10) == 0 {
				v = 0
			}
			s.Add(i, j, v)
			d.Add(i, j, v)
		}
		b := make([]float64, n)
		for i := range b {
			b[i] = rng.NormFloat64()
		}
		xs, xd := make([]float64, n), make([]float64, n)
		errS, errD := s.Solve(b, xs), d.Solve(b, xd)
		if errS != nil {
			if !errors.Is(errS, ErrSingularMatrix) && !errors.Is(errS, ErrNonFinite) {
				t.Fatalf("unexpected error %v", errS)
			}
			return
		}
		// The sparse solve succeeded: its residual must be small even if
		// dense judged the matrix singular (different pivot order).
		res := make([]float64, n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				res[i] += d.At(i, j) * xs[j]
			}
		}
		_ = errD
		bn := 0.0
		for _, v := range b {
			bn = math.Max(bn, math.Abs(v))
		}
		an, xn := 0.0, 0.0
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				an = math.Max(an, math.Abs(d.At(i, j)))
			}
			xn = math.Max(xn, math.Abs(xs[i]))
		}
		for i := range res {
			if math.Abs(res[i]-b[i]) > 1e-8*(bn+an*xn*float64(n)) {
				t.Fatalf("residual %g at row %d (n %d)", res[i]-b[i], i, n)
			}
		}
	})
}

func benchGrid(b *testing.B, mk func(int) stamper, w, h int) {
	n := w*h + 2
	m := mk(n)
	rhs := gridMNA(m, w, h, 2, rand.New(rand.NewSource(1)))
	x := make([]float64, n)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		gridMNA2(m, w, h, 2)
		if err := m.Solve(rhs, x); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkGridSolve(b *testing.B) {
	dense := func(n int) stamper { d, _ := NewDense(n); return d }
	sparse := func(n int) stamper { s, _ := NewSparse(n); return s }
	for _, wh := range [][2]int{{5, 5}, {8, 8}, {16, 16}, {32, 32}} {
		n := wh[0]*wh[1] + 2
		b.Run(fmt.Sprintf("dense/n=%d", n), func(b *testing.B) { benchGrid(b, dense, wh[0], wh[1]) })
		b.Run(fmt.Sprintf("sparse/n=%d", n), func(b *testing.B) { benchGrid(b, sparse, wh[0], wh[1]) })
	}
}
