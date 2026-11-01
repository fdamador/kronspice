// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"errors"
	"math"
	"math/rand"
	"testing"
)

func TestDenseKnownSystem(t *testing.T) {
	d, _ := NewDense(3)
	A := [][]float64{{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}}
	for i := range A {
		for j := range A[i] {
			d.Add(i, j, A[i][j])
		}
	}
	x := make([]float64, 3)
	if err := d.Solve([]float64{8, -11, -3}, x); err != nil {
		t.Fatal(err)
	}
	want := []float64{2, 3, -1}
	for i := range x {
		if math.Abs(x[i]-want[i]) > 1e-12 {
			t.Fatalf("x = %v, want %v", x, want)
		}
	}
}

func TestDenseNeedsPivoting(t *testing.T) {
	// Zero on the diagonal: fails without row exchanges.
	d, _ := NewDense(2)
	d.Add(0, 1, 1)
	d.Add(1, 0, 1)
	x := make([]float64, 2)
	if err := d.Solve([]float64{3, 4}, x); err != nil {
		t.Fatal(err)
	}
	if x[0] != 4 || x[1] != 3 {
		t.Fatalf("x = %v", x)
	}
}

func TestDenseSingular(t *testing.T) {
	d, _ := NewDense(2)
	d.Add(0, 0, 1)
	d.Add(0, 1, 2)
	d.Add(1, 0, 2)
	d.Add(1, 1, 4)
	if err := d.Solve([]float64{1, 2}, make([]float64, 2)); !errors.Is(err, ErrSingularMatrix) {
		t.Fatalf("err = %v, want ErrSingularMatrix", err)
	}
	z, _ := NewDense(2)
	if err := z.Solve([]float64{1, 2}, make([]float64, 2)); !errors.Is(err, ErrSingularMatrix) {
		t.Fatalf("zero matrix err = %v", err)
	}
}

func TestDenseBadInput(t *testing.T) {
	if _, err := NewDense(0); !errors.Is(err, ErrBadSize) {
		t.Fatal("size 0 accepted")
	}
	d, _ := NewDense(1)
	d.Add(0, 0, math.NaN())
	if err := d.Solve([]float64{1}, make([]float64, 1)); !errors.Is(err, ErrNonFinite) {
		t.Fatalf("NaN err = %v", err)
	}
	if err := d.Solve([]float64{1, 2}, make([]float64, 1)); !errors.Is(err, ErrBadSize) {
		t.Fatalf("size mismatch err = %v", err)
	}
}

func TestDenseRandomResidual(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for n := 1; n <= 40; n += 7 {
		d, _ := NewDense(n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				d.Add(i, j, rng.NormFloat64())
			}
			d.Add(i, i, float64(n)) // keep it well-conditioned
		}
		b := make([]float64, n)
		for i := range b {
			b[i] = rng.NormFloat64()
		}
		x := make([]float64, n)
		if err := d.Solve(b, x); err != nil {
			t.Fatal(err)
		}
		for i := 0; i < n; i++ {
			s := 0.0
			for j := 0; j < n; j++ {
				s += d.At(i, j) * x[j]
			}
			if math.Abs(s-b[i]) > 1e-10 {
				t.Fatalf("n=%d residual %g", n, s-b[i])
			}
		}
	}
}

func TestDenseSolveDoesNotAllocate(t *testing.T) {
	d, _ := NewDense(20)
	for i := 0; i < 20; i++ {
		d.Add(i, i, 2)
		if i > 0 {
			d.Add(i, i-1, -1)
			d.Add(i-1, i, -1)
		}
	}
	b, x := make([]float64, 20), make([]float64, 20)
	b[0] = 1
	allocs := testing.AllocsPerRun(100, func() {
		d.Clear()
		for i := 0; i < 20; i++ {
			d.Add(i, i, 2)
		}
		_ = d.Solve(b, x)
	})
	if allocs != 0 {
		t.Fatalf("Solve allocates %v times per run", allocs)
	}
}

func BenchmarkDenseSolve30(b *testing.B) {
	n := 30
	d, _ := NewDense(n)
	for i := 0; i < n; i++ {
		d.Add(i, i, 4)
		if i > 0 {
			d.Add(i, i-1, -1)
			d.Add(i-1, i, -1)
		}
	}
	rhs, x := make([]float64, n), make([]float64, n)
	rhs[0] = 1
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = d.Solve(rhs, x)
	}
}

// The singularity test is relative to each column: a well-conditioned
// system whose unknowns are in wildly different units (one column scaled by
// 1e−12, another by 1e+9) solves on every backend, and a singular one is
// still reported.
func TestColumnScaledSystems(t *testing.T) {
	base := [][]float64{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}
	sc := []float64{1, 1e-12, 1e9}
	for _, be := range []Backend{DenseLU, SparseLU} {
		m, err := NewMatrix(3, be)
		if err != nil {
			t.Fatal(err)
		}
		for i := range base {
			for j := range base[i] {
				if base[i][j] != 0 {
					m.Add(i, j, base[i][j]*sc[j])
				}
			}
		}
		x := make([]float64, 3)
		if err := m.Solve([]float64{1, 2, 3}, x); err != nil {
			t.Fatalf("%v: %v", be, err)
		}
		// Unscaled solution of base·y = b, y_j = x_j·sc_j.
		y := []float64{x[0] * sc[0], x[1] * sc[1], x[2] * sc[2]}
		for i := range base {
			s := 0.0
			for j := range base[i] {
				s += base[i][j] * y[j]
			}
			if math.Abs(s-float64(i+1)) > 1e-12 {
				t.Fatalf("%v: residual row %d: %g", be, i, s)
			}
		}
		sing, _ := NewMatrix(2, be)
		sing.Add(0, 0, 1)
		sing.Add(0, 1, 1e9)
		sing.Add(1, 0, 2)
		sing.Add(1, 1, 2e9)
		if err := sing.Solve([]float64{1, 1}, make([]float64, 2)); !errors.Is(err, ErrSingularMatrix) {
			t.Errorf("%v: singular system: %v", be, err)
		}
	}
}
