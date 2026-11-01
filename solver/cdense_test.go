// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"errors"
	"math/cmplx"
	"math/rand"
	"testing"
)

func TestCDenseRandomResidual(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	for n := 1; n <= 30; n += 5 {
		d, _ := NewCDense(n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				d.Add(i, j, complex(rng.NormFloat64(), rng.NormFloat64()))
			}
			d.Add(i, i, complex(float64(n), 0))
		}
		b := make([]complex128, n)
		for i := range b {
			b[i] = complex(rng.NormFloat64(), rng.NormFloat64())
		}
		x := make([]complex128, n)
		if err := d.Solve(b, x); err != nil {
			t.Fatal(err)
		}
		for i := 0; i < n; i++ {
			var s complex128
			for j := 0; j < n; j++ {
				s += d.At(i, j) * x[j]
			}
			if cmplx.Abs(s-b[i]) > 1e-10 {
				t.Fatalf("n=%d residual %g", n, cmplx.Abs(s-b[i]))
			}
		}
	}
}

func TestCDensePivotAndSingular(t *testing.T) {
	d, _ := NewCDense(2)
	d.Add(0, 1, 1i)
	d.Add(1, 0, 2)
	x := make([]complex128, 2)
	if err := d.Solve([]complex128{1, 4}, x); err != nil {
		t.Fatal(err)
	}
	if x[0] != 2 || cmplx.Abs(x[1]-(-1i)) > 1e-15 {
		t.Fatalf("x = %v", x)
	}
	s, _ := NewCDense(2)
	s.Add(0, 0, 1+1i)
	s.Add(0, 1, 2+2i)
	s.Add(1, 0, 1)
	s.Add(1, 1, 2)
	if err := s.Solve([]complex128{1, 1}, x); !errors.Is(err, ErrSingularMatrix) {
		t.Fatalf("err = %v", err)
	}
}

func TestCDenseSolveDoesNotAllocate(t *testing.T) {
	d, _ := NewCDense(10)
	for i := 0; i < 10; i++ {
		d.Add(i, i, 2+1i)
	}
	b, x := make([]complex128, 10), make([]complex128, 10)
	b[0] = 1
	if a := testing.AllocsPerRun(50, func() { _ = d.Solve(b, x) }); a != 0 {
		t.Fatalf("allocs = %v", a)
	}
}
