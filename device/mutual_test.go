// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

type cStub struct {
	n int
	a []complex128
}

func (c *cStub) Add(i, j int, v complex128) {
	if i >= 0 && j >= 0 {
		c.a[i*c.n+j] += v
	}
}

func mutualPair(t testing.TB) (*Inductor, *Inductor, *Mutual, int) {
	l1 := NewInductor(circuit.Element{Name: "L1", Kind: circuit.Inductor, Nodes: [2]int{1, 0}, Value: 1e-3})
	l2 := NewInductor(circuit.Element{Name: "L2", Kind: circuit.Inductor, Nodes: [2]int{2, 0}, Value: 4e-3})
	k := NewMutual(circuit.Element{Name: "K1", Kind: circuit.Mutual, Value: 0.5, Couples: [2]string{"L1", "L2"}}, 1e-3, 4e-3)
	ctx := NewSetupContext(3)
	for _, d := range []Device{l1, l2, k} {
		if err := d.Setup(ctx); err != nil {
			t.Fatal(err)
		}
	}
	if err := k.Link(ctx); err != nil {
		t.Fatal(err)
	}
	return l1, l2, k, ctx.Size()
}

// The K stamps: −jωM between the branch rows in AC; −(2M/h) with the
// previous partner current on the right-hand side for a trapezoidal step,
// −M/h for backward Euler; nothing at DC. M = k·√(L1·L2) = 1 mH.
func TestMutualStamps(t *testing.T) {
	_, _, k, n := mutualPair(t)
	if math.Abs(k.M()-1e-3) > 1e-18 {
		t.Fatalf("M = %g", k.M())
	}
	c := &cStub{n: n, a: make([]complex128, n*n)}
	_ = k.LoadAC(&ACContext{M: c, Omega: 1e4})
	if z := c.a[k.k1*n+k.k2]; math.Abs(imag(z)+10) > 1e-12 || real(z) != 0 || c.a[k.k2*n+k.k1] != z {
		t.Fatalf("AC stamps %v %v", c.a[k.k1*n+k.k2], c.a[k.k2*n+k.k1])
	}
	x := make([]float64, n)
	x[k.k1], x[k.k2] = 0.2, -0.1
	k.Accept(x, 0, ModeTran)
	for _, tc := range []struct {
		method Method
		f      float64
	}{{Trapezoidal, 2e-3 / 1e-6}, {BackwardEuler, 1e-3 / 1e-6}} {
		m := &denseStub{n: n, a: make([]float64, n*n)}
		rhs := make([]float64, n)
		_ = k.Load(&LoadContext{M: m, RHS: rhs, X: x, Mode: ModeTran, H: 1e-6, Method: tc.method})
		if math.Abs(m.a[k.k1*n+k.k2]+tc.f) > 1e-12*tc.f || m.a[k.k2*n+k.k1] != m.a[k.k1*n+k.k2] {
			t.Errorf("method %d: matrix %g, %g; want %g", tc.method, m.a[k.k1*n+k.k2], m.a[k.k2*n+k.k1], -tc.f)
		}
		if math.Abs(rhs[k.k1]-(-tc.f*-0.1)) > 1e-9 || math.Abs(rhs[k.k2]-(-tc.f*0.2)) > 1e-9 {
			t.Errorf("method %d: rhs %g, %g", tc.method, rhs[k.k1], rhs[k.k2])
		}
	}
	m := &denseStub{n: n, a: make([]float64, n*n)}
	_ = k.Load(&LoadContext{M: m, RHS: make([]float64, n), X: x, Mode: ModeOP})
	for _, v := range m.a {
		if v != 0 {
			t.Fatal("K stamps something at DC")
		}
	}
	k.Operating(x)
	if e := k.StoredEnergy(OperatingPoint{}); math.Abs(e-(1e-3*0.2*-0.1)) > 1e-18 {
		t.Errorf("stored M·i1·i2 = %g", e)
	}
}

func BenchmarkMutualLoad(b *testing.B) {
	_, _, k, n := mutualPair(b)
	m := &denseStub{n: n, a: make([]float64, n*n)}
	ctx := &LoadContext{M: m, RHS: make([]float64, n), X: make([]float64, n), Mode: ModeTran, H: 1e-9, Method: Trapezoidal}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = k.Load(ctx)
	}
}
