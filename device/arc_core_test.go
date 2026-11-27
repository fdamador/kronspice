// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

func testArcElement() circuit.Element {
	return circuit.Element{Name: "A1", Kind: circuit.Arc, Nodes: [2]int{1, 0}, Model: circuit.ArcHabedank,
		Params: map[string]float64{"tauc": 1e-6, "u0": 40, "taum": 2e-6, "p0": 3e3, "g0": 0.5}}
}

func testCoreElement() circuit.Element {
	return circuit.Element{Name: "U1", Kind: circuit.Core, Nodes: [2]int{1, 0}, Extra: []int{2, 0},
		Params: map[string]float64{"ms": 1.6e6, "a": 1100, "k": 400, "c": 0.2, "alpha": 1.6e-3,
			"area": 1e-4, "len": 0.1, "n1": 100, "n2": 25}}
}

// At an operating point the arc is its frozen conductance G0 (Habedank:
// both parts at 2·G0 in series), in the branch row g·v − I = 0.
func TestArcStamp(t *testing.T) {
	a := NewArc(testArcElement())
	ctx := NewSetupContext(2)
	_ = a.Setup(ctx)
	n := ctx.Size()
	st := &denseStub{n: n, a: make([]float64, n*n)}
	lc := &LoadContext{M: st, RHS: make([]float64, n), X: make([]float64, n), Mode: ModeOP}
	if err := a.Load(lc); err != nil {
		t.Fatal(err)
	}
	br := a.br
	if got := st.a[br*n+0]; math.Abs(got-0.5) > 1e-15 {
		t.Errorf("g stamp %g, want 0.5", got)
	}
	if st.a[0*n+br] != 1 || st.a[br*n+br] != -1 {
		t.Errorf("KCL / branch stamps %v", st.a)
	}
	op := OperatingPoint{Voltage: 10, Current: 5, PowerAbsorbed: 50}
	if got := a.ExteriorPower(exchange.Thermal, op); got != -50 {
		t.Errorf("thermal port %g, want −50 W", got)
	}
}

// The core's transient stamp is the derivative of its Faraday residual:
// compare with finite differences of the induced voltage at an iterate.
func TestCoreStampFD(t *testing.T) {
	c := NewCore(testCoreElement())
	ctx := NewSetupContext(3)
	_ = c.Setup(ctx)
	n := ctx.Size()
	x := make([]float64, n)
	x[c.br[0]], x[c.br[1]] = 2, -1 // H = 1750 A/m
	st := &denseStub{n: n, a: make([]float64, n*n)}
	const h = 1e-5
	lc := &LoadContext{M: st, RHS: make([]float64, n), X: x, Mode: ModeTran, H: h, Method: BackwardEuler}
	if err := c.Load(lc); err != nil {
		t.Fatal(err)
	}
	flux := func(i1, i2 float64) float64 { // N1·Ae·B
		hh := (100*i1 + 25*i2) / 0.1
		m, _, _ := c.integrate(0, 0, 0, hh)
		return 100 * 1e-4 * mu0 * (hh + m)
	}
	const d = 1e-6
	for j, br := range c.br {
		xp, xm := [2]float64{2, -1}, [2]float64{2, -1}
		xp[j] += d
		xm[j] -= d
		fd := (flux(xp[0], xp[1]) - flux(xm[0], xm[1])) / (2 * d) / h
		if got := -st.a[c.br[0]*n+br]; math.Abs(got-fd) > 2e-3*math.Abs(fd) { // end-point dM/dH vs the sub-stepped path
			t.Errorf("∂v1/∂i%d: stamp %g, FD %g", j+1, got, fd)
		}
	}
}

// Load must not allocate (CLAUDE.md §8 rule 7).
func BenchmarkArcLoad(b *testing.B) {
	a := NewArc(testArcElement())
	ctx := NewSetupContext(2)
	_ = a.Setup(ctx)
	n := ctx.Size()
	st := &denseStub{n: n, a: make([]float64, n*n)}
	lc := &LoadContext{M: st, RHS: make([]float64, n), X: []float64{30, 2}, Mode: ModeTran, H: 1e-8, Method: Trapezoidal}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = a.Load(lc)
	}
}

func BenchmarkCoreLoad(b *testing.B) {
	c := NewCore(testCoreElement())
	ctx := NewSetupContext(3)
	_ = c.Setup(ctx)
	n := ctx.Size()
	st := &denseStub{n: n, a: make([]float64, n*n)}
	x := make([]float64, n)
	x[c.br[0]] = 0.3
	lc := &LoadContext{M: st, RHS: make([]float64, n), X: x, Mode: ModeTran, H: 1e-6, Method: Trapezoidal}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = c.Load(lc)
	}
}
