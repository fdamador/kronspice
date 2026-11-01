// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// testMachineSpec has two windings, a rotating mutual, a second-harmonic
// saliency and a speed matrix, so every term of the element is exercised.
func testMachineSpec() *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0], s.L0[1][1] = 10e-3, 20e-3
	h1 := s.Harmonic(1)
	h1.Cos[0][1], h1.Cos[1][0] = 5e-3, 5e-3
	h2 := s.Harmonic(2)
	h2.Sin[0][0] = 2e-3
	s.G[0][1], s.G[1][0] = 1e-3, -3e-3
	s.Theta0 = 0.4
	return s
}

// The linearization stamped by Load is the Jacobian of the element's
// equations: the linear model at x0 predicts the model at x0 + δ to
// second order in δ, in a trapezoidal transient step and at an operating
// point.
func TestMachineJacobian(t *testing.T) {
	// Nodes: 1 = winding 1 +, 2 = winding 2 +, 3 = port +; grounds.
	e := circuit.Element{Name: "N1", Kind: circuit.Machine, Nodes: [2]int{1, 0}, Extra: []int{2, 0, 3, 0}}
	r := rand.New(rand.NewPCG(1, 2))
	for _, mode := range []Mode{ModeTran, ModeOP} {
		m := NewMachine(e, testMachineSpec())
		ctx := NewSetupContext(4)
		if err := m.Setup(ctx); err != nil {
			t.Fatal(err)
		}
		n := ctx.Size() // 3 nodes + 2 currents + θ
		x0 := make([]float64, n)
		for i := range x0 {
			x0[i] = r.NormFloat64()
		}
		m.Accept(x0, 0, ModeOP)
		m.thPrev, m.omPrev = 0.3, 2
		lin := func(x []float64) ([]float64, []float64) {
			st := &denseStub{n: n, a: make([]float64, n*n)}
			rhs := make([]float64, n)
			lc := &LoadContext{M: st, RHS: rhs, X: x, Mode: mode, H: 1e-3, Method: Trapezoidal}
			if err := m.Load(lc); err != nil {
				t.Fatal(err)
			}
			return st.a, rhs
		}
		apply := func(a, rhs, y []float64) []float64 {
			f := make([]float64, n)
			for i := 0; i < n; i++ {
				s := -rhs[i]
				for j := 0; j < n; j++ {
					s += a[i*n+j] * y[j]
				}
				f[i] = s
			}
			return f
		}
		a0, r0 := lin(x0)
		for _, eps := range []float64{1e-3, 1e-4} {
			dx := make([]float64, n)
			x1 := make([]float64, n)
			for i := range dx {
				dx[i] = eps * r.NormFloat64()
				x1[i] = x0[i] + dx[i]
			}
			a1, r1 := lin(x1)
			pred, exact := apply(a0, r0, x1), apply(a1, r1, x1)
			worst := 0.0
			for i := range pred {
				worst = math.Max(worst, math.Abs(pred[i]-exact[i]))
			}
			if worst > 50*eps*eps {
				t.Errorf("mode %d, δ = %g: linear model off by %g (not second order)", mode, eps, worst)
			}
		}
	}
}
