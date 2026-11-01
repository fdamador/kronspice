// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// The stamp is the Jacobian of the element's equations: with F(x) = A(x)·x −
// b(x) the nonlinear residual, central differences of F match A(x0)·δ.
func TestBearingJacobian(t *testing.T) {
	for _, p := range []map[string]float64{
		{"kappa": 3e-6, "w0": 15},
		{"kappa": 3e-6, "w0": 15, "thermal": 1, "etamax": 0.1},
		{"kappa": 3e-6, "w0": 15, "vf": 0.7},
		{"kappa": 3e-6, "w0": 15, "vf": 0.7, "thermal": 1, "etamax": 0.1},
	} {
		for _, x0 := range [][]float64{{0.4, 23, 80}, {0.9, -7, 4e-4}, {0.2, 3, -2.5}} {
			e := circuit.Element{Name: "P1", Kind: circuit.BearingMotor, Nodes: [2]int{1, 0}, Extra: []int{2, 0}, Value: 2e-3, Params: p}
			b := NewBearing(e)
			ctx := NewSetupContext(3)
			_ = b.Setup(ctx)
			n := ctx.Size()
			lin := func(x []float64) ([]float64, []float64) {
				st := &denseStub{n: n, a: make([]float64, n*n)}
				rhs := make([]float64, n)
				_ = b.Load(&LoadContext{M: st, RHS: rhs, X: x, Mode: ModeOP})
				return st.a, rhs
			}
			F := func(x []float64) []float64 {
				a, rhs := lin(x)
				f := make([]float64, n)
				for i := 0; i < n; i++ {
					for j := 0; j < n; j++ {
						f[i] += a[i*n+j] * x[j]
					}
					f[i] -= rhs[i]
				}
				return f
			}
			// x0 = V(a), ω, I; small steps so the 1 mA drop transition is resolved
			d := []float64{0.3, -0.7, 0.5 * math.Min(1, math.Abs(x0[2]))}
			a0, _ := lin(x0)
			const h = 1e-6
			xp, xm := make([]float64, n), make([]float64, n)
			for i := range x0 {
				xp[i], xm[i] = x0[i]+h*d[i], x0[i]-h*d[i]
			}
			fp, fm := F(xp), F(xm)
			for i := 0; i < n; i++ {
				var ad float64
				for j := 0; j < n; j++ {
					ad += a0[i*n+j] * d[j]
				}
				fd := (fp[i] - fm[i]) / (2 * h)
				if math.Abs(fd-ad) > 1e-6*(1+math.Abs(ad)) {
					t.Errorf("%v x0=%v row %d: FD %g, stamp %g", p, x0, i, fd, ad)
				}
			}
		}
	}
}

// Load must not allocate (CLAUDE.md §8 rule 7).
func BenchmarkBearingLoad(b *testing.B) {
	e := circuit.Element{Name: "P1", Kind: circuit.BearingMotor, Nodes: [2]int{1, 0}, Extra: []int{2, 0}, Value: 2e-3,
		Params: map[string]float64{"kappa": 3e-6, "w0": 15}}
	d := NewBearing(e)
	ctx := NewSetupContext(3)
	_ = d.Setup(ctx)
	n := ctx.Size()
	st := &denseStub{n: n, a: make([]float64, n*n)}
	lc := &LoadContext{M: st, RHS: make([]float64, n), X: []float64{0.4, 23, 80}, Mode: ModeOP}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = d.Load(lc)
	}
}
