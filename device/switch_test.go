// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// newTestSwitch builds a voltage-controlled switch between nodes 1 and 0,
// controlled by node 2 against ground (unknowns: V1, V2).
func newTestSwitch(t testing.TB, p map[string]float64) (*Switch, int) {
	t.Helper()
	e := circuit.Element{Name: "S1", Kind: circuit.Switch, Nodes: [2]int{1, 0}, CtrlNodes: [2]int{2, 0}, Model: "sw"}
	s := NewSwitch(e, circuit.Model{Name: "sw", Type: "sw", Params: p})
	ctx := NewSetupContext(3)
	if err := s.Setup(ctx); err != nil {
		t.Fatal(err)
	}
	return s, ctx.Size()
}

func loadSwitch(s *Switch, n int, x []float64, mode Mode) ([]float64, []float64) {
	st := &denseStub{n: n, a: make([]float64, n*n)}
	rhs := make([]float64, n)
	_ = s.Load(&LoadContext{M: st, RHS: rhs, X: x, Mode: mode})
	return st.a, rhs
}

// Hysteresis: on above VT + VH, off below VT − VH, and the accepted state
// held in between (ngspice "Switches"; SPICE3 SWload).
func TestSwitchHysteresis(t *testing.T) {
	s, n := newTestSwitch(t, map[string]float64{"vt": 1, "vh": 0.2, "ron": 1, "roff": 1e6})
	for _, c := range []struct {
		ctl  float64
		want bool
	}{{0, false}, {1.1, false}, {1.3, true}, {1, true}, {0.85, true}, {0.7, false}, {1.15, false}} {
		x := []float64{1, c.ctl}
		a, _ := loadSwitch(s, n, x, ModeOP)
		g := 1e-6
		if c.want {
			g = 1
		}
		if math.Abs(a[0]-g) > 1e-12*g || !s.ConvTest(x, 1e-3, 1e-12) {
			t.Errorf("ctl %g: stamp %g, want %g", c.ctl, a[0], g)
		}
		s.Accept(x, 0, ModeOP)
		if s.On() != c.want {
			t.Errorf("ctl %g: on = %v, want %v", c.ctl, s.On(), c.want)
		}
	}
}

// A transient step holds the accepted state; the step is cut just past the
// threshold crossing, unless the crossing is in its last hundredth.
func TestSwitchLimitStep(t *testing.T) {
	s, n := newTestSwitch(t, map[string]float64{"vt": 1, "vh": 0.2})
	s.Accept([]float64{0, 0}, 0, ModeOP)
	x := []float64{0, 2.4} // control ramps 0 → 2.4 over the step: crossing at 1.2, frac 0.5
	a, _ := loadSwitch(s, n, x, ModeTran)
	if a[0] != 1/circuit.SwitchRoffDefault || !s.ConvTest(x, 1e-3, 1e-12) {
		t.Errorf("transient stamp %g: state not held", a[0])
	}
	tc := &TruncContext{H: 1e-6}
	if h := s.LimitStep(x, tc); math.Abs(h-0.505e-6) > 1e-15 {
		t.Errorf("LimitStep = %g, want 0.505e-6", h)
	}
	if h := s.LimitStep([]float64{0, 1.2 / 0.995}, tc); !math.IsInf(h, 1) {
		t.Errorf("crossing in the last hundredth: LimitStep = %g, want +Inf", h)
	}
	s.Accept(x, 0, ModeTran)
	if !s.On() || !s.Restart() {
		t.Error("accepted crossing did not turn the switch on and restart integration")
	}
}

// The smooth switch's stamp is the Jacobian of i = g(c)·v (central
// differences).
func TestSmoothSwitchJacobian(t *testing.T) {
	s, n := newTestSwitch(t, map[string]float64{"vt": 1, "vh": -0.5, "ron": 10, "roff": 1e5})
	for _, x0 := range [][]float64{{2, 0.8}, {-3, 1.3}, {0.5, 1.0}, {1, 2}} {
		F := func(x []float64) []float64 {
			a, rhs := loadSwitch(s, n, x, ModeOP)
			f := make([]float64, n)
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					f[i] += a[i*n+j] * x[j]
				}
				f[i] -= rhs[i]
			}
			return f
		}
		d := []float64{0.3, -0.7}
		a0, _ := loadSwitch(s, n, x0, ModeOP)
		const h = 1e-6
		fp := F([]float64{x0[0] + h*d[0], x0[1] + h*d[1]})
		fm := F([]float64{x0[0] - h*d[0], x0[1] - h*d[1]})
		for i := 0; i < n; i++ {
			var ad float64
			for j := 0; j < n; j++ {
				ad += a0[i*n+j] * d[j]
			}
			if fd := (fp[i] - fm[i]) / (2 * h); math.Abs(fd-ad) > 1e-6*(1+math.Abs(ad)) {
				t.Errorf("x0=%v row %d: FD %g, stamp %g", x0, i, fd, ad)
			}
		}
		// The curve's ends are RON and ROFF.
		if g, _ := s.smoothG(1.5); g != 0.1 {
			t.Errorf("g(VT+|VH|) = %g", g)
		}
		if g, _ := s.smoothG(0.5); g != 1e-5 {
			t.Errorf("g(VT−|VH|) = %g", g)
		}
	}
}

// Load must not allocate (CLAUDE.md §8 rule 7).
func BenchmarkSwitchLoad(b *testing.B) {
	for _, c := range []struct {
		name string
		p    map[string]float64
	}{{"abrupt", map[string]float64{"vt": 1, "vh": 0.2}}, {"smooth", map[string]float64{"vt": 1, "vh": -0.2}}} {
		b.Run(c.name, func(b *testing.B) {
			s, n := newTestSwitch(b, c.p)
			st := &denseStub{n: n, a: make([]float64, n*n)}
			lc := &LoadContext{M: st, RHS: make([]float64, n), X: []float64{1, 1.05}, Mode: ModeOP}
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				_ = s.Load(lc)
			}
		})
	}
}
