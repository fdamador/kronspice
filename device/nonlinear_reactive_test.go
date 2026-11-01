// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// quad integrates f over [a, b] with the composite Simpson rule.
func quad(f func(float64) float64, a, b float64) float64 {
	const n = 20000
	h := (b - a) / n
	s := f(a) + f(b)
	for k := 1; k < n; k++ {
		w := 2.0
		if k%2 == 1 {
			w = 4
		}
		s += w * f(a+float64(k)*h)
	}
	return s * h / 3
}

// The co-energy is ∫₀ˣ f and the energy ∫₀^f(x) x dy = ∫₀ˣ u·f′(u) du, on
// and beyond the curve, on both sides of the origin; a linear law gives
// ½Cv² for both.
func TestStateLawEnergies(t *testing.T) {
	l := newStateLaw(circuit.Element{Curve: []float64{-2, -3, -0.5, -1, 0, 0, 1, 2, 3, 2.5}})
	fx := func(x float64) float64 { y, _ := l.f(x); return y }
	dfx := func(x float64) float64 { _, d := l.f(x); return d }
	for _, x := range []float64{-4, -2, -1.2, -0.25, 0, 0.3, 1, 2.2, 3, 5} {
		if got, want := l.coenergy(x), quad(fx, 0, x); math.Abs(got-want) > 1e-6 {
			t.Errorf("coenergy(%g) = %.9g, want %.9g", x, got, want)
		}
		// dE/dx = x·f′(x) (away from knots), E ≥ 0.
		if e := l.energy(x); e < 0 {
			t.Errorf("energy(%g) = %g < 0", x, e)
		}
		if x == -2 || x == 0 || x == 1 || x == 3 {
			continue
		}
		h := 1e-6
		if got, want := (l.energy(x+h)-l.energy(x-h))/(2*h), x*dfx(x); math.Abs(got-want) > 1e-6 {
			t.Errorf("dE/dx(%g) = %.9g, want x·f′ = %.9g", x, got, want)
		}
	}
	if l.energy(0) != 0 || l.coenergy(0) != 0 {
		t.Errorf("energy(0) = %g, co-energy(0) = %g, want 0", l.energy(0), l.coenergy(0))
	}
	lin := newStateLaw(circuit.Element{Value: 2e-6})
	if e, f := lin.energy(3), lin.coenergy(3); math.Abs(e-9e-6) > 1e-18 || math.Abs(f-9e-6) > 1e-18 {
		t.Errorf("linear energy %g, co-energy %g, want ½Cv² = 9e-6", e, f)
	}
}

// The waveform slope matches a central difference away from corners and
// the left difference at a corner.
func TestWaveSlope(t *testing.T) {
	cases := []struct {
		w  circuit.Waveform
		ts []float64
	}{
		{circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{1, 3, 1e-3, 2e-4, 3e-4, 1e-3, 3e-3}}, []float64{1.1e-3, 1.5e-3, 1.35e-3, 4.25e-3, 2.45e-3}},
		{circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, 1e3, 1e-4, 200, 30}}, []float64{2e-4, 7.7e-4, 1.3e-3}},
		{circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{1, 2, 1e-4, 2e-4, 6e-4, 1e-4}}, []float64{3e-4, 5e-4, 7e-4, 1e-3}},
		{circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 1, 1e-3, 2, 3e-3, 1.5}}, []float64{5e-4, 2e-3}},
	}
	for _, c := range cases {
		w := resolveWave(c.w, 1e-5, 1e-2)
		for _, tt := range c.ts {
			h := 1e-9
			want := (w.value(tt+h) - w.value(tt-h)) / (2 * h)
			if got := w.slope(tt); math.Abs(got-want) > 1e-5*math.Max(1, math.Abs(want)) {
				t.Errorf("%v slope(%g) = %g, want %g", c.w.Kind, tt, got, want)
			}
		}
	}
	// Corner of a PWL: the slope of the segment ending there.
	w := resolveWave(circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 1, 1e-3, 2, 3e-3, 1.5}}, 1e-5, 1e-2)
	if got := w.slope(1e-3); math.Abs(got-1e3) > 1e-9 {
		t.Errorf("corner slope %g, want 1000", got)
	}
}

// A transient Load of the nonlinear elements does not allocate.
func TestNonlinearReactiveLoadDoesNotAllocate(t *testing.T) {
	var devs []Device
	for _, d := range allDevices() {
		switch d.Kind() {
		case circuit.NonlinearCapacitor, circuit.NonlinearInductor:
			devs = append(devs, d)
		}
	}
	ctx := NewSetupContext(2)
	ctx.Step, ctx.Stop = 1e-6, 1e-3
	for _, d := range devs {
		if err := d.Setup(ctx); err != nil {
			t.Fatal(err)
		}
	}
	n := ctx.Size()
	x := make([]float64, n)
	x[0] = 0.4
	lc := &LoadContext{M: &denseStub{n: n, a: make([]float64, n*n)}, RHS: make([]float64, n), X: x,
		Mode: ModeTran, Time: 2e-4, H: 1e-6, Method: Trapezoidal}
	allocs := testing.AllocsPerRun(100, func() {
		for _, d := range devs {
			_ = d.Load(lc)
		}
	})
	if allocs != 0 {
		t.Fatalf("Load allocates %v times", allocs)
	}
}

// Only a varying element models its scale port; a constant one leaves it
// "not modeled".
func TestNonlinearReactivePorts(t *testing.T) {
	still := NewNonlinearCapacitor(circuit.Element{Name: "C", Kind: circuit.NonlinearCapacitor, Nodes: [2]int{1, 0},
		Curve: []float64{-1, -1, 1, 1}})
	if p, _ := still.(*curvedCapacitor).Port(exchange.Mechanical); p.Status != exchange.NotModeled {
		t.Errorf("constant capacitor: mechanical port %v", p.Status)
	}
	moving := NewNonlinearInductor(circuit.Element{Name: "L", Kind: circuit.NonlinearInductor, Nodes: [2]int{1, 0},
		Value: 1, Scale: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.1, 50}}, Channel: exchange.Thermal})
	if p, _ := moving.(*NonlinearInductor).Port(exchange.Thermal); p.Status != exchange.Modeled {
		t.Errorf("varying inductor: thermal port %v", p.Status)
	}
	if p, _ := moving.(*NonlinearInductor).Port(exchange.Mechanical); p.Status != exchange.NotModeled {
		t.Errorf("varying inductor on the thermal channel: mechanical port %v", p.Status)
	}
	if _, ok := still.(Nonlinear); !ok {
		t.Error("a charge curve needs Newton iteration")
	}
	if _, ok := moving.(Nonlinear); ok {
		t.Error("a linear time-varying inductor needs no Newton iteration")
	}
}

// An EXP scale whose defaulted time constants let it reach 0 is refused at
// setup.
func TestExpScaleBound(t *testing.T) {
	d := NewNonlinearCapacitor(circuit.Element{Name: "C", Kind: circuit.NonlinearCapacitor, Nodes: [2]int{1, 0}, Value: 1,
		Scale: circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{0.1, 5, 0, 1, 1e-3, 1e-6}}})
	ctx := NewSetupContext(2)
	ctx.Step, ctx.Stop = 1e-6, 1e-2
	if err := d.Setup(ctx); err == nil {
		t.Fatal("want an error for a scale that reaches 0")
	}
}

// A step crossing a knot is cut to land on it; a value already within
// 10⁻³ of a segment width of the knot, or a crossing in the last 10% of the
// step, is accepted (otherwise steps shrink towards the knot forever).
func TestKnotStep(t *testing.T) {
	l := newStateLaw(circuit.Element{Curve: []float64{-1, -1, 0, 0, 0.1, 0.5, 10, 1}})
	cases := []struct {
		x0, x1, want float64
	}{
		{0.05, 0.25, 0.25},    // up through 0.1 a quarter of the way
		{0.3, -0.7, 0.2},      // down through 0.1 at 20%
		{0.0999999, 0.5, inf}, // on the knot already
		{0.05, 0.0999, inf},   // no knot crossed
		{0.05, 0.1049, inf},   // crossed at 91%
		{-0.5, 0.5, 0.5},      // the first knot crossed (0) decides
		{20, 30, inf},         // beyond the curve
	}
	for _, c := range cases {
		got := l.knotStep(c.x0, c.x1, 1)
		if (c.want == inf) != (got == inf) || (c.want != inf && math.Abs(got-c.want) > 1e-12) {
			t.Errorf("knotStep(%g → %g) = %g, want %g", c.x0, c.x1, got, c.want)
		}
	}
	lin := newStateLaw(circuit.Element{Value: 1e-6})
	if got := lin.knotStep(-1, 1, 1); got != inf {
		t.Errorf("linear law cut a step: %g", got)
	}
}
