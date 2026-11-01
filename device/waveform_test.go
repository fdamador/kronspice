// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func TestPulseValuesAndBreakpoints(t *testing.T) {
	w := resolveWave(circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 5, 1, 0.5, 0.5, 2, 10}}, 0.1, 100)
	cases := map[float64]float64{0: 0, 1: 0, 1.25: 2.5, 1.5: 5, 3.5: 5, 3.75: 2.5, 4: 0, 9: 0, 11.25: 2.5, 12: 5}
	for tt, want := range cases {
		if got := w.value(tt); math.Abs(got-want) > 1e-12 {
			t.Errorf("pulse(%g) = %g, want %g", tt, got, want)
		}
	}
	var got []float64
	for b := w.nextBreakpoint(0); b < 15; b = w.nextBreakpoint(b) {
		got = append(got, b)
	}
	want := []float64{1, 1.5, 3.5, 4, 11, 11.5, 13.5, 14}
	if len(got) != len(want) {
		t.Fatalf("breakpoints %v, want %v", got, want)
	}
	for i := range want {
		if math.Abs(got[i]-want[i]) > 1e-9 {
			t.Fatalf("breakpoints %v, want %v", got, want)
		}
	}
}

func TestPulseDefaults(t *testing.T) {
	// PULSE(0 1) with tstep 1e-3, tstop 1: tr = tf = tstep, pw = per = tstop.
	w := resolveWave(circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1}}, 1e-3, 1)
	if w.a[3] != 1e-3 || w.a[4] != 1e-3 || w.a[5] != 1 || w.a[6] != 1 {
		t.Fatalf("defaults = %v", w.a)
	}
	if got := w.value(0.5e-3); math.Abs(got-0.5) > 1e-12 {
		t.Fatalf("mid-rise = %g", got)
	}
}

func TestSinExpPWL(t *testing.T) {
	s := resolveWave(circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 2, 50}}, 1e-4, 0.1)
	if got := s.value(1.0 / 200); math.Abs(got-3) > 1e-12 { // quarter period
		t.Errorf("sin = %g", got)
	}
	e := resolveWave(circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{0, 1, 0, 1, 5, 1}}, 1e-3, 10)
	if got := e.value(1); math.Abs(got-(1-math.Exp(-1))) > 1e-12 {
		t.Errorf("exp = %g", got)
	}
	p := resolveWave(circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{1, 0, 2, 4, 2, 6}}, 1, 10)
	for tt, want := range map[float64]float64{0: 0, 1.5: 2, 2: 6, 5: 6} {
		if got := p.value(tt); got != want {
			t.Errorf("pwl(%g) = %g, want %g", tt, got, want)
		}
	}
	if b := p.nextBreakpoint(1); b != 2 {
		t.Errorf("pwl breakpoint = %g", b)
	}
}

func TestCapacitorCompanion(t *testing.T) {
	c := NewCapacitor(circuit.Element{Name: "C1", Kind: circuit.Capacitor, Nodes: [2]int{1, 0}, Value: 2e-6})
	ctx := NewSetupContext(2)
	_ = c.Setup(ctx)
	c.vPrev, c.iPrev = 1, 0.5
	m := &denseStub{n: 1, a: make([]float64, 1)}
	rhs := make([]float64, 1)
	_ = c.Load(&LoadContext{M: m, RHS: rhs, Mode: ModeTran, H: 1e-3, Method: Trapezoidal})
	// geq = 2C/h = 4e-3; ihist = −(geq·vPrev + iPrev) = −0.504.
	if math.Abs(m.a[0]-4e-3) > 1e-15 || math.Abs(rhs[0]-0.504) > 1e-15 {
		t.Fatalf("stamp %v rhs %v", m.a, rhs)
	}
	op := c.Operating([]float64{1.1})
	if want := 4e-3*1.1 - 0.504; math.Abs(op.Current-want) > 1e-15 {
		t.Fatalf("current %g, want %g", op.Current, want)
	}
	if e := c.StoredEnergy(op); math.Abs(e-0.5*2e-6*1.21) > 1e-18 {
		t.Fatalf("energy %g", e)
	}
}
