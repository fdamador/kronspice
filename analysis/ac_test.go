// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func TestFrequencies(t *testing.T) {
	f, _ := ACSpec{Variation: "dec", Points: 10, FStart: 1, FStop: 1000}.Frequencies()
	if len(f) != 31 || math.Abs(f[30]-1000) > 1e-9 {
		t.Fatalf("dec: %d points, last %g", len(f), f[len(f)-1])
	}
	f, _ = ACSpec{Variation: "oct", Points: 2, FStart: 1, FStop: 8}.Frequencies()
	if len(f) != 7 || math.Abs(f[2]-2) > 1e-12 {
		t.Fatalf("oct: %v", f)
	}
	f, _ = ACSpec{Variation: "lin", Points: 5, FStart: 0, FStop: 100}.Frequencies()
	if len(f) != 5 || f[4] != 100 || f[1] != 25 {
		t.Fatalf("lin: %v", f)
	}
	if _, err := (ACSpec{Variation: "dec", Points: 10, FStart: 0, FStop: 10}).Frequencies(); !errors.Is(err, ErrBadACSpec) {
		t.Fatalf("err = %v", err)
	}
}

func rcLowpass(t *testing.T) *circuit.Circuit {
	t.Helper()
	c, err := circuit.NewBuilder("rc").
		Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "in", "0").
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1/(2*math.Pi*1e3*1e3)).Build()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func TestACLowpassAndCutoff(t *testing.T) {
	r, err := AC(rcLowpass(t), ACSpec{Variation: "dec", Points: 50, FStart: 1, FStop: 1e6}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, err := Transfer(r, "out", "in")
	if err != nil {
		t.Fatal(err)
	}
	for k, f := range r.Freq {
		want := 1 / complex(1, f/1e3)
		if cmplx.Abs(h[k]-want) > 1e-12 {
			t.Fatalf("H(%g) = %v, want %v", f, h[k], want)
		}
	}
	fl, fh, err := Cutoff3dB(r.Freq, h)
	if err != nil {
		t.Fatal(err)
	}
	// Interpolation between 50-per-decade points: within 0.1 % of 1 kHz.
	if fl != 0 || math.Abs(fh-1e3) > 1 {
		t.Fatalf("cutoff = %g, %g", fl, fh)
	}
	if g, _ := GainAt(r.Freq, h, 1e3); math.Abs(g-1/math.Sqrt2) > 1e-4 {
		t.Fatalf("gain at fc = %g", g)
	}
}

func TestACLedgerBalances(t *testing.T) {
	r, err := AC(rcLowpass(t), ACSpec{Variation: "dec", Points: 5, FStart: 10, FStop: 1e5}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	l := r.Ledger
	for k, f := range r.Freq {
		// Average power from the 1 V source equals ½|I|²R in the resistor.
		ir, _ := r.ElementCurrent("R1")
		pr := 0.5 * cmplx.Abs(ir[k]) * cmplx.Abs(ir[k]) * 1e3
		if math.Abs(l.Inflow[k]-pr) > 1e-15 || math.Abs(l.Dissipation[k]-pr) > 1e-15 {
			t.Fatalf("f=%g: inflow %g, dissipation %g, want %g", f, l.Inflow[k], l.Dissipation[k], pr)
		}
		if l.Reactive[k] >= 0 {
			t.Fatalf("capacitor reactive power should be negative (absorbs leading current), got %g", l.Reactive[k])
		}
	}
	if l.MaxRelativeResidual > 1e-12 {
		t.Fatalf("residual %g", l.MaxRelativeResidual)
	}
}

func TestACBandpassPeak(t *testing.T) {
	L, C, R := 1e-3, 1e-9, 100.0
	c, _ := circuit.NewBuilder("bp").
		Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "in", "0").
		Inductor("L1", "in", "a", L).Capacitor("C1", "a", "out", C).Resistor("R1", "out", "0", R).Build()
	// 5000 points per decade: interpolation error is then far below the
	// 2e-3 bandwidth tolerance (at 200 per decade, points are 1.2 % apart).
	r, err := AC(c, ACSpec{Variation: "dec", Points: 5000, FStart: 1e5, FStop: 2.5e5}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := Transfer(r, "out", "0")
	fl, fh, err := Cutoff3dB(r.Freq, h)
	if err != nil {
		t.Fatal(err)
	}
	// Series RLC: bandwidth R/(2πL), centre 1/(2π√LC) = geometric mean of the edges.
	f0 := 1 / (2 * math.Pi * math.Sqrt(L*C))
	bw := R / (2 * math.Pi * L)
	if math.Abs((fh-fl)-bw) > 2e-3*bw || math.Abs(math.Sqrt(fl*fh)-f0) > 1e-3*f0 {
		t.Fatalf("band %g–%g Hz: bandwidth %g (want %g), centre %g (want %g)", fl, fh, fh-fl, bw, math.Sqrt(fl*fh), f0)
	}
}

func TestDCSweep(t *testing.T) {
	c, _ := circuit.NewBuilder("div").VoltageSource("V1", "in", "0", 1).
		Resistor("R1", "in", "out", 1e3).Resistor("R2", "out", "0", 3e3).
		CurrentSource("I1", "0", "out", 0).Build()
	r, err := DC(c, []DCSweep{{Source: "V1", Start: 0, Stop: 4, Step: 0.5}, {Source: "i1", Start: 0, Stop: 2e-3, Step: 1e-3}}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	vout, _ := r.NodeVoltage("out")
	if len(vout) != 9*3 {
		t.Fatalf("points = %d", len(vout))
	}
	for k := range vout {
		v, i := r.Values[0][k], r.Values[1][k]
		want := 0.75*v + i*750 // Thevenin: 3/4 of V, plus I through 1k∥3k
		if math.Abs(vout[k]-want) > 1e-12 {
			t.Fatalf("V1=%g I1=%g: v(out) = %g, want %g", v, i, vout[k], want)
		}
	}
	if r.MaxRelativeResidual > 1e-12 {
		t.Fatalf("residual %g", r.MaxRelativeResidual)
	}
	// The source's own value is restored afterwards.
	op, _ := OP(c, Options{})
	if v, _ := op.Voltage("out"); math.Abs(v-0.75) > 1e-12 {
		t.Fatalf("override not cleared: v(out) = %g", v)
	}
	if _, err := DC(c, []DCSweep{{Source: "V9", Start: 0, Stop: 1, Step: 1}}, Options{}); !errors.Is(err, ErrUnknownSource) {
		t.Fatalf("err = %v", err)
	}
	if _, err := DC(c, []DCSweep{{Source: "V1", Start: 0, Stop: 1, Step: -1}}, Options{}); !errors.Is(err, ErrBadDCSpec) {
		t.Fatalf("err = %v", err)
	}
}

func TestTransientMeasurements(t *testing.T) {
	// RC step: 10–90 % rise time is τ·ln 9.
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-9, 1e-9, 1, 2}}
	c, _ := circuit.NewBuilder("rc").VoltageWave("V1", "in", "0", w).
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1e-6).Build()
	wf, err := Tran(c, TranSpec{Step: 5e-6, Stop: 10e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("out")
	tr, err := RiseTime(wf.Time, v, 0.1, 0.9)
	if err != nil {
		t.Fatal(err)
	}
	if want := 1e-3 * math.Log(9); math.Abs(tr-want) > 5e-3*want {
		t.Fatalf("rise time %g, want %g", tr, want)
	}
	if ov, _ := Overshoot(v); ov > 1e-6 {
		t.Fatalf("RC overshoot %g", ov)
	}

	// Series RLC: overshoot e^(−απ/ωd).
	R, L, C := 10.0, 1e-3, 1e-6
	alpha := R / (2 * L)
	wd := math.Sqrt(1/(L*C) - alpha*alpha)
	w2 := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-10, 1e-10, 1, 2}}
	c2, _ := circuit.NewBuilder("rlc").VoltageWave("V1", "in", "0", w2).
		Resistor("R1", "in", "a", R).Inductor("L1", "a", "b", L).Capacitor("C1", "b", "0", C).Build()
	wf2, err := Tran(c2, TranSpec{Step: 1e-6, Stop: 3e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	vb, _ := wf2.NodeVoltage("b")
	ov, _ := Overshoot(vb)
	if want := math.Exp(-alpha * math.Pi / wd); math.Abs(ov-want) > 5e-3 {
		t.Fatalf("overshoot %g, want %g", ov, want)
	}
	if st, err := SettlingTime(wf2.Time, vb, 0.02); err != nil || st <= 0 || st > 3e-3 {
		t.Fatalf("settling %g, %v", st, err)
	}
	if _, err := Crossing(wf.Time, v, 2, 1, 0); !errors.Is(err, ErrNoMeasurement) {
		t.Fatalf("crossing of an unreached level: %v", err)
	}
}

func BenchmarkACLowpass100(b *testing.B) {
	c, _ := circuit.NewBuilder("rc").
		Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "in", "0").
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1e-7).Build()
	spec := ACSpec{Variation: "dec", Points: 20, FStart: 1, FStop: 1e5}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := AC(c, spec, Options{}); err != nil {
			b.Fatal(err)
		}
	}
}
