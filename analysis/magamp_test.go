// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// magampCore is a square-loop Jiles–Atherton core (μ0·Ms = 1.51 T,
// 1 cm², 10 cm path); with NG = 200 it can block 2·NG·Ae·Bs = 0.06 V·s,
// more than the 0.040 V·s of a 50 V, 400 Hz half cycle.
var magampCore = map[string]float64{"ms": 1.2e6, "a": 20, "k": 30, "c": 0.05, "alpha": 2e-5,
	"area": 1e-4, "len": 0.1}

const magampVp, magampF, magampRL = 50.0, 400.0, 100.0

// magampOut runs a 400 Hz magnetic amplifier into a 100 Ω load for the
// given control current and returns the mean output over the last
// `average` seconds of `stop`.
func magampOut(t *testing.T, ic, nf, cf, stop, average float64) float64 {
	t.Helper()
	c, err := circuit.NewBuilder("magamp").
		VoltageWave("VAC", "ac", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, magampVp, magampF}}).
		MagAmp("M", "ac", "out", "ctl", "0", "0", circuit.MagAmpSpec{Core: magampCore,
			NG: 200, NC: 100, NF: nf, CF: cf}).
		Resistor("RL", "out", "0", magampRL).
		CurrentSource("IC", "0", "ctl", ic).Build()
	if err != nil {
		t.Fatal(err)
	}
	w, err := Tran(c, TranSpec{Step: 2e-6, Stop: stop, MaxStep: 2e-6}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	l := w.Ledger
	if l.MaxRelativeResidual > 1e-9 {
		t.Errorf("Ic %g: ledger residual %.3g", ic, l.MaxRelativeResidual)
	}
	if math.Abs(l.IntegrationMismatch) > 1e-4*l.EnergyIn {
		t.Errorf("Ic %g: integration mismatch %.3g of %.3g J", ic, l.IntegrationMismatch, l.EnergyIn)
	}
	v, _ := w.NodeVoltage("out")
	s := 0.0
	for k := 1; k < len(w.Time); k++ {
		if w.Time[k-1] >= stop-average-1e-12 {
			s += 0.5 * (v[k] + v[k-1]) * (w.Time[k] - w.Time[k-1])
		}
	}
	return s / average
}

// The mean output rises monotonically with the control current from the
// reset end (the gate blocks most of each half cycle; what remains is the
// magnetizing current NC·|Ic|/NG + Hc·le/NG through the load) to full
// conduction, where the output is the half-wave rectified mean Vp/π less
// the rectifier drop. Every point's ledger closes: the load power comes
// from the AC supply.
func TestMagAmpTransfer(t *testing.T) {
	period := 1 / magampF
	full := magampVp / math.Pi
	prev := math.Inf(-1)
	var out []float64
	for _, ic := range []float64{-0.05, -0.025, 0, 0.025, 0.05} {
		v := magampOut(t, ic, 0, 0, 4*period, 2*period)
		t.Logf("Ic %+.3f A: mean output %.3f V", ic, v)
		if !(v > prev) {
			t.Errorf("Ic %g: output %.4g V not above %.4g V", ic, v, prev)
		}
		prev = v
		out = append(out, v)
	}
	if hi := out[len(out)-1]; !(hi > 0.85*full && hi < full) {
		t.Errorf("full conduction %.4g V, want just below Vp/π = %.4g V", hi, full)
	}
	if lo := out[0]; !(lo < 0.3*full) {
		t.Errorf("reset end %.4g V, want well below Vp/π = %.4g V", lo, full)
	}
}

// Positive feedback (Mali pp. 54–55): with a filter capacitor the load
// current keeps flowing through the feedback winding while the rectifier
// blocks and opposes the reset, so the same step of control current moves
// the output several times as far.
func TestMagAmpFeedbackGain(t *testing.T) {
	period := 1 / magampF
	gain := func(nf float64) float64 {
		a := magampOut(t, -0.05, nf, 20e-6, 8*period, 2*period)
		b := magampOut(t, -0.03, nf, 20e-6, 8*period, 2*period)
		t.Logf("NF %g: %.3f V → %.3f V", nf, a, b)
		return (b - a) / 0.02
	}
	g0, g1 := gain(0), gain(20)
	t.Logf("transfer slope: %.0f V/A without feedback, %.0f V/A with 20 feedback turns", g0, g1)
	if !(g0 > 0) || !(g1 > 3*g0) {
		t.Errorf("feedback slope %.4g V/A, want > 3 × %.4g V/A", g1, g0)
	}
}

func TestMagAmpValidation(t *testing.T) {
	if err := circuit.ValidateMagAmp(circuit.MagAmpSpec{Core: magampCore, NG: 200, NC: 100}); err != nil {
		t.Fatal(err)
	}
	for _, s := range []circuit.MagAmpSpec{
		{Core: magampCore, NG: 0, NC: 100},
		{Core: magampCore, NG: 200, NC: -1},
		{Core: magampCore, NG: 200, NC: 100, NF: -1},
		{Core: magampCore, NG: 200, NC: 100, CF: math.Inf(1)},
	} {
		if err := circuit.ValidateMagAmp(s); !errors.Is(err, circuit.ErrBadValue) {
			t.Errorf("%+v: got %v", s, err)
		}
	}
}
