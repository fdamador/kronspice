// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// royerCore is a small square-loop ferrite-like Jiles–Atherton core:
// μ0·Ms = 0.377 T, 10 mm², 30 mm path.
var royerCore = map[string]float64{"ms": 300e3, "a": 10, "k": 15, "c": 0.1, "alpha": 80e-6,
	"area": 10e-6, "len": 30e-3}

const royerNP, royerNS = 20.0, 10.0

func royerStage(t *testing.T, vin float64) *circuit.Circuit {
	t.Helper()
	c, err := circuit.NewBuilder("royer").
		VoltageWave("VCC", "vcc", "0", circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-6, vin}}).
		Royer("X", "vcc", "0", "s", "0", circuit.RoyerSpec{Core: royerCore,
			NP: royerNP, NB: 4, NS: royerNS, RB: 220, RBias: 4.7e3}).
		Resistor("RL", "s", "0", 100).Build()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// royerCycle returns the mean period of V(c1) − V(c2) over the rising
// crossings after the first, the collector saturation voltage and the
// secondary's peak voltage.
func royerCycle(t *testing.T, w *output.Waveforms) (period, vsat, vs float64) {
	t.Helper()
	c1, _ := w.NodeVoltage("x_c1")
	c2, _ := w.NodeVoltage("x_c2")
	s, _ := w.NodeVoltage("s")
	var cross []float64
	for k := 1; k < len(w.Time); k++ {
		if d0, d1 := c1[k-1]-c2[k-1], c1[k]-c2[k]; d0 < 0 && d1 >= 0 {
			f := -d0 / (d1 - d0)
			cross = append(cross, w.Time[k-1]+f*(w.Time[k]-w.Time[k-1]))
		}
	}
	if len(cross) < 4 {
		t.Fatalf("no sustained oscillation: %d crossings", len(cross))
	}
	period = (cross[len(cross)-1] - cross[1]) / float64(len(cross)-2)
	// The collector voltage in the middle of its low half cycle, over the
	// last full period, and the secondary there.
	t0 := cross[len(cross)-2]
	for _, tm := range []float64{t0 + period/4, t0 + 3*period/4} {
		if v := sampleAt(w.Time, c1, tm); v < sampleAt(w.Time, c2, tm) {
			vsat = v
			vs = math.Abs(sampleAt(w.Time, s, tm))
		}
	}
	return period, vsat, vs
}

// The stage starts from its 1 % base mismatch, switches each time the core
// saturates, and runs at Royer's frequency V/(4·NP·Bsat·Ae) with V the
// supply less the collector saturation voltage. The switching flux density
// that the period implies, B = V/(4·NP·Ae·f), lies between 0.9 and 1.0 of
// μ0·Ms and rises with V: the transistor drops out of saturation when its
// collector current reaches β·Ib, and the base drive grows with the supply,
// so the core is driven a little deeper at 24 V than at 12 V (the formula's
// Bsat is that switching flux, not a material constant). The secondary
// carries the square wave scaled by NS/NP, and the ledger closes.
func TestRoyerFrequency(t *testing.T) {
	type res struct{ f, v float64 }
	var got []res
	for _, vin := range []float64{12, 24} {
		c := royerStage(t, vin)
		w, err := Tran(c, TranSpec{Step: 5e-8, Stop: 1200e-6 / vin, MaxStep: 5e-8}, TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		period, vsat, vs := royerCycle(t, w)
		v := vin - vsat
		f := 1 / period
		bs := 4e-7 * math.Pi * royerCore["ms"]
		want := circuit.RoyerFrequency(v, royerNP, bs, royerCore["area"])
		t.Logf("Vin %g V: f %.4g kHz (Royer %.4g kHz), Vce(sat) %.3f V, secondary %.3f V, steps %d",
			vin, f/1e3, want/1e3, vsat, vs, w.Stats.Accepted)
		if b := bs * want / f; !(b > 0.9*bs && b < 1.001*bs) {
			t.Errorf("Vin %g: f %.4g Hz implies switching flux %.4g T, μ0·Ms = %.4g T", vin, f, b, bs)
		}
		if ws := royerNS / royerNP * v; math.Abs(vs/ws-1) > 0.02 {
			t.Errorf("Vin %g: secondary %.4g V, want NS/NP·V = %.4g V", vin, vs, ws)
		}
		l := w.Ledger
		if l.MaxRelativeResidual > 1e-9 {
			t.Errorf("Vin %g: ledger residual %.3g", vin, l.MaxRelativeResidual)
		}
		if math.Abs(l.IntegrationMismatch) > 1e-4*l.EnergyIn {
			t.Errorf("Vin %g: integration mismatch %.3g of %.3g J", vin, l.IntegrationMismatch, l.EnergyIn)
		}
		got = append(got, res{f, v})
	}
	// f ∝ V/B: twice the supply gives a little less than twice the frequency.
	if r := (got[1].f / got[0].f) / (got[1].v / got[0].v); !(r > 0.9 && r < 1) {
		t.Errorf("(f24/f12)/(V24/V12) = %.4f, want in (0.9, 1)", r)
	}
}

func TestRoyerValidation(t *testing.T) {
	good := circuit.RoyerSpec{Core: royerCore, NP: 20, NB: 4, RB: 220, RBias: 4.7e3}
	if err := circuit.ValidateRoyer(good); err != nil {
		t.Fatal(err)
	}
	for _, s := range []circuit.RoyerSpec{
		{Core: royerCore, NP: 0, NB: 4, RB: 220, RBias: 4.7e3},
		{Core: royerCore, NP: 20, NB: 4, RB: -1, RBias: 4.7e3},
		{Core: royerCore, NP: 20, NB: 4, RB: 220, RBias: 4.7e3, Mismatch: 0.7},
		{Core: royerCore, NP: 20, NB: 4, NS: -1, RB: 220, RBias: 4.7e3},
	} {
		if err := circuit.ValidateRoyer(s); !errors.Is(err, circuit.ErrBadValue) {
			t.Errorf("%+v: got %v", s, err)
		}
	}
}
