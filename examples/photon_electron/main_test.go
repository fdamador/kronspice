// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The paper's numbers (Williamson & van der Mark 1997, eqs. 5, 8, 12).
func TestPaperNumbers(t *testing.T) {
	if q := Charge() / E; math.Abs(q-0.910) > 0.001 {
		t.Errorf("charge %.4f e, paper 0.91 e", q)
	}
	if a := (GFactor() - 2) / 2; math.Abs(a-0.000962) > 1e-6 {
		t.Errorf("anomaly %.6g, expected 0.000962", a)
	}
	if l := Compton(); math.Abs(l-2.42631023538e-12)/l > 1e-9 {
		t.Errorf("Compton wavelength %.12g", l)
	}
	// With q = e the model's g reduces to QED's first-order 2(1 + α/2π).
	if d := math.Abs(AlphaPrime()*(E/Charge())*(E/Charge()) - AlphaFS); d > 1e-15 {
		t.Errorf("α′ does not scale as (q/e)²α: %g", d)
	}
}

// Plain ring: whole-wave modes; crossed ring: half-wave modes, fundamental
// at half the plain one. AC peaks twin the ladder dispersion.
func TestRingModes(t *testing.T) {
	for _, crossed := range []bool{false, true} {
		r := DefaultRing(crossed)
		got, err := r.Modes(3, 1.05*r.ModeFrequency(3))
		if err != nil {
			t.Fatal(err)
		}
		for m, f := range got {
			mm := m + 1
			if crossed {
				mm = m
			}
			if th := r.ModeFrequency(mm); math.Abs(f-th)/th > 2e-4 {
				t.Errorf("crossed=%v mode %d: AC %.6g, theory %.6g", crossed, m+1, f, th)
			}
		}
	}
	p, c := DefaultRing(false), DefaultRing(true)
	if ratio := c.ModeFrequency(0) / p.ModeFrequency(1); math.Abs(ratio-0.5) > 1e-4 {
		t.Errorf("crossed fundamental / plain = %.5f, want 0.5", ratio)
	}
}

// Round trips: the plain ring returns the pulse upright every trip; the
// crossed ring inverts it on odd trips (720° to return). Ledger closes.
func TestRoundTrips(t *testing.T) {
	for _, crossed := range []bool{false, true} {
		run, err := DefaultRing(crossed).Pulse(2)
		if err != nil {
			t.Fatal(err)
		}
		for _, tr := range run.Trips {
			want := 1.0
			if crossed && tr.N%2 == 1 {
				want = -1
			}
			if tr.Relative*want < 0.5 {
				t.Errorf("crossed=%v trip %d: %.2f of launch, want sign %+.0f", crossed, tr.N, tr.Relative, want)
			}
		}
		if res := run.EIn - run.EDiss - run.EStore; math.Abs(res) > 1e-3*run.EIn {
			t.Errorf("crossed=%v ledger residual %.3g of %.3g J", crossed, res, run.EIn)
		}
	}
}

// Part 3 (coarse grid): light set on the double loop leaves under linear
// Maxwell equations; the field starts divergence-free, stays so, and the
// ledger closes.
func TestKnotLightLeaves(t *testing.T) {
	if testing.Short() {
		t.Skip("field solve")
	}
	k := DefaultKnot()
	k.CellsR, k.HalfR, k.Periods = 8, 4, 1
	r, err := k.Run()
	if err != nil {
		t.Fatal(err)
	}
	last := r.Samples[len(r.Samples)-1]
	if last.T != 1 || last.Inside > 0.01 {
		t.Errorf("after %.2g periods %.3g of the energy is still within 2R", last.T, last.Inside)
	}
	if math.Abs(r.Residual) > 1e-10 || r.DivB > 1e-10 {
		t.Errorf("ledger residual %.3g, div B %.3g", r.Residual, r.DivB)
	}
}

// Part 4, pump 1: one electron per cycle at low drive, uphill too, with
// the gates paying for the lift and the heat.
func TestSingleElectronPump(t *testing.T) {
	p := DefaultSEPump()
	r := p.Pump(1e6, 2, 4)
	if math.Abs(r.PerCycle-1) > 1e-6 || math.Abs(r.Residual) > 1e-6 {
		t.Errorf("1 MHz: %.9f electrons per cycle, ledger %.2g", r.PerCycle, r.Residual)
	}
	if hi := p.Pump(3e8, 2, 4); hi.PerCycle > 0.9 {
		t.Errorf("300 MHz should outrun the tunnelling: %.3f per cycle", hi.PerCycle)
	}
	p.V = 20e-6
	u := p.Pump(1e7, 2, 4)
	if math.Abs(u.PerCycle-1) > 1e-4 || u.BiasWork >= 0 || math.Abs(u.Residual) > 1e-6 {
		t.Errorf("uphill: %.6f per cycle, bias work %.3g J, ledger %.2g", u.PerCycle, u.BiasWork, u.Residual)
	}
}

// Part 4, pump 2: the resonator alone matches kronspice; in step the
// electrons multiply, out of step they die; the ledger closes.
func TestMultipactorPump(t *testing.T) {
	m := DefaultMultipactor()
	if tw, err := m.TankTwin(5); err != nil || tw > 0.01 {
		t.Fatalf("resonator twin %.3g (%v)", tw, err)
	}
	r := m.Run(300, true)
	peak := 0.0
	for _, s := range r.Samples {
		peak = math.Max(peak, s.Electrons)
	}
	if peak < 1e4*m.Seed || r.Collected <= 0 || math.Abs(r.Residual) > 1e-5 {
		t.Errorf("peak %.3g electrons, collected %.3g, ledger %.2g", peak, r.Collected, r.Residual)
	}
	off := m
	off.V0 = 300
	if ro := off.Run(40, true); ro.Samples[len(ro.Samples)-1].Electrons > 1 {
		t.Errorf("out of step still has %.3g electrons", ro.Samples[len(ro.Samples)-1].Electrons)
	}
}

// Part 4, pump 3: the bank reaches the ladder's −2·n·peak (to 1 %) and
// the ledger closes.
func TestGroundPump(t *testing.T) {
	g := DefaultGroundPump()
	r, err := g.Run()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.Bank/r.Ideal-1) > 0.01 || r.Bank >= 0 {
		t.Errorf("bank %.1f V, ideal %.1f V", r.Bank, r.Ideal)
	}
	if math.Abs(r.Mismatch) > 1e-4*r.EnergyIn {
		t.Errorf("ledger mismatch %.3g of %.3g J", r.Mismatch, r.EnergyIn)
	}
}

// Part 5: the magnetron cavity matches kronspice, a short evolution
// improves on its start, the winner stays inside the limits and its
// ledger closes.
func TestMagnetronEvolution(t *testing.T) {
	p := magDesign([]float64{1.3e-3, 1.5e-4, 1e-3, 0.4, 1})
	if tw, err := p.Multipactor(800).TankTwin(5); err != nil || tw > 0.01 {
		t.Fatalf("cavity twin %.3g (%v)", tw, err)
	}
	o := DefaultEvolve()
	o.Population, o.Generations, o.Periods = 16, 8, 800
	ev, err := Evolve(o)
	if err != nil {
		t.Fatal(err)
	}
	if ev.Best.Current <= ev.First.Current || ev.Best.Current < 1e-6 {
		t.Errorf("evolved %.3g A vs first generation %.3g A", ev.Best.Current, ev.First.Current)
	}
	b := ev.Best
	if b.PeakField > b.Kilpatrick || b.V0 > MaxGapV || math.Abs(b.Residual) > 1e-4 {
		t.Errorf("field %.3g of Kilpatrick, V0 %.0f V, ledger %.2g", b.PeakField/b.Kilpatrick, b.V0, b.Residual)
	}
	// Without electrons the cavity takes 4β/(1+β)² of the available
	// power (coupling theory); the rest goes to the water load.
	for _, beta := range []float64{1e-3, 0.3, 1, 3} {
		q := p
		q.Beta = beta
		m := q.Multipactor(800)
		r := m.Run(400, false)
		got := (r.SourceWork - r.PortOut) / (400 / q.Mag.F) / q.Mag.Pout
		if want := 4 * beta / sq(1+beta); math.Abs(got-want) > 1e-3*want+1e-6 {
			t.Errorf("β = %g: cavity takes %.5f of the power, coupling theory %.5f", beta, got, want)
		}
	}
}

// Part 6: N gaps sharing the magnetron are N copies of one gap fed Pout/N
// (bookkeeping check), and the 100 kW tube's numbers are those declared.
func TestGapArray(t *testing.T) {
	p := magDesign([]float64{0.6e-3, 0.4e-4, 0.7, 0.44, 1})
	p.N = 1000
	a := p.Evaluate(600, 200)
	q := p
	q.N = 1
	q.Mag.Pout /= 1000
	b := q.Evaluate(600, 200)
	if math.Abs(a.Current-1000*b.Current) > 1e-9*a.Current || math.Abs(a.Delivered-1000*b.Delivered) > 1e-9*math.Abs(a.Delivered) {
		t.Errorf("array %.4g A, %.4g W vs 1000 × single %.4g A, %.4g W", a.Current, a.Delivered, 1000*b.Current, 1000*b.Delivered)
	}
	if m := StellantCWM100; math.Abs(m.Mains()-(100e3/0.85+12.6*69)) > 1e-6 {
		t.Errorf("CWM-100kW mains %.1f W", m.Mains())
	}
}
