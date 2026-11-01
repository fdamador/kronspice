// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func rel(a, b float64) float64 { return math.Abs(a-b) / math.Abs(b) }

// TestModelsMatch: the copied cell and gap models give the numbers of
// examples/sparkgap_diode.
func TestModelsMatch(t *testing.T) {
	c := DefaultCell(Pb, Al)
	if i := c.Current(100); rel(i, 980) > 1e-12 {
		t.Errorf("lead–aluminium at +100 V: %g A, want 980", i)
	}
	if cf := c.FilmC(); rel(cf, 4.378444522813186e-06) > 1e-9 {
		t.Errorf("film %g F", cf)
	}
	if v := DefaultGap().Vpos(); math.Abs(v-2415) > 1 {
		t.Errorf("gap breakdown %g V, want 2415", v)
	}
}

// TestSolverRC: a 1 µF capacitor at 10 V discharging into 1 kΩ follows
// exp(−t/RC) to backward Euler's first-order accuracy, and the ledger
// closes.
func TestSolverRC(t *testing.T) {
	n := NewNet()
	c := n.C("C", "a", "0", 1e-6, 10)
	n.R("R", "a", "0", 1e3)
	if err := run(n, 1e-3, 1e-7, nil); err != nil {
		t.Fatal(err)
	}
	if want := 10 * math.Exp(-1); rel(c.v, want) > 1e-3 {
		t.Errorf("v(1 ms) = %g, want %g", c.v, want)
	}
	l := n.Ledger()
	if math.Abs(l.Residual) > 1e-12 || rel(l.Heat["R"]+l.NumLoss, l.Stored0-l.Stored) > 1e-9 {
		t.Errorf("ledger %+v", l)
	}
}

// TestSupplyTwin: the solver and kronspice agree on the cell bridge.
func TestSupplyTwin(t *testing.T) {
	s, l, err := RunSupply()
	if err != nil {
		t.Fatal(err)
	}
	if rel(s.Vdc, s.KronVdc) > 1e-3 {
		t.Errorf("solver %g V, kronspice %g V", s.Vdc, s.KronVdc)
	}
	if math.Abs(l.Residual) > 1e-6*l.Supplied {
		t.Errorf("ledger residual %g J", l.Residual)
	}
}

// TestSurge: protection keeps the bridge cells below their sparking voltage.
func TestSurge(t *testing.T) {
	u, err := RunSurge(false, 4000)
	if err != nil {
		t.Fatal(err)
	}
	p, err := RunSurge(true, 4000)
	if err != nil {
		t.Fatal(err)
	}
	if u.PeakArm < 140 || p.PeakArm > 140 || p.BridgeSpark > 0.1*u.BridgeSpark || p.Strikes == 0 {
		t.Errorf("unprotected %g V %g J; protected %g V %g J, %d strikes", u.PeakArm, u.BridgeSpark, p.PeakArm, p.BridgeSpark, p.Strikes)
	}
	for _, c := range []SurgeCase{u, p} {
		if math.Abs(c.Ledger.Residual) > 1e-6*c.Ledger.Supplied {
			t.Errorf("%s: ledger residual %g", c.Name, c.Ledger.Residual)
		}
	}
}

// TestPulser: only the zero-crossing-fired output gap is efficient.
func TestPulser(t *testing.T) {
	for k, d := range PulserDesigns {
		r, err := RunPulser(d, 6)
		if err != nil {
			t.Fatal(err)
		}
		eff := r.LoadPower / r.InPower
		if k < 2 && eff > 0.05 {
			t.Errorf("%s: efficiency %g, want it to fail", d.Name, eff)
		}
		if k == 2 && (eff < 0.4 || r.BackEnergy > 0.1 || r.Pulses != 12) {
			t.Errorf("%s: efficiency %g, backward %g J, %d pulses", d.Name, eff, r.BackEnergy, r.Pulses)
		}
		if math.Abs(r.Ledger.Residual) > 1e-6*r.Ledger.Supplied {
			t.Errorf("%s: ledger residual %g", d.Name, r.Ledger.Residual)
		}
	}
}

// TestDoubler: one cell per position sparks over; two hold.
func TestDoubler(t *testing.T) {
	d1, err := RunDoubler(1)
	if err != nil {
		t.Fatal(err)
	}
	d2, err := RunDoubler(2)
	if err != nil {
		t.Fatal(err)
	}
	if d1.Efficiency > 0.3 || d2.Efficiency < 0.8 || d2.Vdc < 170 {
		t.Errorf("1 cell %g, 2 cells %g at %g V", d1.Efficiency, d2.Efficiency, d2.Vdc)
	}
}

// TestVaractor: kronspice's resonance matches 1/(2π√LC) and scales as
// √V(form).
func TestVaractor(t *testing.T) {
	lo, err := RunVaractor(10, 5)
	if err != nil {
		t.Fatal(err)
	}
	hi, err := RunVaractor(140, 70)
	if err != nil {
		t.Fatal(err)
	}
	for _, v := range []Varactor{lo, hi} {
		if rel(v.F0, v.F0Formula) > 2e-3 {
			t.Errorf("formed at %g V: %g Hz, formula %g", v.Vform, v.F0, v.F0Formula)
		}
	}
	if rel(hi.F0/lo.F0, math.Sqrt(14)) > 3e-3 {
		t.Errorf("ratio %g, want √14", hi.F0/lo.F0)
	}
	if f, _ := RunVaractor(140, -5); f.F0 != 0 {
		t.Errorf("forward bias still resonates at %g Hz", f.F0)
	}
}

// TestDoubling: the cell bridge's 2f content is near a perfect bridge's.
func TestDoubling(t *testing.T) {
	hs, err := RunDoubling()
	if err != nil {
		t.Fatal(err)
	}
	id := IdealDoubling()
	if rel(hs[0].H2, id.H2) > 0.05 || rel(id.H2Fraction, 0.951) > 1e-3 {
		t.Errorf("cell bridge 2f %g (ideal %g); ideal share %g", hs[0].H2, id.H2, id.H2Fraction)
	}
}
