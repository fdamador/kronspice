// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// A valve metal against an ordinary one rectifies; two ordinary metals do
// not; two valve metals block both ways; the curves match kronspice's DC
// sweep of the same element.
func TestPairs(t *testing.T) {
	al := DefaultCell(Pb, Al)
	if r := al.Current(100) / -al.Current(-100); r < 1e3 {
		t.Errorf("Pb|Al forward/reverse at 100 V only %.3g", r)
	}
	if fc := DefaultCell(Fe, Cu); math.Abs(fc.Current(100)+fc.Current(-100)) > 1e-12 {
		t.Error("Fe|Cu is not symmetric")
	}
	if aa := DefaultCell(Al, Al); math.Abs(aa.Current(100)) > 1 || math.Abs(aa.Current(-100)) > 1 {
		t.Error("Al|Al does not block")
	}
	if al.Current(-200) > -1 {
		t.Error("no sparking past the oxide's breakdown")
	}
	for _, c := range []Cell{al, DefaultCell(Fe, Cu), DefaultCell(Al, Al)} {
		if tw, err := IVTwin(c, 200); err != nil || tw > 1e-9 {
			t.Errorf("%s: kronspice twin %.2g (%v)", c.Name(), tw, err)
		}
	}
	if g := DefaultCell(Fe, Cu).Galvanic(); math.Abs(g-0.789) > 1e-9 {
		t.Errorf("Fe|Cu galvanic %.4f V", g)
	}
}

// The bridge of aluminium–lead cells rectifies at the efficiency the 1911
// account gives (about 75–80 %; within a few points), one cell sparks over
// in a filtered half-wave rectifier, two hold; ledgers close.
func TestCellRectifiers(t *testing.T) {
	c := DefaultCell(Pb, Al)
	base := Rectifier{Vp: 100, F: 60, Rs: 0.5, Cload: 2200e-6, Rload: 50, Cycles: 30, Skip: 20}
	br := base
	br.Bridge, br.Arm = true, Arm{c, 1}
	r, err := br.Run()
	if err != nil {
		t.Fatal(err)
	}
	if r.Efficiency < 0.7 || r.Efficiency > 0.9 || r.Vdc < 80 {
		t.Errorf("bridge: %.1f V at %.1f %%", r.Vdc, 100*r.Efficiency)
	}
	h1, h2 := base, base
	h1.Arm, h2.Arm = Arm{c, 1}, Arm{c, 2}
	a, err := h1.Run()
	if err != nil {
		t.Fatal(err)
	}
	b, err := h2.Run()
	if err != nil {
		t.Fatal(err)
	}
	if !(a.Efficiency < 0.4 && b.Efficiency > 0.75) {
		t.Errorf("half-wave efficiency: 1 cell %.2f, 2 cells %.2f", a.Efficiency, b.Efficiency)
	}
	for _, x := range []RectResult{r, a, b} {
		if math.Abs(x.Mismatch) > 1e-5*x.EnergyIn {
			t.Errorf("ledger mismatch %.3g of %.3g J", x.Mismatch, x.EnergyIn)
		}
	}
}

// Spark gaps: a resistive half-wave rectifies with ratio 2; a bridge needs
// a ratio above 2 (breaks backwards at 2, rectifies at 3).
func TestGapRectifiers(t *testing.T) {
	g := DefaultGap()
	hw := Rectifier{Vp: 4000, F: 60, Rs: 1000, Rload: 1e6, Cycles: 20, Skip: 10}
	if r := hw.RunGaps(g); r.Backward != 0 || r.Forward == 0 || r.Vdc < 500 {
		t.Errorf("half-wave: %d forward, %d backward, %.0f V", r.Forward, r.Backward, r.Vdc)
	}
	br := Rectifier{Bridge: true, Vp: 6000, F: 60, Rs: 1000, Cload: 0.1e-6, Rload: 1e6, Cycles: 20, Skip: 10}
	if r := br.RunGaps(g); r.Backward == 0 {
		t.Error("ratio 2 bridge did not break down backwards")
	}
	g.NegRatio = 3
	if r := br.RunGaps(g); r.Backward != 0 || r.Vdc < 2000 {
		t.Errorf("ratio 3 bridge: %d backward, %.0f V", r.Backward, r.Vdc)
	}
}

// The gap solver against kronspice below breakdown (every gap is its leak
// resistance): same output voltage.
func TestGapSolverTwin(t *testing.T) {
	g := DefaultGap()
	r := Rectifier{Bridge: true, Vp: 1000, F: 60, Rs: 1000, Cload: 1e-9, Rload: 1e8, Cycles: 3, Skip: 1}
	mine := r.RunGaps(g)
	r.Arm = Arm{Leak{g.Rleak()}, 1}
	ks, err := r.Run()
	if err != nil {
		t.Fatal(err)
	}
	if d := math.Abs(mine.Vdc-ks.Vdc) / math.Max(math.Abs(ks.Vdc), 1e-3); d > 0.02 && math.Abs(mine.Vdc-ks.Vdc) > 0.05 {
		t.Errorf("below breakdown: own solver %.4g V, kronspice %.4g V", mine.Vdc, ks.Vdc)
	}
	if mine.Forward+mine.Backward != 0 {
		t.Error("struck below breakdown")
	}
}

// Silver counter-electrodes rectify against aluminium like lead does;
// the silver plating's life follows Faraday's law.
func TestSilver(t *testing.T) {
	for _, e := range []Electrode{CuAgTips, FeAgCoat} {
		c := DefaultCell(e, Al)
		if r := c.Current(100) / -c.Current(-100); r < 1e3 {
			t.Errorf("%s: forward/reverse %.3g", c.Name(), r)
		}
	}
	c := DefaultCell(FeAgCoat, Al)
	m, q, s := c.SilverLife(FeAgCoat, 1)
	// 10 µm × 100 cm² × 10.49 g/cm³ = 1.049 g; × F / 107.87 g/mol = 938 C
	if math.Abs(m-1.049e-3) > 1e-6 || math.Abs(q-938.3) > 1 || math.Abs(s-q) > 1e-9 {
		t.Errorf("silver plating: %.4g kg, %.4g C, %.4g s at 1 A", m, q, s)
	}
}

// TestPulses: same 10 J for each device; the ledger closes; the gaps are
// the fastest and highest.
func TestPulses(t *testing.T) {
	ps, err := Pulses()
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range ps {
		if e := 0.5 * p.C * p.V0 * p.V0; math.Abs(e-pulseJ) > 1e-9 {
			t.Errorf("%s: stored %g J, want %g", p.Name, e, pulseJ)
		}
		if math.Abs(p.Mismatch) > 1e-4 {
			t.Errorf("%s: ledger mismatch %g", p.Name, p.Mismatch)
		}
	}
	cell, water := ps[0], ps[3]
	if !(water.Peak > 5*cell.Peak && water.Rise < cell.Rise/10) {
		t.Errorf("water gap %g A / %g s should beat the cell %g A / %g s", water.Peak, water.Rise, cell.Peak, cell.Rise)
	}
	// overdamped cell: peak about (100 − 2 V)/(0.1 + 0.1 Ω)
	if math.Abs(cell.Peak-490)/490 > 0.06 {
		t.Errorf("cell peak %g A, want about 490", cell.Peak)
	}
}

// TestSweep: cells hold their amps per volt to 3 kHz and lose efficiency
// above; the water-gap bridge stops above 100 Hz.
func TestSweep(t *testing.T) {
	ds, err := Sweep()
	if err != nil {
		t.Fatal(err)
	}
	pb := ds[0].Rows
	if math.Abs(pb[4].AperV-pb[0].AperV)/pb[0].AperV > 0.01 || pb[0].AperV > IdealAperV() {
		t.Errorf("lead–aluminium bridge: %g A/V at 50 Hz, %g at 3 kHz", pb[0].AperV, pb[4].AperV)
	}
	if pb[len(pb)-1].Efficiency > 0.5*pb[0].Efficiency {
		t.Errorf("efficiency at 100 kHz %g should collapse (50 Hz: %g)", pb[len(pb)-1].Efficiency, pb[0].Efficiency)
	}
	raised := ds[4].Rows
	if raised[0].AperV < 0.1 || raised[0].Backward != 0 {
		t.Errorf("gap bridge at 50 Hz: %g A/V, %d backward strikes", raised[0].AperV, raised[0].Backward)
	}
	if raised[3].AperV > 0.01 || raised[3].Backward == 0 {
		t.Errorf("gap bridge at 1 kHz should short backwards: %g A/V, %d backward", raised[3].AperV, raised[3].Backward)
	}
}
