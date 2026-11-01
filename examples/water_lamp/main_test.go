// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Measured anchors: the water-vapour Paschen curve has its 480 V minimum
// at 0.6 Torr·cm (Škoro 2011); air at 1 atm, 1 mm breaks down near 4.5–5 kV.
func TestPaschen(t *testing.T) {
	if v := PaschenWater(0.6); math.Abs(v-480) > 1e-9 {
		t.Errorf("water vapour minimum %.3f V", v)
	}
	if PaschenWater(0.5) <= 480 || PaschenWater(0.8) <= 480 {
		t.Error("0.6 Torr·cm is not the minimum")
	}
	if v := PaschenAir(760, 0.1); v < 4.2e3 || v > 5.5e3 {
		t.Errorf("air 1 mm: %.0f V", v)
	}
}

// A primed lamp with no statistical wait follows the textbook relaxation
// period T = R·C·ln((Vs − Vm)/(Vs − Vb)) (to the lamp's short on-time), and
// every ledger closes.
func TestRelaxationPeriod(t *testing.T) {
	c := Cases()[0]
	c.Lamp.Seeds = 0
	r, err := c.Run(2, 1)
	if err != nil {
		t.Fatal(err)
	}
	want := c.R * c.C * math.Log((c.Vs-c.Lamp.Vm)/(c.Vs-c.Lamp.Vb))
	if math.Abs(r.MeanInterval/want-1) > 0.02 {
		t.Errorf("period %.4g s, textbook %.4g s", r.MeanInterval, want)
	}
	for _, c := range Cases() {
		sec := math.Min(c.Seconds, 2)
		r, err := c.Run(sec, 7)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.Residual) > 1e-6 {
			t.Errorf("%s: ledger residual %.2g", c.Name, r.Residual)
		}
	}
}

// Without a starter, in the dark, the mean wait matches 1/(S·(1 − Vb/Vs)).
func TestDarkWait(t *testing.T) {
	c := Cases()[1]
	r, err := c.Run(36000, 7)
	if err != nil {
		t.Fatal(err)
	}
	want := 1 / (c.Lamp.Seeds * (1 - c.Lamp.Vb/c.Vs))
	if r.Flashes < 30 || math.Abs(r.MeanInterval/want-1) > 0.3 {
		t.Errorf("%d flashes, mean wait %.3g s, expected %.3g s", r.Flashes, r.MeanInterval, want)
	}
}

// The reduced coil matches kronspice, and every lamp's energy into the
// coil is accounted for.
func TestCoil(t *testing.T) {
	coil := DefaultCoil()
	if tw, err := coil.Twin(); err != nil || tw > 0.01 {
		t.Fatalf("coil twin %.3g (%v)", tw, err)
	}
	c := Cases()[6] // mist gap
	r, err := c.RunCoil(coil, 0.1, 7)
	if err != nil {
		t.Fatal(err)
	}
	if r.Flashes == 0 || r.CoilPeak < 10 || math.Abs(r.Residual) > 1e-4 {
		t.Errorf("%d flashes, coil peak %.3g V, ledger %.2g", r.Flashes, r.CoilPeak, r.Residual)
	}
}

// Part 3: tuning the flash circuit to the coil moves far more of each
// flash into the coil, a quenched gap more than an ordinary arc (mist gap),
// and the ledgers close.
func TestTuning(t *testing.T) {
	rows, err := TunedComparison()
	if err != nil {
		t.Fatal(err)
	}
	mist := rows[3:]
	if !(mist[1].Share > 3*mist[0].Share && mist[2].Share > mist[1].Share) {
		t.Errorf("mist gap shares: untuned %.4f, arc %.4f, quenched %.4f", mist[0].Share, mist[1].Share, mist[2].Share)
	}
	for _, r := range rows {
		if math.Abs(r.Run.Residual) > 1e-3 {
			t.Errorf("%s: ledger %.2g", r.Label, r.Run.Residual)
		}
	}
}

// Part 4: parallel gaps on one capacitor do not add spark current;
// parallel capacitor plates do.
func TestParallel(t *testing.T) {
	g, err := ParallelGaps(0.249e-6, 20, 0.65)
	if err != nil {
		t.Fatal(err)
	}
	if r := g[3].PeakIp / g[0].PeakIp; math.Abs(r-1) > 0.05 {
		t.Errorf("8 gaps: spark current ×%.3f", r)
	}
	p, err := ParallelPlates(0.249e-6, 0.65)
	if err != nil {
		t.Fatal(err)
	}
	if p[3].PeakIp < 4*p[0].PeakIp || p[3].PerFlash < 7*p[0].PerFlash {
		t.Errorf("8 plates: spark ×%.2f, flash energy ×%.2f", p[3].PeakIp/p[0].PeakIp, p[3].PerFlash/p[0].PerFlash)
	}
}
