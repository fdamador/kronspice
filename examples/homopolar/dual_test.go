// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// With both shafts held at one speed the field's symmetric mode grows at
// exactly (−R + (1+K)·M·ω)/(l11 + l12), side by side and stacked.
func TestDualThresholdGrowth(t *testing.T) {
	for _, stacked := range []bool{false, true} {
		for _, x := range []float64{0.8, 1.2} {
			m := DefaultDual()
			m.Stacked = stacked
			w := x * m.Threshold()
			m.Speed = &[2]float64{w, w}
			m.I0 = [2]float64{100, 100}
			tr, err := m.Run(0.05, 1e-5)
			if err != nil {
				t.Fatal(err)
			}
			l11, l12, ks, r := m.loop()
			want := (-r + (1+ks)*m.Coupling()*w) / (l11 + l12)
			if g := growth(tr); math.Abs(g-want) > 1e-3*math.Abs(want)+1e-3 {
				t.Errorf("stacked %v, %.1f ω*: growth %g, want %g", stacked, x, g, want)
			}
		}
	}
}

// The engine and the RK4 equations agree, with feedback coils and stacked,
// and the ledger closes.
func TestDualCircuitMatchesRK4(t *testing.T) {
	for _, m := range []Dual{
		func() Dual { m := DefaultDual(); m.K = 0.5; return m }(),
		func() Dual { m := DefaultDual(); m.K = -0.25; return m }(),
		func() Dual { m := DefaultDual(); m.Stacked = true; return m }(),
	} {
		tr, err := m.Run(0.4, 1e-4)
		if err != nil {
			t.Fatal(err)
		}
		ref := m.RunRK4(0.4, 1e-5)
		scale := 0.0
		for _, v := range ref.I1 {
			scale = math.Max(scale, math.Abs(v))
		}
		for k := 1; k < len(tr.T); k++ { // k = 0 is the imposed UIC point
			j := int(math.Round(tr.T[k] / 1e-5))
			if d := math.Abs(tr.I1[k] - ref.I1[j]); d > 2e-3*scale {
				t.Fatalf("K %g stacked %v: t = %g, I1 %g vs %g", m.K, m.Stacked, tr.T[k], tr.I1[k], ref.I1[j])
			}
			if d := math.Abs(tr.W1[k] - ref.W1[j]); d > 2e-3*math.Abs(ref.W1[j]) {
				t.Fatalf("K %g stacked %v: t = %g, ω1 %g vs %g", m.K, m.Stacked, tr.T[k], tr.W1[k], ref.W1[j])
			}
		}
		if math.Abs(tr.Mismatch) > 1e-6*tr.In {
			t.Errorf("K %g stacked %v: ledger mismatch %g of %g J", m.K, m.Stacked, tr.Mismatch, tr.In)
		}
	}
}

// Stacking lowers the threshold below the copper burst speed; side by side
// it lies above it.
func TestStackedBelowBurst(t *testing.T) {
	side, st := DefaultDual(), DefaultDual()
	st.Stacked = true
	if !(st.Threshold()*st.Disk.Ro < burst() && side.Threshold()*side.Disk.Ro > burst()) {
		t.Errorf("rim at threshold: side %g, stacked %g, burst %g m/s", side.Threshold()*side.Disk.Ro, st.Threshold()*st.Disk.Ro, burst())
	}
	if k := st.MFar() / st.Coupling(); k <= 0 || k >= 1 {
		t.Errorf("stacked self-coupling %g", k)
	}
}
