// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Oxidation heats from the enthalpies of formation, and the vaporization
// cost from its parts.
func TestChemistry(t *testing.T) {
	si := byName("Si")
	if v := si.OxidationHeat() / 1e3; math.Abs(v-32.43) > 0.01 {
		t.Errorf("Si burns %.3f kJ/g, want 32.43 (910.7 kJ/mol ÷ 28.086)", v)
	}
	// Si: 25·1389 + 50210 + 27.2·1851 + 359000 J/mol
	want := (25*(1687-298) + 50210 + 27.2*(3538-1687) + 359000) / 28.086
	if math.Abs(si.VaporizationCost()-want) > 1e-9*want {
		t.Errorf("Si vaporization %.1f J/g, want %.1f", si.VaporizationCost(), want)
	}
	if fe := byName("Fe"); fe.Return() > 1 || si.Return() < 1.5 {
		t.Errorf("returns: Fe %.2f (want < 1), Si %.2f (want > 1.5)", fe.Return(), si.Return())
	}
}

// Deal–Grove dry-O2 rate constants at 1200 °C (1473 K), 1 atm: B ≈ 0.045
// µm²/h (textbook value, checked to 5 %).
func TestDealGrove(t *testing.T) {
	b := DealGroveB(1473) * 3600 * 1e12
	if math.Abs(b/0.0477-1) > 0.05 {
		t.Errorf("B(1473 K) = %.4f µm²/h", b)
	}
	// burn time falls with temperature and rises with radius
	if !(BurnTime(10e-9, 1500, 0.21) < BurnTime(10e-9, 1300, 0.21) && BurnTime(10e-9, 1300, 0.21) < BurnTime(50e-9, 1300, 0.21)) {
		t.Error("burn time not monotonic")
	}
	// seconds need roughly 2000–3000 K for 10 nm particles
	if T := tempFor(10e-9, 3.5, 0.21); T < 2000 || T > 3200 {
		t.Errorf("10 nm in 3.5 s at %.0f K", T)
	}
}

// kronspice welder: 120 A at the arc voltage, about 2.5 kW, ledger closed.
func TestWeld(t *testing.T) {
	w, err := TuneWeld(70, 22.5, 60, 120)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(w.Irms-120) > 0.5 || w.Vrms < 20 || w.Vrms > 25 || w.Power < 2000 || w.Power > 3000 || math.Abs(w.LedgerMismatch) > 1e-3 {
		t.Errorf("weld: %.1f A, %.1f V, %.0f W, ledger %.1e", w.Irms, w.Vrms, w.Power, w.LedgerMismatch)
	}
}

// kronspice power line: power pulses at twice the line frequency, and the
// most the line can drive into a 0.15 m ball 20 m away is far below 100 W;
// the AC maximum equals the matched-load formula ¼·ω·Cm²·V²/(Cm + Cb).
func TestLine(t *testing.T) {
	l, err := RunLine(35e3/math.Sqrt(3), 50, 0.15, 20)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(l.RippleHz-100) > 0.5 || l.BestPower > 1e-3 {
		t.Errorf("ripple %.2f Hz, power %.3g W", l.RippleHz, l.BestPower)
	}
	w := 2 * math.Pi * 50
	ct := l.Cm + l.Cb
	pmax := 0.25 * w * l.Cm * l.Cm * l.Vphase * l.Vphase * 2 / ct
	if math.Abs(l.BestPower/pmax-1) > 0.02 {
		t.Errorf("best power %.4g W, matched-load formula %.4g W", l.BestPower, pmax)
	}
}
