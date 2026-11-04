// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// tankRun is the Q = 3 output design in the given mode, both shafts held at
// 1768 rpm with the full 7.8 Hz sine, and an RK4 step that resolves the
// transformer's leakage time.
func tankRun(mode OutputMode, i0 [2]float64) (Dual, float64) {
	m, _ := outputMachine(10e-6, 3, 10, 0.98)
	m.Output = mode
	w0 := rpmBase * 2 * math.Pi / 60
	m.Speed = m.held(w0)
	m.SpeedAmp, m.SpeedFreq = 1, fSchu
	m.I0 = i0
	leak := (1 - m.KT*m.KT) * m.LSec / m.RAC // the same in every mode (L_s and R_ac scale together)
	return m, math.Min(math.Min(2e-4, m.RLoad*m.CParallel/5), leak/5)
}

// checkTwin compares the circuit with the RK4 twin and checks the ledger.
func checkTwin(t *testing.T, m Dual, h float64) {
	t.Helper()
	ct, err := m.Run(1, 2e-4)
	if err != nil {
		t.Fatal(err)
	}
	ref := m.RunRK4(1, h)
	for k := 1; k < len(ct.T); k += 50 {
		j := int(math.Round(ct.T[k] / h))
		if d := math.Max(math.Abs(ct.I1[k]-ref.I1[j]), math.Abs(ct.I2[k]-ref.I2[j])); d > 0.01*10e3 {
			t.Fatalf("t = %g: circuit I = %g, %g; RK4 %g, %g", ct.T[k], ct.I1[k], ct.I2[k], ref.I1[j], ref.I2[j])
		}
	}
	if math.Abs(ct.Mismatch) > 1e-5*ct.In {
		t.Errorf("ledger mismatch %g of %g J", ct.Mismatch, ct.In)
	}
}

// One tank between the coils, shafts common, shaft 2 counter-rotating, the
// secondary's common end on the shafts: circuit = RK4, ledger closed; and
// the single current is the two per-loop tanks' symmetric mode exactly.
func TestSeriesTank(t *testing.T) {
	m, h := tankRun(OutSeries, [2]float64{10e3, 10e3})
	checkTwin(t, m, h)
	mp, _ := tankRun(OutPerLoop, [2]float64{10e3, 10e3})
	a, b := m.RunRK4(2, h), mp.RunRK4(2, h)
	for k := range a.T {
		if d := math.Abs(a.I1[k] - b.I1[k]); d > 1e-9*10e3 {
			t.Fatalf("t = %g: series %g, per-loop %g", a.T[k], a.I1[k], b.I1[k])
		}
		if a.I2[k] != -a.I1[k] {
			t.Fatalf("t = %g: one current, but I2 = %g against I1 = %g", a.T[k], a.I2[k], a.I1[k])
		}
	}
	// Co-rotating, the disks' EMFs cancel round the loop and the current
	// decays as a plain R–L loop: the DC passes the coils and the primary
	// (R_L and the secondary barely load it at this slow rate), so
	// g ≈ −2·(R_int + R_p)/(2·(L + M_c) + 2·L_p).
	q := m
	w0 := rpmBase * 2 * math.Pi / 60
	q.Speed, q.SpeedAmp = &[2]float64{w0, w0}, 0
	l11, l12, _, r := q.loop()
	want := -2 * (r - q.RLoad + q.RPrim) / (2*(l11+l12) + 2*q.LPrim)
	tr := q.RunRK4(2, h)
	if g := rmsGrowth(tr.T, tr.I1, 0.5); math.Abs(g/want-1) > 0.05 {
		t.Errorf("co-rotating series loop: growth %g 1/s, plain R–L decay %g", g, want)
	}
}

// Both machines across one tank in parallel: circuit = RK4 with unequal
// loop currents; at equal currents it is the two per-loop tanks in parallel;
// the difference current circulates rim to rim through the shafts, never
// sees the tank, and decays at (−R_int − M·(1−K_S)·ω)/(L − M_c).
func TestSharedTank(t *testing.T) {
	m, h := tankRun(OutShared, [2]float64{10e3, 6e3})
	checkTwin(t, m, h)
	ms, _ := tankRun(OutShared, [2]float64{10e3, 10e3})
	mp, _ := tankRun(OutPerLoop, [2]float64{10e3, 10e3})
	a, b := ms.RunRK4(2, h), mp.RunRK4(2, h)
	for k := range a.T {
		if d := math.Abs(a.I1[k] - b.I1[k]); d > 1e-9*10e3 {
			t.Fatalf("t = %g: shared %g, per-loop %g", a.T[k], a.I1[k], b.I1[k])
		}
	}
	m.SpeedAmp = 0
	tr := m.RunRK4(0.5, h)
	l11, l12, ks, r := m.loop()
	want := (-(r - m.RLoad) - (1-ks)*m.Coupling()*m.Speed[0]) / (l11 - l12)
	n := len(tr.T) - 1
	if g := math.Log((tr.I1[n]-tr.I2[n])/(tr.I1[0]-tr.I2[0])) / tr.T[n]; math.Abs(g/want-1) > 1e-3 {
		t.Errorf("difference mode decays at %g 1/s, want %g", g, want)
	}
}

// The rectifier in the series design: the bank charges and the engine's
// ledger closes with the secondary's common on the shafts.
func TestSeriesRectifier(t *testing.T) {
	if testing.Short() {
		t.Skip("rectifier transient")
	}
	m, _ := tankRun(OutSeries, [2]float64{10e3, 10e3})
	z, _ := m.TankScale()
	m.CBank, m.RDC, m.BankV0 = 0.05, 20, 200*z
	c, err := m.Circuit()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := analysisTran(c, 1, 5e-5)
	if err != nil {
		t.Fatal(err)
	}
	vp, _ := wf.NodeVoltage("dp")
	vn, _ := wf.NodeVoltage("dn")
	if v := vp[len(vp)-1] - vn[len(vn)-1]; v < 100 {
		t.Errorf("bank at %g V", v)
	}
	if r := wf.Ledger.IntegrationMismatch / wf.Ledger.EnergyIn; math.Abs(r) > 1e-6 {
		t.Errorf("ledger mismatch %g of the drive work", r)
	}
}
