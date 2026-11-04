// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// The torque-driven stacked machine resonates at f_em = √(2·M_s·G/(L_s·J))/2π,
// and the circuit reproduces the RK4 amplitude there.
func TestElectromechanicalResonance(t *testing.T) {
	m := stackedLoad(50e-6, 5e-3)
	fem, q := m.EMResonance()
	if q < 10 {
		t.Fatalf("Q %g", q)
	}
	amp := func(x float64) float64 {
		d := m
		d.TorqueAmp, d.TorqueFreq = 0.5, x*fem
		tr := d.RunRK4(60, 2e-3)
		return cmplx.Abs(amplitude(tr.T, tr.W1, 40, x*fem))
	}
	a0, lo, hi := amp(1), amp(0.9), amp(1.1)
	if !(a0 > 3*lo && a0 > 3*hi) {
		t.Errorf("no peak at f_em: %g at 0.9, %g at 1, %g at 1.1", lo, a0, hi)
	}
	d := m
	d.TorqueAmp, d.TorqueFreq = 0.5, fem
	tr, err := d.Run(60, 2e-3)
	if err != nil {
		t.Fatal(err)
	}
	if c := cmplx.Abs(amplitude(tr.T, tr.W1, 40, fem)); math.Abs(c/a0-1) > 0.01 {
		t.Errorf("circuit amplitude %g, RK4 %g", c, a0)
	}
	if math.Abs(tr.Mismatch) > 1e-6*tr.In {
		t.Errorf("ledger mismatch %g of %g J", tr.Mismatch, tr.In)
	}
}

// A capacitor across the load keeps DC in the load when small and locks the
// machine into AC when large; the circuit (load and C as their own elements,
// with the ESR) equals the RK4 equations and the ledger closes.
func TestParallelTankLockIn(t *testing.T) {
	m := stackedLoad(50e-6, 5e-3)
	const esr = 0.29e-3
	if lockIn(m, 300, esr) {
		t.Errorf("300 F across the load should keep DC")
	}
	if !lockIn(m, 3000, esr) {
		t.Errorf("3000 F across the load should lock into AC")
	}
	// Small C: the load still carries the DC loop current.
	p, err := tankCase(m, 300, esr).loadPower(30, 2e-3, 20)
	if err != nil {
		t.Fatal(err)
	}
	if want := 50e-6 * m.I0[0] * m.I0[0]; math.Abs(p/want-1) > 0.01 {
		t.Errorf("load power with 300 F: %g W, DC %g W", p, want)
	}
	q := tankCase(m, 3000, esr)
	ct, err := q.Run(20, 2e-3)
	if err != nil {
		t.Fatal(err)
	}
	ref := q.RunRK4(20, 2e-3)
	for k := 1; k < len(ct.T); k++ {
		j := int(math.Round(ct.T[k] / 2e-3))
		if d := math.Abs(ct.I1[k] - ref.I1[j]); d > 0.01*m.I0[0] {
			t.Fatalf("t = %g: circuit I1 %g, RK4 %g", ct.T[k], ct.I1[k], ref.I1[j])
		}
	}
	if math.Abs(ct.Mismatch) > 1e-6*ct.In {
		t.Errorf("ledger mismatch %g of %g J", ct.Mismatch, ct.In)
	}
}

// A smaller capacitor is not rescued by a faster or larger speed swing; and
// a 7.8 Hz shaft oscillation is a forced response, δI = M_s·I·δω/(2πf·L_s).
func TestFastPumpAndSchumann(t *testing.T) {
	if g, _ := pumpBest(100, 0.5, 0.87); g >= 0 {
		t.Errorf("100 F across the load grew at %g 1/s", g)
	}
	m := stackedLoad(50e-6, 5e-3)
	l11, l12, ks, _ := m.loop()
	ls, ms := l11+l12, (1+ks)*m.Coupling()
	q := m
	q.TorqueAmp, q.TorqueFreq = 274, 7.8
	tr := q.RunRK4(10, 2e-4)
	aw := cmplx.Abs(amplitude(tr.T, tr.W1, 5, 7.8))
	ai := cmplx.Abs(amplitude(tr.T, tr.I1, 5, 7.8))
	if want := ms * m.I0[0] * aw / (2 * math.Pi * 7.8 * ls); math.Abs(ai/want-1) > 0.05 {
		t.Errorf("δI %g A, forced estimate %g A", ai, want)
	}
	if (m.W0[0]+aw)*m.Disk.Ro > burst()*1.01 {
		t.Errorf("rim %g m/s past the limit", (m.W0[0]+aw)*m.Disk.Ro)
	}
}

// rpm = 1768·(1 + h·sin 2π·7.8t): the air machine decays at the mean speed's
// rate whatever h; at the iron machine's onset the 7.8 Hz ripple is the exact
// 2·I1(β)/I0(β).
func TestRPMSchumann(t *testing.T) {
	w0 := rpmBase * 2 * math.Pi / 60
	a := rpmMachine(1)
	l11, l12, ks, r := a.loop()
	want := (-r + (1+ks)*a.Coupling()*w0) / (l11 + l12)
	for _, h := range []float64{0, 1} {
		if g, _, _, _ := ripple(rpmRun(a, h, 0, 3)); math.Abs(g-want) > 0.01*math.Abs(want) {
			t.Errorf("air, h = %g: growth %g, mean-speed rate %g", h, g, want)
		}
	}
	k := r / ((1 + ks) * a.Coupling() * w0)
	b := rpmMachine(k)
	l11, l12, _, r = b.loop()
	beta := r / (l11 + l12) / (2 * math.Pi * fSchu)
	g, h1, _, _ := ripple(rpmRun(b, 1, 0, 6))
	if math.Abs(g) > 5e-3 || math.Abs(h1/(2*besselI(1, beta)/besselI(0, beta))-1) > 1e-3 {
		t.Errorf("iron at onset: growth %g, ripple %g vs exact %g", g, h1, 2*besselI(1, beta)/besselI(0, beta))
	}
}
