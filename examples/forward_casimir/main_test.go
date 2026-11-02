// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestCasimirLimit: the Lifshitz energy with a perfect mirror is Forward's
// eq. 1; with aluminium it is smaller and rises toward it at large gaps.
func TestCasimirLimit(t *testing.T) {
	for _, a := range []float64{10e-9, 100e-9, 1e-6} {
		if r := LifshitzEnergy(a, 0) / PerfectEnergy(a); math.Abs(r-1) > 1e-4 {
			t.Errorf("perfect mirror at %g: ratio %g", a, r)
		}
	}
	prev := 0.0
	for _, a := range []float64{10e-9, 30e-9, 100e-9, 300e-9, 1e-6} {
		eta := LifshitzEnergy(a, 107e-9) / PerfectEnergy(a)
		if eta <= prev || eta >= 1 {
			t.Errorf("aluminium η at %g = %g (previous %g)", a, eta, prev)
		}
		prev = eta
	}
	tab := NewTable(Perfect, 5e-9, 10e-6, 400)
	if r := tab.Pressure(50e-9) / PerfectPressure(50e-9); math.Abs(r-1) > 1e-3 {
		t.Errorf("table pressure ratio %g", r)
	}
}

// TestStroke: closing, the supply ends up with the Casimir work less small
// losses, the ledger closes, and the kronspice twin agrees.
func TestStroke(t *testing.T) {
	l := DefaultLeaf(Aluminium)
	c := l.Run(l.A0, l.AMin)
	if math.Abs(c.Residual) > 0.02*c.Casimir || c.Net < 0.95*c.Casimir || c.Net > c.Casimir {
		t.Errorf("Casimir %g, net %g, residual %g", c.Casimir, c.Net, c.Residual)
	}
	tw, err := KronTwin(l, c)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(tw.ToSupply/c.ToSupply-1) > 2e-3 || math.Abs(-tw.Charge/c.ChargeIn-1) > 2e-3 || math.Abs(tw.PumpWork/c.Casimir-1) > 2e-3 {
		t.Errorf("twin: supply %g vs %g, charge %g vs %g, work %g vs %g", tw.ToSupply, c.ToSupply, -tw.Charge, c.ChargeIn, tw.PumpWork, c.Casimir)
	}
}

// TestCycle: closing then reopening returns nothing net (a small loss).
func TestCycle(t *testing.T) {
	l := DefaultLeaf(Aluminium)
	c := l.Run(l.A0, l.AMin)
	o := l.Run(l.AMin, l.A0)
	if net := c.Net + o.Net; net > 0 || net < -0.1*c.Casimir {
		t.Errorf("cycle net %g (Casimir %g)", net, c.Casimir)
	}
	if math.Abs(c.Casimir+o.Casimir) > 1e-3*c.Casimir {
		t.Errorf("Casimir work not conservative: %g, %g", c.Casimir, o.Casimir)
	}
}
