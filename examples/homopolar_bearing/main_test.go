// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The Hertz scaling reproduces the paper's second measured-calculated point
// (2.5 mm² at 1600 N) from the first (1.0 mm² at 400 N).
func TestHertzArea(t *testing.T) {
	d := Defaults()
	d.Bore, d.Preload = DRef, 1600
	if a := d.HertzArea() * 1e6; math.Abs(a-2.5)/2.5 > 0.02 {
		t.Errorf("A_Hz(1600 N) = %.3f mm², paper 2.5", a)
	}
}

// The ledger closes, and the drive power is never below the load power
// plus the losses fixed by measured data: bearing contact drops and
// friction. In EM mode the Huber torque cannot lower that floor; in thermal
// mode it lowers the drive power by at most η × the bearing heat.
func TestNoFreeTorque(t *testing.T) {
	for _, kp := range []float64{0, 1e-3, 1e-2} {
		for _, th := range []bool{false, true} {
			d := Defaults()
			d.KappaPsi, d.Thermal = kp, th
			r, err := d.Solve()
			if err != nil {
				t.Fatal(err)
			}
			if r.Residual > 1e-9 {
				t.Errorf("κψ=%g thermal=%v: residual %g", kp, th, r.Residual)
			}
			floor := r.PLoad + 2*d.FrictionTorque()*d.Omega
			if !th {
				floor += 2 * d.VBearing * math.Abs(r.ILoad) * 0.999
			} else {
				floor += 2 * (1 - d.Eta) * d.VBearing * math.Abs(r.ILoad) * 0.999
			}
			if r.PShaft < floor {
				t.Errorf("κψ=%g thermal=%v: shaft %.4g W below load + measured losses %.4g W", kp, th, r.PShaft, floor)
			}
			if r.Efficiency >= 1 {
				t.Errorf("efficiency %g", r.Efficiency)
			}
		}
	}
	base, _ := Defaults().Solve()
	th := Defaults()
	th.KappaPsi, th.Thermal = 1e-1, true
	r, _ := th.Solve()
	saved := base.PShaft - r.PShaft
	if saved < 0 || saved > th.Eta*base.PBearingHeat*1.0001 {
		t.Errorf("thermal Huber saved %.4g W of drive; ceiling η·H = %.4g W", saved, th.Eta*base.PBearingHeat)
	}
}

// Bigger bearings: less magnetic area, more current capacity, more friction.
func TestBoreTradeOff(t *testing.T) {
	var prevK, prevI, prevF = math.Inf(1), 0.0, 0.0
	for _, bore := range []float64{0.01, 0.04, 0.08, 0.12} {
		d := Defaults()
		d.Bore = bore
		if !(d.K() < prevK && d.CurrentLimit(JLimit) > prevI && d.FrictionTorque() > prevF) {
			t.Errorf("bore %g: K %g, I_max %g, friction %g not monotone", bore, d.K(), d.CurrentLimit(JLimit), d.FrictionTorque())
		}
		prevK, prevI, prevF = d.K(), d.CurrentLimit(JLimit), d.FrictionTorque()
	}
}

// Charging the capacitor: the transient ledger closes and the electrical
// conversion stays below the ceiling the measured bearing drops set.
func TestCapacitorCharge(t *testing.T) {
	d := Defaults()
	vc := 2*d.K()*d.Omega - 2*d.VBearing
	ceiling := vc / (vc + 4*d.VBearing)
	for _, ramp := range []float64{0, 40} {
		tEnd := ramp + 5*d.RL
		st, sh, mm, err := d.Charge(1, ramp, tEnd)
		if err != nil {
			t.Fatal(err)
		}
		conv := st / (sh - 2*d.FrictionTorque()*d.Omega*tEnd)
		if conv >= ceiling || conv <= 0 {
			t.Errorf("ramp %g: conversion %.3f, ceiling %.3f", ramp, conv, ceiling)
		}
		if math.Abs(mm) > 1e-4*sh {
			t.Errorf("ramp %g: integration mismatch %g J of %g J", ramp, mm, sh)
		}
	}
}
