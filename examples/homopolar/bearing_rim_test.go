// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

var rimW0 = rpmBase * 2 * math.Pi / 60

// The Hertz scaling returns the measured point on the measured bearing.
func TestRimHertzArea(t *testing.T) {
	r := &RimBearing{Races: 1, Bore: rimDRef, Preload: rimFRef}
	if a := r.HertzArea(); math.Abs(a/rimAHzRef-1) > 1e-12 {
		t.Errorf("area %g m², measured %g", a, rimAHzRef)
	}
	if i := defaultRim(rimW0).ILimit(); !(i > 1 && i < 100) {
		t.Errorf("rim bank limit %g A", i)
	}
}

// The copper machine with rim bearings: the engine (P elements) equals RK4
// either side of the equilibrium current, the current falls below it and
// rises above it, and the ledger closes.
func TestRimCopperTwin(t *testing.T) {
	m, _ := outputMachine(10e-6, 3, 10, 0.98)
	m.Rim = defaultRim(rimW0)
	w := 2.5 * rimW0
	m.Speed = m.held(w)
	is := m.equilibrium(w)
	if !(is > 0) {
		t.Fatalf("no equilibrium at 2.5× onset: %g", is)
	}
	for _, f := range []float64{0.9, 1.1} {
		q := m
		q.I0 = [2]float64{f * is, f * is}
		ct, err := q.Run(0.5, 2e-4)
		if err != nil {
			t.Fatal(err)
		}
		rk := q.RunRK4(0.5, 1e-4)
		n := len(ct.T) - 1
		a, b := ct.I1[n]/ct.I1[0], rk.I1[len(rk.I1)-1]/rk.I1[0]
		if math.Abs(a/b-1) > 1e-3 {
			t.Errorf("from %g I*: engine %g, RK4 %g", f, a, b)
		}
		if (f < 1) != (a < 1) {
			t.Errorf("from %g I*: ratio %g on the wrong side", f, a)
		}
		if math.Abs(ct.Mismatch) > 1e-6*ct.In {
			t.Errorf("ledger mismatch %g of %g J", ct.Mismatch, ct.In)
		}
	}
}

// The bearing-matched machine: its equilibrium is the set current and is
// unstable at fixed feedback; the PI regulator starts it from a pulse and
// holds the set current at 1.6 T; the books close.
func TestRimRegulated(t *testing.T) {
	rim := defaultRim(rimW0)
	iSet := 0.8 * rim.ILimit()
	s := superMachine(rim, iSet, rimW0, 4.45)
	if math.Abs(s.equilibrium(rimW0)/iSet-1) > 1e-9 {
		t.Fatalf("equilibrium %g, set %g", s.equilibrium(rimW0), iSet)
	}
	for _, f := range []float64{0.95, 1.05} {
		q := s
		q.I0 = [2]float64{f * iSet, f * iSet}
		tr := q.RunRK4(300, 0.05)
		if end := tr.I1[len(tr.I1)-1] / iSet; (f < 1) != (end < f) {
			t.Errorf("from %g I*: %g I* after 300 s, should move away", f, end)
		}
	}
	q := s
	q.Reg = &Regulator{Target: iSet, Kmin: -0.5, Kmax: 1, Kp: 2, Tau: 100}
	q.I0 = [2]float64{7.5, 7.5}
	tr := q.RunRK4(3000, 0.05)
	n := len(tr.T) - 1
	if math.Abs(tr.I1[n]/iSet-1) > 1e-3 {
		t.Errorf("regulated current %g A, set %g A", tr.I1[n], iSet)
	}
	if b := s.fieldB(tr.I1[n]); math.Abs(b/bSatIron-1) > 0.01 {
		t.Errorf("field %g T", b)
	}
	gen, _, _, _, _, res := q.energy(tr, 0)
	if math.Abs(res) > 1e-6*gen {
		t.Errorf("books: residual %g of %g J", res, gen)
	}
}

// Huber torque: EM mode leaves the regulated drive torque unchanged (its
// shaft help is paid for by the back-EMF), thermal mode returns at most
// η of the bearing heat; the engine's ledger closes with the P elements.
func TestRimHuberNoGain(t *testing.T) {
	rim := defaultRim(rimW0)
	iSet := 0.8 * rim.ILimit()
	drive := func(kp float64, thermal bool) (td, th, i float64) {
		r := *rim
		r.KappaPsi, r.Thermal = kp, thermal
		p := superMachine(&r, iSet, rimW0, 4.45)
		p.Reg = &Regulator{Target: iSet, Kmin: -0.5, Kmax: 1, Kp: 2, Tau: 100}
		p.I0 = [2]float64{iSet, iSet}
		tr := p.RunRK4(3000, 0.05)
		j := len(tr.T) - 1
		i = tr.I1[j]
		_, th, _ = r.loop(i, rimW0, -rimW0)
		_, _, ks, _ := p.loop()
		return (1+ks-p.FeedbackK+tr.K[j])*p.Coupling()*i*i - th, th, i
	}
	t0, _, _ := drive(0, false)
	t1, h1, i1 := drive(1e-3, false)
	if h1 <= 0 || math.Abs(i1/iSet-1) > 1e-3 || math.Abs(t1/t0-1) > 1e-3 {
		t.Errorf("EM: drive torque %g against %g without Huber (Huber %g, I %g)", t1, t0, h1, i1)
	}
	_, h2, i2 := drive(1e-2, true)
	if heat := rim.VF*i2 + rim.R/float64(rim.Races)*i2*i2; h2*rimW0 > rimEta*heat*(1+1e-9) {
		t.Errorf("thermal: Huber power %g W above η·H = %g W", h2*rimW0, rimEta*heat)
	}
	r := *rim
	r.KappaPsi = 1e-3
	p := superMachine(&r, iSet, rimW0, 4.45)
	p.I0 = [2]float64{iSet, iSet}
	ct, err := p.Run(50, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ct.Mismatch) > 1e-6*ct.In {
		t.Errorf("engine ledger mismatch %g of %g J", ct.Mismatch, ct.In)
	}
}
