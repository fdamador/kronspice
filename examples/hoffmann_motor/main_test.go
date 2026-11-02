// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestEq5Twin: eq. 5′ gives the frame-I forces at random states; leaving
// out the torsion does not.
func TestEq5Twin(t *testing.T) {
	worst, wrong := twin(Default(), 200)
	if worst > 1e-12 || wrong < 0.5 {
		t.Errorf("with torsion %g, without %g", worst, wrong)
	}
}

// TestTorsionConstant: in the stationary-coil frame Ω does not depend on θ.
func TestTorsionConstant(t *testing.T) {
	q := Default()
	q.LD = 0.02
	for _, p := range []Params{Default(), q} {
		if d := torsionSpread(p); d > 1e-12 {
			t.Errorf("L_D %g: Ω varies by %g", p.LD, d)
		}
	}
}

// TestStart: no start from rest; pushed, the three solutions agree near
// synchronous speed; without torsion the rotor only coasts.
func TestStart(t *testing.T) {
	p := Default()
	s2, err := p.FrameIIb()
	if err != nil {
		t.Fatal(err)
	}
	rest, err := p.Start("rest", p.FrameI(), 0, 1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(rest.Final) > 1e-9 {
		t.Errorf("started from rest: %g rad/s", rest.Final)
	}
	r1, err := p.Start("I", p.FrameI(), 40, 5)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := p.Start("II", s2, 40, 5)
	if err != nil {
		t.Fatal(err)
	}
	r5, err := p.Eq5Run(40, 5)
	if err != nil {
		t.Fatal(err)
	}
	if r1.Final < 0.99*p.Sync() || r1.Final > p.Sync() || math.Abs(r1.Final-r2.Final) > 1e-3*p.Sync() || math.Abs(r1.Final-r5.Final) > 1e-5*p.Sync() {
		t.Errorf("frame I %g, II(b) %g, eq. 5′ %g", r1.Final, r2.Final, r5.Final)
	}
	if math.Abs(r1.Mismatch) > 1e-4 {
		t.Errorf("ledger mismatch %g", r1.Mismatch)
	}
	naive, err := p.Naive()
	if err != nil {
		t.Fatal(err)
	}
	rn, err := p.Start("naive", naive, 40, 5)
	if err != nil {
		t.Fatal(err)
	}
	if want := 40 * math.Exp(-5*p.Friction/p.J); math.Abs(rn.Final-want) > 1e-3*want {
		t.Errorf("without torsion %g rad/s, want friction coast-down %g", rn.Final, want)
	}
}

// TestTorqueCurve: zero torque at standstill, odd in speed, positive
// (motoring) below synchronous speed, frames agree.
func TestTorqueCurve(t *testing.T) {
	p := Default()
	s2, err := p.FrameIIb()
	if err != nil {
		t.Fatal(err)
	}
	for _, w := range []float64{0, 150, 340} {
		a, err := p.Torque(p.FrameI(), w)
		if err != nil {
			t.Fatal(err)
		}
		b, _ := p.Torque(p.FrameI(), -w)
		c, _ := p.Torque(s2, w)
		if w == 0 && math.Abs(a) > 1e-9 {
			t.Errorf("standstill torque %g", a)
		}
		if w > 0 && (a <= 0 || math.Abs(a+b) > 1e-6 || math.Abs(a-c) > 5e-3*math.Abs(a)) {
			t.Errorf("ω %g: %g, at −ω %g, frame II(b) %g", w, a, b, c)
		}
	}
}
