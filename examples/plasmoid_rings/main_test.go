// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestBeltrami: ∇×F = αF, ∇·F = 0, field lines on tori.
func TestBeltrami(t *testing.T) {
	c, d := CurlCheck(44.93)
	if c > 1e-8 || d > 1e-8 {
		t.Errorf("curl %g, div %g", c, d)
	}
	ax := MagneticAxis(44.93)
	if _, dr := Trace(44.93, [3]float64{0.75 * ax, 0, 0}, 5000, 0.1/2000); dr > 1e-10 {
		t.Errorf("ψ drift %g", dr)
	}
}

// TestAppendixB: his charge density, electron count and field; the cold
// plasma dispersion relation; the superluminal quiver of his example.
func TestAppendixB(t *testing.T) {
	e := AppendixB()
	if math.Abs(e.RhoE+32.63) > 0.01 || math.Abs(e.Electrons/8.532e17-1) > 1e-3 || math.Abs(e.DispersionMismatch) > 1e-12 {
		t.Errorf("%+v", e)
	}
	w := 7.920e11
	if E0 := e.C * e.R * w * w / EoM; math.Abs(E0/3.566e9-1) > 1e-3 {
		t.Errorf("E0 %g", E0)
	}
	if v := w * e.C * e.R; v < C0 {
		t.Errorf("quiver %g should exceed c", v)
	}
}

// TestModes: both l = 1 resonances exist, radiate on microsecond scales and
// need different αR.
func TestModes(t *testing.T) {
	e := AppendixB()
	s := Sphere{R: e.R, Wp: e.Wp}
	te, tm := s.FindMode(true), s.FindMode(false)
	if !te.Found || !tm.Found {
		t.Fatalf("modes not found: %+v %+v", te, tm)
	}
	for _, m := range []Mode{te, tm} {
		if m.Lifetime <= 0 || m.Lifetime > 1e-5 {
			t.Errorf("lifetime %g", m.Lifetime)
		}
		if math.Abs(cmplxAbs(s.scaled(m.Omega, m.TE))) > 1e-8 {
			t.Errorf("residual at root")
		}
	}
	if math.Abs(TMZero()-2.7437) > 1e-3 {
		t.Errorf("TM zero %g", TMZero())
	}
}

// TestSpheromak: kronspice's ring circuit decays with the Beltrami
// resistive lifetime and turns all the energy into heat.
func TestSpheromak(t *testing.T) {
	s := Spheromak{R: 0.1, W: 1e6}
	d, err := s.RunDecay(10)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(d.Fitted/(d.Tau/2)-1) > 2e-3 || math.Abs(d.Heat/1e6-1) > 2e-3 || math.Abs(d.Mismatch) > 1e-4 {
		t.Errorf("%+v", d)
	}
	if T := s.TempFor(s.TauB(10) / 2); math.Abs(T/10-1) > 1e-9 {
		t.Errorf("TempFor inverse %g", T)
	}
}

// TestRing: the ring follows R ∝ B^-½, p ∝ B^½, and the induced field's
// work equals the energy gained.
func TestRing(t *testing.T) {
	r, err := DefaultCompressor().Run()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.R1/(r.R0*math.Sqrt(r.B0/r.B1))-1) > 2e-3 || math.Abs(r.P1/r.P0/math.Sqrt(r.B1/r.B0)-1) > 2e-3 {
		t.Errorf("R %g→%g, p ratio %g, B %g→%g", r.R0, r.R1, r.P1/r.P0, r.B0, r.B1)
	}
	if math.Abs(r.WorkE/r.DeltaKE-1) > 1e-3 || r.Invariant > 5e-3 {
		t.Errorf("work %g vs %g, invariant %g", r.WorkE, r.DeltaKE, r.Invariant)
	}
}
