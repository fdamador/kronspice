// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestEq3: Kaluza's closed forms (eq. 3) follow from his metric (eq. 2).
func TestEq3(t *testing.T) {
	if d := Eq3Check(0.7); d > 1e-10 {
		t.Errorf("eq. 3 off by %g", d)
	}
	if a := KaluzaAlpha(); math.Abs(a/3.06e-14-1) > 0.005 {
		t.Errorf("α = %g, Kaluza 3.06e-14", a)
	}
}

// TestGeodesicIsLorentz: the 5-D geodesic conserves p_w and its 4-D path
// is the Lorentz-force orbit with q/m = k·p_w/√(1+p_w²).
func TestGeodesicIsLorentz(t *testing.T) {
	for _, pw := range []float64{0.1, 1, 10} {
		o := RunOrbit(0.05, 2, 0.3, pw)
		if math.Abs(o.Radius/o.RadiusLF-1) > 1e-5 || math.Abs(o.Omega/o.OmegaLF-1) > 1e-5 || o.PwDrift > 1e-9 || o.NormDrift > 1e-9 {
			t.Errorf("p_w %g: radius %g vs %g, ω %g vs %g, drifts %g %g", pw, o.Radius, o.RadiusLF, o.Omega, o.OmegaLF, o.PwDrift, o.NormDrift)
		}
	}
}

// TestKlein: R = 2·l_P/√α, and the electron's charge-to-mass is far above
// the classical bound 2.
func TestKlein(t *testing.T) {
	if r := KleinRadius(); math.Abs(r/(2*LPlanckM/math.Sqrt(AlphaFS))-1) > 1e-4 {
		t.Errorf("Klein radius %g", r)
	}
	if q := ChargeToMass(Ecgs, MeCgs); math.Abs(q/2.041e21-1) > 1e-3 {
		t.Errorf("electron e/(m√G) %g", q)
	}
}

// TestMachine: Kron's geodesic, the reduced rotor and kronspice agree;
// flux linkages and energy are conserved.
func TestMachine(t *testing.T) {
	m := DefaultMachine()
	g := m.Geodesic(1, 1e-4)
	r := m.Reduced(1, 1e-4)
	k, mis, err := m.KronSpice(1, 1e-4, 1e-7)
	if err != nil {
		t.Fatal(err)
	}
	if g.FluxDrift > 1e-12 || g.EnergyDrift > 1e-12 {
		t.Errorf("drifts %g %g", g.FluxDrift, g.EnergyDrift)
	}
	for i := range g.Time {
		if math.Abs(g.Omega[i]-at(r, g.Time[i])) > 1e-8 || math.Abs(g.Omega[i]-at(k, g.Time[i])) > 2e-3 {
			t.Fatalf("t %g: geodesic %g, reduced %g, kronspice %g", g.Time[i], g.Omega[i], at(r, g.Time[i]), at(k, g.Time[i]))
		}
	}
	if math.Abs(mis) > 1e-6 {
		t.Errorf("ledger mismatch %g", mis)
	}
}
