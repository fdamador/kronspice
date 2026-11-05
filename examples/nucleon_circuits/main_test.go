// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// A light ring of radius ħ/(mc) carrying charge e at c has the Dirac
// moment eħ/(2m): for the proton, exactly one nuclear magneton, from the
// kronspice circuit's own current.
func TestLightRingDirac(t *testing.T) {
	p := Ring{Name: "p", Q: 1, R: LightRadius(MpMeV), Up: 1}
	run, err := Circuit([]Ring{p})
	if err != nil {
		t.Fatal(err)
	}
	if m := run.Moment / MuN0; math.Abs(m-1) > 1e-9 {
		t.Errorf("moment %.12f μN, want 1", m)
	}
	if r := p.R * 1e15; math.Abs(r-0.21030891) > 1e-6 {
		t.Errorf("radius %.8f fm, want ħ/(m_p c) = 0.21030891", r)
	}
}

// The circuit's stored magnetic energy equals ½ΣLI² + ΣMI1I2 for three
// coupled rings with currents of both signs, and its ledger closes.
func TestCircuitTwin(t *testing.T) {
	R := LightRadius(300)
	rings := []Ring{{"a", 2.0 / 3, R, 1, -R}, {"b", -1.0 / 3, 0.8 * R, -1, 0}, {"c", 2.0 / 3, 1.3 * R, 1, 0.7 * R}}
	run, err := Circuit(rings)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(run.Stored/run.Analytic-1) > 1e-6 {
		t.Errorf("stored %.9g J, analytic %.9g J", run.Stored, run.Analytic)
	}
	for i, r := range rings {
		if math.Abs(run.Currents[i]/r.Current()-1) > 1e-9 {
			t.Errorf("ring %s current %.6g A, set %.6g A", r.Name, run.Currents[i], r.Current())
		}
	}
	if run.Residual > 1e-9 {
		t.Errorf("ledger residual %.2g", run.Residual)
	}
}

// Three quark rings, spin-weighted: with the ring energy fixed by the
// proton's moment, the neutron's comes out −2/3 of it (−1.862 μN against
// −1.913), and electromagnetism makes the neutron lighter than the proton
// by about 1 MeV at the declared spacing.
func TestQuarkRings(t *testing.T) {
	mq := MpMeV / MuP
	R := LightRadius(mq)
	muP, emP, _, err := Weighted(QuarkNucleon(2.0/3, -1.0/3, R, R))
	if err != nil {
		t.Fatal(err)
	}
	muN, emN, _, err := Weighted(QuarkNucleon(-1.0/3, 2.0/3, R, R))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(muP-MuP) > 1e-6 || math.Abs(muN/muP+2.0/3) > 1e-9 {
		t.Errorf("μp %.6f, μn %.6f (ratio %.6f)", muP, muN, muN/muP)
	}
	if d := emN - emP; d > -0.5 || d < -1.5 {
		t.Errorf("electromagnetic m_n − m_p %.3f MeV, want about −1", d)
	}
}

// The neutron as p + e + ν̄ light rings: the electron ring that gives the
// measured moment is 0.99 fm across, and light on it carries about 200 MeV,
// two orders of magnitude more than the 1.293 MeV the neutron has over the
// proton.
func TestNeutronFromDecayRings(t *testing.T) {
	need := (MuP - MuN) * MuN0
	r := 2 * need / (E * C0)
	if math.Abs(r*1e15-0.9897) > 1e-3 {
		t.Errorf("electron ring radius %.4f fm", r*1e15)
	}
	if e := Hbar * C0 / r / MeV; e < 100*(MnMeV-MpMeV) {
		t.Errorf("ring energy %.1f MeV", e)
	}
}

// The standard-model lifetime from V_ud and λ is within 0.5 % of the
// measured 878.4 s.
func TestLifetime(t *testing.T) {
	if tau := NeutronLifetime(Vud, Lambda); math.Abs(tau/TauN-1) > 0.005 {
		t.Errorf("τ %.1f s", tau)
	}
}
