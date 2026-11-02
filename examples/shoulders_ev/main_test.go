// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(a, b, rel float64) bool { return math.Abs(a-b) <= rel*math.Abs(b) }

// The bead's electrostatics against closed forms, and his density figure.
func TestElectrostatics(t *testing.T) {
	e := Book()
	q := 1e11 * Qe
	if u := 0.6 * q * q / (4 * math.Pi * Eps0 * 0.5e-6); !near(e.SelfEnergy(), u, 1e-12) || !near(u, 2.77, 0.01) {
		t.Errorf("self-energy %g", e.SelfEnergy())
	}
	if !near(e.Pressure(), 1.47e18, 0.01) {
		t.Errorf("pressure %g", e.Pressure())
	}
	big := EV{N: 3.5e14, D: 10e-6}
	if !near(big.Density()/1e6, 6.6e23, 0.03) {
		t.Errorf("his density %g /cm³", big.Density()/1e6)
	}
	if !near(e.MagneticFraction(), 0.01, 1e-12) {
		t.Error("β²")
	}
}

// The explosion: energy conservation for the outer electron (it ends with
// k·Q·e·(1/R − 1/10R)) and the doubling time far below his lifetime.
func TestExplosion(t *testing.T) {
	e := Book()
	t2, _, kev := e.Explode()
	want := e.EscapeEnergy() * 0.9
	if !near(kev, want, 1e-3) {
		t.Errorf("final energy %g eV, want %g", kev, want)
	}
	if Lifetime()/t2 < 1e4 {
		t.Errorf("doubling %g s vs lifetime %g s", t2, Lifetime())
	}
	// non-relativistic small cluster: doubling time from the exact
	// R̈ = a/R² solution, t = √(R³/2a)·(√2 + ln(1+√2)).
	s := EV{N: 100, D: 1e-6}
	a := Ke * s.Charge() * Qe / Me
	r := s.R()
	exact := math.Sqrt(r*r*r/(2*a)) * (math.Sqrt2 + math.Log(1+math.Sqrt2))
	if got, _, _ := s.Explode(); !near(got, exact, 1e-3) {
		t.Errorf("doubling %g, exact %g", got, exact)
	}
}

// Containment numbers and the Gauss radius.
func TestContainment(t *testing.T) {
	if !near(CasimirPlates(1e-6), 1.3e-3, 0.01) {
		t.Errorf("Casimir %g", CasimirPlates(1e-6))
	}
	e := Book()
	if BoyerSphere(e.R()) > 1e-2 || BrillouinField(e.Density()) < 1e5 {
		t.Error("containment")
	}
	if !near(GaussRadius(e.Charge(), 1e3), 0.144, 0.01) {
		t.Errorf("Gauss radius %g", GaussRadius(e.Charge(), 1e3))
	}
}

// kronspice: the source passes at most C·V, energies balance; a strike's
// charge arrives in full and the ledger closes.
func TestKronspice(t *testing.T) {
	b := DefaultBench()
	s, err := b.Source(1e-12, 1414, 20, 5e-9, 3e-9)
	if err != nil {
		t.Fatal(err)
	}
	if q := 1e-12 * 1414; s.Charge > q*1.0001 || s.Charge < 0.99*q {
		t.Errorf("source charge %g of %g", s.Charge, q)
	}
	if s.Mismatch > 1e-3 || !near(s.Stored0, 0.5e-12*1414*1414, 1e-9) {
		t.Errorf("ledger %g, stored %g", s.Mismatch, s.Stored0)
	}
	k, err := b.Strike(2e10, 20e-12, 1e-9)
	if err != nil {
		t.Fatal(err)
	}
	if !near(k.Charge, 2e10*Qe, 1e-4) || k.Mismatch > 1e-3 || !near(k.ToScope, k.EnergyIn, 1e-3) {
		t.Errorf("strike charge %g, ledger %g, scope %g in %g", k.Charge, k.Mismatch, k.ToScope, k.EnergyIn)
	}
	// most of the charge lands on the anode first: energy ≈ Q²/2C
	if q := 2e10 * Qe; k.EnergyIn < 0.8*q*q/(2*b.CA) {
		t.Errorf("strike energy %g", k.EnergyIn)
	}
}
