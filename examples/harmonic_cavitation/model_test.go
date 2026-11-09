// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// FitDelta makes rpm·g(σ(rpm)) peak at Xue's best speed (3800 rpm).
func TestDeltaFitPeaksAtBestSpeed(t *testing.T) {
	p := DefaultParams()
	best, at := 0.0, 0.0
	for rpm := 3000.0; rpm <= 6000; rpm += 1 {
		u := omega(rpm) * Xue2024.R
		if v := rpm * p.G(Sigma(Xue2024.PRef, u)); v > best {
			best, at = v, rpm
		}
	}
	if math.Abs(at-3800) > 2 {
		t.Fatalf("collapse rate peaks at %.0f rpm, want 3800", at)
	}
	if p.G(p.SigmaI+0.01) != 0 || math.Abs(p.G(p.SigmaI-p.Delta)-1) > 1e-12 {
		t.Fatal("g(σ) shape")
	}
}

// A hemisphere: volume 2πa³/3, projected (half-disc) area πa²/2.
func TestHemicap(t *testing.T) {
	a := 0.004
	c := Hemicap(2*a, a, 1)
	if math.Abs(c.V/(2*math.Pi*a*a*a/3)-1) > 1e-12 || math.Abs(c.A/(math.Pi*a*a/2)-1) > 1e-12 {
		t.Fatalf("V %.4g A %.4g", c.V, c.A)
	}
}

// The kronspice shaft circuit returns the model's shaft power with a
// closed ledger, and collapse energy never exceeds the cavity drag on the
// reference band.
func TestShaftLedgerTwin(t *testing.T) {
	p := DefaultParams()
	m := DimpleMachine(0.05, 10, 0.010, 0.5, 1.5, 2.2, 1.5e-3, 3050, pAtm)
	e, err := p.Evaluate(m)
	if err != nil {
		t.Fatal(err)
	}
	pin, res, err := LedgerTwin(e, m.RPM)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(pin/e.PShaft-1) > 1e-9 || res > 1e-9 {
		t.Fatalf("twin %.6g vs %.6g, residual %.1e", pin, e.PShaft, res)
	}
	if !e.Consistent || e.Eta <= 0 || e.Eta >= 1 {
		t.Fatalf("η %.4f consistent %v", e.Eta, e.Consistent)
	}
	// housing cavities raise η here (pulse) — the plain-housing control is lower
	m.Housing = nil
	e2, _ := p.Evaluate(m)
	if !(e2.Eta < e.Eta) {
		t.Fatalf("plain housing η %.4f not below dimpled %.4f", e2.Eta, e.Eta)
	}
}

// Rayleigh: an empty cavity with no gas, viscosity or surface tension
// collapses in 0.915·R0·√(ρ/Δp); with gas the ledger closes.
func TestRayleighCollapse(t *testing.T) {
	r0, dp := 1e-3, 1e5
	c := RunCollapse(r0, dp, 0, 0, 0)
	want := 0.915 * r0 * math.Sqrt(rho/dp)
	if math.Abs(c.T/want-1) > 2e-3 {
		t.Fatalf("collapse %.4g s against Rayleigh %.4g s", c.T, want)
	}
	g := RunCollapse(r0, dp, 1e-3, 0.0717, 853.8e-6)
	if g.LedgerRelResidual > 1e-4 || g.Rmin <= 0 {
		t.Fatalf("gas collapse: R_min %.3g, ledger %.1e", g.Rmin, g.LedgerRelResidual)
	}
	if Simultaneous(35, 35) != 35 || Simultaneous(35, 37) != 1 || Simultaneous(12, 18) != 6 {
		t.Fatal("gcd")
	}
}

// The gap ring in kronspice: peaks at the exact lattice modes, which tend
// to m·c/2πR.
func TestRingModes(t *testing.T) {
	pk, ex, err := RingModes(0.05, 1.5e-3, 0.016, 100, 64, 4)
	if err != nil {
		t.Fatal(err)
	}
	if len(pk) < 4 {
		t.Fatalf("%d peaks", len(pk))
	}
	for m := 0; m < 4; m++ {
		if math.Abs(pk[m]/ex[m]-1) > 1e-4 {
			t.Fatalf("mode %d: %.3f against %.3f", m+1, pk[m], ex[m])
		}
		// lattice dispersion: sin(x)/x with x = π m/64 lowers the mode by ≈ x²/6
		cont := float64(m+1) * 100 / (2 * math.Pi * 0.05)
		x := math.Pi * float64(m+1) / 64
		if math.Abs(ex[m]/cont-math.Sin(x)/x) > 1e-12 || math.Abs(ex[m]/cont-1) > 1.1*x*x/6 {
			t.Fatalf("lattice mode %d %.3f far from continuum %.3f", m+1, ex[m], cont)
		}
	}
}

// Mandelbulb: the latitude convention's equator is the Multibrot set; both
// conventions have (n − 1)-fold symmetry about z; n = 8 shows 7 lobes.
func TestMandelbulb(t *testing.T) {
	for _, n := range []float64{3, 5, 8} {
		b := Bulb{N: n, Latitude: true, MaxIt: 12}
		diff := 0
		for i := 0; i < 2000; i++ {
			x := -1.4 + 2.8*float64(i%50)/49
			y := -1.4 + 2.8*float64(i/50)/39
			if b.In(x, y, 0) != MultibrotIn(n, x, y, 12) {
				diff++
			}
		}
		if diff != 0 {
			t.Fatalf("n %g: equator differs from Multibrot at %d of 2000 points", n, diff)
		}
	}
	for _, lat := range []bool{false, true} {
		b := Bulb{N: 8, Latitude: lat, MaxIt: 10}
		a := 2 * math.Pi / 7
		diff, tot := 0, 0
		for i := 0; i < 3000; i++ {
			x := -1.2 + 2.4*float64(i%15)/14
			y := -1.2 + 2.4*float64((i/15)%15)/14
			z := -1.0 + 2.0*float64(i/225)/13
			xr, yr := x*math.Cos(a)-y*math.Sin(a), x*math.Sin(a)+y*math.Cos(a)
			tot++
			if b.In(x, y, z) != b.In(xr, yr, z) {
				diff++
			}
		}
		if float64(diff)/float64(tot) > 0.01 { // round-off at the boundary only
			t.Fatalf("latitude %v: 7-fold symmetry broken at %d of %d points", lat, diff, tot)
		}
	}
	if n := Lobes(Bulb{N: 8, MaxIt: 10}.Profile(0, 720), 0.01); n != 7 {
		t.Fatalf("n = 8 equator shows %d lobes, want 7", n)
	}
}

// Valleys on a synthetic profile: four equal notches.
func TestValleys(t *testing.T) {
	p := make([]float64, 400)
	for i := range p {
		p[i] = 1
		if (i%100) >= 40 && (i%100) < 60 {
			p[i] = 0.8
		}
	}
	v := Valleys(p, 0.01)
	if len(v) != 4 {
		t.Fatalf("%d valleys", len(v))
	}
	if math.Abs(v[0].Width-20*2*math.Pi/400) > 1e-12 || math.Abs(v[0].Depth-0.2) > 1e-12 {
		t.Fatalf("valley %+v", v[0])
	}
}

// Energy books: heat ≤ electricity, the turbine drive needs more
// electricity for the same shaft power.
func TestDrives(t *testing.T) {
	p := DefaultParams()
	e, _ := p.Evaluate(DimpleMachine(0.036, 3, 0.0043, 0.5, 1.47, 1.4, 0.8e-3, 5099, 1.29e5))
	d := Drives(e)
	for _, x := range d {
		if x.PWaterHeat+x.PRoomHeat > x.PElec*(1+1e-12) || x.CavPerJoule > e.Eta {
			t.Fatalf("%s: heat %.1f + %.1f against %.1f", x.Name, x.PWaterHeat, x.PRoomHeat, x.PElec)
		}
	}
	if !(d[1].PElec > d[0].PElec) {
		t.Fatal("turbine drive should need more electricity")
	}
}

// The GA is deterministic, finds designs within the limits, and does not
// do worse than the reference band's η.
func TestEvolveSmall(t *testing.T) {
	if testing.Short() {
		t.Skip("GA")
	}
	p := DefaultParams()
	a, err := Evolve(p, dimpleSpace, DimpleFrom, dimpleLand, 30, 20, 7, nil)
	if err != nil || len(a) == 0 {
		t.Fatal(err)
	}
	b, _ := Evolve(p, dimpleSpace, DimpleFrom, dimpleLand, 30, 20, 7, nil)
	if a[0].E.Eta != b[0].E.Eta {
		t.Fatal("not deterministic")
	}
	for _, s := range a {
		if s.E.PShaft > PShaftMax*1.0001 || !s.E.Consistent {
			t.Fatalf("infeasible design on the front: %+v", s.E)
		}
	}
	if a[0].E.Eta < 0.12 {
		t.Fatalf("best η %.4f", a[0].E.Eta)
	}
}
