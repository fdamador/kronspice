// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// Kron (Happ two-level) against MNA AC. Off resonance the two solutions
// agree to 1e−9 relative (CLAUDE.md §3a rule 4). On a resonance of
// Q ≈ 1.5e4 the system is ill-conditioned (the forward error is the
// backward error times a condition number of order Q·n): both solutions
// are checked by their backward error instead, ≤ 1e−12, and their
// difference stays below 1e−6.
func TestTwinMNA(t *testing.T) {
	g := Default
	q := g.N / 4
	four := []Emitter{{0, 1, 0}, {q, 1, 30}, {2 * q, 0.5, 0}, {3 * q, 1, -60}}
	diff := func(f float64) (float64, float64, float64) {
		r, err := Solve(g, f, four, 0)
		if err != nil {
			t.Fatal(err)
		}
		vm, err := TwinMNA(g, f, four)
		if err != nil {
			t.Fatal(err)
		}
		var d, vmax float64
		for i := range vm {
			vmax = math.Max(vmax, cmplx.Abs(vm[i]))
		}
		for i := range vm {
			d = math.Max(d, cmplx.Abs(vm[i]-r.V[i])/vmax)
		}
		return d, residual(g, f, four, r.V, 1), residual(g, f, four, vm, 1)
	}
	next := g
	next.Mode++
	fOff := (g.Freq() + next.Freq()) / 2
	if d, _, _ := diff(fOff); d > 1e-9 {
		t.Errorf("off resonance: Kron vs MNA %.2g", d)
	}
	d, bk, bm := diff(g.Freq())
	if bk > 1e-12 || bm > 1e-12 || d > 1e-6 {
		t.Errorf("on resonance: difference %.2g, backward errors Kron %.2g, MNA %.2g", d, bk, bm)
	}
}

// The ring's resonance and Q against the telegrapher analytic values,
// and the energy ledger.
func TestRing(t *testing.T) {
	g := Default
	f := g.Freq()
	cont := float64(g.Mode) * c0 / (2 * math.Pi * g.R)
	// Lumped-ladder dispersion: f/f_cont = sin x / x, x = πn/N (≈ 1e−3 low
	// at 40 sections per wavelength).
	x := math.Pi * float64(g.Mode) / float64(g.N)
	if math.Abs(f/cont-math.Sin(x)/x) > 1e-12 || math.Abs(f/cont-1) > 1.1e-3 {
		t.Errorf("mode %d: discrete %.9g Hz, continuum %.9g Hz", g.Mode, f, cont)
	}
	q := g.N / 4
	r, err := Solve(g, f, []Emitter{{0, 1, 0}, {q, 1, 0}, {2 * q, 1, 0}, {3 * q, 1, 0}}, 1e6)
	if err != nil {
		t.Fatal(err)
	}
	qa := 2 * math.Pi * f * mu0 * g.D / g.H / (2 * g.SurfaceR(f) / g.H) // ωL′/R′
	if math.Abs(r.Q/qa-1) > 2e-3 {
		t.Errorf("Q %.6g, analytic ωL′/R′ %.6g", r.Q, qa)
	}
	if math.Abs(r.Pin-r.Ploss) > 1e-9*r.Pin {
		t.Errorf("ledger: in %g, wall loss %g", r.Pin, r.Ploss)
	}
	for i, x := range r.Quadrant {
		if math.Abs(x-0.25) > 1e-3 {
			t.Errorf("quadrant %d holds %.4f of the energy", i, x)
		}
	}
	// High Q: a single emitter cannot confine energy to its section;
	// only a ring with decay length shorter than the circumference can.
	one, err := Solve(g, f, []Emitter{{0, 1, 0}}, 1e6)
	if err != nil {
		t.Fatal(err)
	}
	lossy := g
	lossy.Rs = 2000 * g.SurfaceR(f)
	loc, err := Solve(lossy, f, []Emitter{{0, 1, 0}}, 1e6)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(one.Quadrant[0]-0.25) > 1e-3 || loc.Quadrant[0] < 0.7 || loc.Q > 10 {
		t.Errorf("localization: high-Q top quadrant %.3f, lossy %.3f at Q %.3g", one.Quadrant[0], loc.Quadrant[0], loc.Q)
	}
}

// The patent's figures and the ceilings they meet.
func TestPatentFigures(t *testing.T) {
	if s := PatentFlux(5e4, 2, 30000*2*math.Pi/60); s < 1e24 || s > 1e25 {
		t.Errorf("rotation flux %.3g, patent ≈ 1e24", s)
	}
	if v := 1e33 * eps0; v < 1e4*c0 {
		t.Errorf("R_v·ν for 1e33 W/m² = %.3g m/s", v)
	}
	r, err := Solve(Default, Default.Freq(), []Emitter{{0, 1, 0}}, 1e6)
	if err != nil {
		t.Fatal(err)
	}
	l := LimitsFor(1, r)
	if math.Abs(l.E/1.1294e11-1) > 1e-3 || l.Pressure < 100*l.CopperUTS || l.DrainTime > 1e-12 {
		t.Errorf("σ = 1 C/m²: %+v", l)
	}
	if math.Abs(l.Schwinger/1.3233e18-1) > 1e-3 {
		t.Errorf("Schwinger field %.5g", l.Schwinger)
	}
	if d := Debye(343); math.Abs(d/7.147e12-1) > 1e-3 {
		t.Errorf("Cu Debye frequency %.4g", d)
	}
}
