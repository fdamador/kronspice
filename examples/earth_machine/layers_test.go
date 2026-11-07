// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Layer 2: the ladder's operating point is the exact ladder solution, the
// Earth's charge is Gauss's law at the surface, and the books balance.
func TestAtmosphereOperatingPoint(t *testing.T) {
	a := StandardAtmosphere()
	r, err := a.Solve()
	if err != nil {
		t.Fatal(err)
	}
	if d := math.Abs(r.VI-r.VIExact) / r.VIExact; d > 1e-4 {
		t.Errorf("V_I circuit %g vs exact %g (rel %g)", r.VI, r.VIExact, d)
	}
	if d := math.Abs(r.VIExact-a.VI) / a.VI; d > 2e-3 { // sphericity only
		t.Errorf("V_I %g vs anchor %g", r.VIExact, a.VI)
	}
	if d := math.Abs(r.QEarth/r.QEarthGeo - 1); d > 1e-3 {
		t.Errorf("Earth charge %g vs −ε0·E0·4πR² %g", r.QEarth, r.QEarthGeo)
	}
	if d := math.Abs(r.R*r.CEff/r.Tau - 1); d > 1e-3 {
		t.Errorf("R·C_eff = %g, ε0/σ0 = %g", r.R*r.CEff, r.Tau)
	}
	if math.Abs(r.QEarth+r.QSpace+r.QIon) > 1e-6*math.Abs(r.QEarth) {
		t.Errorf("charges do not sum to zero")
	}
	if math.Abs(r.LedgerRes) > 1e-4 {
		t.Errorf("ledger residual %g", r.LedgerRes)
	}
}

// Layer 2: the switch-on transient follows Σ I·R_k·(1 − e^(−t/R_kC_k)).
func TestAtmosphereSwitchOn(t *testing.T) {
	a := StandardAtmosphere()
	tt, vi, ex, err := a.SwitchOn(2000, 1)
	if err != nil {
		t.Fatal(err)
	}
	for k := range tt {
		if d := math.Abs(vi[k]-ex[k]) / a.VI; d > 5e-4 {
			t.Fatalf("t = %g: circuit %g, exact %g", tt[k], vi[k], ex[k])
		}
	}
}

// Layer 2: at one cycle per day the ladder is quasi-static, with a lag of
// a few minutes (below the relaxation time ε0/σ0).
func TestAtmosphereCarnegie(t *testing.T) {
	a := StandardAtmosphere()
	ratio, lag, err := a.Carnegie(0.15)
	if err != nil {
		t.Fatal(err)
	}
	s0, _ := a.Calibrate()
	if ratio < 0.999 || ratio > 1 || lag <= 0 || lag > eps0/s0 {
		t.Errorf("ratio %g, lag %g s", ratio, lag)
	}
}

// Layer 3: the perfect-cavity ladder converges to Legendre's modes at
// second order in the band width, and the engine's AC peaks sit on them.
func TestCavityLegendre(t *testing.T) {
	err := func(bands int) float64 {
		c := IdealCavity()
		c.Bands = bands
		return math.Abs(c.Modes(1)[0]/LegendreFreq(1) - 1)
	}
	e1, e2 := err(91), err(181)
	if e2 > 1e-4 || e1/e2 < 3.5 {
		t.Errorf("mode-1 error %g at 91 bands, %g at 181 (ratio %g, want ≈ 4)", e1, e2, e1/e2)
	}
	c := IdealCavity()
	m := c.Modes(5)
	for n := 1; n <= 5; n++ {
		if d := math.Abs(m[n-1]/LegendreFreq(n) - 1); d > 1e-3 {
			t.Errorf("n = %d: ladder %g, Legendre %g", n, m[n-1], LegendreFreq(n))
		}
	}
	c.Alpha = 2 * math.Pi * m[0] / 600
	sp, err2 := c.Sweep(5, 45, 4001)
	if err2 != nil {
		t.Fatal(err2)
	}
	pk := Peaks(sp.F, sp.Pole)
	if len(pk) != 5 {
		t.Fatalf("peaks %v", pk)
	}
	for n := range pk {
		if math.Abs(pk[n]-m[n]) > 0.01 {
			t.Errorf("peak %d at %g Hz, mode at %g", n+1, pk[n], m[n])
		}
	}
	if sp.LedgerRes > 1e-9 {
		t.Errorf("AC ledger residual %g", sp.LedgerRes)
	}
}

// Layer 3: the fit reproduces mode 1 exactly and odd modes vanish at the
// equator.
func TestCavityLossyFit(t *testing.T) {
	c := FitLossy(7.8, 5)
	w0 := 2 * math.Pi * c.Modes(1)[0]
	fd := math.Sqrt(w0*w0-c.Alpha*c.Alpha) / (2 * math.Pi)
	if math.Abs(fd-7.8) > 1e-6 || math.Abs(w0/(2*c.Alpha)-5) > 1e-6 {
		t.Errorf("fit: f %g, Q %g", fd, w0/(2*c.Alpha))
	}
	sp, err := c.Sweep(2, 45, 2151)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range Peaks(sp.F, sp.Equator) { // even modes only: near 13.6, 24.8, 35.9 Hz
		if p < 10 || (p > 16 && p < 22) {
			t.Errorf("equator peak at %g Hz belongs to an odd mode", p)
		}
	}
}

// Layer 7: with a linear ionosphere the two transmitters only add; with
// electron heating the third-order products appear far above the linear
// floor, and the ledger closes in both.
func TestNonlinearIonosphere(t *testing.T) {
	p := DefaultNLIonosphere()
	d, err := p.DriveFor(0.5)
	if err != nil {
		t.Fatal(err)
	}
	p.Drive, p.Linear = d, true
	lin, err := p.Run()
	if err != nil {
		t.Fatal(err)
	}
	p.Linear = false
	nl, err := p.Run()
	if err != nil {
		t.Fatal(err)
	}
	pr := p.Products()
	for _, k := range []string{"2f1-f2", "2f2-f1"} {
		a, b := absBin(nl.Lines, pr[k]), absBin(lin.Lines, pr[k])
		if a < 100*b || a < 1e-3*absBin(nl.Lines, pr["f2"]) {
			t.Errorf("%s: nonlinear %.3g V, linear %.3g V", k, a, b)
		}
	}
	if lin.Ledger > 1e-4 || nl.Ledger > 1e-4 {
		t.Errorf("ledger %.2g, %.2g", lin.Ledger, nl.Ledger)
	}
}
