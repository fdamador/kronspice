// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s = %.6g, want %.6g ± %.2g", what, got, want, tol)
	}
}

// Free precession in the laboratory frame: frequency γB₀ (trapezoidal
// phase error (ωh)²/12 at 400 points per period), T2 decay of |Mxy|, T1
// recovery of Mz, ledger closed.
func TestFreePrecession(t *testing.T) {
	s := Proton()
	s.Minit = [3]float64{1, 0, 0}
	w0 := s.Gamma * s.B0
	r, err := s.Run(2*s.T2, 2*math.Pi/w0/400)
	if err != nil {
		t.Fatal(err)
	}
	near(t, "precession / γB₀", Frequency(r.Time, r.Mx)/w0, 1, 1e-4)
	k := len(r.Time) - 1
	near(t, "|Mxy| at 2T2", math.Hypot(r.Mx[k], r.My[k]), math.Exp(-2), 1e-4)
	near(t, "Mz at 2T2", r.Mz[k], 1-math.Exp(-2*s.T2/s.T1), 1e-4)
	near(t, "ledger mismatch / heat", r.Ledger.IntegrationMismatch/r.Ledger.EnergyDissipated, 0, 1e-3)
	near(t, "source − heat − Δstored", (r.Ledger.EnergyIn-r.Ledger.EnergyOut-r.Ledger.EnergyDissipated-(r.Ledger.StoredEnd-r.Ledger.StoredStart))/r.Ledger.EnergyDissipated, 0, 1e-3)
}

// Rabi nutation at resonance: Mz = M₀ cos(γB₁t); a π/2 pulse tips M onto
// the transverse plane.
func TestRabi(t *testing.T) {
	s := Proton()
	s.B1, s.Wrf, s.T1, s.T2 = 1e-3, s.Gamma*s.B0, 1, 1
	wr := s.Gamma * s.B1
	r, err := s.Run(2*math.Pi/wr, 2*math.Pi/wr/400)
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k < len(r.Time); k += 50 {
		near(t, "Mz(t)", r.Mz[k], math.Cos(wr*r.Time[k]), 1e-3)
	}
	tp := math.Pi / (2 * wr)
	p, err := s.Run(tp, tp/400)
	if err != nil {
		t.Fatal(err)
	}
	m := p.Final()
	near(t, "Mx after π/2", m[0], 0, 1e-3)
	near(t, "My after π/2", m[1], 1, 1e-3)
	near(t, "Mz after π/2", m[2], 0, 1e-3)
}

// The free-induction decay seen by the pickup winding: M_p·dMx/dt at the
// Larmor frequency, peak M_p·γB₀·M₀.
func TestPickup(t *testing.T) {
	s := Proton()
	s.Minit = [3]float64{0, 1, 0}
	w0 := s.Gamma * s.B0
	r, err := s.Run(s.T2, 2*math.Pi/w0/400)
	if err != nil {
		t.Fatal(err)
	}
	vmax := 0.0
	for _, v := range r.Pickup {
		vmax = math.Max(vmax, math.Abs(v))
	}
	near(t, "pickup peak / M_p γB₀ M₀", vmax/(s.Mp*w0*s.M0), 1, 1e-3)
	near(t, "pickup frequency / γB₀", Frequency(r.Time, r.Pickup)/w0, 1, 1e-4)
}

// Continuous-wave steady state against Bloch's solution.
func TestCWLine(t *testing.T) {
	s := Proton()
	s.B1 = 1e-5
	w0 := s.Gamma * s.B0
	for _, dT2 := range []float64{-2, -0.5, 0, 0.5, 2} {
		dw := dT2 / s.T2
		s.Wrf = w0 + dw
		r, err := s.Run(12*s.T1, s.T2/50)
		if err != nil {
			t.Fatal(err)
		}
		bx, by, bz := BlochSteady(s, dw)
		near(t, "absorption", last(r.My), by, 2e-5)
		near(t, "dispersion", last(r.Mx), bx, 2e-5)
		near(t, "Mz", last(r.Mz), bz, 2e-5)
	}
}
