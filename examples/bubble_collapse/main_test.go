// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"strings"
	"testing"
)

const fTest = 25.02e3

// Without losses the circuit's state equation is Rayleigh–Plesset:
// ρ(RR̈ + 3Ṙ²/2) = p_g − 2σ/R − p∞, from the Christoffel term ½M′UU_M.
func TestRayleighPlessetLimit(t *testing.T) {
	b := DefaultBubble(fTest)
	dr := b.drive()
	for _, c := range [][3]float64{{4.5e-6, 0, 0.1e-5}, {20e-6, -30, 0.3e-5}, {2e-6, -400, 0.7e-5}, {30e-6, 5, 1.1e-5}} {
		R, Rd, tt := c[0], c[1], c[2]
		M := b.Rho / (4 * math.Pi * R)
		y := state{R, M * 4 * math.Pi * R * R * Rd}
		d, _ := b.circuitDeriv(tt, y, dr, Losses{})
		// R̈ from dΠ/dt = d/dt(ρ R Ṙ) = ρ(ṘṘ + R R̈)
		Rdd := (d[1]/b.Rho - Rd*Rd) / R
		pg, _ := b.gas(R)
		want := ((pg-2*b.Sigma/R-(b.P0+dr(tt)))/b.Rho - 1.5*Rd*Rd) / R
		if math.Abs(d[0]-Rd) > 1e-12*math.Max(1, math.Abs(Rd)) || math.Abs(Rdd-want) > 1e-9*math.Abs(want) {
			t.Errorf("R=%g Ṙ=%g: Ṙ %g, R̈ %g, Rayleigh–Plesset %g", R, Rd, d[0], Rdd, want)
		}
	}
}

// At small amplitude the bubble is a linear element whose admittance is
// 1/(Z_rad + R_visc + K/jω), Z_rad the pulsating sphere's radiation
// impedance (Kinsler et al. §7.2) and K the gas and surface stiffness.
func TestLinearAdmittance(t *testing.T) {
	b := DefaultBubble(fTest)
	b.Pa = 1e-4 * atm
	cy, _, err := b.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	w := 2 * math.Pi * b.F
	R, h := b.R0, b.h()
	M := b.Rho / (4 * math.Pi * R)
	Rrad := b.Rho * b.C / (4 * math.Pi * R * R)
	zrad := complex(0, w*M) * complex(Rrad, 0) / complex(Rrad, w*M)
	V0, Vh := 4*math.Pi/3*R*R*R, 4*math.Pi/3*h*h*h
	K := b.Gamma*b.pg0()/(V0-Vh) - 2*b.Sigma/(4*math.Pi*R*R*R*R)
	z := zrad + complex(b.Mu/(math.Pi*R*R*R), 0) + complex(0, -K/w)
	want := 1 / z
	if d := cmplx.Abs(cy.Y1-want) / cmplx.Abs(want); d > 1e-4 {
		t.Errorf("Y₁ %v, linear theory %v (%.2g)", cy.Y1, want, d)
	}
	if r := math.Abs(real(cy.Y1)-real(want)) / real(want); r > 1e-3 {
		t.Errorf("Re Y₁ %g, linear theory %g (%.2g)", real(cy.Y1), real(want), r)
	}
}

// The driven cycle: ledger closed, fundamental power = cycle work, and the
// Keller–Miksis twin close on the slow quantities.
func TestCycleLedgerAndTwin(t *testing.T) {
	b := DefaultBubble(fTest)
	cy, s, err := b.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	if r := math.Abs(cy.Residual) / cy.In; r > 1e-8 {
		t.Errorf("ledger residual %.2g of the work in", r)
	}
	if p := 0.5 * real(cy.Y1) * b.Pa * b.Pa; math.Abs(p-cy.In*b.F) > 1e-6*p {
		t.Errorf("½Re(Y₁)Pa² = %g W, cycle work × f = %g W", p, cy.In*b.F)
	}
	kmax, kmin, ktmin, _, err := b.RunKM(6)
	if err != nil {
		t.Fatal(err)
	}
	T := 1 / b.F
	if math.Abs(cy.Rmax-kmax) > 0.01*kmax || math.Abs(cy.Rmin-kmin) > 0.1*kmin ||
		math.Abs(math.Mod(cy.TRmin, T)-math.Mod(ktmin, T)) > 0.005*T {
		t.Errorf("circuit R %g–%g at %g, Keller–Miksis %g–%g at %g", cy.Rmin, cy.Rmax, cy.TRmin, kmin, kmax, ktmin)
	}
	L := b.FlashEstimate(s)
	if L.FWHM < 20e-12 || L.FWHM > 1e-9 || L.Energy/cy.CollapseWork > 1e-4 {
		t.Errorf("flash %g J in %g s", L.Energy, L.FWHM)
	}
}

// The drive chain: the bubble element takes ½Re(Y₁)Pa², the mode's heat
// at 1.3 atm is unchanged, and the ledger closes.
func TestChainWithBubble(t *testing.T) {
	b := DefaultBubble(fTest)
	cy, _, err := b.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	d0, err := readDrive(0, b.Pa)
	if err != nil {
		t.Fatal(err)
	}
	d1, err := readDrive(cy.Y1, b.Pa)
	if err != nil {
		t.Fatal(err)
	}
	want := 0.5 * real(cy.Y1) * b.Pa * b.Pa
	if math.Abs(d1.Bubble-want) > 1e-3*want || math.Abs(d1.Mode-d0.Mode) > 1e-3*d0.Mode {
		t.Errorf("bubble %g W (want %g), mode %g W (without %g)", d1.Bubble, want, d1.Mode, d0.Mode)
	}
	if d1.Residual > 1e-9 || !(d1.Source > d0.Source) {
		t.Errorf("residual %g, source %g vs %g", d1.Residual, d1.Source, d0.Source)
	}
}

// Without losses the collapse runs backwards when the momentum is flipped;
// with either loss it does not.
func TestReversal(t *testing.T) {
	b := DefaultBubble(fTest)
	cy, s, err := b.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	T := 1 / b.F
	var tMax, rm float64
	for _, x := range s {
		if x.T < math.Mod(cy.TRmin, T) && x.R > rm {
			rm, tMax = x.R, x.T
		}
	}
	tStart := 5*T + tMax
	half := cy.TRmin - tStart
	dR, dPi, err := b.Reversal(Losses{}, tStart, half, nil)
	if err != nil {
		t.Fatal(err)
	}
	if dR > 1e-9 || dPi > 1e-9 {
		t.Errorf("lossless: ΔR %g, ΔΠ %g", dR, dPi)
	}
	for _, l := range []Losses{{true, false}, {false, true}} {
		dR, _, err := b.Reversal(l, tStart, half, nil)
		if err != nil {
			t.Fatal(err)
		}
		if dR < 1e-2 {
			t.Errorf("%+v: returned within %g", l, dR)
		}
	}
}

// Saha: n_e² = s·(n − n_e) by construction, vanishing ionization when cold.
func TestIonization(t *testing.T) {
	if x := Ionization(1e28, 300) / 1e28; x > 1e-100 {
		t.Errorf("ionization %g at 300 K", x)
	}
	ne := Ionization(1e28, 2e4)
	if ne <= 0 || ne >= 1e28 {
		t.Errorf("n_e %g", ne)
	}
}

// V > 0 for a field turning from θ̂ to φ̂ (spin along r̂).
func TestStokesSign(t *testing.T) {
	n := 256
	a, b := make([]float64, n), make([]float64, n)
	for i := range a {
		a[i] = math.Cos(2 * math.Pi * 16 * float64(i) / float64(n))
		b[i] = math.Sin(2 * math.Pi * 16 * float64(i) / float64(n))
	}
	V, I := stokes([][]float64{a}, [][]float64{b})
	if math.Abs(V/I-1) > 1e-12 {
		t.Errorf("V/I = %g, want +1", V/I)
	}
}

// Field part on a coarse setup: the measurement sees the hopfion's
// helicity as its exact field does, the chiral pair near 1, the rotating
// dipole and the thermal core near 0, all ledgers closed.
func TestHelicity(t *testing.T) {
	if testing.Short() {
		t.Skip("field solver")
	}
	s := DefaultHelicity()
	s.Probes, s.Dipoles, s.Seeds = 200, 150, 2
	h, err := s.RunAll()
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range h.Runs {
		if r.Residual > 1e-10 {
			t.Errorf("%s: ledger residual %g", r.Name, r.Residual)
		}
	}
	hop, chi, rot := h.Runs[0], h.Runs[1], h.Runs[2]
	if math.Abs(hop.Pc-hop.PcExact) > 0.01 || hop.PcExact < 0.99 {
		t.Errorf("hopfion P_c %g, exact %g", hop.Pc, hop.PcExact)
	}
	if chi.Pc < 0.9 || math.Abs(rot.Pc) > 0.02 {
		t.Errorf("chiral %g, rotating %g", chi.Pc, rot.Pc)
	}
	if math.Abs(h.ThermalMean) > 0.1 {
		t.Errorf("thermal mean %g", h.ThermalMean)
	}
}

// The embedded report template carries one data slot and the three.js
// script the 3D scenes need.
func TestReportTemplate(t *testing.T) {
	s := string(reportTemplate)
	if strings.Count(s, "__DATA__") != 1 || !strings.Contains(s, "three.js/r128/three.min.js") || !strings.HasPrefix(s, "<!doctype html>") {
		t.Error("report template is missing its data slot, three.js or its doctype")
	}
}
