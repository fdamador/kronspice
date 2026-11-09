// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Eq. 8 is force-free and its field integral is the paper's B3.
func TestFieldAndB3(t *testing.T) {
	curl, div := CurlCheck(44.93)
	if curl > 1e-8 || div > 1e-8 {
		t.Errorf("curl %.1e div %.1e", curl, div)
	}
	al := 44.934
	if v := IM(al) * al * al * al; math.Abs(v-35.87) > 0.01 {
		t.Errorf("∫|M|²dV = %.4f/α³, paper 35.87/α³", v)
	}
}

// Appendix B as printed: E0 3.566e9 (from 35.87/α²), ω 7.920e11, k 2642,
// and with his k = 2686: ρ_e −32.63 C/m³, 8.532e17 electrons. With B3 as
// printed (α³) E0 is √α larger.
func TestAppendixB(t *testing.T) {
	al := J1Zero / 0.1
	p := AppendixB("printed", 35.87/(al*al), 2686)
	for _, c := range []struct {
		got, want, tol float64
		name           string
	}{
		{p.E0, 3.566e9, 1e-3, "E0"}, {p.Omega, 7.920e11, 1e-3, "ω"}, {p.K0, 2642, 1e-3, "k0"},
		{p.RhoE, -32.63, 1e-3, "ρ_e"}, {p.Electrons, 8.532e17, 1e-3, "electrons"},
	} {
		if math.Abs(c.got/c.want-1) > c.tol {
			t.Errorf("%s = %.5g, paper %.5g", c.name, c.got, c.want)
		}
	}
	if math.Abs(p.DispersionResidue) > 1e-12 {
		t.Errorf("eq. 49 vs cold-plasma dispersion %.1e", p.DispersionResidue)
	}
	q := AppendixB("B3", IM(al), 0)
	if r := q.E0 / p.E0; math.Abs(r-math.Sqrt(al)) > 0.01*r {
		t.Errorf("E0 ratio %.4f, want √α = %.4f", r, math.Sqrt(al))
	}
}

// The harmonic comb: zeros of j1, and v/c = ξ·kR.
func TestHarmonics(t *testing.T) {
	z := J1Zeros(15)
	want := []float64{4.493409, 7.725252, 10.904122, 14.066194}
	for i, w := range want {
		if math.Abs(z[i]-w) > 1e-5 {
			t.Errorf("x_%d = %.6f, want %.6f", i+1, z[i], w)
		}
	}
	for l, w := range map[int]float64{1: 4.4934, 2: 5.7635, 3: 6.9879} {
		if math.Abs(FirstZero(l)-w) > 1e-4 {
			t.Errorf("first zero of j_%d = %.5f", l, FirstZero(l))
		}
	}
	d := DesignFromC(0.12, 1e5, 0.01)
	if math.Abs(d.Beta-d.C*d.X) > 1e-12 || d.Beta < 1 {
		t.Errorf("v/c = %.3f (C·kR = %.3f); the paper's recipe should exceed c at ball-lightning size", d.Beta, d.C*d.X)
	}
}

// The large-kR lifetimes equal the complex root finder.
func TestLifetimes(t *testing.T) {
	for _, y := range []float64{268.6, 1000} {
		ne := Eps0 * Me * (y*y - J1Zero*J1Zero) * C0 * C0 / (0.01 * Qe * Qe)
		s := Sphere{R: 0.1, Wp: math.Sqrt(ne * Qe * Qe / (Me * Eps0))}
		te, tm := s.FindMode(true), s.FindMode(false)
		if !te.Found || !tm.Found || math.Abs(te.Lifetime/TauTE(0.1, y)-1) > 1e-3 || math.Abs(tm.Lifetime/TauTM(0.1, y)-1) > 5e-3 {
			t.Errorf("kR %.0f: TE %.4g vs %.4g, TM %.4g vs %.4g", y, te.Lifetime, TauTE(0.1, y), tm.Lifetime, TauTM(0.1, y))
		}
	}
	// a 12 cm plasmoid lasting 3 s needs more free electrons than fully stripped air has
	if _, _, ne := ForLifetime(0.12, 3); ne/Loschmidt < AirElectronsPerMolecule {
		t.Errorf("12 cm, 3 s: %.3g × Loschmidt", ne/Loschmidt)
	}
}

// kronspice: the Foster chain returns each mode's ω and Q; the ring-down
// gives Q and closes its ledger; the plasma-loaded line resonates at eq. 49.
func TestKronspice(t *testing.T) {
	al := J1Zero / 0.1
	e := AppendixB("printed", 35.87/(al*al), 2686)
	peaks, err := SpectrumAC(Spectrum(e.R, e.Wp, 2))
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range peaks {
		if p.FErr > 1e-9 || p.QErr > 1e-4 {
			t.Errorf("mode %.6g Hz: f err %.1e, Q err %.1e", p.FMode, p.FErr, p.QErr)
		}
	}
	rd, err := RunRingDown(peaks[0].Mode, 100)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(rd.QFit/rd.Q-1) > 2e-3 || rd.Mismatch > 1e-4 {
		t.Errorf("ring-down Q %.5g vs %.5g, ledger %.1e", rd.QFit, rd.Q, rd.Mismatch)
	}
	ld, err := RunLadder(1.5, 100)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ld.WKron/ld.WLattice-1) > 1e-7 || math.Abs(ld.WLattice/ld.WPaper-1) > 1e-4 {
		t.Errorf("ladder %.8f, lattice %.8f, eq. 49 %.8f", ld.WKron, ld.WLattice, ld.WPaper)
	}
}
