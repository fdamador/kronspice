// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func atomFields() []AtomField {
	var fs []AtomField
	for _, s := range States() {
		fs = append(fs, NewAtomField(s))
	}
	return fs
}

// Every level is normalized and its dipole moment is the Landé value
// −g_j m_j μB (projection theorem).
func TestNormAndDipole(t *testing.T) {
	for _, f := range atomFields() {
		if n := Norm(f.S); math.Abs(n-1) > 1e-6 {
			t.Errorf("%s: ∫|ψ|² = %.8f", f.S.Name, n)
		}
		if m, want := f.Ser.Moment(1), -f.S.GJ()*f.S.M; math.Abs(m-want) > 1e-6 {
			t.Errorf("%s: dipole %.8f μB, Landé %.8f", f.S.Name, m, want)
		}
	}
}

// Selection rules: no even-l moments (parity), no octupole for j = ½
// (2j ≥ k), and the j = 3⁄2 octupoles in the Wigner–Eckart ratio
// m(5m² − 3j(j+1) + 1): (½)/(3⁄2) = −3.
func TestOctupoleSelection(t *testing.T) {
	fs := atomFields()
	for _, f := range fs {
		for _, l := range []int{2, 4} {
			if m := f.Ser.Moment(l); math.Abs(m) > 1e-8 {
				t.Errorf("%s: m%d = %.3g, want 0 (parity)", f.S.Name, l, m)
			}
		}
		if f.S.Jq == 0.5 {
			if m := f.Ser.Moment(3); math.Abs(m) > 1e-8 {
				t.Errorf("%s: octupole %.3g, want 0 (j = ½)", f.S.Name, m)
			}
		}
	}
	if r := fs[3].Ser.Moment(3) / fs[2].Ser.Moment(3); math.Abs(r+3) > 1e-6 {
		t.Errorf("octupole ratio %.8f, want −3", r)
	}
}

// The series field equals a direct current-loop sum outside the cloud.
func TestSeriesTwin(t *testing.T) {
	if e := TwinAtom(NewAtomField(States()[3])).Err; e > 1e-3 {
		t.Errorf("series vs loops: %.2g of |B|", e)
	}
}

// Pole counts: the ground state shows two zones (one pair) at every radius;
// 2p3⁄2 m = ½ shows four out to r_c = √(3 m3/m1) and two beyond.
func TestZones(t *testing.T) {
	fs := atomFields()
	for _, r := range []float64{1, 3, 10, 40} {
		if z := Zones(func(th float64) float64 { return fs[0].Br(r, th) }, 720); z != 2 {
			t.Errorf("1s at %g a0: %d zones", r, z)
		}
	}
	f := fs[3]
	m1, m3 := f.Moments()
	rc := OctupoleRadius(m1, m3)
	if math.Abs(rc-15.61) > 0.05 {
		t.Errorf("r_c = %.3f a0", rc)
	}
	for _, c := range []struct {
		r float64
		z int
	}{{2, 4}, {8, 4}, {14, 4}, {18, 2}, {40, 2}} {
		if z := Zones(func(th float64) float64 { return f.Br(c.r, th) }, 720); z != c.z {
			t.Errorf("2p3/2 m=1/2 at %g a0: %d zones, want %d", c.r, z, c.z)
		}
	}
}

// E1 against tabulated values (Abramowitz & Stegun Table 5.1).
func TestE1(t *testing.T) {
	for _, c := range [][2]float64{{0.5, 0.5597735947761608}, {1, 0.21938393439552029}, {3, 0.013048381094197037}} {
		if v := E1(c[0]); math.Abs(v/c[1]-1) > 1e-12 {
			t.Errorf("E1(%g) = %.16g, want %.16g", c[0], v, c[1])
		}
	}
}

// Heitler–London: the singlet binds at R_e = 1.64 a0 with D_e = 0.1160 Eh
// (3.16 eV); the triplet lies above the separated atoms at every R.
func TestHeitlerLondon(t *testing.T) {
	re, de := HLMinimum()
	if math.Abs(re-1.6425) > 0.002 || math.Abs(de-0.1160) > 0.0005 {
		t.Errorf("R_e %.4f a0, D_e %.5f Eh", re, de)
	}
	for r := 0.8; r < 8; r += 0.2 {
		if e := HLEnergy(r, -1) + 1; e <= 0 {
			t.Errorf("triplet binds at %.1f a0 (%.4g Eh)", r, e)
		}
	}
}

// Spin density: zero for the singlet, two electrons' worth for the triplet.
func TestSpinDensity(t *testing.T) {
	ms := Molecules(1.4011)
	if s := SpinIntegral(ms[0]); s != 0 {
		t.Errorf("singlet ∫ spin density %.3g", s)
	}
	if s := SpinIntegral(ms[2]); math.Abs(s-2) > 1e-3 {
		t.Errorf("triplet ∫ spin density %.6f", s)
	}
}

// The triplet's far field is the dipole of its −g_s μB moment (corrections
// fall as 1/r²); ortho protons alone show four zones at 1 a0; para has no
// moment.
func TestMoleculeField(t *testing.T) {
	ms := Molecules(1.4011)
	if ms[0].Moment() != 0 {
		t.Error("singlet + para has a moment")
	}
	f := NewMolField(ms[2], 0.1, 10, 3)
	for _, z := range []float64{20, 40} {
		_, bz, _ := f.Field(0, z)
		want := 2 * ms[2].Moment() / (z * z * z)
		if d := math.Abs(bz/want - 1); d > 10/(z*z) {
			t.Errorf("triplet axis field at %g a0 off the dipole by %.3g", z, d)
		}
	}
	o := NewMolField(ms[1], 0.1, 10, 3)
	br := func(th float64) float64 {
		bp, bz, _ := o.Field(math.Sin(th), math.Cos(th))
		return bp*math.Sin(th) + bz*math.Cos(th)
	}
	if z := Zones(br, 720); z != 4 {
		t.Errorf("ortho protons at 1 a0: %d zones, want 4", z)
	}
}

// Part 4: kronspice resonances of Kron's radial network give Lyman α and the
// 2p fine structure; the model's spin density and moment give 21 cm and ESR;
// the LC circuit gives the Heitler–London vibration with the measured
// isotope ratio.
func TestSpectrum(t *testing.T) {
	var lv [3]Level
	specs := []struct {
		l  int
		j  float64
		so bool
	}{{0, 0.5, false}, {1, 0.5, true}, {1, 1.5, true}}
	errs := make([]error, 3)
	parallel(3, func(i int) {
		lv[i], errs[i] = FindLevel("t", specs[i].l, specs[i].j, specs[i].so, 0.1, 80)
	})
	for _, e := range errs {
		if e != nil {
			t.Fatal(e)
		}
	}
	for _, l := range lv {
		if d := math.Abs(l.Fine - l.Matrix); d > 1e-12 {
			t.Errorf("l = %d: kronspice resonance vs matrix %.2g Ha", l.L, d)
		}
	}
	lya := ((lv[1].Rich+lv[2].Rich)/2 - lv[0].Rich) * MuRedH * HaHz
	if d := math.Abs(lya/measLya - 1); d > 3e-5 {
		t.Errorf("Lyman α off by %.2g", d)
	}
	fs := (lv[2].Rich - lv[1].Rich) * MuRedH * HaHz
	if d := math.Abs(fs/measFS - 1); d > 1e-4 {
		t.Errorf("fine structure %.1f MHz", fs/1e6)
	}
	s := States()[0]
	if d := math.Abs(HyperfineHz(s)/meas21 - 1); d > 1e-4 {
		t.Errorf("21 cm off by %.2g", d)
	}
	if d := math.Abs(ESRHz(NewAtomField(s))/measESR - 1); d > 3e-5 {
		t.Errorf("ESR off by %.2g", d)
	}
	_, k := HLCurvature()
	h2, _, err := H2Vib(k, MpMe, MpMe)
	if err != nil {
		t.Fatal(err)
	}
	d2, _, err := H2Vib(k, MdMe, MdMe)
	if err != nil {
		t.Fatal(err)
	}
	if want := math.Sqrt(2*k/MpMe) * HaCm; math.Abs(h2/want-1) > 1e-9 {
		t.Errorf("H₂ resonance %.4f cm⁻¹, √(k/μ) %.4f", h2, want)
	}
	if r := d2 / h2; math.Abs(r/(measWe[2]/measWe[0])-1) > 2e-3 {
		t.Errorf("D₂/H₂ %.5f", r)
	}
}
