// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Kron's radial network gives the Coulomb s levels −Z²/2 and −Z²/8 Ha
// (h → 0), for Z = 1 and for the toroidal-photon charge.
func TestNetworkCoulomb(t *testing.T) {
	for _, z := range []float64{1, WvdM().Q} {
		e1, e2 := Extrapolated(z)
		if math.Abs(e1/(-z*z/2)-1) > 1e-7 || math.Abs(e2/(-z*z/8)-1) > 1e-7 {
			t.Errorf("Z %.6f: E1 %.10f, E2 %.10f", z, e1, e2)
		}
	}
}

// The network built in kronspice resonates at the matrix's ground state.
func TestNetworkTwin(t *testing.T) {
	tw, err := Twin(1, 0.1, 30)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(tw.Circuit-tw.Matrix) > 1e-9 {
		t.Errorf("circuit %.12f Ha, matrix %.12f Ha", tw.Circuit, tw.Matrix)
	}
}

// The H23 electron: charge 0.910318 e and g = 2(1 + α′/2π) (Williamson &
// van der Mark eqs. 5, 8, 12; examples/photon_electron).
func TestToroidalElectron(t *testing.T) {
	e := WvdM()
	if math.Abs(e.Q-0.910318) > 1e-6 {
		t.Errorf("charge %.7f e", e.Q)
	}
	if ap := e.Q * e.Q * AlphaFS; math.Abs(e.G-2*(1+ap/(2*math.Pi))) > 1e-15 {
		t.Errorf("g %.12f", e.G)
	}
	if r := e.R * 1e15; math.Abs(r-193.08) > 0.01 {
		t.Errorf("ring radius %.3f fm", r)
	}
}

// The standard atom (measured electron and proton) reproduces measured
// hydrogen: ionization and Lyman-α to the fine structure (1e-4), the
// 21 cm line through Fermi's contact term to 1e-4.
func TestStandardHydrogen(t *testing.T) {
	a := Atom{Measured(), Protons()[3]}
	if d := -a.Level(1)/E/IonEV - 1; math.Abs(d) > 1e-4 {
		t.Errorf("ionization off by %.2g", d)
	}
	if d := a.LymanAlpha()/LymanANm - 1; math.Abs(d) > 1e-4 {
		t.Errorf("Lyman-α off by %.2g", d)
	}
	hf, err := HyperfineCircuit(a, 8)
	if err != nil {
		t.Fatal(err)
	}
	if d := hf.Nu/HFMeasHz - 1; math.Abs(d) > 1e-4 {
		t.Errorf("21 cm off by %.2g", d)
	}
}

// The hyperfine circuit: its interaction energy is Fermi's contact energy
// μ_p·(2/3)μ0 M(0) (to the proton ring's finite size), and the ledger's
// stored energy equals ½ΣLI² + ΣMII.
func TestHyperfineCircuit(t *testing.T) {
	for _, a := range []Atom{{WvdM(), Protons()[0]}, {Measured(), Protons()[3]}} {
		hf, err := HyperfineCircuit(a, 8)
		if err != nil {
			t.Fatal(err)
		}
		if d := hf.Interaction/hf.Analytic - 1; math.Abs(d) > 1e-4 {
			t.Errorf("%s: interaction off by %.2g", a.E.Name, d)
		}
		if math.Abs(hf.Stored) > 1e-12 || hf.Residual > 1e-9 {
			t.Errorf("%s: stored %.2g, residual %.2g", a.E.Name, hf.Stored, hf.Residual)
		}
		if d := hf.Analytic*4/H/a.Hyperfine() - 1; math.Abs(d) > 1e-12 {
			t.Errorf("%s: 4U/h against Fermi's formula %.2g", a.E.Name, d)
		}
	}
}

// The light atom with the toroidal-photon charge is not neutral, by about
// 0.09 e, against a measured bound of 10⁻²¹ e; its lines move by Z².
func TestLightAtom(t *testing.T) {
	a := Atom{WvdM(), Protons()[0]}
	if q := a.Charge(); q < 0.08 || q < 1e19*Neutrality {
		t.Errorf("charge %.4f e", q)
	}
	ref := Atom{Measured(), Protons()[3]}
	if r := a.Level(1) / ref.Level(1); math.Abs(r-a.Z()*a.Z()) > 1e-12 {
		t.Errorf("level ratio %.9f, Z² %.9f", r, a.Z()*a.Z())
	}
}

// The ether knot sized to the electron ring rings at about 1.74 MeV.
func TestKnotExcitation(t *testing.T) {
	k := Knot(WvdM().R)
	if math.Abs(k.ExciteMeV-1.737) > 0.002 {
		t.Errorf("excitation %.4f MeV", k.ExciteMeV)
	}
	if e := 547.5 * math.Sqrt(k.K2*k.K4) / MeV; math.Abs(e/MeMeV-1) > 1e-12 {
		t.Errorf("knot energy %.6f MeV", e)
	}
}

// The frequency response of the whole network: its poles equal the network
// matrix's eigenvalues, extrapolate to −Z²/2n², and their differences give
// the Lyman and Balmer lines within the 0.01 nm the measured values carry
// (q = e) and 1/Z² longer for the toroidal-photon charge.
func TestSpectrum(t *testing.T) {
	hyd, err := RunSpectrum(1, 2000)
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k < NLevels; k++ {
		n := float64(k + 1)
		if math.Abs(hyd.Circuit[k]-hyd.Matrix[k]) > 1e-9 {
			t.Errorf("n %d: pole %.12f, matrix %.12f", k+1, hyd.Circuit[k], hyd.Matrix[k])
		}
		if d := hyd.Extrap[k]/(-1/(2*n*n)) - 1; math.Abs(d) > 5e-6 {
			t.Errorf("n %d: extrapolated level off by %.2g", k+1, d)
		}
	}
	for _, l := range hyd.Lines() {
		if math.Abs(l.Circuit-l.MeasuredNm) > 0.0105 {
			t.Errorf("%s: %.4f nm against %.2f", l.Name, l.Circuit, l.MeasuredNm)
		}
	}
	light, err := RunSpectrum(WvdM().Q, 2000)
	if err != nil {
		t.Fatal(err)
	}
	z := WvdM().Q
	for i, l := range light.Lines() {
		if r := l.Circuit / hyd.Lines()[i].Circuit; math.Abs(r*z*z-1) > 1e-5 {
			t.Errorf("%s: ratio %.6f, 1/Z² %.6f", l.Name, r, 1/(z*z))
		}
	}
}

// Hydrogen under pressure (Michels, de Boer & Bijl 1937): at a wall of
// 2 a0 the ground state is exactly −1/8 Ha; it reaches zero between 1.83 and
// 1.84 a0; compression splits 2s and 2p with 2p lower; the free atom has
// 2s = 2p = −1/8.
func TestConfined(t *testing.T) {
	if e := ConfinedLevel(1, 0, 0, 2); math.Abs(e+0.125) > 1e-6 {
		t.Errorf("E_1s at 2 a0: %.9f", e)
	}
	if lo, hi := ConfinedLevel(1, 0, 0, 1.84), ConfinedLevel(1, 0, 0, 1.83); !(lo < 0 && hi > 0) {
		t.Errorf("zero crossing not in (1.83, 1.84): %.6f, %.6f", hi, lo)
	}
	c := Compress(4)
	if !(c.E2p < c.E2s) || c.PressureGPa <= 0 || c.LymanNm >= 121.57 {
		t.Errorf("r_c = 4: 2s %.4f, 2p %.4f, P %.3g GPa, Lyα %.2f nm", c.E2s, c.E2p, c.PressureGPa, c.LymanNm)
	}
	f := Compress(60)
	if math.Abs(f.E2s-f.E2p) > 1e-6 || math.Abs(f.E2p+0.125) > 1e-6 {
		t.Errorf("free atom: 2s %.8f, 2p %.8f", f.E2s, f.E2p)
	}
}

// The p network (l = 1) in kronspice resonates at its matrix's ground
// state, the 2p level.
func TestTwin2p(t *testing.T) {
	tw, err := TwinL(1, 1, 0.1, 40)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(tw.Circuit-tw.Matrix) > 1e-9 || math.Abs(tw.Matrix+0.125) > 1e-3 {
		t.Errorf("circuit %.12f, matrix %.12f", tw.Circuit, tw.Matrix)
	}
}
