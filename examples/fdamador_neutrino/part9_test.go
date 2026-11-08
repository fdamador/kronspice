// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// besselJ2 by its power series (|y| ≤ 4 here).
func besselJ2(y float64) float64 {
	var s, term float64
	term = y * y / 8 // (y/2)²/2!
	for k := 0; k < 40; k++ {
		s += term
		term *= -(y * y / 4) / (float64(k+1) * float64(k+3))
	}
	return s
}

// The retarded kernel's self term is the exact radiation resistance of a
// uniform-current loop, η0(π/2)·ka·∫0^{2ka} J2(y) dy (Balanis eq. 5-63).
func TestRetardedSelfIsBalanis(t *testing.T) {
	e := DefaultERing()
	a := NewArray(e, Coaxial, e.A, 1)
	for _, ka := range []float64{0.05, 0.432, 1, 2} {
		w := ka * C0 / e.A
		const n = 4000
		h := 2 * ka / n
		var in float64
		for k := 0; k < n; k++ {
			in += besselJ2((float64(k)+0.5)*h) * h
		}
		want := Eta0 * math.Pi / 2 * ka * in
		if got := a.R11(w); math.Abs(got/want-1) > 1e-6 {
			t.Errorf("ka %g: R11 %.9g, Balanis %.9g", ka, got, want)
		}
	}
}

// Far apart, the exact mutual radiation resistance tends to the dipole
// form (Lehmberg), both geometries.
func TestMutualRadiationDipoleLimit(t *testing.T) {
	e := DefaultERing()
	w := e.Omega0()
	for _, g := range []Geom{Coaxial, Coplanar} {
		for _, dd := range []float64{3, 5} {
			a := NewArray(e, g, dd*e.A, 2)
			_, rm := a.Coupling(0, 1, w)
			dip := real(DipoleZ12(g, w*dd*e.A/C0))
			if math.Abs(rm/a.R11(w)-dip) > 0.01*math.Abs(dip) {
				t.Errorf("%s d = %ga: R_m/R11 %.6f, dipole %.6f", g, dd, rm/a.R11(w), dip)
			}
		}
	}
}

// The extrapolated Neumann polygons give the coaxial elliptic formula,
// and the side-by-side static M equals the trapezoidal double integral.
func TestStaticMutual(t *testing.T) {
	a := 0.01
	for _, d := range []float64{0.025, 0.04} {
		got := StaticM(a, Coplanar, d)
		const n = 1024
		h := 2 * math.Pi / n
		var s float64
		for i := 0; i < n; i++ {
			s1, c1 := math.Sincos(float64(i) * h)
			for j := 0; j < n; j++ {
				s2, c2 := math.Sincos(float64(j) * h)
				s += (c1*c2 + s1*s2) / math.Hypot(a*c1-d-a*c2, a*s1-a*s2)
			}
		}
		want := Mu0 / (4 * math.Pi) * a * a * h * h * s
		if math.Abs(got/want-1) > 1e-6 {
			t.Errorf("side by side d = %g: em.Neumann (extrapolated) %.9g, trapezoid %.9g", d, got, want)
		}
	}
}

// A pair's modes are (1, ±1) with the closed-form frequency and Q, and
// kronspice rings each down at that Q with a closed ledger; side by side
// (negative M, made with helper inductors) as well.
func TestPairModesKronspice(t *testing.T) {
	e := DefaultERing()
	for _, c := range []struct {
		g  Geom
		dd float64
	}{{Coaxial, 1}, {Coplanar, 2.5}} {
		a := NewArray(e, c.g, c.dd*e.A, 2)
		in, an := a.PairModes()
		for _, x := range []struct {
			m  AMode
			sg float64
		}{{in, 1}, {an, -1}} {
			f, q := a.PairAnalytic(x.sg, x.m.W)
			// the closed form ignores the damping's own frequency shift (1/8Q²)
			if math.Abs(f/x.m.F-1) > 1/(4*q*q) || math.Abs(q/x.m.Q-1) > 1e-3 {
				t.Errorf("%s sign %+g: f %.9g vs %.9g, Q %.6g vs %.6g", c.g, x.sg, x.m.F, f, x.m.Q, q)
			}
			r, err := a.RunMode(x.m, 30, 200)
			if err != nil {
				t.Fatal(err)
			}
			if math.Abs(r.Q/x.m.Q-1) > 1e-2 || math.Abs(r.Mismatch) > 1e-3 || r.MinRadP < 0 || r.LedgerRel > 1e-6 {
				t.Errorf("%s sign %+g: kronspice Q %.6g vs %.6g, books %.1e, min radiated %.1e, ledger %.1e",
					c.g, x.sg, r.Q, x.m.Q, r.Mismatch, r.MinRadP, r.LedgerRel)
			}
		}
		if c.g == Coaxial && !(an.Q > 10*in.Q) {
			t.Errorf("coaxial antiphase not subradiant: Q %.4g vs %.4g", an.Q, in.Q)
		}
	}
}

// Beats run at the splitting of the circuit as integrated.
func TestPairBeats(t *testing.T) {
	e := DefaultERing()
	bt, err := RunBeats(e, Coaxial, e.A, 30)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(bt.DfMeasured/bt.DfFrozen-1) > 5e-3 || bt.MaxTransfer < 0.9 {
		t.Errorf("Δf %.6g vs %.6g, transfer %.3f", bt.DfMeasured, bt.DfFrozen, bt.MaxTransfer)
	}
}

// Coaxial: in phase attracts, antiphase repels; kronspice's ⟨I1I2⟩ gives
// the analytic force.
func TestPairForceSigns(t *testing.T) {
	e := DefaultERing()
	f, err := Force(e, Coaxial, e.A, 4, 1)
	if err != nil {
		t.Fatal(err)
	}
	if !(f.FInA < 0 && f.FAntiA > 0) || math.Abs(f.FInK/f.FInA-1) > 0.02 || math.Abs(f.FAntiK/f.FAntiA-1) > 0.02 {
		t.Errorf("forces: in %.5f (kron %.5f), anti %.5f (kron %.5f)", f.FInA, f.FInK, f.FAntiA, f.FAntiK)
	}
}

// The chain's darkest mode gets darker with N; kronspice agrees at N = 8.
func TestChainDarkMode(t *testing.T) {
	e := DefaultERing()
	e.TanD = 0
	prev := 0.0
	for _, n := range []int{2, 4, 8} {
		a := NewArray(e, Coaxial, e.A, n)
		_, dk, _ := a.ChainModes()
		if dk.Q <= prev || dk.Dipole > 1e-9 {
			t.Errorf("N = %d: dark Q %.5g (previous %.5g), net moment %.1e", n, dk.Q, prev, dk.Dipole)
		}
		prev = dk.Q
		if n == 8 {
			r, err := a.RunMode(dk, 100, 100)
			if err != nil {
				t.Fatal(err)
			}
			if math.Abs(r.Q/dk.Q-1) > 1e-2 {
				t.Errorf("N = 8: kronspice Q %.6g vs %.6g", r.Q, dk.Q)
			}
		}
	}
}
