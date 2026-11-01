// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func TestHullCircuitMatchesWaveguideFormulas(t *testing.T) {
	h := Example()
	// Resonance = n·c/circumference times the ladder's lattice factor
	// sin(πn/N)/(πn/N) (0.997 at 24 sections per wavelength).
	x := math.Pi * float64(h.Mode) / float64(h.N)
	if f0 := float64(h.Mode) * C0 / h.Circumference() * math.Sin(x) / x; math.Abs(h.Freq()/f0-1) > 1e-9 {
		t.Errorf("f = %g, n·c/ℓ = %g", h.Freq(), f0)
	}
	for _, st := range []bool{false, true} {
		h.Standing, h.Vdc = st, 2500
		r, err := h.Solve()
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.Pin-h.Pin)/h.Pin > 1e-6 {
			t.Errorf("standing=%v: radiators %g W, want %g", st, r.Pin, h.Pin)
		}
		if math.Abs(r.PLoss+r.PFeed-r.Pin)/r.Pin > 2e-3 {
			t.Errorf("standing=%v: loss %g W + feed %g W ≠ input %g W", st, r.PLoss, r.PFeed, r.Pin)
		}
		// The solver's gmin conductances leak about 1 nA through the feed resistor
		// feed, so allow 1e-4 relative.
		if math.Abs(r.Vdc-2500) > 2500e-4 {
			t.Errorf("bias at the hull %g V, want 2500", r.Vdc)
		}
		want := r.Predicted
		if st {
			want *= math.Sqrt2 // power split into two counter-running waves, antinodes add
		}
		if math.Abs(r.VrfMax/want-1) > 0.01 {
			t.Errorf("standing=%v: RF peak %g V, formula %g V", st, r.VrfMax, want)
		}
		if math.Abs(r.PeakGap-(r.Vdc+r.VrfMax)) > 1e-9 {
			t.Errorf("peak %g ≠ bias + RF", r.PeakGap)
		}
		// Stored energy = Q·P/ω.
		if w := r.Q * r.Pin / (2 * math.Pi * r.F); math.Abs(r.Stored/w-1) > 0.01 {
			t.Errorf("stored %g J, Q·P/ω = %g J", r.Stored, w)
		}
	}
}

func TestFurmanYield(t *testing.T) {
	s := CopperAsReceived
	if y := s.Yield(s.Emax); math.Abs(y-s.DeltaMax) > 1e-12 {
		t.Errorf("δ(Emax) = %g, want δmax %g", y, s.DeltaMax)
	}
	if s.Yield(s.Emax*0.9) > s.DeltaMax || s.Yield(s.Emax*1.1) > s.DeltaMax {
		t.Error("yield not maximal at Emax")
	}
}

func TestMultipactorAndBias(t *testing.T) {
	h := Example()
	if g, _ := Multipactor(h.Gap, h.Freq(), 1000, 0, CopperAsReceived); g <= 1 {
		t.Errorf("1 kV unbiased: growth %g, expected a resonant zone", g)
	}
	if g, _ := Multipactor(h.Gap, h.Freq(), 1000, 500, CopperAsReceived); g >= 1 {
		t.Errorf("1 kV with 500 V bias: growth %g, expected suppression", g)
	}
	if g, _ := Multipactor(h.Gap, h.Freq(), 10, 0, CopperAsReceived); g >= 1 {
		t.Errorf("10 V: growth %g, below first crossover", g)
	}
}

func TestKilpatrick(t *testing.T) {
	for _, f := range []float64{10e6, 352e6, 3e9} {
		e := KilpatrickField(f) / 1e6
		if got := 1.64 * e * e * math.Exp(-8.5/e); math.Abs(got-f/1e6)/(f/1e6) > 1e-9 {
			t.Errorf("f = %g: inverse gives %g MHz", f, got)
		}
	}
}
