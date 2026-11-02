// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// TestBohr: Puthoff's eq. 20 with the Coulomb orbit gives the Bohr radius
// and −13.6057 eV, where eqs. 18 and 19 are equal.
func TestBohr(t *testing.T) {
	r, w, ev := Orbit()
	if math.Abs(r/A0-1) > 1e-8 || math.Abs(-ev/Ry-1) > 1e-8 {
		t.Errorf("r = %g a0, E = %g eV", r/A0, ev)
	}
	if math.Abs(Pabs(w)/Prad(r, w)-1) > 1e-12 {
		t.Errorf("absorbed %g, radiated %g", Pabs(w), Prad(r, w))
	}
	// absorbed/radiated = ħ/(mωr²) ∝ r^(−½) on Coulomb orbits
	for _, f := range []float64{0.5, 2} {
		rr := f * A0
		ww := Kepler(rr)
		if got := Pabs(ww) / Prad(rr, ww); math.Abs(got-1/math.Sqrt(f)) > 1e-9 {
			t.Errorf("ratio at %g a0: %g", f, got)
		}
	}
	if tau := Relax(); tau <= 0 || math.Abs(tau/9.337e-11-1) > 1e-3 {
		t.Errorf("relaxation %g s", tau)
	}
	if c := Collapse(); math.Abs(c/1.556e-11-1) > 1e-3 {
		t.Errorf("collapse %g s", c)
	}
}

// TestOscillator: in the simulated zero-point field the oscillator holds
// the expected energy (≈ ½ħω0) and radiates what it absorbs.
func TestOscillator(t *testing.T) {
	o := Osc{Gamma: 0.01, Band: 2, N: 1 << 18, Dt: 0.05, Runs: 8, Seed: 3}
	r := o.Run()
	if math.Abs(r.Energy-r.ExpEnergy) > 4*r.Spread || math.Abs(r.ExpEnergy-0.5) > 0.01 {
		t.Errorf("energy %g ± %g, expected %g", r.Energy, r.Spread, r.ExpEnergy)
	}
	if math.Abs(r.Absorbed-r.Radiated) > 0.02*r.Radiated || math.Abs(r.Absorbed-r.ExpPower) > 0.1*r.ExpPower {
		t.Errorf("absorbed %g, radiated %g, expected %g", r.Absorbed, r.Radiated, r.ExpPower)
	}
}

// TestTank: kronspice's LC tank in its resistor's zero-point noise holds
// the expected energy and its ledger closes.
func TestTank(t *testing.T) {
	k := Tank{R: 0.01, N: 1 << 14, Dt: 0.25, Step: 0.02, Band: 2, Runs: 4, Seed: 5}
	r, err := k.Run()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.Energy-r.ExpEnergy) > 5*r.Spread+0.01 {
		t.Errorf("energy %g ± %g, expected %g", r.Energy, r.Spread, r.ExpEnergy)
	}
	if r.Mismatch > 0.01 || math.Abs(r.Supplied-r.Heat-r.Stored) > 0.01*r.Supplied {
		t.Errorf("ledger: supplied %g, heat %g, stored %g, mismatch %g", r.Supplied, r.Heat, r.Stored, r.Mismatch)
	}
}

// TestFFT: the inverse FFT of one mode is a cosine of the right amplitude.
func TestFFT(t *testing.T) {
	n := 64
	a := make([]complex128, n)
	a[3], a[n-3] = 0.5, 0.5
	fft(a, true)
	for i := range a {
		if want := math.Cos(2 * math.Pi * 3 * float64(i) / float64(n)); math.Abs(real(a[i])-want) > 1e-12 {
			t.Fatalf("sample %d: %g, want %g", i, real(a[i]), want)
		}
	}
}
