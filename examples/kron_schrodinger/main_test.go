// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The H atom and H₂⁺ on Kron's network converge to the exact energies
// (−0.5 Ha; −1.1026342 Ha at R = 2 a₀) as the grid is refined. The
// finite-difference energies lie below the exact ones (the grid samples
// the Coulomb wells too deeply; a difference scheme is not variational) and
// approach them with the order set by the cusps, about 1.5 on fine grids
// (main, h down to 0.05) and nearer 1 on these coarse test grids, so the
// extrapolation is checked to 5 mHa here and to 0.3 mHa in run_output.txt.
func TestKronNetworkEnergies(t *testing.T) {
	var eh, e2 []float64
	for _, h := range []float64{0.4, 0.2, 0.1} {
		e, _, err := Hydrogen(h).Ground(-0.7)
		if err != nil {
			t.Fatal(err)
		}
		eh = append(eh, e)
		e, _, err = H2Plus(h, 2).Ground(-1.4)
		if err != nil {
			t.Fatal(err)
		}
		e2 = append(e2, e)
	}
	x, p := Extrapolated(eh[0], eh[1], eh[2])
	if math.Abs(x+0.5) > 5e-3 || p < 0.8 || p > 2.2 {
		t.Errorf("H atom: extrapolated %.6f (order %.2f), exact −0.5; runs %v", x, p, eh)
	}
	x, p = Extrapolated(e2[0], e2[1], e2[2])
	if math.Abs(x+1.1026342) > 5e-3 || p < 0.8 || p > 2.2 {
		t.Errorf("H2+: extrapolated %.6f (order %.2f), exact −1.1026342; runs %v", x, p, e2)
	}
	for i := 1; i < 3; i++ {
		if math.Abs(eh[i]+0.5) >= math.Abs(eh[i-1]+0.5) || math.Abs(e2[i]+1.1026342) >= math.Abs(e2[i-1]+1.1026342) {
			t.Errorf("error does not fall with h: H %v, H2+ %v", eh, e2)
		}
	}
}

// The network as a kronspice LC circuit resonates exactly at the driven
// network's eigenvalue (same equations, same grid).
func TestKronCircuitTwin(t *testing.T) {
	g := H2Plus(0.4, 2)
	eg, _, err := g.Ground(-1.4)
	if err != nil {
		t.Fatal(err)
	}
	v0 := 0.0
	for i := 0; i < g.Nr; i++ {
		for j := 0; j < g.Nz; j++ {
			v0 = math.Max(v0, -g.potential(i, j))
		}
	}
	ec, err := g.CircuitGround(v0+1, eg-0.3, eg+0.3)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ec-eg) > 1e-7 {
		t.Errorf("circuit resonance %.10f, eigenvalue %.10f", ec, eg)
	}
}

// CubicMin finds the minimum of a known cubic.
func TestCubicMin(t *testing.T) {
	f := func(r float64) float64 { return 0.3*(r-2)*(r-2) - 0.05*(r-2)*(r-2)*(r-2) - 0.6 }
	var x, y []float64
	for r := 1.5; r <= 2.6; r += 0.1 {
		x, y = append(x, r), append(y, f(r))
	}
	r, e := CubicMin(x, y)
	if math.Abs(r-2) > 1e-9 || math.Abs(e+0.6) > 1e-12 {
		t.Errorf("min at %.12f, %.12f", r, e)
	}
}
