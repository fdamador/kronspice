// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func TestJacobi(t *testing.T) {
	a := [][]float64{{4, 1, 0.5}, {1, 3, 0.2}, {0.5, 0.2, 1}}
	vals, vecs := jacobiEigen(a)
	for c := 0; c < 3; c++ {
		for r := 0; r < 3; r++ {
			s := 0.0
			for k := 0; k < 3; k++ {
				s += a[r][k] * vecs[k][c]
			}
			if math.Abs(s-vals[c]*vecs[r][c]) > 1e-12 {
				t.Fatalf("eigenpair %d off by %g", c, s-vals[c]*vecs[r][c])
			}
		}
	}
	if !(vals[0] < vals[1] && vals[1] < vals[2]) {
		t.Errorf("not ascending: %v", vals)
	}
}

// The round bubble matches the independent Python reference
// (ref/void_shells_ref.py: ether 106.48697, light 285.33622 at ε = 300),
// lies between the pressure and wall limits of g, and the prepared shell
// relaxes to the very same state.
func TestBubbleAndShell(t *testing.T) {
	b, _ := relaxed("ball", EpsRef, []int{0}, Ball(5), Relax)
	if math.Abs(b.Ee/106.48697-1) > 1e-6 || math.Abs(b.EL/285.33622-1) > 1e-6 {
		t.Errorf("ball: E_e %g, E_L %g against the reference", b.Ee, b.EL)
	}
	if b.G0 <= 4.0/3 || b.G0 >= 1.5 || b.L.Voids != 1 || b.L.LightIn < 0.999 {
		t.Errorf("ball: g %g, voids %d, light in %g", b.G0, b.L.Voids, b.L.LightIn)
	}
	s, _ := relaxed("shell", EpsRef, []int{0}, ShellAmp(3, 8), Relax)
	if math.Abs((s.Ee+s.EL)/(b.Ee+b.EL)-1) > 1e-8 || s.L.Voids != 1 || s.L.CoreDens > 1e-6 {
		t.Errorf("shell did not dissolve into the bubble: M %g vs %g, voids %d, core %g", s.Ee+s.EL, b.Ee+b.EL, s.L.Voids, s.L.CoreDens)
	}
}

// Two photons: the onion's ether layer is squeezed out, one void left.
func TestOnionMerges(t *testing.T) {
	on, _ := relaxed("onion", EpsRef, []int{0, 1}, OnionAmp(3, 5.5, 9), Relax)
	if on.L.Voids != 1 || on.G0 <= 4.0/3 || on.G0 >= 1.5 {
		t.Errorf("onion: voids %d, g %g", on.L.Voids, on.G0)
	}
}

// Real time: energy is conserved and the bubble keeps its light.
func TestRealTimeBall(t *testing.T) {
	s := NewState(grid(), EpsRef, V0, []int{0}, Ball(5))
	e0, l0 := s.Energies()
	for k := 0; k < 300; k++ {
		s.Step(0.01)
	}
	ee, el := s.Energies()
	if d := (ee + el - e0 - l0) / (e0 + l0); math.Abs(d) > 1e-3 {
		t.Errorf("energy drift %g", d)
	}
	if l := s.Measure(); l.CoreDens > 0.01 || l.LightIn < 0.999 { // sound dips can count as extra small voids; the centre and the light are what matter
		t.Errorf("after t = 3: %+v", l)
	}
}
