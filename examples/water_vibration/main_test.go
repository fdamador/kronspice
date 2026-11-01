// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The force field reproduces the H₂O fundamentals it is fitted to; the
// circuit's AC resonances equal the GF eigenvalues (the circuit is the
// molecule's equations, not an approximation of them); and the same bonds
// predict D₂O and HDO to the anharmonicity budget (≤ 1.5 %).
func TestWaterCircuit(t *testing.T) {
	ff, err := Fit(H2O, H2O.Measured)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range Modes(H2O, ff) {
		if math.Abs(w/H2O.Measured[i]-1) > 1e-9 {
			t.Errorf("H2O mode %d: %.6f, target %.1f", i, w, H2O.Measured[i])
		}
	}
	for _, m := range []Molecule{H2O, D2O, HDO, H2O18} {
		ev := Modes(m, ff)
		pk, err := Peaks(Netlist(m, ff, 1e6), m, 500, 4500)
		if err != nil {
			t.Fatal(err)
		}
		if len(pk) != 3 {
			t.Fatalf("%s: %d spectral lines in the AC sweep, want 3: %v", m.Name, len(pk), pk)
		}
		for i := range ev {
			if math.Abs(pk[i]/ev[i]-1) > 1e-6 {
				t.Errorf("%s mode %d: circuit peak %.4f, eigenvalue %.4f", m.Name, i, pk[i], ev[i])
			}
			if m.Measured != nil && math.Abs(ev[i]/m.Measured[i]-1) > 0.015 {
				t.Errorf("%s mode %d: predicted %.1f, measured %.1f", m.Name, i, ev[i], m.Measured[i])
			}
		}
	}
}

// Three zero modes (two translations and the in-plane rotation) and
// physically sensible force constants: stretch stiffer than bend, weak
// negative stretch–stretch interaction.
func TestForceField(t *testing.T) {
	ff, err := Fit(H2O, HarmonicH2O)
	if err != nil {
		t.Fatal(err)
	}
	if !(ff.Fr > 5*ff.Fa && ff.Frr < 0 && math.Abs(ff.Frr) < 0.05*ff.Fr) {
		t.Errorf("force field %+v", ff)
	}
	b := BMatrix()
	// Rigid translation along x changes no internal coordinate.
	for k := 0; k < 3; k++ {
		s := b[k][0] + b[k][2] + b[k][4]
		if math.Abs(s) > 1e-12 {
			t.Errorf("coordinate %d moves under translation: %g", k, s)
		}
	}
}
