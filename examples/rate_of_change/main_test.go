// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// TestCircuitTwin: the water cell built as C∞ ∥ (R–C) and solved by
// kronspice AC equals jω(ε0A/d)ε(ω) to 1e-9 (CLAUDE.md §3a rule 4 tolerance).
func TestCircuitTwin(t *testing.T) {
	worst, err := CircuitTwin(em.WaterDebye25, []float64{1e3, 1e6, 1e9, 19.2e9, 1e11})
	if err != nil {
		t.Fatal(err)
	}
	if worst > 1e-9 {
		t.Errorf("circuit vs ε(ω): %.3g", worst)
	}
	// A Lorentz term too (series R–L–C branch).
	m := em.Dispersive{Name: "Lorentz test", EpsInf: 2, Lorentz: []em.LorentzTerm{{DeltaEps: 3, Omega0: 2 * math.Pi * 1e9, Gamma: 1e8}}}
	worst, err = CircuitTwin(m, []float64{1e8, 1e9, 3e9})
	if err != nil {
		t.Fatal(err)
	}
	if worst > 1e-9 {
		t.Errorf("Lorentz circuit vs ε(ω): %.3g", worst)
	}
}

// TestSlap: a slow ramp sees ε(0); a ramp much faster than τ sees ε∞
// plus the small part that relaxes during the ramp (analytically
// ε∞ + Δε·τ/T·(1 − e^(−T/τ)) at the end of a ramp of length T); every
// ramp settles to ε(0).
func TestSlap(t *testing.T) {
	w := em.WaterDebye25
	tau, de := w.Debye[0].Tau, w.Debye[0].DeltaEps
	for _, T := range []float64{1e-6, 1e-10, 1e-12, 1e-13} {
		a, b := Slap(w, T)
		want := w.EpsInf + de*(1-tau/T*(-math.Expm1(-T/tau)))
		if math.Abs(a-want) > 1e-6*want || math.Abs(b-w.Static()) > 1e-6 {
			t.Errorf("rise %g: %.6g (want %.6g), settled %.6g", T, a, want, b)
		}
	}
}

func TestRuns(t *testing.T) { main() }
