// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func below(t *testing.T, what string, got, limit float64) {
	t.Helper()
	if !(got <= limit) {
		t.Errorf("%s = %.3g, want ≤ %.1g", what, got, limit)
	}
}

// Step 1 on the circuit: the gauge-theory identities hold on Kron's
// mesh-frame solution to round-off, and the energy terms are the ledger's.
func TestAbelianCircuit(t *testing.T) {
	ab, err := RunAbelian(2 * math.Pi * 1e3)
	if err != nil {
		t.Fatal(err)
	}
	if ab.Meshes != 3 || ab.Branches != 6 || ab.Nodes != 4 {
		t.Fatalf("counts b/N/m = %d/%d/%d", ab.Branches, ab.Nodes, ab.Meshes)
	}
	below(t, "Curl∘Grad", ab.CurlGrad, 1e-15)
	below(t, "KVL (zero holonomy)", ab.KVL, 1e-14)
	below(t, "pure gauge", ab.Pure, 1e-14)
	below(t, "mesh EMF", ab.MeshEMF, 1e-14)
	below(t, "KCL", ab.KCL, 1e-14)
	below(t, "continuity", ab.Continuity, 1e-14)
	if ab.FluxQuanta[0] < 1e9 || ab.FluxQuanta[2] != 0 {
		t.Errorf("flux quanta %v: mesh 0 should link ~1e10 Φ₀, mesh 2 (no inductor) none", ab.FluxQuanta)
	}
	if ab.Magnetic <= 0 || ab.Electric <= 0 {
		t.Errorf("energy terms %g, %g", ab.Electric, ab.Magnetic)
	}
	below(t, "linear composition defect", ab.LinearCompositionDefect, 1e-14)
}

// Step 1 on the field: the U(1) field on the Yee grid is the FDTD solve.
func TestYeeTwin(t *testing.T) {
	tw, err := RunTwin(1500)
	if err != nil {
		t.Fatal(err)
	}
	below(t, "E edge by edge vs FDTD", tw.EdgeDefect, 1e-11)
	below(t, "energy vs FDTD ledger", tw.EnergyDefect, 1e-12)
	below(t, "energy drift", tw.EnergyDrift, 1e-12)
	below(t, "Gauss", tw.Gauss, 1e-12)
	below(t, "ω vs Yee", math.Abs(tw.Omega/tw.OmegaYee-1), 1e-5)
	below(t, "ω vs continuum", math.Abs(tw.Omega/tw.OmegaContinuum-1), 2e-3)
}

// Step 2: SU(3) loops do not commute, abelian ones do, and nothing in the
// circuit — linear or nonlinear — composes non-abelianly.
func TestNonAbelian(t *testing.T) {
	ab, err := RunAbelian(2 * math.Pi * 1e3)
	if err != nil {
		t.Fatal(err)
	}
	na, err := RunNonAbelian(ab, 50)
	if err != nil {
		t.Fatal(err)
	}
	if na.MeanDefect < 0.3 || na.MinDefect < 0.05 {
		t.Errorf("SU(3) composition defect mean %.3f min %.3f, expected order 1", na.MeanDefect, na.MinDefect)
	}
	below(t, "Cartan defect", na.CartanDefect, 1e-14)
	below(t, "Tr W invariance", na.TraceInvariance, 1e-12)
	below(t, "Gauss covariance", na.GaussCovariance, 1e-12)
	if na.GaussInvarianceBroken < 0.1 {
		t.Errorf("non-abelian Gauss generator changed by only %.3g under a gauge transformation", na.GaussInvarianceBroken)
	}
	below(t, "linear composition defect", na.LinearCompositionDefect, 1e-14)
	below(t, "nonlinear KCL", na.NonlinearKCL, 1e-12)
	below(t, "nonlinear KVL", na.NonlinearKVL, 1e-12)
	below(t, "nonlinear composition defect", na.NonlinearComposition, 1e-14)
	if na.DiodeCurrents[0] <= na.DiodeCurrents[1] || na.DiodeVoltages[0] <= na.DiodeVoltages[1] {
		t.Errorf("diodes: D1 %.3g A %.3f V, D2 %.3g A %.3f V; the lighter load should pass more", na.DiodeCurrents[0], na.DiodeVoltages[0], na.DiodeCurrents[1], na.DiodeVoltages[1])
	}
}

// The 2019 SI constants (CODATA 2018 exact values).
func TestMetrology(t *testing.T) {
	m := RunMetrology()
	if math.Abs(m.KJ/483597.8484e9-1) > 1e-9 {
		t.Errorf("K_J = %.10g Hz/V", m.KJ)
	}
	if math.Abs(m.RK-25812.80745) > 1e-5 {
		t.Errorf("R_K = %.10g Ω", m.RK)
	}
	if m.CopperZ != 29 {
		t.Errorf("copper Z = %d", m.CopperZ)
	}
}
