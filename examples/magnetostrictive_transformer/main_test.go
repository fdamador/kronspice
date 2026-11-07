// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The declared materials land in their measured ranges (Etrema data:
// Terfenol-D μr 3–10; Galfenol d33 15–30 nm/A, μr 75–100; composite d33
// ≈ 5.5 nm/A, k33 ≈ 0.34).
func TestMaterials(t *testing.T) {
	d, mur, _ := Terfenol.Small()
	if mur < 3 || mur > 10 || d < 8e-9 || d > 20e-9 {
		t.Errorf("Terfenol d %.3g, μr %.3g", d, mur)
	}
	d, mur, _ = Galfenol.Small()
	if mur < 75 || mur > 100 || d < 15e-9 || d > 30e-9 {
		t.Errorf("Galfenol d %.3g, μr %.3g", d, mur)
	}
	d, _, k := TerfenolComposite.Small()
	if math.Abs(d-5.5e-9) > 0.5e-9 || math.Abs(k-0.34) > 0.05 {
		t.Errorf("composite d %.3g, k %.3g", d, k)
	}
}

// Small-signal reciprocity of the coupled core (an energy-consistent
// magnetoelastic model must be reciprocal), the ledger, and the products
// that appear only with both the drive and the sound.
func TestRodAndRing(t *testing.T) {
	cs := Cores()
	for _, c := range []Core{cs[0], cs[2]} {
		k, err := RunCase(c)
		if err != nil {
			t.Fatal(err)
		}
		if rc := k.Reciprocity(); rc > 1e-8 {
			t.Errorf("%s %s: reciprocity %.2g", c.Mat.Name, c.Shape, rc)
		}
		for _, r := range []Result{k.Fwd, k.Rev} {
			if math.Abs(r.Residual) > 2e-2*r.In {
				t.Errorf("%s %s: ledger %.2g of %.3g J", c.Mat.Name, c.Shape, r.Residual, r.In)
			}
		}
		// Products: present with sound, absent without it.
		with := k.Fwd.Amp("fa+fe") + k.Fwd.Amp("fb+fe")
		without := k.NoSound.Amp("fa+fe") + k.NoSound.Amp("fb+fe")
		if with < 1e-5 || without > 1e-3*with {
			t.Errorf("%s %s: products with sound %.3g A, without %.3g A", c.Mat.Name, c.Shape, with, without)
		}
	}
}
