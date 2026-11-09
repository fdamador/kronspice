// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The A_{M,K} starts carry degree K (Hopf charge M·K) on the grid.
func TestStartDegree(t *testing.T) {
	for _, c := range []struct{ M, K int }{{1, 1}, {1, 2}, {2, 1}} {
		a := NewAxi(8, 8, 0.1, 1, 1, 1, c.M, c.K)
		if d := a.Degree(); math.Abs(d-float64(c.K)) > 0.08*float64(c.K) {
			t.Errorf("A%d%d degree %.4f, want %d", c.M, c.K, d, c.K)
		}
	}
}

// The energy gradient (with the winding M) equals finite differences.
func TestGradientWinding(t *testing.T) {
	for _, M := range []int{1, 2} {
		a := NewAxi(3, 3, 0.25, 1, 1, 1.3, M, 1)
		var g [3][]float64
		for c := 0; c < 3; c++ {
			g[c] = make([]float64, a.NR*a.NZ)
		}
		a.Energy(&g)
		worst := 0.0
		for _, n := range []int{a.at(0, a.NZ/2), a.at(3, a.NZ/2+1), a.at(5, a.NZ/2-2)} {
			for c := 0; c < 3; c++ {
				const eps = 1e-6
				a.m[c][n] += eps
				ep := a.Total()
				a.m[c][n] -= 2 * eps
				em := a.Total()
				a.m[c][n] += eps
				fd := (ep - em) / (2 * eps)
				worst = math.Max(worst, math.Abs(fd-g[c][n])/(1+math.Abs(fd)))
			}
		}
		if worst > 1e-5 {
			t.Errorf("M = %d: gradient error %.2e", M, worst)
		}
	}
}

// Step 0 and step 3 numbers.
func TestBenchValues(t *testing.T) {
	if g := GCirculation(WvdMRadiusFm()); math.Abs(g-1) > 1e-12 {
		t.Errorf("charge at c on λ/4π: g = %v, want 1", g)
	}
	if q := WvdMCharge(); math.Abs(q-0.910318) > 2e-6 {
		t.Errorf("WvdM charge %.6f", q)
	}
	if q := 1 / math.Sqrt(4*math.Pi*alpha); math.Abs(q-3.30227) > 1e-4 {
		t.Errorf("Rañada unit %.5f e", q)
	}
	// the proton's 0.8409 fm moves 1S–2S by ≈ 1 MHz (7/8 of the ≈ 1.1 MHz 1S shift)
	if s := H1S2SShiftMHz(0.8409); math.Abs(s-0.97) > 0.03 {
		t.Errorf("1S–2S finite-size shift %.3f MHz", s)
	}
}

// Steps 1, 2, 4 on a coarse grid: the single loop holds and keeps its
// charge, its inertia about the axis agrees both ways, its spin-½ size
// bound lies above the Compton length, and the windows do not overlap;
// the double loop unties.
func TestElectronKnot(t *testing.T) {
	if testing.Short() {
		t.Skip("relaxation")
	}
	a := NewAxi(8, 8, 0.15, 1, 1, 1, 1, 1)
	e2, e4, _ := a.Relax(20000, 0.05, 1e-7)
	if d := a.Degree(); math.Abs(d-1) > 0.06 {
		t.Errorf("A11 degree after relaxation %.4f", d)
	}
	mo := a.Inertia(16)
	if math.Abs(mo.IParVia/mo.IParFull-1) > 1e-9 {
		t.Errorf("I∥ %.6f vs %.6f", mo.IParVia, mo.IParFull)
	}
	k := Knot{Q: 1, Fermion: true, E: e2 + e4, Core: a.Core(), RmsRho: math.Sqrt(mo.Rho2), RmsR: math.Sqrt(mo.R2), IPar: mo.IPar, IPer: mo.IPerp}
	if k.MinSizeSpin() < lbarC {
		t.Errorf("spin-½ size bound %.1f fm below λ̄", k.MinSizeSpin())
	}
	if lo, hi := Overlap(k.Windows()); lo <= hi {
		t.Errorf("windows overlap: [%g, %g] fm", lo, hi)
	}
	if _, ok := k.E0At(WvdMRadiusFm() / k.Core); ok {
		t.Error("reading A: a real E0 with the spin energy, want none")
	}
	d := NewAxi(8, 8, 0.15, 1, 1, 1, 1, 2)
	d.Relax(2000, 0.05, 0)
	if deg := d.Degree(); deg > 1.5 {
		t.Errorf("double loop kept degree %.3f, expected it to untie", deg)
	}
	t.Logf("A11 E %.2f core %.3f I∥ %.1f I⊥ %.1f ℓ_min %.1f fm g %.3f", k.E, k.Core, k.IPar, k.IPer, k.MinSizeSpin(), k.GRigid())
}
