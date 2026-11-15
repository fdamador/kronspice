// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func TestFreeDirac(t *testing.T) {
	p := [3]float64{0.3, 0.4, 1.2}
	e, h := HelicityOnPositive(p)
	want := math.Sqrt(1 + 0.09 + 0.16 + 1.44)
	for i, x := range e {
		s := 1.0
		if i < 2 {
			s = -1
		}
		_ = s
		if math.Abs(math.Abs(x)-want) > 1e-12 {
			t.Errorf("eigenvalue %g, want ±%g", x, want)
		}
	}
	if math.Abs(math.Abs(h[0])-0.5) > 1e-12 || math.Abs(math.Abs(h[1])-0.5) > 1e-12 || h[0]*h[1] > 0 {
		t.Errorf("helicities %v", h)
	}
	if math.Abs(SpinSquared()-0.75) > 1e-15 {
		t.Errorf("S² = %g", SpinSquared())
	}
}

func TestLandauG(t *testing.T) {
	for _, b := range []float64{1e-6, 1e-2, 1} {
		g, sp, cy := GFromLandau(LandauDirac(b, 30))
		if math.Abs(g-2) > 1e-9 || math.Abs(sp-(math.Sqrt(1+2*b)-1)) > 1e-12*math.Max(1, sp) || math.Abs(cy-sp) > 1e-12 {
			t.Errorf("b = %g: g %g, spin %g, cyclotron %g", b, g, sp, cy)
		}
	}
}

// The shooting agrees with the exact formula, and a coarse integration does
// not: the agreement is the integration's, not a tautology.
func TestHydrogen(t *testing.T) {
	for _, c := range [][2]int{{1, -1}, {2, -1}, {2, 1}, {2, -2}} {
		e := HydrogenLevel(c[0], c[1], Alpha)
		ex := DiracExact(c[0], c[1], Alpha)
		if math.Abs(e-ex) > 1e-14 {
			t.Errorf("%v: %g vs exact %g", c, e, ex)
		}
	}
	ex := DiracExact(1, -1, Alpha)
	if c := hydrogenLevelSteps(1, -1, Alpha, 50); math.Abs(c-ex) < 1e-13 {
		t.Errorf("a coarse integration should miss the level: %g vs %g", c, ex)
	}
	fs := (DiracExact(2, -2, Alpha) - DiracExact(2, 1, Alpha)) * HzPerMC2
	if math.Abs(fs/10.9497e9-1) > 1e-4 {
		t.Errorf("fine structure %g Hz", fs)
	}
}

func TestQED(t *testing.T) {
	if f := SchwingerF2(0, 400); math.Abs(f-0.5) > 1e-12 {
		t.Errorf("F2(0) = %g", f)
	}
	if a := A4Petermann(); math.Abs(a+0.32847896557919378) > 1e-14 {
		t.Errorf("A4 = %.17g", a)
	}
	// published comparison (arXiv:2311.12112): +34(16)e-14 with Rb, −101(27)e-14 with Cs
	_, rb := AeSeries(AlphaInvRb)
	_, cs := AeSeries(AlphaInvCs)
	if d := (AeMeasured - rb) * 1e14; math.Abs(d-34) > 3 {
		t.Errorf("Rb: measured − theory = %g e-14", d)
	}
	if d := (AeMeasured - cs) * 1e14; math.Abs(d+101) > 3 {
		t.Errorf("Cs: measured − theory = %g e-14", d)
	}
}

func TestTrap(t *testing.T) {
	tr := DefaultTrap()
	fc, w, depth, _, _, _, err := tr.AxialDip(5*tr.DipWidth(), 8001)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(fc-tr.Fz) > 0.01 || depth > 1e-3 || math.Abs(w*math.Sqrt(3)/tr.DipWidth()-1) > 2e-3 {
		t.Errorf("dip at %g, depth %g, width %g (×√3 vs %g)", fc-tr.Fz, depth, w, tr.DipWidth())
	}
	_, ae := AeSeries(AlphaInvRb)
	f := tr.Frequencies(GFromAe(ae))
	if math.Abs(f.NuM/43e3-1) > 0.03 || math.Abs(f.NuA/173e6-1) > 0.01 {
		t.Errorf("magnetron %g, anomaly %g", f.NuM, f.NuA)
	}
	rc, _ := ReadTank(f.NucBar, 1e6)
	ra, _ := ReadTank(f.NuA, 1e5)
	if g := GFromClocks(rc, f.NuZ, f.NuM, ra); math.Abs(g-1-ae) > 1e-13 {
		t.Errorf("g/2 read %.15f vs %.15f", g, 1+ae)
	}
}
