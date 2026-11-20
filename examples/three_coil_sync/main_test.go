// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// The turn-by-turn inductance of a close-wound helix agrees with Lorenz's
// current-sheet formula (Nagaoka) to the round-wire correction, a few
// per cent for these coils (Rosa 1906).
func TestHelixAgainstNagaoka(t *testing.T) {
	for _, h := range []Helix{
		{Turns: 40, R: 25.5e-3, Wire: 1e-3, Pitch: 1.1e-3},
		{Turns: 400, R: 37.7e-3, Wire: 0.4e-3, Pitch: 0.45e-3},
	} {
		got, want := h.SelfL(), em.Nagaoka(h.R, h.Length(), float64(h.Turns))
		if rel := math.Abs(got-want) / want; rel > 0.05 {
			t.Errorf("N %d: L %.4g H, Nagaoka %.4g H (%.1f %%)", h.Turns, got, want, 100*rel)
		}
	}
}

// Mutual inductances are symmetric and every coupling is below 1.
func TestCouplings(t *testing.T) {
	d := Baseline()
	_, _, c2a, _, c3 := d.Coils()
	if a, b := Mutual(c2a, c3), Mutual(c3, c2a); math.Abs(a-b) > 1e-12*a {
		t.Errorf("M not symmetric: %g vs %g", a, b)
	}
	e := d.Electrical(60e3)
	for name, k := range map[string]float64{"K1": e.K1, "K2ab": e.K2ab, "K2a3": e.K2a3, "K2b3": e.K2b3} {
		if !(k > 0 && k < 1) {
			t.Errorf("%s = %g", name, k)
		}
	}
	// The 3×3 matrix of coil 2's strands and coil 3 must be positive definite.
	l := [3][3]float64{{e.L2a, e.M2ab, e.M2a3}, {e.M2ab, e.L2b, e.M2b3}, {e.M2a3, e.M2b3, e.L3}}
	d1 := l[0][0]
	d2 := l[0][0]*l[1][1] - l[0][1]*l[1][0]
	d3 := l[0][0]*(l[1][1]*l[2][2]-l[1][2]*l[2][1]) - l[0][1]*(l[1][0]*l[2][2]-l[1][2]*l[2][0]) + l[0][2]*(l[1][0]*l[2][1]-l[1][1]*l[2][0])
	if d1 <= 0 || d2 <= 0 || d3 <= 0 {
		t.Errorf("inductance matrix not positive definite: %g %g %g", d1, d2, d3)
	}
}

// The Royer alone oscillates near its tank frequency with the ZVS peak
// π·Vcc on each drain (half-sine, choke-fed), and the books close.
func TestRoyerAlone(t *testing.T) {
	d := Baseline()
	d.HVOff = true
	r := Simulate(context.Background(), d, false)
	if r.Err != nil {
		t.Fatal(r.Err)
	}
	if rel := math.Abs(r.FRoyer-r.FR0) / r.FR0; rel > 0.05 {
		t.Errorf("Royer %.5g Hz, tank %.5g Hz", r.FRoyer, r.FR0)
	}
	if rel := math.Abs(r.VdsPk-math.Pi*vcc) / (math.Pi * vcc); rel > 0.1 {
		t.Errorf("drain peak %.1f V, π·Vcc %.1f V", r.VdsPk, math.Pi*vcc)
	}
	if math.Abs(r.Mismatch) > 1e-3 {
		t.Errorf("ledger mismatch %.2g of input", r.Mismatch)
	}
	if r.POut <= 0 || r.POut > r.PInRoy {
		t.Errorf("P_out %.3g W against %.3g W in", r.POut, r.PInRoy)
	}
}

// Output never exceeds the total input, with either or both drives.
func TestNoSurplus(t *testing.T) {
	base := Baseline()
	hv := base
	hv.RoyerOff = true
	for _, r := range parallel(context.Background(), []Design{base, hv}, false) {
		if r.Err != nil {
			t.Fatal(r.Err)
		}
		if r.POut > r.PInHV+r.PInRoy {
			t.Errorf("P_out %.4g W > input %.4g W", r.POut, r.PInHV+r.PInRoy)
		}
		if math.Abs(r.Mismatch) > 1e-3 {
			t.Errorf("ledger mismatch %.2g", r.Mismatch)
		}
	}
}

// Decode follows Space order and keeps integers whole.
func TestDecode(t *testing.T) {
	s := Space()
	if len(s) != 10 {
		t.Fatalf("%d genes", len(s))
	}
	d := Decode([]float64{40, 6, 5e-3, 400, 8, 680e-9, 100e-9, 50, 1, 90})
	if d.N3 != 40 || d.N1p != 8 || d.PhiHV != 90 || math.Abs(d.FHV-d.RoyerFreq0()) > 1e-9 {
		t.Errorf("decoded %+v", d)
	}
}

// The two evolved designs reproduce their recorded steady state
// (run_output.txt, run_free.txt) and never deliver more than they take in.
func TestEvolvedDesigns(t *testing.T) {
	for _, c := range []struct {
		name string
		d    Design
		pOut float64
	}{{"synchronized", EvolvedSync(), 374.587}, {"free-running", EvolvedFree(), 577.443}} {
		r := Simulate(context.Background(), c.d, false)
		if r.Err != nil {
			t.Fatalf("%s: %v", c.name, r.Err)
		}
		if rel := math.Abs(r.POut-c.pOut) / c.pOut; rel > 1e-4 {
			t.Errorf("%s: P_out %.4f W, recorded %.4f W", c.name, r.POut, c.pOut)
		}
		if r.POut >= r.PInHV+r.PInRoy {
			t.Errorf("%s: P_out %.4g W ≥ input %.4g W", c.name, r.POut, r.PInHV+r.PInRoy)
		}
		if math.Abs(r.Mismatch) > 1e-3 {
			t.Errorf("%s: ledger mismatch %.2g", c.name, r.Mismatch)
		}
	}
}
