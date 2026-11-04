// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The scaled geometry gives the copper machine's couplings
// (examples/homopolar: 0.519 and 0.387), falling with the gap.
func TestStackedCouplings(t *testing.T) {
	ks, km := StackedCouplings(StackedGapHomopolar)
	if math.Abs(ks-0.519) > 2e-3 || math.Abs(km-0.387) > 2e-3 {
		t.Errorf("KS %g, KM %g", ks, km)
	}
	prev := 1.0
	for _, g := range []float64{0.6, 1, 2, 4} {
		k, _ := StackedCouplings(g)
		if k >= prev || k <= 0 {
			t.Errorf("KS not falling with the gap: %g at %g", k, g)
		}
		prev = k
	}
}

// The stacked circuit (coil inductors with their Mutual, sense windings)
// equals the RK4 equations, sign included, and its ledger closes.
func TestStackedCircuitMatchesRK4(t *testing.T) {
	for _, g := range []float64{StackedGapHomopolar, 2} {
		tr, mis, err := stackedRunEngine(g, 10, 0.01)
		if err != nil {
			t.Fatal(err)
		}
		ref := StackedRikitake(g).RunRK4(10, 0.01)
		for k := 1; k < len(tr.T); k++ {
			j := int(math.Round(tr.T[k] / 0.01))
			if d := math.Abs(tr.I1[k] - ref.I1[j]); d > 0.02 {
				t.Fatalf("gap %g, t = %g: I1 %g vs RK4 %g", g, tr.T[k], tr.I1[k], ref.I1[j])
			}
		}
		if math.Abs(mis) > 1e-4 {
			t.Errorf("gap %g: ledger mismatch %g J", g, mis)
		}
	}
}

// Any self-coupling drains A = ω1 − ω2, and the field then holds one
// polarity.
func TestStackedLocksPolarity(t *testing.T) {
	p := StackedRikitake(3)
	tr := p.RunRK4(2000, 0.5)
	k := len(tr.T) - 1
	if a := tr.W1[k] - tr.W2[k]; math.Abs(a) > 1e-2 {
		t.Errorf("A = %g after t = 2000", a)
	}
	if r := p.StatRun(5000, 0.002, 200, 1e-2); len(r.Rev) != 0 {
		t.Errorf("%d reversals", len(r.Rev))
	}
}
