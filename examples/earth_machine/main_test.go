// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The circuit and the RK4 reference agree before chaos separates them.
func TestCircuitMatchesRK4(t *testing.T) {
	p := Rikitake()
	tr, _, err := p.RunCircuit(10, 0.01)
	if err != nil {
		t.Fatal(err)
	}
	ref := p.RunRK4(10, 0.01)
	for k := 1; k < len(tr.T); k++ {
		j := int(math.Round(tr.T[k] / 0.01))
		if d := math.Abs(tr.I1[k] - ref.I1[j]); d > 0.02 {
			t.Fatalf("t=%.2f: I1 circuit %.4f vs RK4 %.4f", tr.T[k], tr.I1[k], ref.I1[j])
		}
	}
}

// Rikitake's machine reverses on its own; the energy books close.
func TestReversesAndBalances(t *testing.T) {
	p := Rikitake()
	tr, mis, err := p.RunCircuit(300, 0.01)
	if err != nil {
		t.Fatal(err)
	}
	if n := len(Reversals(tr, 50, 1e-2)); n < 5 {
		t.Fatalf("only %d reversals", n)
	}
	in, _, _, res := p.Books(tr)
	if math.Abs(res/in) > 1e-4 {
		t.Fatalf("energy residual %.3g of %.3g", res, in)
	}
	// The engine's own ledger, UIC start included (analysis/tran.go).
	if math.Abs(mis/in) > 1e-5 {
		t.Fatalf("engine integration mismatch %.3g of %.3g", mis, in)
	}
}

// Counter-rotating with the same wiring: negative feedback, the field dies.
func TestCounterRotatingSameWiringDecays(t *testing.T) {
	p := Rikitake()
	p.G2, p.W2 = -1, -1
	tr, _, err := p.RunCircuit(100, 0.01)
	if err != nil {
		t.Fatal(err)
	}
	if m := absMax(tr.I1, len(tr.I1)-100); m > 1e-6 {
		t.Fatalf("I1 still %.3g", m)
	}
}

// Counter-rotating with coil 2 rewound is Rikitake's machine with I2, ω2
// negated: I1 is the same.
func TestCounterRotatingRewoundMirrors(t *testing.T) {
	a := Rikitake()
	b := Rikitake()
	b.G2, b.W2, b.Wind2, b.I2 = -1, 0, -1, -0.5
	ta, _, _ := a.RunCircuit(50, 0.01)
	tb, _, _ := b.RunCircuit(50, 0.01)
	for k := range ta.T {
		if math.Abs(ta.I1[k]-tb.I1[k]) > 1e-9 || math.Abs(ta.I2[k]+tb.I2[k]) > 1e-9 {
			t.Fatalf("t=%.2f: not a mirror", ta.T[k])
		}
	}
}

func TestEquilibriumGoldenRatio(t *testing.T) {
	k2, per := Equilibrium(1)
	if math.Abs(k2-(1+math.Sqrt(5))/2) > 1e-12 || math.Abs(per-4.2018) > 1e-3 {
		t.Fatalf("k²=%g period=%g", k2, per)
	}
}
