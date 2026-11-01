// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/solver"
)

// A divider with a negative lower leg: v(out) = 10·(−500)/(1000 − 500) =
// −10 V. The negative resistance delivers v²/|R| = 0.2 W, booked as
// model-supplied; the source delivers 0.2 W; R1 dissipates 0.4 W.
func TestNegativeResistorDivider(t *testing.T) {
	c := mustDeck(t, "neg\nV1 in 0 10\nR1 in out 1k\nR2 out 0 -500\n.op\n")
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("out"); !near(v, -10, 1e-12) {
		t.Fatalf("v(out) = %g, want -10", v)
	}
	l := r.Ledger
	if !near(l.ModelSupplied, 0.2, 1e-15) || !near(l.Inflow, 0.2, 1e-15) || !near(l.Dissipation, 0.4, 1e-15) {
		t.Fatalf("ledger: inflow %g, model-supplied %g, dissipation %g", l.Inflow, l.ModelSupplied, l.Dissipation)
	}
	if l.RelativeResidual() > 1e-12 {
		t.Fatalf("residual %g", l.RelativeResidual())
	}
	if len(l.Supplied) != 1 || l.Supplied[0].Element != "R2" {
		t.Fatalf("supplied rows = %+v", l.Supplied)
	}
	// The negative resistance takes nothing in through its thermal port.
	r2, _ := r.Element("R2")
	if r2.ExteriorPower[exchange.Thermal] != 0 || r2.Ports[exchange.Thermal].Status != exchange.NotModeled {
		t.Fatalf("R2 thermal port: %+v, %g W", r2.Ports[exchange.Thermal], r2.ExteriorPower[exchange.Thermal])
	}
}

// Conductances that cancel exactly leave a node with no conductance to
// ground: a typed singular-matrix error, never a panic.
func TestNegativeResistorCancellationIsSingular(t *testing.T) {
	c := mustDeck(t, "s\nV1 a 0 1\nR1 a b 1k\nR2 b 0 1k\nR3 b 0 -500\n.op\n")
	if _, err := OP(c, Options{}); !errors.Is(err, solver.ErrSingularMatrix) {
		t.Fatalf("err = %v, want ErrSingularMatrix", err)
	}
}

// AC gain above 1 from a negative leg: H = −2k/(1k − 2k) = 2 at every
// frequency, and the average-power ledger balances with model-supplied power.
func TestNegativeResistorACGain(t *testing.T) {
	c := mustDeck(t, "g\nV1 in 0 AC 1\nR1 in out 1k\nR2 out 0 -2k\n.ac dec 2 1 1k\n")
	r, err := AC(c, ACSpec{Variation: "dec", Points: 2, FStart: 1, FStop: 1e3}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := r.NodeVoltage("out")
	for k := range v {
		if cmplx.Abs(v[k]-2) > 1e-12 {
			t.Fatalf("H = %v, want 2", v[k])
		}
		// |V|²/(2|R|) = 4/4000 W.
		if !near(r.Ledger.ModelSupplied[k], 1e-3, 1e-15) {
			t.Fatalf("model-supplied %g", r.Ledger.ModelSupplied[k])
		}
	}
}

// A parallel LC tank with a negative resistance across it oscillates with a
// growing envelope e^(αt), α = 1/(2|R|C) (the circuit behind tunnel-diode
// oscillators). The ledger shows the stored energy growing by exactly what
// the negative resistance supplies.
func TestNegativeResistorGrowingOscillation(t *testing.T) {
	c := mustDeck(t, "osc\nC1 a 0 1u IC=1\nL1 a 0 1m\nRN a 0 -10k\n.tran 1u 2m uic\n")
	wf, err := Tran(c, TranSpec{Step: 1e-6, Stop: 2e-3, UIC: true, MaxStep: 2e-7}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	checkLedger(t, wf)
	v, _ := wf.NodeVoltage("a")
	// Exact solution with v(0) = 1, i_L(0) = 0: v = e^(αt)(cos ωt + (α/ω) sin ωt).
	alpha := 1 / (2 * 1e4 * 1e-6)
	w := math.Sqrt(1/(1e-3*1e-6) - alpha*alpha)
	worst := 0.0
	for k, tt := range wf.Time {
		exact := math.Exp(alpha*tt) * (math.Cos(w*tt) + alpha/w*math.Sin(w*tt))
		worst = math.Max(worst, math.Abs(v[k]-exact))
	}
	// Peak |v| is e^(0.1) ≈ 1.105 V; ~10 periods of trapezoidal phase error.
	if worst > 1e-2 {
		t.Fatalf("worst error %g V", worst)
	}
	l := wf.Ledger
	gain := l.StoredEnd - l.StoredStart
	if gain <= 0 || math.Abs(l.EnergyModelSupplied-gain) > 1e-3*gain {
		t.Fatalf("stored energy grew %g J, model-supplied %g J", gain, l.EnergyModelSupplied)
	}
	if l.EnergyIn != 0 {
		t.Fatalf("energy in from the exterior = %g, want 0", l.EnergyIn)
	}
}
