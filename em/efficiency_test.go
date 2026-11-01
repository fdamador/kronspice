// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em_test

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/netlist"
)

func parse(t *testing.T, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// DC efficiency from the ledger: a source with 50 Ω internal resistance
// into R_L gives η = R_L/(R_L + 50) — ½ at the matched load.
func TestOPEfficiency(t *testing.T) {
	for _, rl := range []float64{50, 150, 10} {
		c := parse(t, fmt.Sprintf("divider\nV1 s 0 10\nRS s a 50\nRL a 0 %g\n.end\n", rl))
		r, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		eta, err := em.OPEfficiency(r, "RL")
		if err != nil || math.Abs(eta-rl/(rl+50)) > 1e-12 {
			t.Errorf("R_L = %g: η = %g (%v), want %g", rl, eta, err, rl/(rl+50))
		}
	}
	c := parse(t, "nothing\nI1 0 a 0\nR1 a 0 1\n.end\n")
	r, _ := analysis.OP(c, analysis.Options{})
	if _, err := em.OPEfficiency(r, "R1"); err == nil {
		t.Error("efficiency without inflow")
	}
}

// Transient efficiency over whole periods of a sinusoidal steady state
// equals the AC phasor efficiency Re(V_L·I_L*)/Re(V_S·I_S*) of the same
// circuit (R_S, then R_L ∥ C), from the engine's own AC analysis.
func TestTranEfficiencyMatchesAC(t *testing.T) {
	const f = 1e3
	src := "rc load\nV1 s 0 SIN(0 1 1k) AC 1\nRS s a 50\nRL a 0 200\nCL a 0 1u\n.end\n"
	c := parse(t, src)
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: 1e-6, Stop: 20e-3, MaxStep: 1e-6}, analysis.TranOptions{Options: analysis.Options{RelTol: 1e-7}})
	if err != nil {
		t.Fatal(err)
	}
	eta, eLoad, eIn, err := em.TranEfficiency(wf, 16e-3, 20e-3, "RL")
	if err != nil {
		t.Fatal(err)
	}
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	va, _ := ac.NodeVoltage("a")
	vs := complex(1, 0)
	is := (vs - va[0]) / 50
	pin := real(vs * cmplx.Conj(is))
	pl := real(va[0]*cmplx.Conj(va[0])) / 200
	if math.Abs(eta-pl/pin) > 1e-4 {
		t.Errorf("transient η = %.6f (%.4g J of %.4g J), AC η = %.6f", eta, eLoad, eIn, pl/pin)
	}
}

// Tellegen's theorem in its general form: the node voltages of one
// circuit times the branch currents of another with the same graph sum to
// zero, and so do the voltages and currents of one circuit at two
// different instants.
func TestTellegenGeneral(t *testing.T) {
	r := rand.New(rand.NewSource(41))
	topo := func(vals []float64) string {
		var b strings.Builder
		b.WriteString("tellegen\n")
		fmt.Fprintf(&b, "V1 1 0 %g\nR1 1 2 %g\nR2 2 0 %g\nR3 2 3 %g\nI1 3 0 %g\nR4 3 1 %g\nR5 3 4 %g\nV2 4 0 %g\n.end\n",
			vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7])
		return b.String()
	}
	for trial := 0; trial < 50; trial++ {
		var a, b [8]float64
		for k := range a {
			a[k], b[k] = math.Pow(10, 2*r.Float64()), math.Pow(10, 2*r.Float64())
		}
		ca, cb := parse(t, topo(a[:])), parse(t, topo(b[:]))
		ra, err := analysis.OP(ca, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		rb, err := analysis.OP(cb, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		ib := make([]float64, len(cb.Elements))
		for k, e := range cb.Elements {
			el, _ := rb.Element(e.Name)
			ib[k] = el.Current
		}
		// Node voltages of A in A's node order; B has the same names and
		// order (same netlist text apart from values).
		sum, scale, err := em.Tellegen(ca, ra.NodeVoltages, ib)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(sum) > 1e-12*scale {
			t.Fatalf("Σ v_A·i_B = %g (scale %g)", sum, scale)
		}
	}
	// Two instants of one transient.
	c := parse(t, "rlc\nV1 1 0 PULSE(0 1 0 1u 1u 5u 10u)\nR1 1 2 10\nL1 2 3 1m\nC1 3 0 1u\nR2 3 0 100\n.end\n")
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: 1e-7, Stop: 50e-6}, analysis.TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	k1, k2 := len(wf.Time)/3, 2*len(wf.Time)/3
	v := make([]float64, c.NumNodes())
	for n := 1; n < c.NumNodes(); n++ {
		v[n] = wf.V[n][k1]
	}
	i := make([]float64, len(c.Elements))
	for e := range c.Elements {
		i[e] = wf.I[e][k2]
	}
	sum, scale, err := em.Tellegen(c, v, i)
	if err != nil || math.Abs(sum) > 1e-9*scale {
		t.Errorf("Σ v(t₁)·i(t₂) = %g (scale %g, %v)", sum, scale, err)
	}
}
