// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/fluid"
)

// piezoInWater is an electrical → mechanical → acoustic chain: a piezo
// disc (C0, force factor 0.8 N/V) on a mechanical resonator, a piston of
// radius 1 cm into water, and the piston's radiation load.
func piezoInWater(t testing.TB) *circuit.Circuit {
	a := 0.01
	Ma, Rpar, _, err := fluid.PistonRadiation(a, fluid.Water)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.NewBuilder("piezo in water").Domain(exchange.Mechanical, "v").Domain(exchange.Acoustic, "p").
		Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, Value: 1, ACMag: 10}, "e", "0").
		Resistor("RS", "e", "x", 50).Capacitor("C0", "x", "0", 2e-9).Gyrator("Y1", "x", "0", "v", "0", 0.8).
		Capacitor("CM", "v", "0", 2e-3).Inductor("LK", "v", "0", 1/5e7).Resistor("RB", "v", "0", 1/2.0).
		Gyrator("Y2", "v", "0", "p", "0", math.Pi*a*a).
		Inductor("LMA", "p", "0", Ma).ResistorTo("RA", "p", "0", Rpar, exchange.Acoustic).Build()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// The multi-domain chain in Kron's frames: transducers are antisymmetric
// mutual admittances of the primitive network; mesh and node-pair solutions
// equal MNA to 1e−9 at DC and across the resonance, and the open paths to
// the electrical, thermal and acoustic exterior junctions carry the MNA
// ledger's powers (the electrical source's power in, the damper's and
// resistor's heat and the radiated sound out).
func TestKronMultiDomainTwin(t *testing.T) {
	c := piezoInWater(t)
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	spec := analysis.ACSpec{Variation: "dec", Points: 20, FStart: 1e3, FStop: 2e4}
	mna, err := analysis.AC(c, spec, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	for _, form := range []Form{MeshForm, NodeForm} {
		k, err := OP(c, form)
		if err != nil {
			t.Fatal(err)
		}
		if w, _ := worstRel(real2c(k.NodeVoltages), real2c(op.NodeVoltages)); w > twinTol {
			t.Errorf("%s OP: %g", form, w)
		}
		ka, err := AC(c, mna.Freq, form)
		if err != nil {
			t.Fatal(err)
		}
		n, _ := Build(c)
		for fi, f := range mna.Freq {
			var got, want []complex128
			for i := 1; i < len(mna.V); i++ {
				got, want = append(got, ka.V[i][fi]), append(want, mna.V[i][fi])
			}
			// Velocities and pressures differ by orders of magnitude from
			// volts: compare each domain against its own scale.
			for i := range got {
				if d := cmplx.Abs(got[i] - want[i]); d > twinTol*cmplx.Abs(want[i])+1e-300 {
					t.Errorf("%s at %g Hz: node %s %v, MNA %v", form, f, mna.NodeNames[i+1], got[i], want[i])
				}
			}
			ext, err := ExteriorAC(n, form, f)
			if err != nil {
				t.Fatal(err)
			}
			kl, ml := ext.Ledger, mna.Ledger
			for _, pr := range [][2]float64{{ml.Inflow[fi], kl.Inflow}, {ml.Outflow[fi], kl.Outflow}, {ml.Dissipation[fi], kl.Dissipation}} {
				if math.Abs(pr[0]-pr[1]) > twinTol*kl.Scale {
					t.Errorf("%s at %g Hz: ledger MNA %g, Kron %g", form, f, pr[0], pr[1])
				}
			}
			_, acoustic := ext.Junction(exchange.Acoustic)
			if math.Abs(-acoustic-kl.Outflow) > 1e-12*kl.Scale || !(acoustic < 0) {
				t.Errorf("%s at %g Hz: acoustic junction %g W, outflow %g W", form, f, acoustic, kl.Outflow)
			}
		}
	}
}
