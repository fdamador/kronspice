// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// tetherCircuit builds the circuit Solve solves (same elements and nodes).
func tetherCircuit(t Tether) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("electrodynamic tether")
	b.VoltageSource("VEMF", "top", "m", t.EMF()).
		Resistor("RT", "m", "bot", t.Resistance()).
		Resistor("RCT", "top", "x", t.RContactor).
		Resistor("RCB", "bot", "y", t.RContactor)
	if t.Load > 0 {
		b.Resistor("RL", "x", "y", t.Load)
	} else {
		b.VoltageSource("VS", "x", "y", t.Supply)
	}
	b.Resistor("RREF", "y", "0", 1e9)
	return b.Build()
}

// The generator (500 Ω load) and thruster (3 kV supply) configurations of
// main, at their operating points. Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	gen := Defaults()
	th := Defaults()
	th.Load, th.Supply = 0, 3000
	export.Test(t, "export",
		export.Item{Name: "tether_generator", Build: func() (*circuit.Circuit, error) { return tetherCircuit(gen) }, Control: []string{".op"}},
		export.Item{Name: "tether_thruster", Build: func() (*circuit.Circuit, error) { return tetherCircuit(th) }, Control: []string{".op"}},
	)
}
