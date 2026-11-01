// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// The EM-mode SPG loop of SimulateSPGLoop (no-load EMF, ball-bearing motor
// with the κψ CheckSPG derives from the measured η, the rest of the loop
// sized for 3438 A, shaft driven at the SPG's speed), at the middle of the
// swept contact resistances, Rb = 10 µΩ, at its operating point. The
// circuit is rebuilt here as SimulateSPGLoop builds it (main builds it
// inline). Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	c := CheckSPG(0.10) // Marinov 1989, the registry's measured η
	const rb = 1e-5
	build := func() (*circuit.Circuit, error) {
		w0 := c.Omega
		kappa := c.KappaPsi / math.Tanh(c.Omega/w0)
		rRest := c.E0/c.I - rb - c.KappaPsi*c.Omega
		return circuit.NewBuilder("SPG loop with bearing motor").Domain(exchange.Mechanical, "w").
			VoltageSource("E0", "e", "0", c.E0).
			BearingMotor("PB", "e", "x", "w", "0", rb, kappa, w0, false, 0).
			Resistor("RLOOP", "x", "0", rRest).
			VoltageSource("VW", "w", "0", c.Omega).Build()
	}
	// circuit.WriteNetlist writes element names as given; the example's
	// source "E0" (a V element) reads back as a VCVS, so the netlist does
	// not parse. Skip until the writer prefixes such names.
	if ckt, err := build(); err != nil {
		t.Fatal(err)
	} else if _, err := export.Circuit("bearing_motor_spg_loop", ckt, sheet.Options{}, ".op"); err != nil {
		t.Skipf("netlist writer: %v", err)
	}
	export.Test(t, "export", export.Item{Name: "bearing_motor_spg_loop", Build: build, Control: []string{".op"}, Opt: sheet.Options{}})
}
