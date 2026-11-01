// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/kron"
)

// The motor in Kron's frame I and rotor frame II(b) with main's 1.5 s
// start-up transient, and the parametric tank (damping 0.5 × pump rate).
// Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	s1 := Primitive(0.12, 0.12, 0.11)
	tran := []string{".tran 0.5m 1.5 uic"}
	export.Test(t, "export",
		export.Item{Name: "kron_machine_motor_frame1", Build: func() (*circuit.Circuit, error) { return Motor(s1) }, Control: tran},
		export.Item{Name: "kron_machine_motor_frame2", Build: func() (*circuit.Circuit, error) {
			s2, err := kron.RotatingFrame(s1, kron.RotorFrame(4, [][2]int{{2, 3}}, 1))
			if err != nil {
				return nil, err
			}
			return Motor(s2)
		}, Control: tran},
		export.Item{Name: "kron_machine_parametric", Build: parametricTank, Control: []string{".tran 10u 60m uic"}},
	)
}

// parametricTank is the circuit of parametric(1e-3, 0.05, 1e3, 0.5).
func parametricTank() (*circuit.Circuit, error) {
	const l0, m, f0, damp = 1e-3, 0.05, 1e3, 0.5
	w0 := 2 * math.Pi * f0
	s := circuit.NewMachineSpec(1)
	s.L0[0][0] = l0
	s.Harmonic(2).Cos[0][0] = m * l0
	r := damp * m * w0 / 4 * 2 * l0
	b := circuit.NewBuilder("parametric").Domain(exchange.Mechanical, "w").MachineModel("pump", s).
		Machine("N1", "pump", [][2]string{{"a", "b"}}, "w", "0").
		CapacitorIC("C1", "a", "0", 1/(w0*w0*l0), 1).VoltageSource("VW", "w", "0", w0)
	b.Resistor("R1", "b", "0", math.Max(r, 1e-9))
	return b.BuildTransient()
}
