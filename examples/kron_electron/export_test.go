// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// Kelvin's ring as a Kron machine (main's coil switch-on) and the ring
// tilting in a field. Netlist, KiCad schematic and results sheet go to
// export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	export.Test(t, "export",
		export.Item{Name: "kron_electron_kelvin", Build: func() (*circuit.Circuit, error) {
			return KelvinCircuit(Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.01}, Coil{R: 0.1, Z: 0.05, NI: 2e4}, 1, false)
		}, Control: []string{".tran 5m 2 uic"}},
		export.Item{Name: "kron_electron_tilt", Build: func() (*circuit.Circuit, error) {
			return TiltCircuit(Tilt{A: 0.05, E: 1e-6, Mass: 0.01, OmegaSpin: 1000, B: 1, Theta0: 0.05})
		}, Control: []string{".tran 50m 60 uic"}},
	)
}
