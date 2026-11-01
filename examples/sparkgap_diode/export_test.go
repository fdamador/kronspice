// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The aluminium–lead water-cell bridge rectifier (Part 2). Netlist, KiCad
// schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	r := Rectifier{Name: "water-cell bridge rectifier", Bridge: true, Arm: Arm{DefaultCell(Pb, Al), 1}, Vp: 100, F: 60, Rs: 0.5, Cload: 2200e-6, Rload: 50}
	export.Test(t, "export", export.Item{Name: "water_cell_bridge", Build: func() (*circuit.Circuit, error) { return r.Build() },
		Control: []string{".tran 41.7u 0.5"}})
}
