// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/internal/export"
)

// The two test circuits of H22: the linear two-mesh RLC whose meshes are
// the Wilson loops of step 1, and the diode pair of step 2. Netlist, KiCad
// schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	export.Test(t, "export",
		export.Item{Name: "quantumtricity_two_mesh", Src: TwoMesh + ".ac dec 20 10 1meg\n"},
		export.Item{Name: "quantumtricity_diode_pair", Src: DiodePair + ".op\n"},
	)
}
