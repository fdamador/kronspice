// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// Tesla's two-disk, belt-joined machine and the single disk with a rim
// brush (the two baselines main reports), at their operating points.
// Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	one := Defaults()
	one.TwoDisk = false
	export.Test(t, "export",
		export.Item{Name: "homopolar_two_disk_belt", Build: Defaults().Circuit, Control: []string{".op"}, Opt: sheet.Options{}},
		export.Item{Name: "homopolar_single_disk_rim", Build: one.Circuit, Control: []string{".op"}, Opt: sheet.Options{}},
	)
}
