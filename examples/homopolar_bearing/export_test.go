// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// The reference unit (load current through the shaft bearings) at its
// operating point, and main's capacitor-load case: 1 F in series with R_L,
// field stepped (1 ms ramp), run for 5·R_L·C as Charge does. The sheet
// is at the operating point, where that run starts: field at zero, C
// uncharged, so only the driven shafts are non-zero there. Netlist,
// KiCad schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	d := Defaults()
	tEnd := 5 * d.RL * 1
	export.Test(t, "export",
		export.Item{Name: "homopolar_bearing_reference", Build: d.Circuit, Control: []string{".op"}, Opt: sheet.Options{}},
		export.Item{Name: "homopolar_bearing_capacitor", Build: func() (*circuit.Circuit, error) { return d.build(true, 1, 0) },
			Control: []string{fmt.Sprintf(".tran %.17g %.17g", tEnd/400, tEnd)}, Opt: sheet.Options{}},
	)
}
