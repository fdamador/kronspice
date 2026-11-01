// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"strings"
	"testing"

	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/report"
)

// The GA problem circuit (default 40-section LC line plant, starting
// pre-filter values) with the GA's AC sweep, torn at RT as the diakoptic
// variant of main tears it. Netlist, KiCad schematic and results sheet go
// to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	src := strings.Replace(Netlist("line", 40), ".end\n", ".ac dec 50 100 1meg\n.end\n", 1)
	export.Test(t, "export", export.Item{
		Name: "torn_ga_line40", Src: src,
		Opt: report.Options{Tear: []string{"RT"}},
	})
}
