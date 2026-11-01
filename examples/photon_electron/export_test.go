// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The plain and crossed 200-section ring lines with a 1 A AC drive at
// (a0, b0). Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	item := func(name string, crossed bool) export.Item {
		r := DefaultRing(crossed)
		ctl := ".ac lin 2001 1k 450k"
		return export.Item{Name: name, Control: []string{ctl}, Build: func() (*circuit.Circuit, error) {
			return r.Build(circuit.Element{Kind: circuit.CurrentSource, ACMag: 1})
		}}
	}
	export.Test(t, "export", item("ring_plain", false), item("ring_crossed", true))
}
