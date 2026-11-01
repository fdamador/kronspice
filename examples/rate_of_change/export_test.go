// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/netlist"
)

// The water plate cell (1 cm², 1 mm) as C∞ ∥ (R–C), swept 1 MHz–100 GHz.
// Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	build := func() (*circuit.Circuit, error) {
		d, err := netlist.ParseString(em.WaterDebye25.CellNetlist(CellArea, CellGap, ".ac dec 10 1meg 100g"))
		if err != nil {
			return nil, err
		}
		return circuit.FromDeck(d)
	}
	export.Test(t, "export", export.Item{Name: "water_cell_debye", Build: build, Control: []string{".ac dec 10 1meg 100g"}})
}
