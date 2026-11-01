// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/fluid"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// The drive chain in water without and with the matching inductor, and in
// heavy water with it, swept over ±3 % around the breathing mode as main
// does; the sheet is taken at the mode frequency f0. Netlist, KiCad
// schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	var items []export.Item
	for _, x := range []struct {
		name  string
		m     fluid.Medium
		withL bool
	}{
		{"acoustic_drive_water", fluid.Water, false},
		{"acoustic_drive_water_matched", fluid.Water, true},
		{"acoustic_drive_d2o_matched", fluid.HeavyWater, true},
	} {
		c, f0, err := chain(x.m, x.withL)
		if err != nil {
			t.Fatal(err)
		}
		items = append(items, export.Item{
			Name:    x.name,
			Build:   func() (*circuit.Circuit, error) { return c, nil },
			Control: []string{fmt.Sprintf(".ac lin 2001 %.17g %.17g", 0.97*f0, 1.03*f0)},
			Opt:     sheet.Options{Freq: f0},
		})
	}
	export.Test(t, "export", items...)
}
