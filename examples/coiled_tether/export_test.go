// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The solar-driven resonant loops of main (ordinary coil with tank
// capacitor, bifilar coil, bifilar with matched load). Netlist, KiCad
// schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	ord, bif := Coils()
	cb := bif.SelfCapacitance()
	item := func(name string, c float64, load bool) export.Item {
		l := bif.Inductance()
		f := 1 / (2 * math.Pi * math.Sqrt(l*c))
		rw := bif.Resistance(f)
		rl := 0.0
		if load {
			rl = rw
		}
		v := math.Sqrt(2 * (rw + rl) * SolarPower())
		stop := 7 * 2 * l / (rw + rl)
		return export.Item{Name: name, Build: func() (*circuit.Circuit, error) { return Loop(name, rw, l, c, rl, v, f) },
			Control: []string{fmt.Sprintf(".tran %g %g 0 %g", 1/f/200, stop, 1/f/200)}}
	}
	_ = ord
	export.Test(t, "export",
		item("coiled_tether_ordinary_tank", cb, false),
		item("coiled_tether_bifilar", cb, false),
		item("coiled_tether_bifilar_load", cb, true),
	)
}
