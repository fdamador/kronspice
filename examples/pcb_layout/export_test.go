// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/report"
)

// The coupled pair w = 3 mm, s = 2 mm (the middle geometry of
// TestCoupledPair), L and C from the cross-section field solve, at the
// crosstalk frequency 100 MHz. Netlist, KiCad schematic and results sheet
// go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	export.Test(t, "export", export.Item{
		Name: "pcb_layout_pair_w3_s2",
		Build: func() (*circuit.Circuit, error) {
			p := Pair{W: 3e-3, S: 2e-3}
			r, err := p.Extract()
			if err != nil {
				return nil, err
			}
			return p.Circuit(r)
		},
		Control: []string{fmt.Sprintf(".ac lin 1 %g %g", fNext, fNext)},
		Opt:     report.Options{Freq: fNext},
	})
}
