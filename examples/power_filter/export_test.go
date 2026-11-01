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

// The unfiltered rectifier and one with 5th/7th tuned branches (1 mH with
// the E12 capacitors nearest to resonance at 300 Hz and 420 Hz: 270 µF
// and 150 µF), with main's 14-cycle transient at 600 steps per cycle.
// Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	h := 1 / (f0 * 600)
	tran := []string{fmt.Sprintf(".tran %.17g %.17g 0 %.17g", h, cycles/f0, h)}
	tuned := Filter{L5: 1e-3, C5: 270e-6, L7: 1e-3, C7: 150e-6}
	for _, x := range []struct{ n, l, c float64 }{{5, tuned.L5, tuned.C5}, {7, tuned.L7, tuned.C7}} {
		t.Logf("branch %g: tuned to %.0f Hz", x.n, 1/(2*math.Pi*math.Sqrt(x.l*x.c)))
	}
	export.Test(t, "export",
		export.Item{Name: "power_filter_unfiltered", Build: func() (*circuit.Circuit, error) { return Circuit(Filter{}) }, Control: tran},
		export.Item{Name: "power_filter_tuned_5_7", Build: func() (*circuit.Circuit, error) { return Circuit(tuned) }, Control: tran},
	)
}
