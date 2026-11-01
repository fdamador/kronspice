// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// radiatorAmps is the radiator amplitude Solve uses: the four radiators
// scaled so they deliver h.Pin at the ring resonance.
func radiatorAmps(t *testing.T, h Hull) float64 {
	ckt, err := h.circuit(1e-3)
	if err != nil {
		t.Fatal(err)
	}
	f := h.Freq()
	ac, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	var pin float64
	for _, e := range ckt.Elements {
		if e.Kind != circuit.CurrentSource {
			continue
		}
		is, _ := ac.ElementCurrent(e.Name)
		v, _ := ac.NodeVoltage(ckt.NodeNames[e.Nodes[1]])
		pin += 0.5 * real(v[0]*cmplx.Conj(is[0]))
	}
	return 1e-3 * math.Sqrt(h.Pin/pin)
}

// The example hull with 2.5 kV bias, traveling and standing wave, at the
// ring resonance (the AC point Solve uses), radiators at 1 kW. Netlist,
// KiCad schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	var items []export.Item
	for _, st := range []bool{false, true} {
		h := Example()
		h.Standing, h.Vdc = st, 2500
		amps := radiatorAmps(t, h)
		name := "double_hull_traveling"
		if st {
			name = "double_hull_standing"
		}
		f := h.Freq()
		items = append(items, export.Item{
			Name:    name,
			Build:   func() (*circuit.Circuit, error) { return h.circuit(amps) },
			Control: []string{fmt.Sprintf(".ac lin 1 %.17g %.17g", f, f)},
			Opt:     sheet.Options{Freq: f},
			// Q ≈ 9,400 and the sheet is on resonance: the forward error is
			// the backward error (~1e-15) times a condition number of order
			// Q·n (n = 96 sections), so torn and untorn agree to ~6e-9. Off
			// resonance the twin is exact (1e-12 at 1.01·f0, 6e-14 at 0.9·f0;
			// checked when this test was written).
			Tol: 1e-7,
		})
	}
	export.Test(t, "export", items...)
}
