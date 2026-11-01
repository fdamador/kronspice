// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// loop builds the series loop Run builds for a case with native drive
// and scale waveforms (scale zero: constant C). The source is named VU,
// not U as in Run: circuit.WriteNetlist writes element names as given,
// and SPICE reads the first letter as the element type. The source also
// carries "AC U0" (no effect on the transient) so the results sheet can
// solve the loop at resonance: at DC the capacitor is open and node c has
// no path to ground in any torn piece.
func loop(drive, scale circuit.Waveform) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("parametric LCR (Amador & Eckardt 2012)")
	b.Add(circuit.Element{Name: "VU", Kind: circuit.VoltageSource, Wave: drive, ACMag: U0}, "s", "0")
	b.InductorIC("L", "s", "a", L, 0).Resistor("R", "a", "c", R)
	if scale.Kind == circuit.WaveNone {
		b.CapacitorIC("C", "c", "0", C, 0)
	} else {
		b.VaryingCapacitor("C", "c", "0", C, scale)
	}
	return b.BuildTransient()
}

// The paper's Fig. 4 (constant C driven at resonance) and Fig. 5
// (C(1 + 0.2 sin 2ω0t) driven at f0), with main's 30 ms transient at 100
// steps per period; the results sheet is the AC solution at f0. Netlist, KiCad schematic and results sheet go to
// export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	f0 := W0 / (2 * math.Pi)
	h := 1 / f0 / 100
	tran := []string{fmt.Sprintf(".tran %.17g 0.03 0 %.17g uic", h, h)}
	export.Test(t, "export",
		export.Item{Name: "parametric_lcr_fig4_constant_c", Build: func() (*circuit.Circuit, error) {
			return loop(sinWave(0, U0, W0, 0), circuit.Waveform{})
		}, Control: tran, Opt: sheet.Options{Freq: f0}},
		export.Item{Name: "parametric_lcr_fig5_pumped_c", Build: func() (*circuit.Circuit, error) {
			return loop(sinWave(0, U0, W0, 0), sinWave(1, 0.2, 2*W0, 0))
		}, Control: tran, Opt: sheet.Options{Freq: f0}},
	)
}
