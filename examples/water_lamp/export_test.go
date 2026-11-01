// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The reduced bifilar coil with its primary, driven by a 2 µs current
// pulse (the twin-test circuit). Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	c := DefaultCoil()
	build := func() (*circuit.Circuit, error) {
		b := circuit.NewBuilder("bifilar coil, first mode, pulsed primary")
		b.CurrentWave("IP", "0", "p", circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 0.2e-6, 0.2e-6, 2e-6, 1}})
		b.Inductor("LP", "p", "0", c.Lp)
		b.Inductor("LS", "s", "0", c.Ls)
		b.Capacitor("CS", "s", "0", c.Cs)
		b.Resistor("RS", "s", "0", c.Rs)
		b.Mutual("K1", "LP", "LS", c.M/math.Sqrt(c.Lp*c.Ls))
		return b.BuildTransient()
	}
	export.Test(t, "export", export.Item{Name: "bifilar_first_mode", Build: build, Control: []string{".tran 0.1u 240u"}})
}
