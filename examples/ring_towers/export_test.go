// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/report"
)

// withAC gives every tower radiator an AC specification equal to its SIN
// waveform (amplitude Amp, phase φ − 90°, since sin(ωt + φ) is the phasor
// at φ − 90° in SPICE's cosine convention). The transient is unchanged;
// the results sheet can then solve the ring at its resonance.
func withAC(c *circuit.Circuit, s State, err error) (*circuit.Circuit, error) {
	if err != nil {
		return nil, err
	}
	for _, r := range s.Rad {
		k := c.Index("I" + r.Name)
		if k < 0 {
			return nil, fmt.Errorf("radiator I%s not found", r.Name)
		}
		c.Elements[k].ACMag, c.Elements[k].ACPhase = Amp, r.Phase-90
	}
	return c, nil
}

// The four towers synchronized on the sectored ring biased 100 V (A–B)
// and 50 V (C–D), with main's transient (T/64 steps to 8 amplitude-decay
// times), torn at its four DC-blocking capacitors (45°, 135°, 225°, 315°)
// into the four tower sectors: once at the ring frequency (the RF passes
// the blocking capacitors, so the links carry the waves between sectors)
// and once at the DC operating point (each sector holds its own bias, the
// links carry nothing). The plain, unbiased ring is not exported: it has
// no DC path to ground (only capacitors to ground and current sources;
// main runs it with UIC), so neither the operating point nor the AC
// analysis can be computed for it. Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	f := Freq()
	T := 1 / f
	lap := N * math.Sqrt(Lsec*Csec)
	tEnd := math.Ceil(8*lap/Loss/T) * T
	tran := []string{fmt.Sprintf(".tran %.17g %.17g 0 %.17g uic", T/64, tEnd, T/64)}
	four := FourStates()[2]
	var blocks []string
	for k := 0; k < N; k++ {
		if sectorEdge(k) {
			blocks = append(blocks, fmt.Sprintf("CB%d", k))
		}
	}
	sectored := func() (*circuit.Circuit, error) {
		c, err := build(four, DefaultBias(100, 50, true), true)
		return withAC(c, four, err)
	}
	export.Test(t, "export",
		export.Item{Name: "ring_towers_four_sectored_rf", Build: sectored, Control: tran, Opt: report.Options{Freq: f, Tear: blocks}},
		export.Item{Name: "ring_towers_four_sectored_bias", Build: sectored, Control: tran, Opt: report.Options{Tear: blocks}},
	)
}
