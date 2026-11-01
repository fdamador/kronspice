// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// The base design of the tests (200-turn bifilar coil, 40 mm plates) as the
// dry (Seebeck) and wet (galvanic) element, and the dry one wound
// ordinarily for comparison, each built as Measure builds it (coil
// resistance at the resonance) and swept like Resonance; the sheet is at
// the tank's parallel resonance f0 (the wet tank has none; see below). (main's GA fronts come from a long run,
// so the fixed base design stands in.) Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	var items []export.Item
	f0 := 0.0 // the dry bifilar tank's resonance, also used for the wet one
	for _, x := range []struct {
		name    string
		v       Variant
		bifilar bool
	}{
		{"bial_bifilar_dry", Dry, true},
		{"bial_bifilar_wet", Wet, true},
		{"bial_bifilar_dry_ordinary", Dry, false},
	} {
		d := designFrom(x.v, DefaultAssumptions(), x.bifilar, baseX())
		r, err := Measure(d, analysis.Options{})
		switch {
		case errors.Is(err, ErrNoResonance) && f0 > 0:
			// The wet element's electrolyte resistance and double layer
			// damp the tank: no parallel resonance. Take its DC coil
			// resistance and the dry tank's f0 as the analysis point.
			r.F0, r.Rcoil = f0, d.Coil.Resistance(0)
		case err != nil:
			t.Fatalf("%s: %v", x.name, err)
		}
		if f0 == 0 {
			f0 = r.F0
		}
		items = append(items, export.Item{
			Name:    x.name,
			Build:   func() (*circuit.Circuit, error) { return d.Build(r.Rcoil, false) },
			Control: []string{".ac dec 40 1 10g"},
			Opt:     sheet.Options{Freq: r.F0},
			// On the parallel resonance (Q ≈ 190–240) the 1 A probe sees
			// |Z| ≈ 1e5 Ω while the element's series resistance is ~1 mΩ:
			// the torn solve is backward stable (componentwise backward
			// error ~1e-16, TearStats), but the forward difference to MNA is
			// that times a condition number of order Q·1e8 value spread, so
			// ~3e-9 on resonance against 1e-11 at 1.1·f0 and 4e-13 at
			// 0.01·f0 (checked when this test was written).
			Tol: 1e-7,
		})
	}
	export.Test(t, "export", items...)
}
