// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The four coils with their primaries, driven by an AC source, swept
// around their resonances. Netlist, KiCad schematic and results sheet go
// to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	item := func(name string, mk func() *Coil, sweep string) export.Item {
		return export.Item{Name: name, Build: func() (*circuit.Circuit, error) {
			c := mk()
			return c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, c.SectionResistance(1e6))
		}, Control: []string{sweep}}
	}
	export.Test(t, "export",
		item("pulsed_ordinary_coil", func() *Coil { return NewCoil(Solenoid(), false) }, ".ac dec 50 100k 10meg"),
		item("pulsed_bifilar_coil", func() *Coil { return NewCoil(Solenoid(), true) }, ".ac dec 50 10k 1meg"),
		item("pulsed_conical_bifilar_coil", func() *Coil { return NewCoil(Cone(), true) }, ".ac dec 50 1k 10meg"),
		item("pulsed_flat_bifilar_coil", func() *Coil { return NewCoil(Flat(), true) }, ".ac dec 50 1k 10meg"),
		export.Item{Name: "pulsed_ordinary_coil_tuned_pickup", Build: func() (*circuit.Circuit, error) {
			c := NewCoil(Solenoid(), false)
			const f1 = 1.1666e6 // the coil's first resonance (Resonance)
			return c.BuildWithPickup(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, f1,
				Pickup{Centre: 0.5 * c.Geom.Z[c.Geom.N()-1], Tuned: true, CFactor: 1, RLoad: 0.2})
		}, Control: []string{".ac lin 1401 816.6k 1.633meg"}},
	)
}
