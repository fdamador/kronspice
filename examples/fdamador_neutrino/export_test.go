// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The E-ring resonator with a 1 V series AC drive, and the plain and
// crossed inverted (CL) ring lines with a 1 A AC drive at (a0, b0).
// Netlists, KiCad schematics and results sheets go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	e := DefaultERing()
	// Tol 1e-6: the 1 GΩ DC bleeders sit beside 0.29 Ω section resistors
	// (ratio ≈ 3×10⁹), so the torn and untorn solves differ at ~4e-8 from
	// conditioning alone.
	ering := export.Item{Name: "e_ring", Control: []string{".ac lin 2001 1.9g 2.2g"}, Tol: 1e-6, Build: func() (*circuit.Circuit, error) {
		return e.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, 0)
	}}
	dual := func(name string, crossed bool) export.Item {
		d := DefaultDualRing(crossed)
		return export.Item{Name: name, Control: []string{".ac lin 2001 10meg 200meg"}, Build: d.Build}
	}
	export.Test(t, "export", ering, dual("dual_ring_plain", false), dual("dual_ring_crossed", true))
}
