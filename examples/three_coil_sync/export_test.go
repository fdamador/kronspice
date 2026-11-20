// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The hand-picked baseline and the two evolved designs (seed 1) as netlists,
// KiCad schematics and results sheets in export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	item := func(name string, d Design) export.Item {
		return export.Item{Name: name, Build: func() (*circuit.Circuit, error) {
			c, _, err := d.Build()
			return c, err
		}, Control: []string{fmt.Sprintf(".tran %g %g 0 %g", maxStep, tStop, maxStep)}}
	}
	export.Test(t, "export", item("three_coil_baseline", Baseline()), item("three_coil_evolved_free", EvolvedFree()), item("three_coil_evolved_sync", EvolvedSync()))
}
