// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
	sheet "github.com/fdamador/kronspice/report"
)

// The matched water chain with the bubble's fundamental admittance at the
// centre node (RBUBBLE, CBUBBLE), swept ±0.5 % around the mode; the sheet
// is taken at the mode frequency. Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	d0, err := readDrive(0, 1.3*atm)
	if err != nil {
		t.Fatal(err)
	}
	b := DefaultBubble(d0.FPeak) // as main does
	cy, _, err := b.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	c, f0, err := chain(cy.Y1)
	if err != nil {
		t.Fatal(err)
	}
	export.Test(t, "export", export.Item{
		Name:    "bubble_on_flask",
		Build:   func() (*circuit.Circuit, error) { return c, nil },
		Control: []string{fmt.Sprintf(".ac lin 4001 %.17g %.17g", 0.995*f0, 1.005*f0)},
		Opt:     sheet.Options{Freq: f0},
	})
}
