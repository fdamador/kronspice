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

// The circuits: Kron's radial network for Z = 1 (dual form, AC at the
// frequency of E = −0.6 Ha, off resonance) and the 21 cm loops of the
// standard atom. Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	l := NewLadder(1, 0.1, 30)
	_, _, v0, err := l.Network()
	if err != nil {
		t.Fatal(err)
	}
	f := 1 / math.Sqrt(-0.6+v0) / (2 * math.Pi)
	export.Test(t, "export",
		export.Item{Name: "kron_radial_network", Control: []string{fmt.Sprintf(".ac lin 1 %.12g %.12g", f, f)},
			Build: func() (*circuit.Circuit, error) { c, _, _, err := l.Network(); return c, err }},
		export.Item{Name: "hyperfine_loops", Control: []string{".tran 10p 2n"},
			Build: func() (*circuit.Circuit, error) {
				c, _, _, err := HyperfineBuild(Atom{Measured(), Protons()[3]}, 8)
				return c, err
			}},
	)
}
