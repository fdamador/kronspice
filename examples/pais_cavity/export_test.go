// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/internal/export"
	"github.com/fdamador/kronspice/report"
)

// The ring of case A (four emitters at 0/90/180/270°, in phase) at the
// mode frequency, torn into the example's 16 Kron sectors. Netlist,
// KiCad schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	g := Default
	f := g.Freq()
	q := g.N / 4
	four := []Emitter{{0, 1, 0}, {q, 1, 0}, {2 * q, 1, 0}, {3 * q, 1, 0}}
	src := strings.Replace(Netlist(g, f, four), ".end\n", fmt.Sprintf(".ac lin 1 %.17g %.17g\n.end\n", f, f), 1)
	// The example's Kron sectors: cut the inductor that closes each sector.
	var cut []string
	for k := 1; k <= g.S; k++ {
		cut = append(cut, fmt.Sprintf("L%d", k*g.N/g.S-1))
	}
	export.Test(t, "export", export.Item{
		Name: "pais_ring_4emitters", Src: src,
		Opt: report.Options{Freq: f, Tear: cut},
		// On resonance Q ≈ 1.5e4: the forward error is the backward error
		// (~1e-15) times a condition number of order Q·n ≈ 1e7, so the twin
		// agrees to ~1e-8, not 1e-9 (main_test.go checks backward errors).
		Tol: 1e-6,
	})
}
