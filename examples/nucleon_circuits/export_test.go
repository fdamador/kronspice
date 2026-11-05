// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The ring circuits: the proton as one light ring, the neutron as proton
// and electron rings (light model), and the leading spin arrangement of
// the uud and udd quark-ring nucleons. Netlist, KiCad schematic and results
// sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	item := func(name string, rings []Ring) export.Item {
		return export.Item{Name: name, Control: []string{".tran 10p 2n"}, Build: func() (*circuit.Circuit, error) { return Build(rings) }}
	}
	rq := 2 * MuP * MuN0 / (E * C0)
	reN := 2 * (MuP - MuN) * MuN0 / (E * C0)
	R := LightRadius(MpMeV / MuP)
	export.Test(t, "export",
		item("proton_one_ring", []Ring{{Name: "p", Q: 1, R: LightRadius(MpMeV), Up: 1}}),
		item("neutron_p_e_rings", []Ring{{Name: "p", Q: 1, R: rq, Up: 1}, {Name: "e", Q: -1, R: reN, Up: 1}}),
		item("proton_uud", QuarkNucleon(2.0/3, -1.0/3, R, R)[0].Rings),
		item("neutron_udd", QuarkNucleon(-1.0/3, 2.0/3, R, R)[0].Rings),
	)
}
