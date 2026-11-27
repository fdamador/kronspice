// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/internal/export"
)

// The argon Marx generator (charging and firing) and the argon varactor
// pump. Netlist, KiCad schematic and results sheet go to export/
// (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	m, p := DefaultMarx(Argon), DefaultPump(Argon)
	ece := MarxECE(Argon)[1]
	export.Test(t, "export",
		export.Item{Name: "marx_argon_charge", Src: m.ChargeNetlist()},
		export.Item{Name: "marx_argon_fire", Src: m.FireNetlist()},
		export.Item{Name: "marx_argon_fire_ece2_torsion", Src: ece.FireNetlist()},
		export.Item{Name: "varactor_pump_argon", Src: p.Netlist()})
}
