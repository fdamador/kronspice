// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/export"
)

// The proton's free precession in the laboratory frame and the
// continuous-wave run at resonance in the rotating frame. Netlist, KiCad
// schematic and results sheet go to export/ (KRONSPICE_EXPORT=1).
func TestExport(t *testing.T) {
	export.Test(t, "export",
		export.Item{Name: "kron_nmr_precession", Build: func() (*circuit.Circuit, error) {
			s := Proton()
			s.Minit = [3]float64{1, 0, 0}
			return s.Circuit()
		}, Control: []string{".tran 50p 8u uic"}},
		export.Item{Name: "kron_nmr_cw_resonance", Build: func() (*circuit.Circuit, error) {
			s := Proton()
			s.B1, s.Wrf = 1e-5, s.Gamma*s.B0
			return s.Circuit()
		}, Control: []string{".tran 40n 60u uic"}},
	)
}
