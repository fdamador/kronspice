// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/spacetime"
)

// TestSelfField: a 1 kg point mass at 1 m shifts clocks by G/c² =
// 7.4261×10⁻²⁸ (CODATA G and c), and a disk's J is ½MR²ω.
func TestSelfField(t *testing.T) {
	m := spacetime.Machine{Radius: 2, Omega: 3}
	sf := spacetime.MachineSelfField(m, 1, 1)
	if want := spacetime.G / (spacetime.C * spacetime.C); math.Abs(sf.PhiOverC2-want)/want > 1e-15 || math.Abs(sf.PhiOverC2-7.4261e-28) > 1e-31 {
		t.Errorf("GM/dc² %.6g", sf.PhiOverC2)
	}
	if sf.J != 6 {
		t.Errorf("J %g, want 6", sf.J)
	}
}

func TestRuns(t *testing.T) { main() }
