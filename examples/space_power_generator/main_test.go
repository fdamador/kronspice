// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol*math.Abs(want) {
		t.Errorf("%s = %.6g, want %.6g", what, got, want)
	}
}

// The book's own numbers (pp. 159–160) and the ledger's bound.
func TestPublishedTest(t *testing.T) {
	l := Published().Analyse()
	near(t, "W_r (6.20)", l.Wr, 5677, 1e-9)
	near(t, "W_E (6.21)", l.WE, 11724, 1e-3)
	near(t, "η incremental (6.22)", l.EtaInc, 2.066, 1e-3)
	near(t, "η SPG (6.23)", l.EtaSPG, 2.49, 2e-3)
	near(t, "ratio to all inputs", l.TotalRatio, 1.265, 1e-3)
	near(t, "largest loaded EMF", l.ELoadMax, 1.3706, 1e-3)
	near(t, "largest loaded EMF, lossless DM", l.ELoadMaxIdeal, 1.6513, 1e-3)
	// Consistency: at the bound, E_load·I equals the shaft power and the
	// torque k_load·I equals the measured torque.
	near(t, "E_load·I = shaft power", l.ELoadMax*Published().I, l.ShaftIncrease, 1e-12)
	near(t, "torque at the bound", l.ELoadMax*Published().I/l.Omega, l.TorqueMeasured, 1e-12)
}
