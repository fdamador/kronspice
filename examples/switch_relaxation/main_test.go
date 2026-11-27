// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Each period agrees with the exact piecewise-exponential value. The
// tolerance is the time-step resolution of the turn-on instants (local
// maxima of the recorded points, steps of a few µs against periods of
// 4.5 and 17.9 ms); the crossings themselves are located by the switch.
func TestPeriods(t *testing.T) {
	r, err := Run(".")
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range r.Cells {
		if c.Cycles < 3 || math.Abs(c.Simulated/c.Exact-1) > 5e-4 {
			t.Errorf("cell %s: period %g s over %d cycles, exact %g s", c.Node, c.Simulated, c.Cycles, c.Exact)
		}
	}
	if r.Mismatch > 1e-3 {
		t.Errorf("ledger integration mismatch %g of the energy in", r.Mismatch)
	}
}
