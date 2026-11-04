// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// A UIC start can impose inductor currents that disagree with each other:
// here a 1 nH inductor in series with the 1 H one is given the opposite
// current, as a winding wired against its loop is. At t = 0 the solve routes
// the difference through the only other path, the 10 MΩ reference
// resistors, and books about 3×10⁷ W of heat that never happens: within the
// first step the small inductor follows the large one. The integrated ledger
// must leave that imposed point out (examples/earth_machine found it as a
// −3.75×10⁴ J mismatch in a run whose whole budget was a few kJ).
func TestUICImposedPointNotIntegrated(t *testing.T) {
	// Floating loop a →L1→ b →L2→ c →R1→ a: L1 = 1 H carries 1 A; L2 =
	// 1 nH is wound c → b, so its 1 A initial current opposes the loop.
	// τ = 1 s; the stored ½L1·I² = 0.5 J all goes to heat.
	c, err := circuit.NewBuilder("uic loop").
		InductorIC("L1", "a", "b", 1, 1).
		InductorIC("L2", "c", "b", 1e-9, 1).
		Resistor("R1", "c", "a", 1).
		Resistor("RA", "a", "0", 1e7).
		Resistor("RB", "b", "0", 1e7).
		Resistor("RC", "c", "0", 1e7).
		BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 0.01, Stop: 10, MaxStep: 0.01, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	l := wf.Ledger
	if d := math.Abs(l.StoredStart - 0.5); d > 1e-8 {
		t.Errorf("stored start %g J, want 0.5", l.StoredStart)
	}
	want := 0.5 * (1 - math.Exp(-20))
	if d := math.Abs(l.EnergyDissipated - want); d > 1e-3*want {
		t.Errorf("heat %g J, want %g", l.EnergyDissipated, want)
	}
	if m := math.Abs(l.IntegrationMismatch); m > 1e-3*0.5 {
		t.Errorf("integration mismatch %g J of 0.5 J", l.IntegrationMismatch)
	}
}
