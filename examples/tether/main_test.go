// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// M9n: the tether's Lorentz-force EMF is vBL; as a generator the drag
// power F·v equals −EMF·I (the orbit supplies the electrical energy) and
// the load power plus the resistive losses equals EMF·I; driven above the
// EMF the force reverses into thrust and F·v = −EMF·I is the share of the
// supply's power that goes into the orbit; the momentum ledger books the
// reaction on the exterior field source.
func TestTether(t *testing.T) {
	g := Defaults()
	if want := g.V * g.B * g.Length; math.Abs(g.EMF()-want) > 1e-12*want {
		t.Errorf("EMF %g, want vBL %g", g.EMF(), want)
	}
	for _, c := range []Tether{g, func() Tether { x := g; x.Load, x.Supply = 0, 3000; return x }()} {
		r, err := c.Solve()
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.PMechanical+r.EMF*r.Current) > 1e-9*math.Abs(r.EMF*r.Current) {
			t.Errorf("F·v = %g W, −EMF·I = %g W", r.PMechanical, -r.EMF*r.Current)
		}
		loss := r.Current * r.Current * (c.Resistance() + 2*c.RContactor)
		if c.Load > 0 {
			if r.Force[0] >= 0 || math.Abs(r.PElectrical+loss-r.EMF*r.Current) > 1e-9*r.EMF*r.Current {
				t.Errorf("generator: force %v, load %g + losses %g vs EMF·I %g", r.Force, r.PElectrical, loss, r.EMF*r.Current)
			}
		} else if r.Force[0] <= 0 || math.Abs(r.PElectrical-loss-r.PMechanical) > 1e-9*r.PElectrical {
			t.Errorf("thruster: force %v, supply %g vs losses %g + orbit %g", r.Force, r.PElectrical, loss, r.PMechanical)
		}
		if n := r.Momentum.Net(); n != r.Force || r.LedgerResidue > 1e-12 {
			t.Errorf("momentum ledger net %v, force %v, energy residual %g", n, r.Force, r.LedgerResidue)
		}
	}
}
