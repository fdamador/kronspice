// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The bench reproduces the paper's own numbers where it states them, then
// records the comparisons with measurement (H13 stays "untested" until a
// Zeeman prediction exists; these are the quantities it must also meet).
func TestPaperNumbers(t *testing.T) {
	m := FromCharge()
	if math.Abs(m.Re-4e-11)/4e-11 > 0.01 {
		t.Errorf("r_e = %.4g cm, paper ≈ 4e-11", m.Re)
	}
	// Paper: r_e = 4e-11 cm in eq. 4 gives 4.7326e-10.
	if q := (Model{Re: 4e-11}).Charge(); math.Abs(q-4.7326e-10)/4.7326e-10 > 1e-3 {
		t.Errorf("eq. 4 at r_e = 4e-11: %.5g, paper 4.7326e-10", q)
	}
	// Paper: void volume 2.67e-31 cm³ at r_e = 4e-11.
	if v := 4 * math.Pi / 3 * math.Pow(4e-11, 3); math.Abs(v-2.68e-31)/2.68e-31 > 0.01 {
		t.Errorf("void volume %.4g", v)
	}
	// Paper eqs. 10–12: 1 g ≈ 8.4, 8.8 and 8.6 × 10⁶ cm⁴/s; eq. 11 is
	// m_e(cm⁴/s)/m_e(g) with the paper's rounded r_e = 4e-11 cm.
	if g := (Model{Re: 4e-11}).MassUnits() / Me; math.Abs(g-8.8e6)/8.8e6 > 0.02 {
		t.Errorf("eq. 11: 1 g = %.4g cm⁴/s, paper 8.8e6", g)
	}
	if g := m.GFactor(); math.Abs(g-15.0/8) > 1e-12 {
		t.Errorf("g = %.15g, want 15/8", g)
	}
	if a := m.Alpha(); math.Abs(a-15*math.Pi*math.Pi/16) > 1e-9 {
		t.Errorf("α = %.10g, UPSM eq. 5.2 gives 15π²/16", a)
	}
	for _, c := range m.Checks() {
		t.Logf("%-48s model/measured = %.4g", c.Name, c.Ratio())
	}
}
