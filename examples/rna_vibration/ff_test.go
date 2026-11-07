// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The analytic gradient of every energy part equals central differences of
// the energy, at the starting (unminimized, so every term is strained)
// geometry of each molecule.
func TestGradient(t *testing.T) {
	ms, err := Molecules()
	if err != nil {
		t.Fatal(err)
	}
	for name, m := range ms {
		x := m.Coords()
		for _, part := range []int{PartBond, PartAngle, PartTorsion, PartIntra, PartInter} {
			g := make([]float64, len(x))
			m.Energy(x, part, g)
			const h = 1e-6
			worst, scale := 0.0, 1e-3
			for i := range x {
				x0 := x[i]
				x[i] = x0 + h
				ep := m.Energy(x, part, nil)
				x[i] = x0 - h
				em := m.Energy(x, part, nil)
				x[i] = x0
				num := (ep - em) / (2 * h)
				worst = math.Max(worst, math.Abs(num-g[i]))
				scale = math.Max(scale, math.Abs(g[i]))
			}
			if worst > 1e-6*scale {
				t.Errorf("%s part %d: gradient error %.3g of %.3g", name, part, worst, scale)
			}
		}
	}
}
