// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"math"
	"testing"
)

// Bewley 1961 Eqs. 20‑46, 20‑47: the phase-frame inductances of a
// salient-pole machine are L₁ + L₂·cos 2θ (self) and −M₁ + L₂·cos 2(θ−60°)
// (mutual a–b; Bewley prints θ+60° with phase b measured the other way), and Park's transformation returns them to the constant
// diag(ld, lq, l0) at every rotor angle.
func TestSalientPoleSpec(t *testing.T) {
	const ld, lq, l0 = 1.1e-3, 0.7e-3, 0.2e-3
	spec, err := NewSalientPoleSpec(ld, lq, l0)
	if err != nil {
		t.Fatal(err)
	}
	l := squareZero(3)
	L1, L2, M1 := (ld+lq+l0)/3, (ld-lq)/3, (ld+lq-2*l0)/6
	s := math.Sqrt(2.0 / 3.0)
	for _, th := range []float64{0, 0.3, 1.1, 2.5, -0.8} {
		spec.Inductance(th, l, nil, nil)
		if w := L1 + L2*math.Cos(2*th); math.Abs(l[0][0]-w) > 1e-15 {
			t.Errorf("θ=%g: Laa = %g, want %g", th, l[0][0], w)
		}
		if w := -M1 + L2*math.Cos(2*(th-math.Pi/3)); math.Abs(l[0][1]-w) > 1e-15 {
			t.Errorf("θ=%g: Lab = %g, want %g", th, l[0][1], w)
		}
		// Park: rows d, q, 0 at θ.
		c1, c2, c3 := math.Cos(th), math.Cos(th-2*math.Pi/3), math.Cos(th+2*math.Pi/3)
		s1, s2, s3 := math.Sin(th), math.Sin(th-2*math.Pi/3), math.Sin(th+2*math.Pi/3)
		p := [3][3]float64{{s * c1, s * c2, s * c3}, {-s * s1, -s * s2, -s * s3}, {s / math.Sqrt2, s / math.Sqrt2, s / math.Sqrt2}}
		want := [3]float64{ld, lq, l0}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				var v float64
				for a := 0; a < 3; a++ {
					for b := 0; b < 3; b++ {
						v += p[i][a] * l[a][b] * p[j][b]
					}
				}
				w := 0.0
				if i == j {
					w = want[i]
				}
				if math.Abs(v-w) > 1e-15 {
					t.Errorf("θ=%g: P·L·Pᵀ[%d][%d] = %g, want %g", th, i, j, v, w)
				}
			}
		}
	}
	if err := spec.CheckPositive(); err != nil {
		t.Errorf("CheckPositive: %v", err)
	}
	if _, err := NewSalientPoleSpec(1, 0, 1); err == nil {
		t.Error("lq = 0 accepted")
	}
}
