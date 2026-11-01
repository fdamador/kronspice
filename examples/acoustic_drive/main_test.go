// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/fluid"
)

// At the pressure peak, the power into the liquid is P²/(2·Rmode) — the
// mode's losses at the center amplitude — and the ledger balances; the
// matching inductor raises the pressure per volt and puts the peak back at
// the mode frequency; D₂O resonates lower than H₂O.
func TestChain(t *testing.T) {
	var f0s []float64
	for _, m := range []fluid.Medium{fluid.Water, fluid.HeavyWater} {
		var perVolt [2]float64
		for i, withL := range []bool{false, true} {
			c, f0, err := chain(m, withL)
			if err != nil {
				t.Fatal(err)
			}
			r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 2001, FStart: 0.97 * f0, FStop: 1.03 * f0}, analysis.Options{})
			if err != nil {
				t.Fatal(err)
			}
			p, _ := r.NodeVoltage("p")
			best := 0
			for k := range p {
				if cmplx.Abs(p[k]) > cmplx.Abs(p[best]) {
					best = k
				}
			}
			_, _, _, rpar, _ := fluid.SphericalResonator(radius, m, q)
			im, _ := r.ElementCurrent("RMODE")
			heat := 0.5 * real(p[best]*cmplx.Conj(im[best]))
			if want := 0.5 * cmplx.Abs(p[best]) * cmplx.Abs(p[best]) / rpar; math.Abs(heat-want) > 1e-9*want {
				t.Errorf("%s: liquid heat %g W, P²/2R = %g W", m.Name, heat, want)
			}
			if r.Ledger.MaxRelativeResidual > 1e-9 {
				t.Errorf("%s: ledger residual %g", m.Name, r.Ledger.MaxRelativeResidual)
			}
			perVolt[i] = cmplx.Abs(p[best])
			if withL && math.Abs(r.Freq[best]-f0) > 1e-3*f0 {
				t.Errorf("%s matched: peak %g Hz, mode %g Hz", m.Name, r.Freq[best], f0)
			}
			if i == 0 {
				f0s = append(f0s, f0)
			}
		}
		if !(perVolt[1] > 10*perVolt[0]) {
			t.Errorf("%s: matching gives %g Pa/V against %g Pa/V", m.Name, perVolt[1], perVolt[0])
		}
	}
	if !(f0s[1] < f0s[0]) {
		t.Errorf("D2O f0 %g, H2O %g", f0s[1], f0s[0])
	}
}
