// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func peaks(src []srcPoint) (sigma, source, q, net, abs float64) {
	for _, p := range src {
		sigma = math.Max(sigma, p.sigma)
		source = math.Max(source, math.Abs(p.source))
		q += p.sigma * p.area
		net += p.source * p.area
		abs += math.Abs(p.source) * p.area
	}
	return
}

// The source is IMV's: total charge Q, peak |Ω|σ_max e^(−½)/σ (eq. B3),
// zero ring integral (eq. 24).
func TestSource(t *testing.T) {
	d := DefaultDisc()
	sigma, source, q, net, abs := peaks(d.Sources(0.3))
	if math.Abs(q/d.Q-1) > 1e-12 {
		t.Errorf("charge %.6g against %.6g", q, d.Q)
	}
	want := math.Abs(d.Omega) * sigma * math.Exp(-0.5) / d.Sigma()
	if math.Abs(source/want-1) > 2e-3 {
		t.Errorf("peak source %.6g against B3 %.6g", source, want)
	}
	if math.Abs(net) > 1e-12*abs {
		t.Errorf("ring integral %.3g of %.3g", net, abs)
	}
}

// Smoothing the sector from σ to √(σ² + λ²) at fixed charge lowers the peak
// source by σ²/(σ² + λ²) (IMV eq. B7).
func TestSmoothing(t *testing.T) {
	a := DefaultDisc()
	b := a
	b.FWHMDeg = 2 * a.FWHMDeg
	_, sa, _, _, _ := peaks(a.Sources(0))
	_, sb, _, _, _ := peaks(b.Sources(0))
	s1, s2 := a.Sigma(), b.Sigma()
	if got, want := sb/sa, s1*s1/(s2*s2); math.Abs(got/want-1) > 5e-3 {
		t.Errorf("smoothing ratio %.5f against %.5f", got, want)
	}
}

// The signatures: S odd in Ω and linear in QΩ, zero for the whole ring lit;
// the S² force even in Ω. Standard prediction: S = 0 throughout.
func TestSignatures(t *testing.T) {
	d := DefaultDisc()
	x := [3]float64{0.03, 0.02, 0.015}
	s0 := SAt(d.Sources(0), x, d.soft())
	if s0 == 0 {
		t.Fatal("no S at the probe")
	}
	r := d
	r.Omega = -d.Omega
	if got := SAt(r.Sources(0), x, d.soft()); math.Abs(got+s0) > 1e-14*math.Abs(s0) {
		t.Errorf("reversed S %.6g against %.6g", got, -s0)
	}
	two := d
	two.Omega, two.Q = 2*d.Omega, 3*d.Q
	if got := SAt(two.Sources(0), x, d.soft()); math.Abs(got/(6*s0)-1) > 1e-12 {
		t.Errorf("S not ∝ QΩ: %.6g against %.6g", got, 6*s0)
	}
	u := d
	u.Uniform = true
	if got := SAt(u.Sources(0), x, d.soft()); got != 0 {
		t.Errorf("whole ring lit: S = %.3g", got)
	}
	cp := MMFit()
	gF := SGrid(d.Sources(0), 0.06, 0.006, d.soft())
	gR := SGrid(r.Sources(0), 0.06, 0.006, d.soft())
	fF, fR := gF.Force(x, cp, 0.003), gR.Force(x, cp, 0.003)
	for a := 0; a < 3; a++ {
		if math.Abs(fF[a]-fR[a]) > 1e-12*math.Abs(fF[2]) {
			t.Errorf("S² force not even in Ω: %v against %v", fF, fR)
		}
	}
}

// The quasi-static kernel has the zero-net-source (dipole) far field: S
// falls as 1/R² well away from the disc.
func TestFarField(t *testing.T) {
	d := DefaultDisc()
	src := d.Sources(0)
	near := SAt(src, [3]float64{0, 2, 0.5}, d.soft())
	far := SAt(src, [3]float64{0, 4, 1}, d.soft())
	if got := near / far; math.Abs(got/4-1) > 0.02 {
		t.Errorf("S(R)/S(2R) = %.4f, dipole 4", got)
	}
}
