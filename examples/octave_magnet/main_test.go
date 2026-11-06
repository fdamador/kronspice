// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// The current sheet's on-axis field equals the exact solenoid formula, for
// the uniform and the four-segment bar.
func TestSheetAxis(t *testing.T) {
	for _, b := range []Bar{Uniform(), FourPairs()} {
		sh := b.Sheet(2000)
		for i := 0; i <= 40; i++ {
			z := -1.5*b.L + 3*b.L*float64(i)/40
			_, bz, _ := b.Field(0, z, sh)
			if d := math.Abs(bz-b.AxisExact(z)) / (Mu0 * b.M0); d > 1e-5 {
				t.Errorf("%s z %.4f: off by %.2g of μ0M0", b.Name, z, d)
			}
		}
	}
}

// The elliptic loop field equals Biot–Savart off the axis.
func TestLoopFieldTwin(t *testing.T) {
	a := 6.35e-3
	l := em.Loop{P: em.Circle([3]float64{0, 0, 0}, [3]float64{0, 0, 1}, a, 4000), I: 1}
	for _, p := range [][2]float64{{1.7, 0.6}, {0.4, 0.3}, {3, -2}, {1.2, 0}} {
		rho, z := p[0]*a, p[1]*a
		bs := l.B([3]float64{rho, 0, z})
		br, bz, _ := LoopField(a, 1, rho, z)
		n := math.Hypot(br, bz)
		if math.Abs(br-bs[0])/n > 1e-5 || math.Abs(bz-bs[2])/n > 1e-5 {
			t.Errorf("(%.1f, %.1f)a: elliptic (%.6g, %.6g), Biot–Savart (%.6g, %.6g)", p[0], p[1], br, bz, bs[0], bs[2])
		}
	}
}

// Poles: a uniform bar has charge on its two end faces only (one pair, net
// zero); four reversed segments have five charge sites (net zero).
func TestPoles(t *testing.T) {
	for _, c := range []struct {
		b     Bar
		sites int
	}{{Uniform(), 2}, {FourPairs(), 5}} {
		ps := c.b.Poles()
		var q float64
		for _, p := range ps {
			q += p.Q
		}
		if len(ps) != c.sites || math.Abs(q) > 1e-9 {
			t.Errorf("%s: %d sites, net charge %.3g", c.b.Name, len(ps), q)
		}
	}
}

// The bench scan beside the bar: one sign change per charge site, 2 for one
// pair and 5 for four pairs end to end; on the axis 0 and 3.
func TestScan(t *testing.T) {
	gap := Uniform().A + 3e-3
	for _, c := range []struct {
		b            Bar
		side, onaxis int
	}{{Uniform(), 2, 0}, {FourPairs(), 5, 3}} {
		sh := c.b.Sheet(1000)
		_, s := c.b.Scan(gap, 601, sh)
		_, a := c.b.Scan(0, 601, sh)
		if SignChanges(s) != c.side || SignChanges(a) != c.onaxis {
			t.Errorf("%s: %d beside, %d on axis", c.b.Name, SignChanges(s), SignChanges(a))
		}
	}
}

// The kronspice search coil reads the flux the field model gives:
// ∫V dt = ΣM·I = 2πψ.
func TestSearchCoil(t *testing.T) {
	b := Uniform()
	sh := b.Sheet(200)
	cf, mf, err := b.SearchCoil(b.A+3e-3, 0.3*b.L, sh)
	if err != nil {
		t.Fatal(err)
	}
	_, _, psi := b.Field(b.A+3e-3, 0.3*b.L, sh)
	if math.Abs(cf/mf-1) > 1e-6 || math.Abs(mf/(2*math.Pi*psi)-1) > 1e-6 {
		t.Errorf("circuit %.8g, mutual %.8g, 2πψ %.8g Wb", cf, mf, 2*math.Pi*psi)
	}
}
