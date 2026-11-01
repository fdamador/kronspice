// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// M9n: B from a changing or rotating E field. Rowland's disk matches its
// analytic axis field; a rotating uniform E field's B satisfies the
// Ampère–Maxwell law with its displacement current; a symmetric charged
// ring makes B only when its charge moves; a rotating charge pair's B is
// consistent with the displacement current of its rotating E.
func TestRotatingE(t *testing.T) {
	w := 20000 * 2 * math.Pi / 60
	for _, z := range []float64{0, 0.05} {
		n, a := RowlandB(1e-5, 0.1, w, z)
		if math.Abs(n-a) > 1e-4*a {
			t.Errorf("Rowland B(%g) = %g, analytic %g", z, n, a)
		}
	}
	b, amp, far := RotatingUniformE(1e6, 2*math.Pi*1e6, 0.05)
	if want := em.Mu0 * em.Eps0 / 2 * 1e6 * 2 * math.Pi * 1e6 * 0.05; math.Abs(b-want) > 1e-12*want || amp > 1e-8 || far > 1e-5 {
		t.Errorf("rotating E: |B| %g (want %g), Ampère residual %g, Faraday %g", b, want, amp, far)
	}
	rest, spin := RingCases(1e-6, 0.05, 2*math.Pi*1000)
	if want := em.Mu0 * (1e-6 * 1000) / (2 * 0.05); rest != 0 || math.Abs(spin-want) > 1e-4*want {
		t.Errorf("ring: at rest %g, spinning %g (want μ₀I/2a = %g)", rest, spin, want)
	}
	if bc, ar := RotatingPair(1e-9, 0.05, 2*math.Pi*1000); !(bc > 0) || ar > 1e-3 {
		t.Errorf("rotating pair: |B| %g, Ampère–Maxwell residual %g", bc, ar)
	}
}

// Röntgen–Eichenwald (electric mirror of the homopolar disk): a dielectric
// spinning in a fixed E field is magnetized, M = P × v, so its mid-plane
// field is −(εr − 1)·E·v/c². The factor 1/c² is μ₀ε₀ and is not inserted
// by hand: the ring sum uses only the bound charge ε₀(εr − 1)E and
// Biot–Savart's μ₀. Tolerances: the ring sum is converged in rings and
// segments to 1e−5 (checked at 4× resolution); what remains is the finite
// thickness, which falls linearly in h/R (≈ 0.9 h/R at ρ = R/2), so the
// test checks 1.5 % at h/R = 0.01 and that halving h halves the deficit.
func TestRontgenEichenwald(t *testing.T) {
	const epsR, e0, R, rho = 3.0, 1e6, 0.1, 0.05
	w := 2 * math.Pi * 50
	c2 := em.SpeedOfLight * em.SpeedOfLight
	n1, thin := RontgenEichenwald(epsR, e0, R, 0.002, w, rho)
	n2, _ := RontgenEichenwald(epsR, e0, R, 0.001, w, rho)
	if want := -(epsR - 1) * e0 * w * rho / c2; math.Abs(thin-want) > 1e-12*math.Abs(want) {
		t.Errorf("thin-disk B %g, want −(εr−1)Ev/c² = %g", thin, want)
	}
	d1, d2 := 1-n1/thin, 1-n2/thin
	if d2 < 0 || d2 > 0.015 || math.Abs(d1/d2-2) > 0.05 {
		t.Errorf("ring-sum deficit %g (h = 2 mm), %g (h = 1 mm): want < 1.5 %% and ratio 2", d1, d2)
	}
}
