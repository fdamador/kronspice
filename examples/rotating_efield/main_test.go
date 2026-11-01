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
