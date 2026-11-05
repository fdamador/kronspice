// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// The Gauss linking integral: two linked circles (a Hopf link) give 1,
// two separated circles give 0.
func TestLinkingKnownCurves(t *testing.T) {
	circle := func(cx, cy, cz float64, plane int) [][3]float64 {
		var c [][3]float64
		for i := 0; i < 400; i++ {
			u := 2 * math.Pi * float64(i) / 400
			p := [3]float64{math.Cos(u), math.Sin(u), 0}
			if plane == 1 { // in the x–z plane
				p = [3]float64{math.Cos(u), 0, math.Sin(u)}
			}
			c = append(c, [3]float64{p[0] + cx, p[1] + cy, p[2] + cz})
		}
		return c
	}
	a := circle(0, 0, 0, 0)
	if l := Linking(a, circle(1, 0, 0, 1)); math.Abs(math.Abs(l)-1) > 1e-3 {
		t.Errorf("Hopf link: linking %.4f, want ±1", l)
	}
	if l := Linking(a, circle(5, 0, 0, 1)); math.Abs(l) > 1e-3 {
		t.Errorf("separated circles: linking %.4f, want 0", l)
	}
}

// The hopfion is a null field (E·B = 0 and |E| = c|B|, i.e. F·F = 0) that
// satisfies Maxwell's equations with the sign the code picks, and its two
// traced B lines link once, at t = 0 and after it has spread (t = 3).
func TestHopfionExact(t *testing.T) {
	h := Hopfion{}
	if h.Sign() != 1 && h.Sign() != -1 {
		t.Fatalf("sign %g", h.Sign())
	}
	for _, p := range [][4]float64{{0.3, -0.2, 0.5, 0}, {-0.7, 0.1, 0.2, 1.3}, {1.5, 1, -0.4, 2.5}} {
		F := h.F(p[0], p[1], p[2], p[3])
		ff := F[0]*F[0] + F[1]*F[1] + F[2]*F[2]
		norm := real(F[0]*cmplx.Conj(F[0]) + F[1]*cmplx.Conj(F[1]) + F[2]*cmplx.Conj(F[2]))
		if cmplx.Abs(ff) > 1e-12*norm {
			t.Errorf("not null at %v: |F·F| = %.2e of |F|² %.2e", p, cmplx.Abs(ff), norm)
		}
	}
	for _, tt := range []float64{0, 3} {
		b := func(q [3]float64) [3]float64 {
			F := h.F(q[0], q[1], q[2], tt)
			return [3]float64{imag(F[0]), imag(F[1]), imag(F[2])}
		}
		l1, c1 := trace(b, Seeds[0], 0.01, 40)
		l2, c2 := trace(b, Seeds[1], 0.01, 40)
		if !c1 || !c2 {
			t.Fatalf("t = %g: exact field lines did not close", tt)
		}
		if l := Linking(l1, l2); math.Abs(math.Abs(l)-1) > 1e-3 {
			t.Errorf("t = %g: exact linking %.4f, want ±1", tt, l)
		}
	}
}

// The knot in the solver (coarse grid): the energy leaves the centre as the
// exact solution says, the helicity (the linking, in integral form) is
// kept as the exact solution keeps it, curl A = B to round-off, and the
// energy ledger closes.
func TestKnotInSolver(t *testing.T) {
	if testing.Short() {
		t.Skip("3-D field solve")
	}
	k := DefaultKnot()
	k.Cells, k.Half, k.Times = 6, 5, []float64{0, 1, 3}
	r, err := k.Run()
	if err != nil {
		t.Fatal(err)
	}
	last := r.Samples[len(r.Samples)-1]
	if last.InsideFDTD > 0.05 || math.Abs(last.InsideFDTD-last.InsideExact) > 0.02 {
		t.Errorf("t = 3: energy within 2 L0 %.4f (exact %.4f), want spread out and matching", last.InsideFDTD, last.InsideExact)
	}
	if math.Abs(last.HelicityFDTD-last.HelicityExact) > 0.02 || last.HelicityFDTD < 0.95 {
		t.Errorf("t = 3: helicity %.4f (exact %.4f), want kept", last.HelicityFDTD, last.HelicityExact)
	}
	if math.Abs(r.H0Ratio-1) > 0.05 {
		t.Errorf("starting helicity %.4f of exact", r.H0Ratio)
	}
	for _, s := range r.Samples {
		if s.CurlAErr > 1e-12 {
			t.Errorf("t = %g: curl A − B = %.1e of max B", s.T, s.CurlAErr)
		}
	}
	if math.Abs(r.Residual) > 1e-10 {
		t.Errorf("ledger residual %.1e", r.Residual)
	}
}
