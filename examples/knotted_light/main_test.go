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

// H30: the knot's energy density at t = 0 is spherically symmetric (so the
// radial integral of knotEnergy is its total), and that total agrees with
// a 3-D sum over a large cube.
func TestKnotEnergy(t *testing.T) {
	h := Hopfion{}
	u := func(x, y, z float64) float64 {
		F := h.F(x, y, z, 0)
		return real(F[0]*cmplx.Conj(F[0]) + F[1]*cmplx.Conj(F[1]) + F[2]*cmplx.Conj(F[2]))
	}
	for _, r := range []float64{0.3, 0.9, 2.5} {
		ref := u(r, 0, 0)
		for _, d := range [][3]float64{{0, 1, 0}, {0, 0, 1}, {1, 1, 1}, {-0.3, 0.8, -0.5}} {
			n := math.Sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
			if v := u(r*d[0]/n, r*d[1]/n, r*d[2]/n); math.Abs(v-ref) > 1e-12*ref {
				t.Errorf("r = %g, direction %v: %.6g against %.6g", r, d, v, ref)
			}
		}
	}
	U := knotEnergy()
	var s float64
	const h3, L = 0.05, 8.0
	for x := -L + h3/2; x < L; x += h3 {
		for y := -L + h3/2; y < L; y += h3 {
			for z := -L + h3/2; z < L; z += h3 {
				s += u(x, y, z) * h3 * h3 * h3
			}
		}
	}
	if math.Abs(s/U-1) > 2e-3 { // the cube misses the tail beyond 8 L0
		t.Errorf("3-D sum %.6g against radial %.6g", s, U)
	}
}

// H30, coarse grid: the incoming half enters through the Huygens surface,
// the knot is rebuilt at the centre as the exact solution says, its lines
// link once at t = 0, it leaves through the same surface, little leaks
// past it, and the ledger closes.
func TestImplosionThroughSurface(t *testing.T) {
	if testing.Short() {
		t.Skip("3-D field solve")
	}
	cy := DefaultImplosion()
	cy.Cells, cy.Half, cy.Start, cy.Times = 6, 5, -8, []float64{0, 8}
	r, err := cy.Run()
	if err != nil {
		t.Fatal(err)
	}
	peak, end := r.Samples[0], r.Samples[1]
	if peak.In < 0.9 || peak.Out > 1e-3 {
		t.Errorf("t = 0: put in %.4f, taken out %.4f knot energies; want almost all in, none out", peak.In, peak.Out)
	}
	if math.Abs(peak.InsideFDTD-peak.InsideExact) > 0.05 {
		t.Errorf("t = 0: within 2 L0 %.4f, exact %.4f", peak.InsideFDTD, peak.InsideExact)
	}
	if math.Abs(peak.Link-1) > 0.02 || peak.Gap > 0.2 {
		t.Errorf("t = 0: linking %.4f with gap %.3f L0, want 1", peak.Link, peak.Gap)
	}
	if math.Abs(end.Out-end.In) > 0.05 || end.Leak > 0.05 {
		t.Errorf("t = 8: in %.4f, out %.4f, leak %.4f", end.In, end.Out, end.Leak)
	}
	if math.Abs(r.Residual) > 1e-10 {
		t.Errorf("ledger residual %.2g", r.Residual)
	}
	t.Logf("in %.4f out %.4f leak %.2g; t = 0: %.4f vs %.4f, linking %.4f (gap %.3f); residual %.1e",
		end.In, end.Out, end.Leak, peak.InsideFDTD, peak.InsideExact, peak.Link, peak.Gap, r.Residual)
}

// H30 Part 3: the Euler–Heisenberg vacuum returns none of an explosion
// and reflects almost nothing of a wave crossing the strongest field it
// allows.
func TestVacuumInward(t *testing.T) {
	v, err := RunVacuum()
	if err != nil {
		t.Fatal(err)
	}
	// The nonlinear share of the field energy is ≈ 6a (u = ε₀E²/2 + 3ξE⁴).
	if r := v.Nonlinear / (6 * v.A); math.Abs(r-1) > 1e-2 {
		t.Errorf("nonlinear energy share %.4g, 6a = %.4g", v.Nonlinear, 6*v.A)
	}
	// The halves did act on each other while crossing …
	if !(v.Altered > 1e-14) {
		t.Errorf("no nonlinear interaction while the halves crossed (%.3g)", v.Altered)
	}
	// … but nothing comes back between them: round-off, below even the
	// linear lattice's dispersion ripple.
	if v.Returned > 1e-24 || v.Returned > v.ReturnedLn {
		t.Errorf("nonlinear energy back at the centre %.3g (linear lattice ripple %.3g)", v.Returned, v.ReturnedLn)
	}
	// The strongest counter-propagating field turns round less than the
	// Fresnel bound for its index step.
	if !(v.Reflected < v.Fresnel) {
		t.Errorf("reflected %.3g, Fresnel bound %.3g", v.Reflected, v.Fresnel)
	}
	t.Logf("a %.3g; returned %.2g (linear lattice %.2g); altered %.2g; reflected %.2g (Fresnel %.2g); nonlinear share %.3g",
		v.A, v.Returned, v.ReturnedLn, v.Altered, v.Reflected, v.Fresnel, v.Nonlinear)
}
