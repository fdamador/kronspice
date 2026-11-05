// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The ether's curvature B_n for Rañada's Hopf map is Rañada's knotted
// magnetic field, B = (2(y − xz), −2(x + yz), x² + y² − z² − 1)/(1 + r²)³,
// times −16, at every point: the magnetic term of the ether energy is the
// knotted light's magnetic energy.
func TestEtherCurvatureIsKnottedLight(t *testing.T) {
	const e = 1e-6
	bn := func(x, y, z float64) [3]float64 {
		var d [3][3]float64
		p := [3]float64{x, y, z}
		for a := 0; a < 3; a++ {
			q, r := p, p
			q[a] += e
			r[a] -= e
			u, v := HopfMap(q[0], q[1], q[2], 1), HopfMap(r[0], r[1], r[2], 1)
			for c := 0; c < 3; c++ {
				d[a][c] = (u[c] - v[c]) / (2 * e)
			}
		}
		n := HopfMap(x, y, z, 1)
		return [3]float64{dot(n, cross(d[1], d[2])), dot(n, cross(d[2], d[0])), dot(n, cross(d[0], d[1]))}
	}
	for _, p := range [][3]float64{{0.3, -0.2, 0.5}, {-0.9, 0.4, 0.1}, {1.4, 1.1, -0.7}, {0.05, 0.02, 2.2}} {
		x, y, z := p[0], p[1], p[2]
		d := math.Pow(1+x*x+y*y+z*z, 3)
		ran := [3]float64{2 * (y - x*z) / d, -2 * (x + y*z) / d, (x*x + y*y - z*z - 1) / d}
		b := bn(x, y, z)
		for c := 0; c < 3; c++ {
			if math.Abs(b[c]+16*ran[c]) > 1e-6*(1+math.Abs(16*ran[c])) {
				t.Errorf("at %v component %d: B_n %.8g, −16·B_Rañada %.8g", p, c, b[c], -16*ran[c])
			}
		}
	}
}

// The analytic gradients of both discrete energies match finite
// differences, so the relaxation follows the energy it reports.
func TestGradient(t *testing.T) {
	g := NewGrid(11, 3, 1, 1.3, 0.7)
	var gr [3][]float64
	for c := 0; c < 3; c++ {
		gr[c] = make([]float64, len(g.n[c]))
	}
	g.Energy(&gr)
	worst := 0.0
	for _, m := range []int{g.idx(5, 5, 5), g.idx(4, 6, 3), g.idx(2, 7, 8)} {
		for c := 0; c < 3; c++ {
			o, e := g.n[c][m], 1e-6
			g.n[c][m] = o + e
			ep := g.Total()
			g.n[c][m] = o - e
			em := g.Total()
			g.n[c][m] = o
			fd := (ep - em) / (2 * e)
			worst = math.Max(worst, math.Abs(fd-gr[c][m])/(math.Abs(fd)+1e-9))
		}
	}
	if worst > 1e-6 {
		t.Errorf("3-D gradient off by %.2g", worst)
	}
	a := NewAxi(3, 3, 0.25, 1, 1.3, 0.7)
	var ga [3][]float64
	for c := 0; c < 3; c++ {
		ga[c] = make([]float64, len(a.m[c]))
	}
	a.Energy(&ga)
	worst = 0
	for _, ij := range [][2]int{{0, 12}, {1, 12}, {3, 10}, {5, 14}, {0, 5}, {9, 20}} {
		m := a.at(ij[0], ij[1])
		for c := 0; c < 3; c++ {
			o, e := a.m[c][m], 1e-6
			a.m[c][m] = o + e
			ep := a.Total()
			a.m[c][m] = o - e
			em := a.Total()
			a.m[c][m] = o
			fd := (ep - em) / (2 * e)
			worst = math.Max(worst, math.Abs(fd-ga[c][m])/(math.Abs(fd)+1e-6))
		}
	}
	if worst > 1e-5 {
		t.Errorf("axisymmetric gradient off by %.2g", worst)
	}
}

// The axisymmetric energy of a field equals the 3-D energy of the same
// field rotated about the axis (two discretizations of one functional):
// E2 to 1 %; E4 to 6 %, because the 3-D grid (h = 0.1) is twice as coarse
// as the axisymmetric one (0.05) and central differences underestimate
// the quartic term by O(h²).
func TestAxiEqualsThreeD(t *testing.T) {
	a := NewAxi(5, 5, 0.05, 1, 1, 1)
	e2a, e4a := a.Energy(nil)
	g := a.To3D(101, 5)
	e2g, e4g := g.Energy(nil)
	if math.Abs(e2g/e2a-1) > 0.01 || math.Abs(e4g/e4a-1) > 0.06 {
		t.Errorf("E2 %.3f vs %.3f, E4 %.3f vs %.3f", e2a, e2g, e4a, e4g)
	}
	t.Logf("axisymmetric E2 %.3f E4 %.3f; 3-D E2 %.3f E4 %.3f", e2a, e4a, e2g, e4g)
}

// Coarse grid: with both terms the knot settles (energy stops falling)
// near the virial balance (E2/E4 = 0.83 at h = 0.15, rising to 0.997 as
// h → 0 in the run's grid study) and keeps a core of order ℓ; with the
// magnetic term alone it spreads, with the stiffness alone it shrinks.
func TestEtherHoldsKnot(t *testing.T) {
	if testing.Short() {
		t.Skip("relaxation")
	}
	a := NewAxi(8, 8, 0.15, 1, 1, 1)
	e2, e4, n := a.Relax(20000, 0.05, 1e-7, 0, nil)
	core, _ := a.Size()
	if n >= 20000 || math.Abs(e2/e4-1) > 0.2 || core < 0.8 || core > 1.6 {
		t.Errorf("both terms: steps %d, E2/E4 %.4f, core %.3f", n, e2/e4, core)
	}
	b := NewAxi(8, 8, 0.15, 1, 0, 1)
	_, r0 := b.Size()
	b.Relax(2000, 0.05, 0, 0, nil)
	_, r1 := b.Size()
	if r1 <= r0*1.1 {
		t.Errorf("magnetic only: rms radius %.3f → %.3f, want spreading", r0, r1)
	}
	c := NewAxi(8, 8, 0.15, 1, 1, 0)
	c0, _ := c.Size()
	s0 := c.Total()
	c.Relax(100, 0.05, 0, 0, nil)
	c1, _ := c.Size()
	if c1 >= c0*0.9 {
		t.Errorf("stiffness only: core %.3f → %.3f, want shrinking", c0, c1)
	}
	c.Relax(1000, 0.05, 0, 0, nil)
	if e := c.Total(); e > 0.05*s0 {
		t.Errorf("stiffness only: energy %.3f of %.3f after 1100 steps, want the knot gone through the lattice", e, s0)
	}
	t.Logf("both: E %.3f E2/E4 %.4f core %.3f; magnetic only rms %.3f → %.3f; stiffness only core %.3f → %.3f",
		e2+e4, e2/e4, core, r0, r1, c0, c1)
}
