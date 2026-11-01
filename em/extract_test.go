// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/material"
)

func within(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol*math.Abs(want) {
		t.Errorf("%s = %.6g, want %.6g (%.2f %% off, limit %.2g %%)", what, got, want, 100*(got-want)/want, 100*tol)
	} else {
		t.Logf("%s = %.6g, analytic %.6g (%+.3f %%)", what, got, want, 100*(got-want)/want)
	}
}

// Coaxial line (CLAUDE.md §3c rule 3): C′ = 2πε/ln(b/a), L′ = μ₀/2π·ln(b/a),
// with a dielectric of εr = 2.1 — within 0.5 % on a grid of 20 cells per
// inner radius (curved boundaries cut at their true position).
func TestCoaxSection(t *testing.T) {
	const a, b, er = 1e-3, 3.5e-3, 2.1
	c0 := 4e-3
	s := Section{W: 2 * c0, H: 2 * c0, Step: a / 20,
		Conductors: []Conductor2D{Disk("in", c0, c0, a), Outside("gnd", c0, c0, b)},
		EpsR:       func(x, y float64) float64 { return er }}
	r, err := s.Solve()
	if err != nil {
		t.Fatal(err)
	}
	within(t, "coax C′", r.C[0][0], 2*math.Pi*Eps0*er/math.Log(b/a), 5e-3)
	within(t, "coax L′", r.L[0][0], Mu0/(2*math.Pi)*math.Log(b/a), 5e-3)
	line := FromSection(r, 0, 1, 0, 0)
	within(t, "coax Z₀", real(line.Z0(0)), math.Sqrt(Mu0/(Eps0*er))/(2*math.Pi)*math.Log(b/a), 5e-3)
}

// Stripline: a zero-thickness strip of width w midway between ground
// planes b apart has the exact Z₀ = (η₀/4)·K(k)/K(k′), k = sech(πw/2b),
// k′ = tanh(πw/2b) (S. B. Cohn, IRE Trans. MTT 2, 1954), so
// C′ = 4ε₀·K(k′)/K(k) in vacuum.
// This is also the parallel-plate case with its fringing, exactly.
func TestStriplineSection(t *testing.T) {
	const bgap, w = 2e-3, 2e-3
	k := 1 / math.Cosh(math.Pi*w/(2*bgap))
	K, _ := EllipticKE(k)
	Kp, _ := EllipticKE(math.Sqrt(1 - k*k))
	want := 4 * Eps0 * Kp / K
	// The field is singular (∝ r^−½) at a zero-thickness edge, so the
	// error falls only linearly with the step (2.5, 1.24, 0.62, 0.31 % at
	// 20 … 160 cells per gap); Richardson extrapolation of two grids
	// removes that first-order term.
	var c [2]float64
	for i, n := range []float64{40, 80} {
		W := 12 * bgap
		s := Section{W: W, H: bgap, Step: bgap / n,
			Conductors: []Conductor2D{Strip("s", W/2-w/2, W/2+w/2, bgap/2)}}
		r, err := s.Solve()
		if err != nil {
			t.Fatal(err)
		}
		c[i] = r.C[0][0]
	}
	within(t, "stripline C′ (b/80)", c[1], want, 1e-2)
	within(t, "stripline C′ (extrapolated)", 2*c[1]-c[0], want, 1e-3)
}

// Two-wire line: C′ = πε₀/acosh(D/2a) and L′ = (μ₀/π)·acosh(D/2a)
// between the wires, from the 2×2 matrices of the two wires in a large
// grounded box.
func TestTwoWireSection(t *testing.T) {
	const a, D = 1e-3, 6e-3
	box := 120e-3 // the grounded box adds ≈ (D/box)² of error
	s := Section{W: box, H: box, Step: a / 4,
		Conductors: []Conductor2D{Disk("p", box/2-D/2, box/2, a), Disk("n", box/2+D/2, box/2, a)}}
	r, err := s.Solve()
	if err != nil {
		t.Fatal(err)
	}
	c := r.C
	cdiff := -c[0][1] + (c[0][0]+c[0][1])*(c[1][1]+c[0][1])/(c[0][0]+c[1][1]+2*c[0][1])
	within(t, "two-wire C′", cdiff, math.Pi*Eps0/math.Acosh(D/(2*a)), 5e-3)
	ldiff := r.L[0][0] + r.L[1][1] - 2*r.L[0][1]
	within(t, "two-wire L′", ldiff, Mu0/math.Pi*math.Acosh(D/(2*a)), 5e-3)
}

// Microstrip on FR-4 (εr = 4.4, h = 1.6 mm, w = 3 mm) against the
// Hammerstad–Jensen formulas (E. Hammerstad, Ø. Jensen, "Accurate models
// for microstrip computer-aided design," IEEE MTT-S Digest 1980): Z₀ and
// ε_eff within 2 % (the box and the zero-thickness strip are the
// approximations on our side).
func TestMicrostripSection(t *testing.T) {
	const h, w, er = 1.6e-3, 3e-3, 4.4
	step := h / 16
	W, H := 40e-3, 20e-3
	y := step * math.Round(h/step)
	s := Section{W: W, H: H, Step: step,
		Conductors: []Conductor2D{Strip("trace", W/2-w/2, W/2+w/2, y)},
		EpsR: func(x, yy float64) float64 {
			if yy < y {
				return er
			}
			return 1
		}}
	r, err := s.Solve()
	if err != nil {
		t.Fatal(err)
	}
	z0 := math.Sqrt(r.L[0][0] / r.C[0][0])
	eeff := r.C[0][0] / r.C0[0][0]
	zhj, ehj := HammerstadJensen(w, h, er)
	within(t, "microstrip Z₀", z0, zhj, 2e-2)
	within(t, "microstrip ε_eff", eeff, ehj, 2e-2)
}

// Coaxial circles: Neumann's double integral over 360-gons equals
// Maxwell's elliptic-integral form; a solenoid current sheet integrated
// from it equals Nagaoka's closed form (1e−6), and Wheeler's formula
// (em.Coil) is within its stated 1 % for l ≥ 0.8r.
func TestCoilInductance(t *testing.T) {
	for _, c := range [][3]float64{{0.05, 0.03, 0.02}, {0.1, 0.1, 0.005}, {0.02, 0.08, 0.1}} {
		l1 := Loop{P: Circle([3]float64{}, [3]float64{0, 0, 1}, c[0], 360), I: 1}
		l2 := Loop{P: Circle([3]float64{0, 0, c[2]}, [3]float64{0, 0, 1}, c[1], 360), I: 1}
		within(t, "Neumann M", Neumann(l1, l2), MutualCircles(c[0], c[1], c[2]), 1e-4)
	}
	for _, g := range [][2]float64{{0.02, 0.01}, {0.02, 0.04}, {0.02, 0.2}} {
		R, l := g[0], g[1]
		within(t, "current sheet vs Nagaoka", SolenoidSheet(R, l, 100), Nagaoka(R, l, 100), 1e-6)
	}
	for _, turns := range []int{20, 60, 200} {
		coil := Coil{Turns: turns, Diameter: 0.04, WireDiameter: 0.5e-3, Pitch: 1e-3, InsulationEps: 3, Conductor: material.Copper}
		within(t, "Wheeler vs Nagaoka", coil.Inductance(), Nagaoka(0.02, coil.Length(), float64(turns)), 1e-2)
	}
}

// HammerstadJensen is Z₀ and ε_eff of a zero-thickness microstrip.
func HammerstadJensen(w, h, er float64) (float64, float64) {
	u := w / h
	eta0 := math.Sqrt(Mu0 / Eps0)
	f := 6 + (2*math.Pi-6)*math.Exp(-math.Pow(30.666/u, 0.7528))
	z0air := eta0 / (2 * math.Pi) * math.Log(f/u+math.Sqrt(1+4/(u*u)))
	a := 1 + math.Log((math.Pow(u, 4)+math.Pow(u/52, 2))/(math.Pow(u, 4)+0.432))/49 + math.Log(1+math.Pow(u/18.1, 3))/18.7
	b := 0.564 * math.Pow((er-0.9)/(er+3), 0.053)
	eeff := (er+1)/2 + (er-1)/2*math.Pow(1+10/u, -a*b)
	return z0air / math.Sqrt(eeff), eeff
}
