// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The straight-wire validation case of
// CLAUDE.md §3g on Kron's three-dimensional network (stage B of
// doc/kron-figure-catalogue.md §10).

// wireH is the magnetic field [A/m] at (x, y) of a line current I along
// +z at (x0, y0) inside a square perfectly conducting tube [0, a]²:
// Ampère's B of the line plus its images (A_z = 0 on the walls, so the
// images alternate in sign across each wall). The four charges of one
// image cell form a quadrupole, so the lattice sum converges absolutely;
// p images per direction leave a remainder of order 1/p.
func wireH(I, x0, y0, a, x, y float64, p int) (hx, hy float64) {
	for i := -p; i <= p; i++ {
		for j := -p; j <= p; j++ {
			for _, s := range [4][3]float64{{1, 1, 1}, {-1, 1, -1}, {1, -1, -1}, {-1, -1, 1}} {
				xs := s[0]*x0 + 2*float64(i)*a
				ys := s[1]*y0 + 2*float64(j)*a
				dx, dy := x-xs, y-ys
				r2 := dx*dx + dy*dy
				hx += -s[2] * I * dy / (2 * math.Pi * r2)
				hy += s[2] * I * dx / (2 * math.Pi * r2)
			}
		}
	}
	return hx, hy
}

// A straight wire: an impressed current I on every z link of the column
// through the centre of a square conducting tube, from end wall to end
// wall, at 10 MHz (the box is 4×10⁻³ λ across, quasi-static). Two checks:
//
//   - Ampère's law. The circulation of H·Δ round every square of dual
//     edges about the wire equals the current through it: exactly (to
//     round-off) when the displacement and conduction currents of the
//     enclosed links are counted, and I itself to the order of those
//     currents (the quasi-static limit).
//   - Biot–Savart. H beside the wire equals the line current's field in
//     the tube (wireH), to the lattice's discretization error, which falls
//     with distance from the wire as the cells become small against it.
func TestMaxwellNetworkWire(t *testing.T) {
	const (
		nc   = 16 // cells across
		dx   = 0.005
		f    = 10e6
		amps = 1.0
	)
	n := [3]int{nc + 1, nc + 1, 3}
	g, err := NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		t.Fatal(err)
	}
	m := &MaxwellNetwork{Grid: g, Eps: 1 / (mu0Test * c0Test * c0Test), Mu: mu0Test, Sigma: 1e-7, Ideal: 1e7, Leak: 1, SeriesR: 1e-9}
	b := circuit.NewBuilder("Kron's Maxwell network: straight wire in a square tube")
	if err := m.Build(b); err != nil {
		t.Fatal(err)
	}
	w := nc / 2
	for k := 0; k < n[2]-1; k++ {
		l := g.Link(2, w, w, k)
		b.Add(circuit.Element{Name: "IW" + string(rune('0'+k)), Kind: circuit.CurrentSource, ACMag: amps}, m.Point(l, 0), m.Point(l, 1))
	}
	c, err := b.Build()
	if err != nil {
		t.Fatal(err)
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	if lg := r.Ledger.MaxRelativeResidual; lg > 1e-9 {
		t.Errorf("AC ledger relative residual %.3g", lg)
	}
	face := func(axis, i, j, k int) complex128 {
		p := g.Face(axis, i, j, k)
		cur, ok := r.ElementCurrent(m.Arc(p, 0))
		if !ok {
			t.Fatalf("no face %d (%d, %d, %d, %d)", p, axis, i, j, k)
		}
		return cur[0]
	}
	capv, _ := m.Constants()
	omega := 2 * math.Pi * f
	const k = 1 // middle plane
	worstExact, worstQS := 0.0, 0.0
	for h := 0; h < w-1; h++ {
		// Counter-clockwise round the nodes |i − w|, |j − w| ≤ h: east
		// side +H_y, north −H_x, west −H_y, south +H_x.
		var circ complex128
		for j := w - h; j <= w+h; j++ {
			circ += face(1, w+h, j, k) - face(1, w-h-1, j, k)
		}
		for i := w - h; i <= w+h; i++ {
			circ += face(0, i, w-h-1, k) - face(0, i, w+h, k)
		}
		// Current through the square: the impressed I plus, on every
		// enclosed z link, conduction and displacement (G + jωC)·V.
		enclosed := complex(amps, 0)
		for i := w - h; i <= w+h; i++ {
			for j := w - h; j <= w+h; j++ {
				l := g.Link(2, i, j, k)
				enclosed += complex(m.Sigma*capv[l]/m.Eps, omega*capv[l]) * linkVoltage(t, m, r, l, 0)
			}
		}
		worstExact = math.Max(worstExact, cmplx.Abs(circ-enclosed)/amps)
		worstQS = math.Max(worstQS, cmplx.Abs(circ-amps)/amps)
	}
	if worstExact > 1e-6 || worstQS > 1e-3 {
		t.Errorf("Ampère: circulation − enclosed current %.2g, − I %.2g (of I)", worstExact, worstQS)
	}
	t.Logf("Ampère round %d squares: circulation = enclosed current to %.2g of I, = I to %.2g (quasi-static)", w-1, worstExact, worstQS)

	// H_y on the wire's row east of it (y-normal faces, centres at
	// x = (i + ½)Δ), against the tube's line-current field.
	a := float64(nc) * dx
	x0 := float64(w) * dx
	var errs []float64
	for i := w; i < nc-1; i++ {
		x := (float64(i) + 0.5) * dx
		_, hy := wireH(amps, x0, x0, a, x, x0, 300)
		got := real(face(1, i, w, k)) / dx
		e := math.Abs(got-hy) / math.Abs(hy)
		errs = append(errs, e)
		t.Logf("r = %.1f cells: H_y circuit %.6g A/m, line current in the tube %.6g (%.2g)", float64(i-w)+0.5, got, hy, e)
	}
	// Next to the wire the lattice spreads the current over a cell (at
	// half a cell H is I/4Δ, not I/πΔ). Beyond it the error is the
	// discrete Laplacian's, second order, e ≈ κ·(Δ/r)² with one κ, and
	// below 1 % from 5.5 cells out.
	for q := 2; q < len(errs); q++ {
		r := float64(q) + 0.5
		if kappa := errs[q] * r * r; kappa < 0.2 || kappa > 0.45 {
			t.Errorf("H at %.1f cells off by %.2g: κ = e·(r/Δ)² = %.3g", r, errs[q], kappa)
		}
		if r >= 5.5 && errs[q] > 1e-2 {
			t.Errorf("H at %.1f cells off by %.2g", r, errs[q])
		}
	}
}
