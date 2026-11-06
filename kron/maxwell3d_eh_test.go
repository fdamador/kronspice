// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// The Euler–Heisenberg vacuum on Kron's three-dimensional network: a
// cavity in a magnetar-strength magnetic background B₀ = 0.1 B_c
// (4.4×10⁸ T), the strongest field the weak-field expansion takes.
//
// About a static background the QED vacuum is a linear anisotropic
// medium for small fields (em.EHLinearize); with B₀ along a grid axis its
// tensors are diagonal, so in Kron's network the background only changes
// element values: the capacitors of links along B₀ by 1 + 10a, across it
// by 1 − 4a; the inductors of faces normal to B₀ by 1/(1 − 12a), the
// others by 1/(1 − 4a), with a = (α/90π)(B₀/B_c)² = 2.6×10⁻⁷. The
// topology, Faraday's and Ampère's laws, is untouched. The background is
// a declared input, as an operating point is for kronspice's AC analysis.
//
// The TM110 mode (E along z) is then detuned by
//
//	ω² = (k̃_x²·ν_yy + k̃_y²·ν_xx)/ε_zz,
//
// k̃ the lattice wavenumbers and ν = 1/μ. With B₀ along z the light's
// electric field lies along B₀: ω = ω₀/n∥, n∥ = 1 + 7a (Toll 1952; Adler
// 1971). With B₀ along x it lies across B₀, and its magnetic field is
// partly along and partly across it. The circuit's measured shifts must
// equal the lattice formula, and the lattice formula the continuum one.
// The circuit resolves ω to about 1e-8 (the transformers' finite L_w at
// Ideal = 1e7, and the pole fit), so the 2 % allowed on shifts of 0.6 and
// 1.8 × 10⁻⁶ is that floor, not a looser physics.
func TestMaxwellNetworkEulerHeisenberg(t *testing.T) {
	n := [3]int{7, 6, 4}
	const dx = 0.01
	eps0 := 1 / (mu0Test * c0Test * c0Test)
	g, err := NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		t.Fatal(err)
	}
	a, b := float64(n[0]-1)*dx, float64(n[1]-1)*dx
	kxl, kyl := 2/dx*math.Sin(math.Pi*dx/(2*a)), 2/dx*math.Sin(math.Pi*dx/(2*b))
	kxc, kyc := math.Pi/a, math.Pi/b
	b0 := 0.1 * em.CriticalB
	aEH := em.EHBirefringenceParameter(b0)

	type medium struct {
		name    string
		eps, mu [3]float64 // relative factors per axis
		w       float64    // measured resonance [rad/s]
	}
	mk := func(name string, bg [3]float64) medium {
		s, err := em.EHLinearize([3]float64{}, bg)
		if err != nil {
			t.Fatal(err)
		}
		md := medium{name: name}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				if i != j && (s.Eps[i][j] != 0 || s.Nu[i][j] != 0) {
					t.Fatalf("%s: tensors not diagonal on the grid axes", name)
				}
			}
			md.eps[i] = s.Eps[i][i] / em.Eps0
			md.mu[i] = 1 / (em.Mu0 * s.Nu[i][i])
		}
		return md
	}
	media := []medium{
		mk("no background", [3]float64{}),
		mk("B₀ = 0.1 B_c along z", [3]float64{0, 0, b0}),
		mk("B₀ = 0.1 B_c along x", [3]float64{b0, 0, 0}),
	}
	// ω² = (k_x²·ν_yy + k_y²·ν_xx)/ε_zz in units of the vacuum's.
	omega := func(md medium, kx, ky float64) float64 {
		return math.Sqrt((kx*kx/md.mu[1] + ky*ky/md.mu[0]) / md.eps[2] / (kx*kx + ky*ky))
	}
	for q := range media {
		md := &media[q]
		m := &MaxwellNetwork{Grid: g, Eps: eps0, Mu: mu0Test, Sigma: 1e-7, EpsAxis: md.eps, MuAxis: md.mu,
			Ideal: 1e7, Leak: 1, SeriesR: 1e-9}
		bl := circuit.NewBuilder("Euler–Heisenberg vacuum in Kron's network: " + md.name)
		if err := m.Build(bl); err != nil {
			t.Fatal(err)
		}
		src := g.Link(2, 3, 2, 1)
		bl.Add(circuit.Element{Name: "IJ", Kind: circuit.CurrentSource, ACMag: 1e-3}, m.Point(src, 0), m.Point(src, 1))
		c, err := bl.Build()
		if err != nil {
			t.Fatal(err)
		}
		wl := c0Test * math.Hypot(kxl, kyl)
		md.w = 2 * math.Pi * maxwellPole(t, m, c, g.Link(2, 3, 3, 2), wl/(2*math.Pi), 6)
	}
	w0 := media[0].w
	for _, md := range media[1:] {
		got := md.w/w0 - 1
		lattice := omega(md, kxl, kyl) - 1
		continuum := omega(md, kxc, kyc) - 1
		t.Logf("%s: ω/ω₀ − 1 circuit %.5g, lattice %.5g, continuum %.5g (a = %.4g; −7a = %.5g)",
			md.name, got, lattice, continuum, aEH, -7*aEH)
		if math.Abs(got-lattice) > 0.02*math.Abs(lattice) {
			t.Errorf("%s: circuit shift %.5g against the lattice's %.5g", md.name, got, lattice)
		}
		if math.Abs(lattice-continuum) > 0.01*math.Abs(continuum) {
			t.Errorf("%s: lattice shift %.5g against the continuum's %.5g", md.name, lattice, continuum)
		}
	}
	// Light polarized along B₀ sees n∥ = 1 + 7a.
	if r := (media[1].w/w0 - 1) / (-7 * aEH); math.Abs(r-1) > 0.02 {
		t.Errorf("B₀ along z: shift is %.4g of −7a", r)
	}
}
