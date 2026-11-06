// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/output"
)

// maxwellCircuit builds Kron's network on an n-node box of cell dx with
// a 1 mA impressed current on the z link at src. σ = 1e-7 S/m (loss
// tangent ≈ 5e-7 at the TM110 frequency) joins the spheres at DC for the
// operating point; the Yee twin carries it exactly.
func maxwellCircuit(t *testing.T, n [3]int, dx float64, ties TieRule, ideal float64, src [3]int) (*MaxwellNetwork, *circuit.Circuit, []float64) {
	t.Helper()
	g, err := NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		t.Fatal(err)
	}
	m := &MaxwellNetwork{Grid: g, Eps: 1 / (mu0Test * c0Test * c0Test), Mu: mu0Test, Sigma: 1e-7,
		Ties: ties, Ideal: ideal, Leak: 1, SeriesR: 1e-9}
	b := circuit.NewBuilder("Kron 1962 Fig. 1: Maxwell network with ideal transformers")
	if err := m.Build(b); err != nil {
		t.Fatal(err)
	}
	k := g.Link(2, src[0], src[1], src[2])
	if k < 0 || !g.free(k) {
		t.Fatalf("source link %v not free", src)
	}
	b.Add(circuit.Element{Name: "IJ", Kind: circuit.CurrentSource, ACMag: 1e-3}, m.Point(k, 0), m.Point(k, 1))
	c, err := b.Build()
	if err != nil {
		t.Fatal(err)
	}
	j := make([]float64, g.Links())
	j[k] = 1e-3
	return m, c, j
}

// linkVoltage is link k's E·Δ in an AC result (frequency i).
func linkVoltage(t *testing.T, m *MaxwellNetwork, r *output.ACResult, k, i int) complex128 {
	t.Helper()
	a, ok1 := r.NodeVoltage(m.Point(k, 0))
	z, ok2 := r.NodeVoltage(m.Point(k, 1))
	if !ok1 || !ok2 {
		t.Fatalf("link %d points missing", k)
	}
	return a[i] - z[i]
}

// maxwellTwin solves the circuit by MNA at f and returns the largest
// departure from the Yee twin of the link voltages and of every arc's
// current (relative to the largest of each), and the ledger residual.
func maxwellTwin(t *testing.T, m *MaxwellNetwork, c *circuit.Circuit, j []float64, f float64) (dv, di, ledger float64) {
	t.Helper()
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	vy, iy, err := m.MaxwellYee(2*math.Pi*f, j)
	if err != nil {
		t.Fatal(err)
	}
	vmax, imax := 0.0, 0.0
	for _, v := range vy {
		vmax = math.Max(vmax, cmplx.Abs(v))
	}
	for _, i := range iy {
		imax = math.Max(imax, cmplx.Abs(i))
	}
	for k := range vy {
		if m.Grid.free(k) {
			dv = math.Max(dv, cmplx.Abs(linkVoltage(t, m, r, k, 0)-vy[k])/vmax)
		}
	}
	for p := range iy {
		if !m.Faces(p) {
			continue
		}
		for cn := 0; cn < 4; cn++ {
			cur, ok := r.ElementCurrent(m.Arc(p, cn))
			if !ok {
				t.Fatalf("no arc %s", m.Arc(p, cn))
			}
			di = math.Max(di, cmplx.Abs(cur[0]-iy[p])/imax)
		}
	}
	return dv, di, r.Ledger.MaxRelativeResidual
}

// The three-dimensional network with its ideal transformers is the
// space-discrete Yee system. Solved by kronspice's MNA, every link voltage
// E·Δ and every quarter inductor's current equal the direct solve of
// Ampère and Faraday on the same grid, the four arcs of a face carrying
// the face current H·Δ; the departure is the transformers' magnetizing
// current and falls as 1/Ideal (Kron's L = M = ∞ as a limit). This holds
// for Kron's own arrangement (two transformers per large mesh, diagonally
// opposite, Balasubramanian p. 153) as well as for three ties per mesh:
// the lattice of spheres fixes the remaining arc currents. Ledger closed.
func TestMaxwellNetworkIsYee(t *testing.T) {
	n := [3]int{6, 5, 4}
	const dx = 0.01
	a, b := float64(n[0]-1)*dx, float64(n[1]-1)*dx
	wl := c0Test * math.Hypot(2/dx*math.Sin(math.Pi*dx/(2*a)), 2/dx*math.Sin(math.Pi*dx/(2*b)))
	fo := 0.83 * wl / (2 * math.Pi)
	for _, tc := range []struct {
		name string
		ties TieRule
	}{{"Kron (two diagonal ties per mesh)", KronTies}, {"chain (three ties per mesh)", ChainTies}} {
		var dv, di [2]float64
		for i, ideal := range []float64{1e5, 1e7} {
			m, c, j := maxwellCircuit(t, n, dx, tc.ties, ideal, [3]int{2, 2, 1})
			var lg float64
			dv[i], di[i], lg = maxwellTwin(t, m, c, j, fo)
			if lg > 1e-9 {
				t.Errorf("%s, Ideal %g: AC ledger relative residual %.3g", tc.name, ideal, lg)
			}
			if i == 1 {
				t.Logf("%s: %d links, %d faces, %d elements", tc.name, m.Grid.Links(), len(m.Grid.Plaquettes), len(c.Elements))
			}
		}
		if dv[1] > 1e-6 || di[1] > 1e-6 {
			t.Errorf("%s, Ideal 1e7: circuit vs Yee %.2g (E·Δ), %.2g (H·Δ)", tc.name, dv[1], di[1])
		}
		// 1/Ideal: a hundredfold larger winding inductance, about a
		// hundredfold smaller departure.
		if dv[0]/dv[1] < 30 || di[0]/di[1] < 30 {
			t.Errorf("%s: departure does not fall as 1/Ideal: E·Δ %.2g → %.2g, H·Δ %.2g → %.2g", tc.name, dv[0], dv[1], di[0], di[1])
		}
		t.Logf("%s: circuit − Yee, E·Δ %.2g → %.2g, H·Δ (every arc) %.2g → %.2g for Ideal 1e5 → 1e7",
			tc.name, dv[0], dv[1], di[0], di[1])
	}
	// Leaky windings (k = 1 − 1e-8) with their leakage taken off the
	// quarter inductors are the same network: still Yee at Ideal = 1e7.
	m, c, j := maxwellCircuit(t, n, dx, KronTies, 1e7, [3]int{2, 2, 1})
	m.Leakage = 1e-8
	bl := circuit.NewBuilder("leaky windings")
	if err := m.Build(bl); err != nil {
		t.Fatal(err)
	}
	src := m.Grid.Link(2, 2, 2, 1)
	bl.Add(c.Elements[len(c.Elements)-1], m.Point(src, 0), m.Point(src, 1))
	cl, err := bl.Build()
	if err != nil {
		t.Fatal(err)
	}
	dvl, dil, _ := maxwellTwin(t, m, cl, j, fo)
	if dvl > 1e-6 || dil > 1e-6 {
		t.Errorf("leaky windings, compensated: circuit vs Yee %.2g (E·Δ), %.2g (H·Δ)", dvl, dil)
	}
	t.Logf("leaky windings (k = 1 − 1e-8, leakage off the arcs): circuit − Yee %.2g (E·Δ), %.2g (H·Δ)", dvl, dil)
}

// maxwellPole finds the resonance near f0 from |V| on a probe link, which
// peaks at the pole: narrowing sweeps of 21 points from ±5 %, each round
// five times finer, then the vertex of the V-shaped 1/|V| (|V| ≈
// A/|ω − ω0| near a nearly lossless pole) from a line fitted on each side.
func maxwellPole(t *testing.T, m *MaxwellNetwork, c *circuit.Circuit, probe int, f0 float64, rounds int) float64 {
	t.Helper()
	lo, hi := 0.95*f0, 1.05*f0
	var fr float64
	for round := 0; round < rounds; round++ {
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 21, FStart: lo, FStop: hi}, analysis.Options{SkipEnergyCheck: true})
		if err != nil {
			t.Fatal(err)
		}
		fs := r.Freq
		inv := make([]float64, len(fs))
		best := 0
		for i := range fs {
			inv[i] = 1 / cmplx.Abs(linkVoltage(t, m, r, probe, i))
			if inv[i] < inv[best] {
				best = i
			}
		}
		if best < 2 || best > len(fs)-3 {
			t.Fatalf("no peak inside [%g, %g]", lo, hi)
		}
		// The pole lies between fs[best−1] and fs[best+1]; 1/|V| falls
		// linearly to it from each side, with slopes that differ when other
		// modes add a background. Intersect the line through the two
		// samples left of the pole with the one through the two right of it.
		sl := (inv[best-1] - inv[best-2]) / (fs[best-1] - fs[best-2])
		sr := (inv[best+2] - inv[best+1]) / (fs[best+2] - fs[best+1])
		fr = (inv[best+1] - inv[best-1] + sl*fs[best-1] - sr*fs[best+1]) / (sl - sr)
		step := fs[1] - fs[0]
		lo, hi = fs[best]-2*step, fs[best]+2*step
	}
	return fr
}

// TM110 resonance of the network, Kron's arrangement: the circuit's pole
// is the lattice's own (Yee spatial dispersion), its continuum error is
// that dispersion's, and em.FDTD on the same grid rings at the same ω
// once its leapfrog time step is accounted for,
// ω_FDTD = (2/dt)·asin(ω·dt/2).
func TestMaxwellNetworkCavity(t *testing.T) {
	n := [3]int{7, 6, 4}
	const dx = 0.01
	m, c, _ := maxwellCircuit(t, n, dx, KronTies, 1e7, [3]int{3, 2, 1})
	g := m.Grid
	a, b := float64(n[0]-1)*dx, float64(n[1]-1)*dx
	wl := c0Test * math.Hypot(2/dx*math.Sin(math.Pi*dx/(2*a)), 2/dx*math.Sin(math.Pi*dx/(2*b)))
	wc := c0Test * math.Pi * math.Hypot(1/a, 1/b)

	fr := maxwellPole(t, m, c, g.Link(2, 3, 3, 2), wl/(2*math.Pi), 4)
	w := 2 * math.Pi * fr
	if e := math.Abs(w-wl) / wl; e > 1e-6 {
		t.Errorf("TM110: circuit ω %.9g, lattice %.9g (%.2g)", w, wl, e)
	}
	if e := math.Abs(w/wc - wl/wc); e > 1e-6 {
		t.Errorf("TM110: continuum error %.4g, the lattice dispersion's %.4g", w/wc-1, wl/wc-1)
	}

	f, err := em.NewFDTD(em.FDTDConfig{N: n, D: [3]float64{dx, dx, dx}})
	if err != nil {
		t.Fatal(err)
	}
	f.SetE(func(x, y, z float64) [3]float64 {
		return [3]float64{0, 0, math.Sin(math.Pi*x/a) * math.Sin(math.Pi*y/b)}
	})
	px, py, pz := 3*dx, 3*dx, 2.5*dx
	prev := 0.0
	var cross []float64
	for f.Steps < 6000 {
		f.Step()
		e, _ := f.FieldAt(px, py, pz)
		if f.Steps > 1 && (prev > 0) != (e[2] > 0) {
			cross = append(cross, f.Time()-f.Dt*e[2]/(e[2]-prev))
		}
		prev = e[2]
	}
	wf := math.Pi * float64(len(cross)-1) / (cross[len(cross)-1] - cross[0])
	wcf := 2 / f.Dt * math.Asin(w*f.Dt/2)
	if e := math.Abs(wf-wcf) / wcf; e > 2e-6 {
		t.Errorf("FDTD ω %.9g, circuit ω through the leapfrog %.9g (%.2g)", wf, wcf, e)
	}
	t.Logf("TM110 (%d×%d×%d cells): circuit ω/ω_lattice − 1 = %.2g; continuum error %.4g (lattice %.4g); FDTD vs circuit through the leapfrog %.2g (uncorrected %.2g)",
		n[0]-1, n[1]-1, n[2]-1, w/wl-1, w/wc-1, wl/wc-1, wf/wcf-1, wf/w-1)
}
