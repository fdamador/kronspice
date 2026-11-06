// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// The coaxial TEM-pulse validation case of CLAUDE.md §3g on Kron's
// three-dimensional network (stage B of doc/kron-figure-catalogue.md §10).
//
// A square line (6 × 6 cells, inner conductor 2 × 2) 24 cells long,
// shorted at both ends. The pulse is launched by an impressed current
// sheet on the cross-section at z = 6 cells, each transverse link
// carrying j = −(2/η)·Δφ·s(t), Δφ the lattice's own TEM shape (discrete
// Laplace, inner 1, outer 0): a modal current −(2/Z₀)·s(t) into the line
// seen both ways (Z₀/2) sets V = s(t) at the sheet, a Gaussian of σ = 2
// cells in time (sampled every half cell), and sends half each way. Run as a kronspice transient,
// at a probe plane 6 cells on:
//
//   - the direct pulse arrives with V ≈ 1 at the lattice's group speed,
//     c less 1/(8σ²) for a Gaussian of σ cells (3 % at σ = 2: the same
//     dispersion FDTD has), its peaks reshaped by it by a few %;
//   - V (∫E·dl from inner to outer) over I (∮H·dl round the inner
//     conductor) is the line's impedance Z₀ = η/Σ(Δφ)², the lattice's own
//     (L′C′ = με, C′ = ε·Σ(Δφ)²);
//   - the pulse sent backwards returns from the near short inverted, still
//     travelling +z (V/I = Z₀); the direct pulse returns from the far short
//     with V inverted and I kept (travelling −z);
//   - the ledger closes;
//   - em.FDTD, with the inner conductor as a perfect Conductor and the same
//     edge sources, gives the same V(t) at the probe.
func TestMaxwellNetworkCoaxPulse(t *testing.T) {
	const (
		nc, nz = 6, 24
		lo, hi = 2, 4 // inner conductor nodes
		dx     = 0.01
		ks     = 6   // source plane
		kp     = 12  // probe plane
		jp     = 3   // probe row (centre)
		width  = 2.0 // Gaussian σ [cells of time Δ/c]
		tc     = 10.0
	)
	if testing.Short() {
		t.Skip("transient of a 50,000-unknown network (minutes)")
	}
	eps := 1 / (mu0Test * c0Test * c0Test)
	eta := math.Sqrt(mu0Test / eps)
	cell := dx / c0Test
	n := [3]int{nc + 1, nc + 1, nz + 1}
	g, err := NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		t.Fatal(err)
	}
	inner := func(i, j int) bool { return i >= lo && i <= hi && j >= lo && j <= hi }
	g.Conductor(func(i, j, k int) bool { return inner(i, j) })

	// Transverse electrostatics on the lattice.
	var phi [nc + 1][nc + 1]float64
	for i := range phi {
		for j := range phi[i] {
			if inner(i, j) {
				phi[i][j] = 1
			}
		}
	}
	for it := 0; it < 20000; it++ {
		d := 0.0
		for i := 1; i < nc; i++ {
			for j := 1; j < nc; j++ {
				if inner(i, j) {
					continue
				}
				v := (phi[i-1][j] + phi[i+1][j] + phi[i][j-1] + phi[i][j+1]) / 4
				d = math.Max(d, math.Abs(v-phi[i][j]))
				phi[i][j] += 1.5 * (v - phi[i][j])
			}
		}
		if d < 1e-15 {
			break
		}
	}
	// dphi is E·Δ of the TEM shape on the transverse link along axis c
	// (0 = x, 1 = y) from node (i, j).
	dphi := func(c, i, j int) float64 {
		if c == 0 {
			return phi[i][j] - phi[i+1][j]
		}
		return phi[i][j] - phi[i][j+1]
	}
	sum := 0.0
	for i := 0; i < nc; i++ {
		for j := 0; j <= nc; j++ {
			sum += dphi(0, i, j)*dphi(0, i, j) + dphi(1, j, i)*dphi(1, j, i)
		}
	}
	z0Lattice := eta / sum
	// s(t): the Gaussian sampled every half cell of time over 0 … 2·tc
	// and zero after, linear between samples (the PWL source; FDTD is
	// given the same function). Few breakpoints keep the step control
	// from restarting every tenth of a cell.
	const dtS = 0.5
	nS := int(2 * tc / dtS)
	samp := make([]float64, nS+1)
	for r := 0; r < nS; r++ {
		u := float64(r)*dtS - tc
		samp[r] = math.Exp(-u * u / (2 * width * width))
	}
	shape := func(tt float64) float64 {
		x := tt / cell / dtS
		r := int(math.Floor(x))
		if r < 0 || r >= nS {
			return 0
		}
		return samp[r] + (x-float64(r))*(samp[r+1]-samp[r])
	}

	// Ideal = 1e3 with leaky windings (k = 1 − 1e-4, leakage taken off the
	// arcs). In a time step the spheres' attachment points are held to one
	// another only through the windings' magnetizing path, h/(2·L_w), so
	// L_w must stay small enough for that to clear the sparse LU's relative
	// pivot threshold (≈ 4e-10 at 40,000 unknowns) even on the short
	// restart steps after a source corner; the price is a departure from
	// Yee of about 1e-3 (TestMaxwellNetworkIsYee scaling).
	m := &MaxwellNetwork{Grid: g, Eps: eps, Mu: mu0Test, Sigma: 1e-7, Ideal: 1e3, Leakage: 1e-4, Leak: 1, SeriesR: 1e-9}
	b := circuit.NewBuilder("Kron's Maxwell network: square coaxial line, TEM pulse")
	if err := m.Build(b); err != nil {
		t.Fatal(err)
	}
	type src struct {
		c, i, j int
		w       float64 // amps per unit s(t)
	}
	var srcs []src
	for c := 0; c < 2; c++ {
		for i := 0; i <= nc; i++ {
			for j := 0; j <= nc; j++ {
				l := g.Link(c, i, j, ks)
				if l < 0 || !g.free(l) {
					continue
				}
				if w := -2 / eta * dphi(c, i, j); w != 0 {
					srcs = append(srcs, src{c, i, j, w})
				}
			}
		}
	}
	for q, s := range srcs {
		l := g.Link(s.c, s.i, s.j, ks)
		var args []float64
		for r := 0; r <= nS; r++ {
			args = append(args, float64(r)*dtS*cell, s.w*samp[r])
		}
		b.CurrentWave(fmt.Sprintf("IS%d", q), m.Point(l, 0), m.Point(l, 1), circuit.Waveform{Kind: circuit.WavePWL, Args: args})
	}
	c, err := b.Build()
	if err != nil {
		t.Fatal(err)
	}
	stop := 46 * cell
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: cell / 4, Stop: stop, MaxStep: cell / 5}, analysis.TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("%d links, %d faces, %d elements (%d sources), %d time points", g.Links(), len(g.Plaquettes), len(c.Elements), len(srcs), len(wf.Time))

	// V along the probe row from the inner conductor to the wall, I round
	// the inner conductor (the dual square just outside it).
	volt := make([]float64, len(wf.Time))
	for i := hi; i < nc; i++ {
		l := g.Link(0, i, jp, kp)
		a, _ := wf.NodeVoltage(m.Point(l, 0))
		z, _ := wf.NodeVoltage(m.Point(l, 1))
		for q := range volt {
			volt[q] += a[q] - z[q]
		}
	}
	curr := make([]float64, len(wf.Time))
	addFace := func(axis, i, j int, s float64) {
		cur, ok := wf.ElementCurrent(m.Arc(g.Face(axis, i, j, kp), 0))
		if !ok {
			t.Fatalf("no face (%d, %d, %d)", axis, i, j)
		}
		for q := range curr {
			curr[q] += s * cur[q]
		}
	}
	for j := lo; j <= hi; j++ {
		addFace(1, hi, j, 1)
		addFace(1, lo-1, j, -1)
	}
	for i := lo; i <= hi; i++ {
		addFace(0, i, lo-1, 1)
		addFace(0, i, hi, -1)
	}
	// Extremum in a window: vertex of the parabola through the three
	// samples round the largest one.
	extremum := func(y []float64, t0, t1, sign float64) (tp, yp float64) {
		best := -1
		for q, tt := range wf.Time {
			if tt < t0 || tt > t1 || q == 0 || q == len(y)-1 {
				continue
			}
			if best < 0 || sign*y[q] > sign*y[best] {
				best = q
			}
		}
		ta, tb, tcc := wf.Time[best-1], wf.Time[best], wf.Time[best+1]
		ya, yb, yc := y[best-1], y[best], y[best+1]
		d1, d2 := (yb-ya)/(tb-ta), (yc-yb)/(tcc-tb)
		a2 := (d2 - d1) / (tcc - ta)
		tp = (ta+tb)/2 - d1/(2*a2)
		yp = yb + d1*(tp-tb) + a2*(tp-ta)*(tp-tb)
		return tp, yp
	}
	type event struct {
		name       string
		at         float64 // expected arrival [cells]
		vs, is     float64 // expected signs of V and I
		tV, v, cur float64
	}
	evs := []event{
		{name: "direct", at: tc + kp - ks, vs: 1, is: 1},
		{name: "from the near short", at: tc + kp + ks, vs: -1, is: -1},
		{name: "from the far short", at: tc + 2*nz - kp - ks, vs: -1, is: 1},
	}
	for q := range evs {
		e := &evs[q]
		e.tV, e.v = extremum(volt, (e.at-5)*cell, (e.at+5)*cell, e.vs)
		_, e.cur = extremum(curr, (e.at-5)*cell, (e.at+5)*cell, e.is)
		t.Logf("%-20s V %+.5f at t = %.3f cells (expected %.0f), I %+.6g A, |V/I| %.4f Ω", e.name, e.v, e.tV/cell, e.at, e.cur, math.Abs(e.v/e.cur))
		// The lattice's own speed: v_g = c·cos(kΔ/2) delays a baseband
		// Gaussian of σ cells by a fraction 1/(8σ²) (3 % at σ = 2).
		want := tc + (e.at-tc)*(1+1/(8*width*width))
		if math.Abs(e.tV/cell-want) > 0.015*(e.at-tc) {
			t.Errorf("%s: arrival %.3f cells, lattice %.3f (c: %.0f)", e.name, e.tV/cell, want, e.at)
		}
		if math.Abs(e.v-e.vs) > 5e-2 {
			t.Errorf("%s: V %.5f, want %+.0f", e.name, e.v, e.vs)
		}
		if r := math.Abs(e.v/e.cur)/z0Lattice - 1; math.Abs(r) > 2e-2 {
			t.Errorf("%s: |V/I| %.4f Ω against the lattice Z₀ %.4f Ω", e.name, math.Abs(e.v/e.cur), z0Lattice)
		}
	}
	t.Logf("lattice Z₀ %.4f Ω; direct pulse speed %.5f c", z0Lattice, float64(kp-ks)/(evs[0].tV/cell-tc))
	lg := wf.Ledger
	if lg.MaxRelativeResidual > 1e-6 {
		t.Errorf("transient ledger: max relative residual %.3g", lg.MaxRelativeResidual)
	}
	t.Logf("ledger: source energy %.6g J, stored at the end %.6g J, dissipated %.3g J, integration mismatch %.2g J, max relative residual %.2g",
		lg.EnergyIn, lg.StoredEnd, lg.EnergyDissipated, lg.IntegrationMismatch, lg.MaxRelativeResidual)

	// em.FDTD with the same grid, conductor and edge sources.
	cond := func(x, y, z float64) bool {
		const tol = 1e-9
		return x > lo*dx-tol && x < hi*dx+tol && y > lo*dx-tol && y < hi*dx+tol
	}
	fd, err := em.NewFDTD(em.FDTDConfig{N: n, D: [3]float64{dx, dx, dx}, Conductor: cond, Courant: 0.25, Ledger: true})
	if err != nil {
		t.Fatal(err)
	}
	for _, s := range srcs {
		w := s.w
		if err := fd.AddSource(em.EdgeSource{Comp: s.c, I: s.i, J: s.j, K: ks, Current: func(tt float64) float64 { return w * shape(tt) }}); err != nil {
			t.Fatal(err)
		}
	}
	vAt := func() float64 {
		v := 0.0
		for i := hi; i < nc; i++ {
			v += fd.E[0][i+n[0]*(jp+n[1]*kp)] * dx
		}
		return v
	}
	worst, q := 0.0, 0
	for fd.Time() < stop {
		tt := fd.Time()
		for q+1 < len(wf.Time) && wf.Time[q+1] < tt {
			q++
		}
		if q+1 < len(wf.Time) {
			s := (tt - wf.Time[q]) / (wf.Time[q+1] - wf.Time[q])
			vc := volt[q] + s*(volt[q+1]-volt[q])
			worst = math.Max(worst, math.Abs(vAt()-vc))
		}
		fd.Step()
	}
	// Same space lattice, two time integrators (leapfrog at Courant 0.25,
	// trapezoidal at ≤ cell/5) and the circuit's 1/Ideal departure.
	if worst > 2e-2 {
		t.Errorf("FDTD vs circuit V(t) at the probe: %.3g of the 1 V pulse", worst)
	}
	fl := fd.Ledger()
	t.Logf("FDTD (Courant 0.25) vs circuit, V(t) at the probe over %d steps: largest difference %.3g of the 1 V pulse; FDTD source work %.6g J (circuit %.6g J), its ledger residual %.2g",
		fd.Steps, worst, fl.SourceWork, lg.EnergyIn, fl.Residual()/fl.SourceWork)
}
