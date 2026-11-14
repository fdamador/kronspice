// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

func rel(a, b float64) float64 { return math.Abs(a-b) / math.Abs(b) }

// Part 1: the n = 1 vortex against the literature — axis slope
// 0.5827811878 (Rorai, Sreenivasan & Fisher, arXiv:1212.6389) and line energy ln(1.46 R/ξ) (Fetter
// 2009) — and an empty axis for every winding.
func TestVortex(t *testing.T) {
	v, err := SolveVortex(1, 40, 8000)
	if err != nil {
		t.Fatal(err)
	}
	if c := v.AxisCoefficient(); rel(c, 0.5827811878) > 1e-3 {
		t.Errorf("axis slope %.6f, literature 0.5827811878", c)
	}
	if k := math.Exp(v.LineEnergy() - math.Log(40)); math.Abs(k-1.46) > 0.01 {
		t.Errorf("line energy ln(%.3f R/ξ), Fetter ln(1.46 R/ξ)", k)
	}
	if f := v.F[len(v.F)-1]; math.Abs(f-1) > 1e-3 {
		t.Errorf("f(R) = %g, want ≈ 1", f)
	}
	for n := 1; n <= 3; n++ {
		v, err := SolveVortex(n, 40, 4000)
		if err != nil {
			t.Fatal(err)
		}
		if v.F[0] != 0 || v.F[1] <= 0 {
			t.Errorf("n = %d: f(0) = %g, f(h) = %g", n, v.F[0], v.F[1])
		}
		// The swirl outruns the local wave speed only inside the core region.
		if e, c := v.ErgoRadius(), v.CoreRadius(); !(e > c && e < 2*c) {
			t.Errorf("n = %d: ergoregion %.3f, core %.3f", n, e, c)
		}
	}
}

// Part 2: draining bathtub closed forms, and the Hawking lifetime scaling
// as r_h² with the electron-sized horizon gone in ~10⁻²⁰ s.
func TestHorizon(t *testing.T) {
	b := Bathtub{A: 3, B: 4, C: 1}
	if b.Horizon() != 3 || b.Ergosurface() != 5 {
		t.Errorf("bathtub r_h %g r_e %g, want 3, 5", b.Horizon(), b.Ergosurface())
	}
	if r := HawkingLifetime(2e-13) / HawkingLifetime(1e-13); rel(r, 4) > 1e-12 {
		t.Errorf("lifetime ratio %g, want 4", r)
	}
	if tau := HawkingLifetime(RingRadius()); tau < 1e-21 || tau > 1e-19 {
		t.Errorf("electron-sized horizon lifetime %g s", tau)
	}
	if r := HawkingRadiusFor(AgeU); rel(HawkingLifetime(r), AgeU) > 1e-12 {
		t.Errorf("HawkingRadiusFor does not invert HawkingLifetime")
	}
}

// Part 3: GR numbers of an electron-mass hole, and no window for an ether
// horizon between absorption (too big) and Hawking emission (too small).
func TestCoupling(t *testing.T) {
	h := ElectronHole()
	if rel(h.Rs, 1.3526e-57) > 1e-3 {
		t.Errorf("r_s = %g", h.Rs)
	}
	if h.SpinA/h.RG < 1e44 || h.ChargeR/h.RG < 1e21 {
		t.Errorf("not naked? a/M %g, r_Q/M %g", h.SpinA/h.RG, h.ChargeR/h.RG)
	}
	if h.HawkingLifetime > 1e-100 {
		t.Errorf("Hawking lifetime %g s", h.HawkingLifetime)
	}
	// Page's σ over the membrane estimate is (k r_h)²/2 (the GR barrier).
	rh, k := 1e-20, 2*math.Pi/Compton()
	if r := PageRate(rh) / MembraneRate(rh); rel(r, k*k*rh*rh/2) > 1e-12 {
		t.Errorf("Page/membrane %g, want (k r_h)²/2 = %g", r, k*k*rh*rh/2)
	}
	if rel(1/MembraneRate(MembraneRadiusFor(AgeU)), AgeU) > 1e-12 || rel(1/PageRate(PageRadiusFor(AgeU)), AgeU) > 1e-12 {
		t.Errorf("radius bounds do not invert the rates")
	}
	if !(MembraneRadiusFor(AgeU) < PageRadiusFor(AgeU) && PageRadiusFor(AgeU) < HawkingRadiusFor(AgeU)/1e20) {
		t.Errorf("window: absorption ≤ %g / %g m, Hawking ≥ %g m", MembraneRadiusFor(AgeU), PageRadiusFor(AgeU), HawkingRadiusFor(AgeU))
	}
	if z := Z0(); math.Abs(z-376.730313) > 1e-5 {
		t.Errorf("Z0 = %g", z)
	}
}

// Part 4: kronspice AC on the crossed ring with horizon resistors twins the
// ladder's complex dispersion (frequency and Q), the horizon's power is
// booked and the ledger closes, and the spin-½ boundary is kept.
func TestRingWithHorizon(t *testing.T) {
	qs := []float64{0, 30}
	if testing.Short() {
		qs = []float64{30}
	}
	for _, qh := range qs {
		r := DefaultRing()
		if qh > 0 {
			r.RH = r.HorizonRFor(0, qh)
		}
		f, q := r.Mode(0)
		m, err := r.Measure(0, 0.25, 3001)
		if err != nil {
			t.Fatal(err)
		}
		if rel(m.F, f) > 1e-4 || rel(m.Q, q) > 5e-3 {
			t.Errorf("Q_h %g: AC %.2f Hz, Q %.3f; ladder %.2f Hz, Q %.3f", qh, m.F, m.Q, f, q)
		}
		if res := math.Abs(m.PIn-m.PWire-m.PHor-m.PBleed) / m.PIn; res > 1e-6 || m.MaxLedgerResidual > 1e-9 {
			t.Errorf("Q_h %g: power split residual %g, kronspice ledger %g", qh, res, m.MaxLedgerResidual)
		}
		if qh > 0 {
			// The horizon's share of the loss is Q/Q_h of the input.
			if share := m.PHor / m.PIn; math.Abs(share-q/qh) > 0.02 {
				t.Errorf("horizon share %.4f, Q/Q_h %.4f", share, q/qh)
			}
		}
		if r.SignAfterTrip(0) != -1 {
			t.Errorf("crossed ring lost its sign flip")
		}
	}
	// Combined Q: 1/Q = 1/Q_wire + 1/Q_h to second order.
	r := DefaultRing()
	_, q0 := r.Mode(0)
	r.RH = r.HorizonRFor(0, 100)
	if _, q := r.Mode(0); rel(1/q, 1/q0+1/100.0) > 1e-3 {
		t.Errorf("1/Q = %g, 1/Q0 + 1/Q_h = %g", 1/q, 1/q0+0.01)
	}
}

// Part 5: the sink's mass ledger closes, weak sinks stay subsonic with the
// flux linear in γ (2γ·πr_s²·ρ∞ for n = 0), strong ones make a horizon
// outside the sink, and the flux falls again past its maximum (Zeno).
func TestSink(t *testing.T) {
	rs := 1 / math.Sqrt2 // 1 ξ
	weak, err := SolveSink(0, 0.01, rs, 60, 0.002)
	if err != nil {
		t.Fatal(err)
	}
	if want := 2 * math.Pi * 0.01 * rs * rs; rel(-weak.FluxFar, want) > 2e-3 {
		t.Errorf("weak sink flux %g, want 2γ·πr_s² = %g", -weak.FluxFar, want)
	}
	if len(weak.Horizons()) != 0 {
		t.Errorf("weak sink has horizons %v", weak.Horizons())
	}
	gc, sc, err := CriticalGamma(0, rs, 60, 0.002)
	if err != nil {
		t.Fatal(err)
	}
	if rel(sc.LossRate, -sc.FluxFar) > 1e-5 {
		t.Errorf("mass ledger: loss %g, inflow %g", sc.LossRate, -sc.FluxFar)
	}
	strong, err := SolveSink(0, 10, rs, 60, 0.002)
	if err != nil {
		t.Fatal(err)
	}
	if h := strong.OuterHorizon(); !(h > rs && h < 2*rs) {
		t.Errorf("strong sink outer horizon %g, sink edge %g", h, rs)
	}
	if !(gc > 0.5 && gc < 3) || -strong.FluxFar >= -sc.FluxFar {
		t.Errorf("γ_c %g; flux at γ_c %g, at γ = 10 %g (expected Zeno decline)", gc, -sc.FluxFar, -strong.FluxFar)
	}
}

// Parts 6–7: heat bound, gravity fit, growth bound, Coulomb ratio.
func TestFarField(t *testing.T) {
	if b := HeatBoundPerElectron(); b < 1e-38 || b > 1e-37 {
		t.Errorf("heat bound %g W", b)
	}
	if rel(HawkingPowerPaidByFlow(HeatHorizonBound()), HeatBoundPerElectron()) > 1e-12 {
		t.Errorf("HeatHorizonBound does not invert the Hawking power")
	}
	e := SinkElectron{RhoE: RhoCrit}
	if rel(e.RhoE*e.Q()*e.Q()/(4*math.Pi), G*Me*Me) > 1e-12 {
		t.Errorf("gravity fit")
	}
	if rel(SinkElectron{RhoE: RhoForGrowth(GMSunRate)}.GrowthRate(), GMSunRate) > 1e-12 {
		t.Errorf("growth bound does not invert")
	}
	if rel(CoulombOverGravity(), 4.17e42) > 2e-3 {
		t.Errorf("Coulomb/gravity %g", CoulombOverGravity())
	}
}

// The kronspice lattice: KCL is continuity (inflow through the faces =
// Q), source power = resistor power, and the far field tends to
// Q/(4πd²).
func TestSinkGrid(t *testing.T) {
	g := &SinkGrid{N: 13, Sinks: [][3]int{{6, 6, 6}}, Q: []float64{1}}
	if err := g.Solve(); err != nil {
		t.Fatal(err)
	}
	if math.Abs(g.Bound-1) > 1e-12 || rel(g.PIn, g.PRes) > 1e-12 || g.Ledger > 1e-9 {
		t.Errorf("faces %g, P_in %g, P_R %g, ledger %g", g.Bound, g.PIn, g.PRes, g.Ledger)
	}
	r3 := -g.Speed(9, 6, 6) * 4 * math.Pi * 9
	r2 := -g.Speed(8, 6, 6) * 4 * math.Pi * 4
	if !(r3 < r2 && r3 > 0.8 && r3 < 1.4) {
		t.Errorf("speed·4πd² at d = 2, 3: %g, %g", r2, r3)
	}
}

// Part 8a: the fountain is an exact steady flow (Bragg–Hawthorne residual
// at finite-difference level), destroys nothing (zero net flux), has no
// 2:1 surface without the axis vortex, and the fitted q = 2 surface has
// tube/ring ½.
func TestFountain(t *testing.T) {
	for _, l := range []int{1, 2} {
		f := NewFountain(l, 1.5*NewFountain(l, 0).K*NewFountain(l, 0).PsiO())
		if r := f.Residual(0.6, 1.0); r > 1e-6 {
			t.Errorf("l = %d: Bragg–Hawthorne residual %g", l, r)
		}
		if q := f.SphereFlux(0.5, 4000); math.Abs(q) > 1e-5 {
			t.Errorf("l = %d: net flux %g", l, q)
		}
		bare := NewFountain(l, 0)
		for _, s := range bare.QProfile([]float64{0.05, 0.5, 0.95}, FountainDT) {
			if s.Q >= 1 {
				t.Errorf("l = %d, no axis vortex: q = %g ≥ 1", l, s.Q)
			}
		}
	}
	// l = 2: out at both poles, in at the equator.
	f := NewFountain(2, 0)
	if ur, _, _ := f.Velocity(0.5, 0.01); ur <= 0 {
		t.Errorf("north axis u_r = %g, want outflow", ur)
	}
	if ur, _, _ := f.Velocity(0.5, math.Pi-0.01); ur <= 0 {
		t.Errorf("south axis u_r = %g, want outflow", ur)
	}
	if ur, _, _ := f.Velocity(0.5, math.Pi/2); ur >= 0 {
		t.Errorf("equator u_r = %g, want inflow", ur)
	}
	if testing.Short() {
		return
	}
	_, s, ok := FitTorus(1, 2, 0.5, FountainDT)
	if !ok || math.Abs(s.Q-2) > 1e-6 || math.Abs(s.Minor/s.Major-0.5) > 1e-6 || s.Closure > 1e-3 {
		t.Errorf("fit: ok %v q %g tube/ring %g closure %g", ok, s.Q, s.Minor/s.Major, s.Closure)
	}
}

// Part 8b: the membrane circuit in kronspice — maximum power at the match
// (Ω_F = Ω_H/2), the ledger closes, the core is +R_H with its sources off
// and −R_L at its operating point (Γ_core·Γ_load = 1), and the spin cannot
// keep the leaky ring going for even a period.
func TestSpinningCore(t *testing.T) {
	pts, err := SpinSweep([]float64{0.5, 1, 2})
	if err != nil {
		t.Fatal(err)
	}
	if !(pts[1].PLoad > pts[0].PLoad && pts[1].PLoad > pts[2].PLoad) || math.Abs(pts[1].OmegaF-0.5) > 1e-12 {
		t.Errorf("match: P_load %g %g %g, Ω_F/Ω_H %g", pts[0].PLoad, pts[1].PLoad, pts[2].PLoad, pts[1].OmegaF)
	}
	for _, p := range pts {
		if math.Abs(p.PSpin-p.PLoad-p.PHorizon) > 1e-12 || p.Residual > 1e-9 {
			t.Errorf("R_L/R_H %g: spin %g ≠ load %g + horizon %g", p.RLoverRH, p.PSpin, p.PLoad, p.PHorizon)
		}
		if p.RLoverRH != 1 && cmplx.Abs(p.GammaCore*p.GammaLoad-1) > 1e-12 {
			t.Errorf("Γ_core·Γ_load = %v", p.GammaCore*p.GammaLoad)
		}
	}
	cs, err := CoreSmith()
	if err != nil {
		t.Fatal(err)
	}
	if cmplx.Abs(cs.Points[0].Z-1) > 1e-6 || cs.Points[0].Region.String() != "passive" {
		t.Errorf("core with sources off: Z = %v, %v", cs.Points[0].Z, cs.Points[0].Region)
	}
	if tau := SpinLifetime(SuperradianceOmega()); tau*OmegaC()/(2*math.Pi) > 0.1 {
		t.Errorf("spin keeps the ring for %g periods", tau*OmegaC()/(2*math.Pi))
	}
	if rel(RimSpeed(RingRadius()), C0) > 1e-9 {
		t.Errorf("rim speed at r_c = R: %g c", RimSpeed(RingRadius())/C0)
	}
}
