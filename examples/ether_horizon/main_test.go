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

// Part 9a: the round void's equilibrium carries m_e c² and is stable, with
// the bag model's virial (pressure: E = (4/3)·light; tension: E = (3/2)·light).
func TestRoundVoid(t *testing.T) {
	for _, x := range []float64{XTM11, XTE11} {
		R, p := PressureForElectron(x)
		if rel(VoidEnergy(R, x, p, 0), Me*C0*C0) > 1e-12 || rel(VoidRadius(x, p, 0), R) > 1e-9 || VoidStiffness(R, x, p, 0) <= 0 {
			t.Errorf("pressure-held void, x = %g: E %g, R %g vs %g", x, VoidEnergy(R, x, p, 0), VoidRadius(x, p, 0), R)
		}
		Rs, sg := TensionForElectron(x)
		if rel(VoidEnergy(Rs, x, 0, sg), Me*C0*C0) > 1e-12 || rel(VoidRadius(x, 0, sg), Rs) > 1e-9 || VoidStiffness(Rs, x, 0, sg) <= 0 {
			t.Errorf("tension-held void, x = %g", x)
		}
		if rel(Hbar*C0*x/R*4/3, Me*C0*C0) > 1e-12 {
			t.Errorf("virial (pressure)")
		}
	}
	// TM11 root: [x·j1(x)]' = 0.
	g := func(x float64) float64 { return x * sphj(1, x) }
	if d := (g(XTM11+1e-6) - g(XTM11-1e-6)) / 2e-6; math.Abs(d) > 1e-8 {
		t.Errorf("TM11 root: derivative %g", d)
	}
}

// Part 9b: the breathing tank in kronspice keeps W·R constant (the photon's
// ħω ∝ 1/R), halves its frequency when the void doubles, and the energy it
// loses is the work on the moving wall.
func TestBreathingVoid(t *testing.T) {
	steps := 800.0
	if testing.Short() {
		steps = 300
	}
	r, err := BreathingTank{L0: 1e-6, C0: 1e-9, I0: 1, SEnd: 2, Periods: 50, StepsPerPeriod: steps}.Run()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.Invariant-1) > 2e-4 || math.Abs(r.BooksResidual) > 2e-4 || math.Abs(r.Freq1/r.Freq0-0.5) > 1e-4 {
		t.Errorf("W·R %g, books %g, frequency ratio %g", r.Invariant, r.BooksResidual, r.Freq1/r.Freq0)
	}
}

// Part 9c: the flow over the wall vanishes at the poles and the equator,
// the void is round with no flow and stretches along the axis with it.
func TestVoidShape(t *testing.T) {
	if WallSpeed(0) != 0 || math.Abs(WallSpeed(math.Pi/2)) > 1e-15 || rel(WallSpeed(math.Pi/4), 2.5) > 1e-12 {
		t.Errorf("wall speed")
	}
	if Elongation(0) != 1 || !(Elongation(0.1) > 1 && Elongation(0.3) > Elongation(0.1)) {
		t.Errorf("elongation %g %g %g", Elongation(0), Elongation(0.1), Elongation(0.3))
	}
	// Mean of ½u² over the sphere has no P2, P4 parts beyond those used: ε6 = 0.
	if e := ShapeModes(1)[6]; math.Abs(e) > 1e-9 {
		t.Errorf("ε6 = %g", e)
	}
	if we := WeberFor(1.1); we < 0.2 || we > 0.35 {
		t.Errorf("We for 10 %% elongation %g", we)
	}
	if rel(HollowTail(2, 8, 1), 0.5) > 1e-12 {
		t.Errorf("hollow tail")
	}
}

// Part 10: only a squeeze that grows with R holds the void at one size;
// Meyl's self-similar 1/R field is neutral or runs away; −1/R² balances
// unstably.
func TestContractionLaws(t *testing.T) {
	want := map[string]string{
		"pressure, +p·R³":                 "stable",
		"tension, +σ·R²":                  "stable",
		"linear, +k·R":                    "stable",
		"logarithmic, +k·ln R":            "stable",
		"Meyl self-similar, −k/R (k = ½)": "none",
		"Meyl self-similar, −k/R (k = 1)": "neutral",
		"steeper, −k/R² (k = ½)":          "unstable",
	}
	for _, c := range Contractions() {
		roots, kinds, neutral := Balance(c)
		got := "none"
		switch {
		case neutral:
			got = "neutral"
		case len(roots) == 1 && kinds[0] > 0:
			got = "stable"
		case len(roots) == 1:
			got = "unstable"
		case len(roots) > 1:
			got = "several"
		}
		if got != want[c.Name] {
			t.Errorf("%s: %s, want %s", c.Name, got, want[c.Name])
		}
		if len(roots) == 1 && math.Abs(roots[0]-1) > 1e-6 {
			t.Errorf("%s: root %g, want 1", c.Name, roots[0])
		}
	}
}

func TestNucleonVoid(t *testing.T) {
	rp := VoidAtMass(Mp, XTM11, false)
	if math.Abs(rp*1e15-0.7694) > 2e-3 {
		t.Errorf("proton TM11 pressure void %g fm", rp*1e15)
	}
	rt := VoidAtMass(Mp, XTM11, true)
	if r := rt / RpCODATA; r < 1 || r > 1.04 {
		t.Errorf("proton TM11 tension void / measured = %g", r)
	}
	if r := VoidAtMass(Mn, XTM11, true) / rt; math.Abs(r-Mp/Mn) > 1e-12 {
		t.Errorf("neutron/proton void ratio %g", r)
	}
	if r := (MuN / MuP) / BagMomentRatio; math.Abs(r-1) > 0.03 {
		t.Errorf("moment ratio off by %g", r-1)
	}
	d, b := ProbeBound()
	if math.Abs(d) > 1e-3 || b > 4.1e-3 || b < 3.9e-3 {
		t.Errorf("probe bound: d %g, bound %g", d, b)
	}
}

func TestOnion(t *testing.T) {
	o := Onion{N: 3, Q: 1 / MuOverE, K: 0.01, L: 1e-3, F0: 1e3}
	p, err := o.Modes()
	if err != nil {
		t.Fatal(err)
	}
	if len(p) != 3 {
		t.Fatalf("%d peaks, want 3: %v", len(p), p)
	}
	for k, f := range o.Freqs() {
		if r := p[k] / f; r < 0.98 || r > 1.001 {
			t.Errorf("shell %d: coupled/uncoupled %g", k, r)
		}
	}
	r1, r2 := p[1]/p[0], p[2]/p[1]
	if math.Abs(r2/r1-1) > 0.01 {
		t.Errorf("self-similar ladder not geometric: %g, %g", r1, r2)
	}
	if math.Abs(r2/TauOverM-1) < 5 {
		t.Errorf("onion step %g should miss the τ/μ ratio %g by >5×", r2, TauOverM)
	}
}

// Part 11: g = 1 + 1/n for a squeeze R^n; H23 as written gives g = 1;
// the crossed ring's lowest note is half the plain ring's; the charge's
// self-energy on a torus overshoots α/2π, a shell at 2πR hits it.
func TestRingG(t *testing.T) {
	if _, g := WvdMMoment(); math.Abs(g-1) > 1e-9 {
		t.Errorf("H23 g = %g, want 1", g)
	}
	for n, want := range map[float64]float64{3: 4.0 / 3, 2: 1.5, 1: 2} {
		if g := SplitG(n, 0); math.Abs(g-want) > 1e-12 {
			t.Errorf("n = %g: g = %g, want %g", n, g, want)
		}
	}
	R, _, ratio := StringRing(Me)
	if math.Abs(R/(Hbar/(Me*C0))-1) > 1e-12 || math.Abs(ratio-2) > 1e-12 {
		t.Errorf("string ring R %g, ratio %g", R, ratio)
	}
	if g := SplitG(1, 2*Alpha*BetaSchwinger); math.Abs(g/2-1-AeSchwinger) > 1e-15 {
		t.Errorf("string ring with β = 1/4π: a = %g", g/2-1)
	}
	if n := NForAnomaly(); math.Abs(SplitG(n, 0)/2-1-AeMeasured) > 1e-12 {
		t.Errorf("n for anomaly: %g", n)
	}
}

func TestTorusSelfEnergy(t *testing.T) {
	for _, a := range []float64{0.01, 0.05} {
		if d := TorusInvD(a, 400); math.Abs(d/ThinTorusInvD(a)-1) > 2e-3 {
			t.Errorf("a = %g: ⟨1/d⟩ %g vs thin %g", a, d, ThinTorusInvD(a))
		}
	}
	prev := math.Inf(1)
	for _, a := range []float64{0.01, 0.1, 0.5, 1} {
		b := Beta(TorusInvD(a, 400))
		if b >= prev || Anomaly(b) < 3*AeMeasured {
			t.Errorf("a = %g: β %g (prev %g), anomaly %g", a, b, prev, Anomaly(b))
		}
		prev = b
	}
	if s := ShellForSchwinger(false); math.Abs(Anomaly(1/(2*s))-AeSchwinger) > 1e-15 {
		t.Errorf("shell %g gives %g", s, Anomaly(1/(2*s)))
	}
}

func TestCrossedHalf(t *testing.T) {
	plain := DefaultRing()
	plain.Crossed = false
	mp, err := plain.Measure(1, 0.02, 2001)
	if err != nil {
		t.Fatal(err)
	}
	mc, err := DefaultRing().Measure(0, 0.02, 2001)
	if err != nil {
		t.Fatal(err)
	}
	if r := mc.F / mp.F; math.Abs(r-0.5) > 1e-4 {
		t.Errorf("crossed/plain = %g, want 0.5", r)
	}
}

// Part 12: the vortex ring is a constant-force string (g = 2 with no
// pressure or charge, virial exact); pressure lowers g, the wall charge
// raises it; the donut never holds a thin shape; the horn-torus anomaly is
// 95.9 % of measured; kronspice's guided ring note is √(f_c² + f_z²).
func TestDonut(t *testing.T) {
	for _, pi := range DonutPressures {
		for _, s := range DonutShapes {
			d := DonutAt(s, pi, WallInvD(s))
			if v := d.EL + d.EU - d.EV - 3*d.EP; math.Abs(v) > 1e-9*d.M {
				t.Errorf("Π %g a/R %g: virial %g", pi, s, v)
			}
			if g := 2*(1+d.EU/d.EL) - 2*d.EP/d.EL; math.Abs(g-d.G) > 1e-9 {
				t.Errorf("Π %g a/R %g: g %g vs virial %g", pi, s, d.G, g)
			}
			if d.D2 <= 0 {
				t.Errorf("Π %g a/R %g: unstable in size", pi, s)
			}
		}
		if _, fat := ShapeScan(ShapeGrid(), pi, WallInvD); !fat {
			t.Errorf("Π %g: a thin donut is held", pi)
		}
	}
	for _, s := range DonutShapes {
		if g := DonutAt(s, 0, 0).G; math.Abs(g-2) > 1e-12 {
			t.Errorf("a/R %g: bare swirl g = %g", s, g)
		}
	}
	h := DonutAt(1, 0, WallInvD(1))
	if r := (h.G/2 - 1) / AeMeasured; math.Abs(r-0.9588) > 2e-3 {
		t.Errorf("horn anomaly ratio %g", r)
	}
	if r := h.RadiusOverLbar(); math.Abs(r-3.820) > 0.01 {
		t.Errorf("horn R/ƛ %g", r)
	}
}

func TestGuideRing(t *testing.T) {
	fz := GuideRing{N: 200, L: 1e-6, Cd: 1e-9, Ld: 1e30}.Analytic()
	g := GuideRing{N: 200, L: 1e-6, Cd: 1e-9}
	g.Ld = 1 / (math.Pow(2*math.Pi*3*fz, 2) * g.Cd)
	pk, err := g.Peak()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(pk/math.Hypot(g.Cutoff(), fz)-1) > 1e-4 {
		t.Errorf("guided note %g vs %g", pk, math.Hypot(g.Cutoff(), fz))
	}
}

// Part 13: the Go solver reproduces the independent Python reference
// (ref/ether13_ref.py) at the start, conserves energy without friction,
// keeps the light in the void and the thin donut's vortex, and the
// dissipative flow only lowers the energy.
func TestEtherSimReference(t *testing.T) {
	s := NewRingState(NewEtherGrid(0.3, 22, 30), 300, 20, 8, 2)
	ee, el, _ := s.Energies()
	for _, c := range []struct{ got, want float64 }{{ee, 731.406}, {el, 262.4841}, {s.Lam, 0.76553}} {
		if math.Abs(c.got/c.want-1) > 2e-5 {
			t.Errorf("start state %g, reference %g", c.got, c.want)
		}
	}
}

func TestEtherSimRealTime(t *testing.T) {
	s := NewRingState(Part13Grid(), 300, Part13V0, 8, 2)
	e0, l0, _ := s.Energies()
	for k := 0; k < 200; k++ {
		s.Step(0.01)
	}
	ee, el, _ := s.Energies()
	sh := s.Shape()
	if d := (ee + el - e0 - l0) / (e0 + l0); math.Abs(d) > 1e-3 {
		t.Errorf("energy drift %g", d)
	}
	if math.Abs(ee/626.5559-1) > 1e-4 || math.Abs(s.Lam/1.3579-1) > 1e-4 {
		t.Errorf("t = 2: E_e %g, λ %g against the reference 626.5559, 1.3579", ee, s.Lam)
	}
	if !sh.Any || sh.OnAxis || sh.Winding != 1 || sh.LightIn < 0.999 || sh.ZC > -0.2 {
		t.Errorf("thin donut at t = 2: %+v", sh)
	}
}

func TestEtherSimRelax(t *testing.T) {
	s := NewRingState(Part13Grid(), 100, Part13V0, 8, 2)
	e0, l0, _ := s.Energies()
	prev := e0 + l0
	for k := 0; k < 4; k++ {
		s.Relax(500, 20)
		ee, el, _ := s.Energies()
		if ee+el > prev*(1+1e-9) {
			t.Errorf("relax step %d raised the energy: %g > %g", k, ee+el, prev)
		}
		prev = ee + el
	}
	// light mode residual
	g := s.G
	h := make([]float64, len(s.Chi))
	g.LightOp(s.Chi, h, s.dens, s.V0, s.Lam)
	r, n := 0.0, 0.0
	for i := 0; i < g.NR; i++ {
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			r += h[k] * h[k] / g.Rho[i]
			n += g.Rho[i] * s.Chi[k] * s.Chi[k] * s.Lam * s.Lam
		}
	}
	if math.Sqrt(r/n) > 1e-3 {
		t.Errorf("light mode residual %g", math.Sqrt(r/n))
	}
}
