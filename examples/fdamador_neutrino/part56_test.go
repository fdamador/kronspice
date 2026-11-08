// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// Part 5a: the time-reversed ring-down is absorbed with no reflection
// beyond its switch-on (e⁻⁸ of the power); the literal reversal of the
// lossy ring-down adds the predicted (R_d/(R_rad+R_d))²; drives that are
// not time-reversed are reflected; every run's books close.
func TestImplosionCPA(t *testing.T) {
	e := DefaultERing()
	set, err := e.Implosions()
	if err != nil {
		t.Fatal(err)
	}
	cpa, exact := set.Runs[0], set.Runs[1]
	if r := 1 - cpa.Absorbed; r > 2*math.Exp(-8) {
		t.Errorf("CPA reflected %.2e of the incident energy, want ≤ 2e⁻⁸", r)
	}
	if d := (cpa.Absorbed - exact.Absorbed) - set.RefCPA*set.RefCPA; math.Abs(d) > 1e-5 {
		t.Errorf("exact-reversal extra reflection %.3e, predicted %.3e", cpa.Absorbed-exact.Absorbed, set.RefCPA*set.RefCPA)
	}
	if p := cpa.PeakStored / cpa.EInc; p < 0.99 {
		t.Errorf("CPA peak stored %.4f of the incident energy", p)
	}
	for _, r := range set.Runs[2:] {
		if r.Absorbed > 0.95 {
			t.Errorf("%s absorbed %.4f: a drive that is not the reversed ring-down should be partly reflected", r.Name, r.Absorbed)
		}
	}
	if a := set.Runs[4].Absorbed; a > 0.1 {
		t.Errorf("the decaying (not reversed) drive absorbed %.3f, want < 0.1", a)
	}
	for _, r := range set.Runs {
		if r.LedgerRel > 1e-4 || r.WaveRel > 1e-4 {
			t.Errorf("%s: ledger %.1e, wave books %.1e", r.Name, r.LedgerRel, r.WaveRel)
		}
	}
}

// Part 5a pulsing ring: in a steady cycle the outside's net supply (in −
// out) is the dielectric heat plus the change of stored energy.
func TestTrainBooks(t *testing.T) {
	if testing.Short() {
		t.Skip("six-cycle transient")
	}
	tr, err := DefaultERing().Train("locked", 40, 6, 1)
	if err != nil {
		t.Fatal(err)
	}
	c := tr.Cycles[len(tr.Cycles)-2]
	if d := math.Abs(c.EInc-c.EOut-c.Heat-(c.StoredEnd-c.StoredStart)) / c.EInc; d > 1e-4 {
		t.Errorf("cycle books: in %.4g, out %.4g, heat %.3g, Δstored %.3g (rel %.1e)", c.EInc, c.EOut, c.Heat, c.StoredEnd-c.StoredStart, d)
	}
	if tr.LedgerRel > 1e-4 {
		t.Errorf("train ledger %.1e", tr.LedgerRel)
	}
}

// Part 5b on a coarse grid (r/2, shortened windows): the torus rings near
// the lumped f0 and Q; the time-reversed outgoing mode is taken in (little
// outside the surface at the peak); the rebuilt field is E round the ring
// and B poloidal; the field ledger closes; a constant burst takes in less.
func TestFieldImplosionCoarse(t *testing.T) {
	if testing.Short() {
		t.Skip("field-solver runs (about 2 minutes)")
	}
	fr := DefaultFieldRing()
	fr.Cells, fr.FilterT, fr.RingDownTE, fr.ImplodeTE, fr.ExplodeTE = 2, 10, 3, 4, 1
	r, err := fr.Run()
	if err != nil {
		t.Fatal(err)
	}
	e := DefaultERing()
	if d := rel(r.F0, e.Omega0()/(2*math.Pi)); d > 0.05 {
		t.Errorf("f0 %.4g Hz, lumped %.4g (rel %.2f)", r.F0, e.Omega0()/(2*math.Pi), d)
	}
	if d := rel(r.Q, e.QRad()); d > 0.1 {
		t.Errorf("Q %.4g, lumped %.4g", r.Q, e.QRad())
	}
	// The shortened growth window (4 τ_E) starts at e⁻⁴ = 0.018 of the
	// final power, and that abrupt start is not absorbed: allow 1.5× it.
	if a, lim := r.Peak/r.Delivered, 1-1.5*math.Exp(-fr.ImplodeTE); a < lim {
		t.Errorf("held at the peak %.4f of the delivered energy, want ≥ %.4f", a, lim)
	}
	if r.EPhiFrac < 0.99 || r.BPhiFrac > 0.01 {
		t.Errorf("pattern: E_φ share %.4f, B_φ share %.4f", r.EPhiFrac, r.BPhiFrac)
	}
	if math.Abs(r.Residual) > 1e-10 || math.Abs(r.RDResidual) > 1e-10 {
		t.Errorf("field ledger residuals %.1e, %.1e", r.Residual, r.RDResidual)
	}
	if b := r.BurstPeak / r.BurstEInc; b > 0.9 {
		t.Errorf("burst peak %.3f of its incident energy, want clearly below the reversed drive", b)
	}
}

// Part 6: a state that cancels everywhere does not exist — the two waves
// of a harmonic are orthogonal on the ring.
func TestNoFullCancellation(t *testing.T) {
	r := DefaultFRing(false)
	for _, m := range []int{1, 2, 5} {
		k := r.Kappa(m)
		var g complex128
		for j := 0; j < r.N; j++ {
			g += cmplx.Exp(complex(0, -2*k*float64(j)))
		}
		if cmplx.Abs(g) > 1e-12 {
			t.Errorf("harmonic %d: ⟨CW|CCW⟩ = %.2e, want 0", m, cmplx.Abs(g))
		}
	}
}

// Test 1: every wave type is stationary on the perfect ring (plain and
// crossed); C turns once a period.
func TestTypesPersist(t *testing.T) {
	ps, err := Persistence(100)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range ps {
		if p.MinOwn < 1-1e-9 || p.MaxOther > 1e-9 || p.MinPurity < 1-1e-9 {
			t.Errorf("%s %s: own min %.12f, other max %.2e, purity %.12f", p.Ring, p.Type, p.MinOwn, p.MaxOther, p.MinPurity)
		}
		if p.Type == "C" && math.Abs(p.Rotation-1) > 2e-3 {
			t.Errorf("%s C: %.5f turns a period", p.Ring, p.Rotation)
		}
		if p.LedgerRel > 1e-3 {
			t.Errorf("%s %s: ledger %.1e", p.Ring, p.Type, p.LedgerRel)
		}
	}
}

// Test 2: a defect makes the two-state beat; Δf from kronspice equals the
// exact eigen-splitting, CW converts fully (θ = 45°), and a standing start
// s sections from the defect gives sin²2θ with θ = κs.
func TestFlavourBeat(t *testing.T) {
	b, err := FlavourBeat(false, 0.1, "C", 0, 1.25)
	if err != nil {
		t.Fatal(err)
	}
	if d := rel(b.DfFit, b.DfExact); d > 0.01 {
		t.Errorf("Δf kronspice %.5g, exact %.5g (rel %.1e)", b.DfFit, b.DfExact, d)
	}
	if d := rel(b.DfExact, b.DfFirst); d > 0.02 {
		t.Errorf("Δf exact %.5g vs first order %.5g", b.DfExact, b.DfFirst)
	}
	if math.Abs(b.AmpFit-1) > 0.01 {
		t.Errorf("CW → CCW amplitude %.4f, want 1", b.AmpFit)
	}
	if b.LedgerRel > 1e-3 {
		t.Errorf("ledger %.1e", b.LedgerRel)
	}
	s, err := FlavourBeat(false, 0.2, "S", 1, 1.25)
	if err != nil {
		t.Fatal(err)
	}
	if want := sq(math.Sin(2 * s.ThetaPred * math.Pi / 180)); math.Abs(s.AmpFit-want) > 0.01 {
		t.Errorf("standing start s = 1: sin²2θ %.4f, want %.4f", s.AmpFit, want)
	}
}

// Test 2, three states: the fit reaches the measured angles and splitting
// ratio (each within 1.5σ, ratio within 5 %), and a kronspice run follows
// the ring's full mode expansion; with no CP phase P(e→μ) = P(μ→e).
func TestThreeFlavour(t *testing.T) {
	if testing.Short() {
		t.Skip("PMNS search and two long transients")
	}
	cl := FitPMNS()
	m := NuFIT6()
	if math.Abs(cl.Th12-m.Th12) > 1.5*m.S12 || math.Abs(cl.Th23-m.Th23) > 1.5*m.S23 || math.Abs(cl.Th13-m.Th13) > 1.5*m.S13 ||
		math.Abs(cl.Ratio/m.Ratio-1) > 0.05 {
		t.Errorf("fit θ12 %.2f, θ23 %.2f, θ13 %.2f, ratio %.2f", cl.Th12, cl.Th23, cl.Th13, cl.Ratio)
	}
	a, err := ThreeFlavour(cl, cl.Assign[0], 1)
	if err != nil {
		t.Fatal(err)
	}
	b, err := ThreeFlavour(cl, cl.Assign[1], 1)
	if err != nil {
		t.Fatal(err)
	}
	if a.MaxDevAll > 1e-3 {
		t.Errorf("kronspice vs full mode expansion: %.2e", a.MaxDevAll)
	}
	var d float64
	for i := range a.TAll {
		d = math.Max(d, math.Abs(a.PAll[cl.Assign[1]][i]-interp(b.TAll, b.PAll[cl.Assign[0]], a.TAll[i])))
	}
	if d > 0.03 {
		t.Errorf("P(e→μ) − P(μ→e) up to %.3f", d)
	}
}

// Test 3: through one port S leaves, D stays (dark), C is left as D with
// half its energy; the reversed explosion rebuilds S (one port) and C (two
// ports in quadrature), but one port gives only half of C.
func TestImplodeTypes(t *testing.T) {
	r := DefaultFRing(false)
	k := r.Kappa(1)
	R := r.PortR(k, 40.5)
	conj := func(f Fractions) float64 { return f.D }
	ex := map[string]Explosion{}
	for _, ty := range []string{"S", "D", "C"} {
		x, _, err := r.Explode(ty, []Comp{r.Type(ty, k, 0)}, []int{0}, R, 80, 0, conj)
		if err != nil {
			t.Fatal(err)
		}
		ex[ty] = x
	}
	if ex["S"].Remaining > 1e-3 || ex["D"].Remaining < 0.99 || math.Abs(ex["C"].Remaining-0.5) > 0.01 || ex["C"].End.D < 0.999 {
		t.Errorf("explosions: S left %.4f, D left %.4f, C left %.4f (D share %.4f)", ex["S"].Remaining, ex["D"].Remaining, ex["C"].Remaining, ex["C"].End.D)
	}
	if testing.Short() {
		return
	}
	s, err := r.RebuildByReversal("S", []Comp{r.Type("S", k, 0)}, []int{0}, R, 80, 10, 0, conj)
	if err != nil {
		t.Fatal(err)
	}
	c2, err := r.RebuildByReversal("C2", []Comp{r.Type("C", k, 0)}, []int{0, 6}, R, 80, 10, 0, conj)
	if err != nil {
		t.Fatal(err)
	}
	c1, err := r.RebuildByReversal("C1", []Comp{r.Type("C", k, 0)}, []int{0}, R, 80, 10, 0, conj)
	if err != nil {
		t.Fatal(err)
	}
	if s.Fidelity < 0.9999 || s.Absorbed < 0.995 || c2.Fidelity < 0.9999 || math.Abs(c1.Fidelity-0.5) > 0.01 {
		t.Errorf("rebuilds: S fidelity %.5f absorbed %.5f; C two ports %.5f; C one port %.5f", s.Fidelity, s.Absorbed, c2.Fidelity, c1.Fidelity)
	}
	for _, x := range []Rebuild{s, c2, c1} {
		if x.LedgerRel > 1e-3 {
			t.Errorf("%s ledger %.1e", x.Case, x.LedgerRel)
		}
	}
}
