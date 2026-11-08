// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// Part 8a: the triplet m = +1, 0, −1 spans three dimensions (Gram matrix
// the identity), Part 6's types C, S, D only two (one zero eigenvalue);
// each triplet state persists on the perfect ring with J·ω/U = +1, 0, −1.
func TestTripletIndependent(t *testing.T) {
	ti, err := TripletIndependence(20)
	if err != nil {
		t.Fatal(err)
	}
	for a := 0; a < 3; a++ {
		if math.Abs(ti.GramTriplet[a]-1) > 1e-12 {
			t.Errorf("triplet Gram eigenvalue %d = %.15f, want 1", a, ti.GramTriplet[a])
		}
	}
	if math.Abs(ti.GramTypes[0]) > 1e-12 || math.Abs(ti.GramTypes[1]-1) > 1e-12 || math.Abs(ti.GramTypes[2]-2) > 1e-12 {
		t.Errorf("types Gram eigenvalues %v, want 0, 1, 2", ti.GramTypes)
	}
	want := [3]float64{1, 0, -1}
	for a, p := range ti.Persist {
		if p.MinOwn < 1-1e-9 || math.Abs(p.Jw0-want[a]) > 1e-9 || math.Abs(p.JwT-want[a]) > 1e-6 || p.LedgerRel > 1e-3 {
			t.Errorf("%s: own share %.12f, J·ω/U %.9f → %.9f, ledger %.1e", p.Name, p.MinOwn, p.Jw0, p.JwT, p.LedgerRel)
		}
	}
}

// Series-capacitor defects never touch the uniform state (Cmat·1 = Cg·1).
func TestSeriesDefectsLeaveUniform(t *testing.T) {
	r := DefaultFRing(false)
	r.Def = map[int]float64{2: 0.2064, 23: 0.1628, 7: -0.3}
	_, v := r.Modes()
	var best float64
	for i := range v[0] {
		var s float64
		for j := range v {
			s += v[j][i] / math.Sqrt(float64(r.N))
		}
		best = math.Max(best, s*s)
	}
	// 1e-9: the Jacobi solver resolves the eigenvectors of a matrix whose
	// entries span C = 1 nF to Cg = 1 pF to about 1e-10
	if 1-best > 1e-9 {
		t.Errorf("no mode is the uniform pattern: best overlap %.15f", best)
	}
}

// Reciprocal ring with e = m = 0: μ–τ reflection symmetry, so θ23 = 45°
// and δ_CP = 90° or 270° exactly; rotation breaks it.
func TestMuTauSymmetry(t *testing.T) {
	k := TripletKnobs{M1: 0.03, M2: 0.136, Psi: 0.09, Cg: 1e-9}
	c := TripletMix(k, [3]int{1, 0, 2})
	d := math.Min(math.Abs(c.DeltaCP-90), math.Abs(c.DeltaCP-270))
	if c.MuTauAsym() > 1e-12 || math.Abs(c.Th23-45) > 1e-9 || d > 1e-6 {
		t.Errorf("reciprocal: asymmetry %.1e, θ23 %.12f, δ_CP %.9f", c.MuTauAsym(), c.Th23, c.DeltaCP)
	}
	k.Nu = 2e-3
	cr := TripletMix(k, [3]int{1, 0, 2})
	if cr.MuTauAsym() < 1e-3 || math.Abs(cr.Th23-45) < 0.01 {
		t.Errorf("rotating: asymmetry %.1e, θ23 %.4f, want both off the symmetric values", cr.MuTauAsym(), cr.Th23)
	}
}

// 8c: the angular momentum the field loses under a defect equals the
// time integral of the torque −∂H/∂α computed independently; +1 → −1 is
// a change of 2 units; without the defect J is kept.
func TestTorqueLedger(t *testing.T) {
	tq, err := Torque(0.1, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(tq.DeltaJ+2) > 0.02 {
		t.Errorf("ΔJ·ω/U = %.5f, want −2 (less the loss)", tq.DeltaJ)
	}
	if tq.MaxDiff > 0.02 || math.Abs(tq.DeltaJ-tq.TorqueAt) > 0.01 {
		t.Errorf("ΔJ %.5f vs ∫τ dt %.5f, largest difference %.1e", tq.DeltaJ, tq.TorqueAt, tq.MaxDiff)
	}
	if tq.LedgerRel > 1e-3 {
		t.Errorf("ledger %.1e", tq.LedgerRel)
	}
	fr, err := FreeRingJ(50)
	if err != nil {
		t.Fatal(err)
	}
	if fr > 1e-6 {
		t.Errorf("free ring: J changes by %.1e", fr)
	}
}

// 8b in kronspice: a reciprocal triplet cluster started in m = 0 follows
// the exact mode expansion (same readout); ledger closes.
func TestTripletRunTwin(t *testing.T) {
	if testing.Short() {
		t.Skip("kronspice transient (about 10 s)")
	}
	c := TripletMix(TripletKnobs{M1: 0.4856, M2: 0.6507, Psi: -0.0132, Cg: 5.73e-9}, [3]int{1, 0, 2})
	r, err := RunTriplet(c, 0.5)
	if err != nil {
		t.Fatal(err)
	}
	if r.MaxDevAll > 2e-3 || r.LedgerRel > 1e-3 {
		t.Errorf("kronspice vs exact expansion %.2e, ledger %.1e", r.MaxDevAll, r.LedgerRel)
	}
	if r.MinM0 > 0.2 {
		t.Errorf("m = 0 share stays above %.3f: no mixing", r.MinM0)
	}
}

// 8e: Coriolis gives m·ħΩ (none for m = 0); on the Part 2 ring Ω_cor is
// its own ω₀; the centrifugal term stays below ½ħω·⟨ρ²⟩/a² at v < c;
// RatioOmega returns the Ω that gives the asked splitting ratio.
func TestOutsideRotation(t *testing.T) {
	b := DefaultBench().Run()
	burst := 1881.0
	rm := OutsideRotation("part 2", 1, b.Omega0, b.RadiusNeeded, 1, 0.1, burst)
	if math.Abs(rm.CorOverOmega-1) > 1e-9 {
		t.Errorf("Ω_cor/ω₀ = %.12f, want 1", rm.CorOverOmega)
	}
	if z := OutsideRotation("part 2", 0, b.Omega0, b.RadiusNeeded, 1, 0.1, burst); !math.IsInf(z.OmCor, 1) || z.CorAtBurst != 0 {
		t.Errorf("m = 0 got a Coriolis shift: %v, %v", z.OmCor, z.CorAtBurst)
	}
	if rm.RimCen < C0 {
		t.Errorf("centrifugal 0.1 eV reached below c (%.3g m/s) with ħω = 0.1 eV", rm.RimCen)
	}
	w0, w1 := 2*math.Pi*2.0571e9, 2*math.Pi*4.4354e9
	om := RatioOmega(w0, w1, 33.55)
	lp, lm := sq(w1+om), sq(w1-om)
	r := math.Min(math.Abs(sq(w0)-lp), math.Abs(sq(w0)-lm)) / (lp - lm)
	if math.Abs(r/33.55-1) > 1e-6 {
		t.Errorf("ratio at Ω %.4g: %.6f", om, r)
	}
}

// 8d on grid r/2 with short windows: the triplet is orthonormal; the
// turning states carry equal and opposite J_z, constant magnetic energy,
// no magnetic dipole and a turning electric dipole; the E-ring's axial
// moment matches the lumped loop and passes through zero.
func TestTripletFieldCoarse(t *testing.T) {
	if testing.Short() {
		t.Skip("field solver (about 40 s)")
	}
	sf := DefaultSpinField()
	sf.ProbeT, sf.FilterT, sf.ProbeLong, sf.FilterLong = 6, 9, 12, 20
	res, err := RunTripletField(sf, 24, nil)
	if err != nil {
		t.Fatal(err)
	}
	for a := 0; a < 3; a++ {
		if math.Abs(res.GramEig[a]-1) > 1e-6 {
			t.Errorf("Gram eigenvalue %d = %.9f", a, res.GramEig[a])
		}
	}
	p, z, m := res.States[0], res.States[1], res.States[2]
	if p.M != 1 || m.M != -1 || math.Abs(p.JzMMean+m.JzMMean) > 1e-9*math.Abs(p.JzMMean) {
		t.Errorf("J_z: +1 %.6g, −1 %.6g", p.JzMMean, m.JzMMean)
	}
	for _, s := range []TState{p, m} {
		ax, pp := s.MomentPerQuantum()
		if s.UBMax-s.UBMin > 0.01 || math.Hypot(ax, pp) > 1e-9*s.DipolePerQuantum()*C0 || s.PRange[0] < 0.95 {
			t.Errorf("%s: U_B [%.4f, %.4f], |m| %.3g, |p|c %.3g, |p(t)| range %v", s.Name, s.UBMin, s.UBMax, math.Hypot(ax, pp), s.DipolePerQuantum()*C0, s.PRange)
		}
	}
	e := DefaultERing()
	lump := math.Sqrt(2*Hbar*e.Omega0()/e.L()) * math.Pi * e.A * e.A
	ax, _ := z.MomentPerQuantum()
	if r := ax / lump; r < 0.8 || r > 1.05 || z.MRange[0] > 0.05 || z.UBMin > 0.05 {
		t.Errorf("E-ring: field/lumped moment %.3f, |m(t)| min %.3f, U_B min %.3f", r, z.MRange[0], z.UBMin)
	}
	if cmplx.Abs(z.MHat[0])+cmplx.Abs(z.MHat[1]) > 1e-9*cmplx.Abs(z.MHat[2]) {
		t.Errorf("E-ring moment not axial: %v", z.MHat)
	}
}
