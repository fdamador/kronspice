// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Total ν–e elastic cross-sections far above m_e against the tabulated
// values σ/E [cm²/MeV] of J. A. Formaggio, G. P. Zeller, Rev. Mod. Phys.
// 84, 1307 (2012), Table: ν_e 9.49e-45, ν_μ 1.58e-45, ν̄_e 3.98e-45,
// ν̄_μ 1.31e-45 (theirs include radiative corrections, ours are tree level
// with sin²θ_W = 0.23122: held to 2.5 %). The numerical T integral equals
// the closed form (G_F²m_eE/2π)[(g_V ± g_A)² + (g_V ∓ g_A)²/3] up to the
// m_e/E term.
func TestWeakKnownValues(t *testing.T) {
	const e = 100.0
	for _, c := range []struct {
		f    Flavour
		anti bool
		want float64
	}{{NuE, false, 9.49e-45}, {NuMu, false, 1.58e-45}, {NuE, true, 3.98e-45}, {NuMu, true, 1.31e-45}} {
		got := SigmaWeak(c.f, c.anti, e, 1e-9) / e
		if math.Abs(got/c.want-1) > 0.025 {
			t.Errorf("flavour %d anti %v: σ/E %.4g, tabulated %.4g", c.f, c.anti, got, c.want)
		}
		gv, ga := WeakCouplings(c.f, 0)
		if c.anti {
			ga = -ga
		}
		pre := gFMeV * gFMeV * meMeV / (2 * math.Pi) * hcMeVcm * hcMeVcm
		closed := pre * ((gv+ga)*(gv+ga) + (gv-ga)*(gv-ga)/3)
		if math.Abs(got/closed-1) > 2*meMeV/e {
			t.Errorf("flavour %d anti %v: numeric %.6g, closed form %.6g", c.f, c.anti, got, closed)
		}
	}
}

// The dipole term is 2.495e-45 cm²·(μ/1e-10 μB)²·(1/T − 1/E) (Vogel &
// Engel 1989, T and E in MeV); the ring form factor is 1 − 0.2x² at small
// x and falls as x⁻³ far out.
func TestDipoleFormula(t *testing.T) {
	for _, x := range [][2]float64{{1, 0.1}, {10, 3}} {
		got := DipoleDSigmaDT(1e-20, x[0], x[1], 0)
		want := 2.495e-45 * (1/x[1] - 1/x[0])
		if math.Abs(got/want-1) > 1e-3 {
			t.Errorf("E %g T %g: %.5g, Vogel–Engel %.5g", x[0], x[1], got, want)
		}
	}
	if f := FormFactor2(0.01); math.Abs(f-(1-0.2e-4)) > 1e-8 {
		t.Errorf("FormFactor2(0.01) = %.12f, want 1 − 2e-5", f)
	}
	r := FormFactor2(200) / FormFactor2(100)
	if math.Abs(r-0.125) > 0.01 {
		t.Errorf("FormFactor2 tail ratio x 100 → 200: %.4f, want 1/8", r)
	}
}

// A toroidal solenoid (n turns of radius b round a torus of major radius
// a) has toroidal dipole ½·n·I·π·b²·a along the axis, exactly for any b/a,
// and no charge or magnetic dipole.
func TestToroidalSolenoid(t *testing.T) {
	a, b, i := 0.1, 0.03, 2.0
	n := 64
	m := Moments(ToroidalSolenoid(a, b, i, n, 512))
	want := 0.5 * float64(n) * i * math.Pi * b * b * a
	if math.Abs(math.Abs(m.T[2])/want-1) > 1e-9 || math.Hypot(m.T[0], m.T[1]) > 1e-9*want {
		t.Errorf("toroidal dipole %v, want ±%.9g ẑ", m.T, want)
	}
	if norm3(m.M) > 1e-12*want/a || m.Q != 0 {
		t.Errorf("magnetic dipole %v, charge %g: want 0", m.M, m.Q)
	}
}

// The ring's states: m = 0 has |μ| = πa²I0 and nothing else; m = ±1 has
// |p| = πaI0/ω, toroidal ω|T|/(c²|p|) = (ka)²/5, no charge, no ⟨r²⟩, no
// magnetic dipole; every time average vanishes. The m = 0 moment with one
// quantum equals Part 2's.
func TestRingMoments(t *testing.T) {
	e := DefaultERing()
	w := e.Omega0()
	s0 := ComputeState(e.A, 0, 1e-3, w)
	if math.Abs(s0.MAmp/(math.Pi*e.A*e.A*1e-3)-1) > 1e-12 || s0.PAmp != 0 || s0.TAmp > 1e-12*1e-3*e.A*e.A*e.A {
		t.Errorf("m = 0: μ %.9g (want %.9g), p %g, T %g", s0.MAmp, math.Pi*e.A*e.A*1e-3, s0.PAmp, s0.TAmp)
	}
	if s0.MAvg > 1e-12*s0.MAmp {
		t.Errorf("m = 0: time-averaged μ %.3g of the amplitude", s0.MAvg/s0.MAmp)
	}
	for _, n := range []int{1, -1} {
		s := ComputeState(e.A, n, 1e-3, 2*w)
		if math.Abs(s.PAmp/(math.Pi*e.A*1e-3/(2*w))-1) > 1e-12 {
			t.Errorf("m = %d: |p| %.9g, want πaI0/ω", n, s.PAmp)
		}
		ka := 2 * w * e.A / C0
		if math.Abs(s.ToroidalOverP/(ka*ka/5)-1) > 1e-9 {
			t.Errorf("m = %d: ω|T|/(c²|p|) %.9g, want (ka)²/5 = %.9g", n, s.ToroidalOverP, ka*ka/5)
		}
		if s.QMax > 1e-14*s.PAmp/e.A || s.R2Max > 1e-14*s.PAmp*e.A || s.MAmp > 1e-12*s.PAmp*C0 || s.PAvg > 1e-12*s.PAmp {
			t.Errorf("m = %d: q %g, Σλr² %g, μ %g, ⟨p⟩ %g", n, s.QMax, s.R2Max, s.MAmp, s.PAvg)
		}
	}
	b := DefaultBench()
	r := b.Run()
	a := r.RadiusNeeded
	got := math.Pi * a * a * I0Quantum0(a, a/b.AoverR, r.Omega0) / MuB
	if math.Abs(got/r.MomentAmplitude-1) > 1e-12 {
		t.Errorf("Part 2 moment %.9g μB, bench %.9g", got, r.MomentAmplitude)
	}
}

// The charge-radius shift (2/3)M_W²⟨r²⟩sin²θ_W equals √2πα⟨r²⟩/(3G_F) at
// tree level (M_W² sin²θ_W = πα/(√2G_F)); with measured M_W and the MS-bar
// sin²θ_W the two differ by the radiative corrections, about 7 %.
func TestChargeRadiusShift(t *testing.T) {
	r2 := 1e-32
	tree := math.Sqrt2 * math.Pi * alphaFS * r2 / (3 * gFMeV * hcMeVcm * hcMeVcm)
	if g := DeltaGV(r2); math.Abs(g/tree-1) > 0.08 {
		t.Errorf("δg_V %.5g, tree level %.5g", g, tree)
	}
	gv0, _ := WeakCouplings(NuMu, 0)
	gv1, _ := WeakCouplings(NuMu, r2)
	if math.Abs(gv1-gv0-DeltaGV(r2)) > 1e-15 {
		t.Errorf("g_V not shifted by δg_V")
	}
}

// Bethe–Peierls with the effective-range factor peaks at ω = 2B with
// (8π/3)α(ħ²/M)/(8B)/(1 − κρ), about 2.4 mb (measured near 2.4–2.5 mb).
func TestSigmaE1Peak(t *testing.T) {
	kappa := math.Sqrt(mNucMeV*bDeut) / hcMeVfm
	want := 8 * math.Pi / 3 * alphaFS * hcMeVfm * hcMeVfm / mNucMeV / (8 * bDeut) / (1 - kappa*rhoTrip) * 1e-26
	pk := SigmaE1(2 * bDeut)
	if math.Abs(pk/want-1) > 1e-12 || SigmaE1(2*bDeut*1.01) > pk || SigmaE1(2*bDeut*0.99) > pk {
		t.Errorf("peak %.6g cm² at 2B, closed form %.6g", pk, want)
	}
	if pk < 2.2e-27 || pk > 2.7e-27 {
		t.Errorf("peak %.3g mb, outside 2.2–2.7 mb", pk*1e27)
	}
}

// The deuteron break-up scales as μ² and tends to the point-dipole value
// as the ring shrinks (the form-factor table path against the point path).
func TestDeutEMLimit(t *testing.T) {
	ff := newFFTable()
	pt := SigmaDeutEM(1, 10, 0, ff)
	small := SigmaDeutEM(1, 10, 1e-18, ff)
	if math.Abs(small/pt-1) > 1e-3 {
		t.Errorf("a → 0: %.6g, point %.6g", small, pt)
	}
	if math.Abs(SigmaDeutEM(4, 10, 0, ff)/pt-4) > 1e-12 {
		t.Errorf("not ∝ μ²")
	}
	if SigmaDeutEM(1, 10, 1.46e-12, ff) > 1e-3*pt {
		t.Errorf("the ring form factor should suppress the break-up")
	}
}

// σ_weak ∝ E far above m_e; the ring's dipole cross-section saturates.
func TestEnergyDependence(t *testing.T) {
	sl := math.Log(SigmaWeak(NuE, false, 30, 1e-3)/SigmaWeak(NuE, false, 10, 1e-3)) / math.Log(3)
	if math.Abs(sl-1) > 0.02 {
		t.Errorf("weak slope %.4f, want 1", sl)
	}
	a := 1.46e-12
	sr := math.Log(SigmaDipole(0.003, 30, 1e-3, a)/SigmaDipole(0.003, 10, 1e-3, a)) / math.Log(3)
	if sr > 0.05 {
		t.Errorf("ring slope %.4f, want ≈ 0", sr)
	}
}

// kronspice twin: a moment held fixed does not change the ring, the work
// done on it is M·I_e·ΔI_ring, its time average vanishes; ledger.
func TestProbeTwin(t *testing.T) {
	tw, err := RunProbeTwin(DefaultERing(), 20)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(tw.QWith/tw.QWithout-1) > 1e-9 || tw.MaxDiff > 1e-9 {
		t.Errorf("ring changed by the probe: Q %.9g vs %.9g, max diff %.2g", tw.QWith, tw.QWithout, tw.MaxDiff)
	}
	if math.Abs(tw.WorkSource-tw.WorkPredicted) > 0.01*math.Abs(tw.WorkPredicted)+1e-6*tw.PeakInteraction {
		t.Errorf("work %.6g J, M·I_e·ΔI %.6g J", tw.WorkSource, tw.WorkPredicted)
	}
	if tw.MeanPowerOverPeak > 1e-3 || tw.LedgerRel > 1e-9 {
		t.Errorf("mean power / peak %.2g, ledger %.2g", tw.MeanPowerOverPeak, tw.LedgerRel)
	}
}
