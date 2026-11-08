// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// Part 7a, test 1: the gyrators split harmonic 1 by exactly 4·Gy·sin κ/c_κ
// (the exact quadratic, the complex mode solver and the first-order value
// agree); kronspice measures the same splitting at its trapezoidal
// frequencies; standing turns into the other standing wave and back
// (amplitude 1); the gyrators do no net work; ledger closes.
func TestSagnacSplit(t *testing.T) {
	for _, crossed := range []bool{false, true} {
		s, err := Sagnac(crossed, 0.01, 120)
		if err != nil {
			t.Fatal(err)
		}
		if d := rel(s.DwModes, s.DwExact); d > 1e-9 {
			t.Errorf("%s: modes %.9g vs quadratic %.9g", s.Ring, s.DwModes, s.DwExact)
		}
		if d := rel(s.DwFirst, s.DwExact); d > 1e-12 {
			t.Errorf("%s: 4·Gy·sin κ/c_κ %.12g vs quadratic %.12g", s.Ring, s.DwFirst, s.DwExact)
		}
		if d := rel(s.DwKron, s.DwWarped); d > 1e-3 {
			t.Errorf("%s: kronspice %.6g vs exact at trapezoidal frequencies %.6g (rel %.1e)", s.Ring, s.DwKron, s.DwWarped, d)
		}
		if math.Abs(s.ConvAmp-1) > 0.01 {
			t.Errorf("%s: S → D amplitude %.4f, want 1", s.Ring, s.ConvAmp)
		}
		if s.GyroPowerMax > 1e-12 || s.GyroWork > 1e-12 || s.GyroPass < 1e-4 {
			t.Errorf("%s: gyrators absorb %.1e, pass %.1e of ω·W, net work %.1e", s.Ring, s.GyroPowerMax, s.GyroPass, s.GyroWork)
		}
		if s.LedgerRel > 1e-3 {
			t.Errorf("%s: ledger %.1e", s.Ring, s.LedgerRel)
		}
	}
}

// The rim-speed conversion aΩ = ν·c²/(2v_g) is spacetime.SagnacSplitting
// times v_g/c; speeds at or above c are flagged impossible.
func TestRimTwin(t *testing.T) {
	for _, vg := range []float64{C0, C0 / 10, 1e5} {
		r := Rim(1e-5, vg, 0.01, 2e7)
		if r.BelowC && r.TwinRel > 1e-12 {
			t.Errorf("v_g %.3g: twin differs by %.1e", vg, r.TwinRel)
		}
	}
	if r := Rim(1e-3, 1e3, 0.01, 2e7); r.BelowC || r.Possible {
		t.Errorf("rim %.3g m/s should be at or above c", r.RimSpeed)
	}
}

// Without rotation U is real: J = 0 and δ_CP is 0 or 180°. With rotation
// J ≠ 0, and reversing the rotation reverses J (δ → 360° − δ).
func TestJarlskog(t *testing.T) {
	c0 := RotClusterMix(0.2064, 0.1628, 2, 23, 0, [3]int{0, 1, 2})
	if math.Abs(c0.J) > 1e-12 || math.Min(math.Abs(c0.DeltaCP-180), math.Min(c0.DeltaCP, 360-c0.DeltaCP)) > 1e-6 {
		t.Errorf("no rotation: J %.2e, δ_CP %.6f", c0.J, c0.DeltaCP)
	}
	// Part 6's angles at ν = 0
	if math.Abs(c0.Th12-33.66) > 0.01 || math.Abs(c0.Th23-42.28) > 0.01 || math.Abs(c0.Th13-8.56) > 0.01 {
		t.Errorf("ν = 0 angles %.3f %.3f %.3f, Part 6 has 33.66 42.28 8.56", c0.Th12, c0.Th23, c0.Th13)
	}
	cp := RotClusterMix(0.2064, 0.1628, 2, 23, 1e-3, [3]int{0, 1, 2})
	cm := RotClusterMix(0.2064, 0.1628, 2, 23, -1e-3, [3]int{0, 1, 2})
	if math.Abs(cm.J) < 0.01 {
		t.Errorf("ν = −1e-3: J %.4f, want clearly non-zero", cm.J)
	}
	if math.Abs(cp.J+cm.J) > 1e-9 || math.Abs(cp.DeltaCP+cm.DeltaCP-360) > 1e-6 {
		t.Errorf("reversal: J %.6f / %.6f, δ %.4f / %.4f", cp.J, cm.J, cp.DeltaCP, cm.DeltaCP)
	}
	if cm.Modes.ResidualA > 1e-10 {
		t.Errorf("mode residual %.1e", cm.Modes.ResidualA)
	}
}

// Test 2: with a defect and an equal rotation splitting, CW → CCW converts
// with sin²2θ = ½ at Δf = √(Δd² + Δr²) (first order).
func TestTwoStateRot(t *testing.T) {
	ts, err := TwoStateRot(0.1, 1, 1.25)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ts.AmpPred-0.5) > 1e-9 || math.Abs(ts.AmpFit-0.5) > 0.01 {
		t.Errorf("sin²2θ predicted %.4f, kronspice %.4f", ts.AmpPred, ts.AmpFit)
	}
	if d := rel(ts.DfFit, ts.DfExact); d > 0.01 {
		t.Errorf("Δf kronspice %.6g, exact %.6g", ts.DfFit, ts.DfExact)
	}
	if d := rel(ts.DfExact, ts.DfQuad); d > 0.02 {
		t.Errorf("Δf exact %.6g, √(Δd² + Δr²) %.6g", ts.DfExact, ts.DfQuad)
	}
	if ts.LedgerRel > 1e-3 {
		t.Errorf("ledger %.1e", ts.LedgerRel)
	}
}

// Test 3 in kronspice: at the pinned CP point the runs follow the exact
// complex mode expansion, and P(e→μ) − P(μ→e) is much larger with the
// rotation than without (where only the readout's asymmetry is left).
func TestTViolation(t *testing.T) {
	if testing.Short() {
		t.Skip("four long transients")
	}
	cl := RotClusterMix(0.2040, 0.1651, 2, 23, -5.567e-4, [3]int{0, 1, 2})
	z := RotClusterMix(0.2040, 0.1651, 2, 23, 0, [3]int{0, 1, 2})
	e, mu := cl.Assign[0], cl.Assign[1]
	a, err := RotThreeFlavour(cl, e, 1)
	if err != nil {
		t.Fatal(err)
	}
	b, err := RotThreeFlavour(cl, mu, 1)
	if err != nil {
		t.Fatal(err)
	}
	a0, err := RotThreeFlavour(z, e, 1)
	if err != nil {
		t.Fatal(err)
	}
	b0, err := RotThreeFlavour(z, mu, 1)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range []RotThree{a, b, a0, b0} {
		if r.MaxDevAll > 2e-3 {
			t.Errorf("kronspice vs exact expansion %.2e", r.MaxDevAll)
		}
		if r.LedgerRel > 1e-3 || r.GyroRel > 1e-9 {
			t.Errorf("ledger %.1e, gyrator work %.1e", r.LedgerRel, r.GyroRel)
		}
	}
	m, p := TAsym(a, b, e, mu)
	m0, p0 := TAsym(a0, b0, e, mu)
	if math.Abs(m-p) > 0.01 || p < 3*p0 {
		t.Errorf("T asymmetry: rotating %.4f (predicted %.4f), still %.4f (%.4f)", m, p, m0, p0)
	}
}

// Helicity normalization: h·ω/u = ±1 for circular light at every instant,
// 0 for linear.
func TestHelicityNorm(t *testing.T) {
	for _, tt := range []float64{0, 0.3e-9, 1.7e-9} {
		if h := HelicityPlaneWave(-1, false, tt); math.Abs(h-1) > 1e-12 {
			t.Errorf("circular (+): %.15f", h)
		}
		if h := HelicityPlaneWave(1, false, tt); math.Abs(h+1) > 1e-12 {
			t.Errorf("circular (−): %.15f", h)
		}
		if h := HelicityPlaneWave(1, true, tt); math.Abs(h) > 1e-12 {
			t.Errorf("linear: %.15f", h)
		}
	}
}

// The grid integrals against a closed form: E = E0·(x, y, 0)/L, B = B0·ẑ
// gives (E×B)_φ = −E0B0ρ/L, so J_z = −ε0E0B0/L·∫(x² + y²)dV =
// −ε0E0B0/L·16b⁵/3 over the cube |x|, |y|, |z| < b; ∫E·B = 0.
func TestFieldIntegrals(t *testing.T) {
	fr := DefaultFieldRing()
	g := spinGrid{fr: fr, n: 41, dx: 1e-3}
	g.half = float64(g.n-1) / 2 * g.dx
	g.b = 15 * g.dx
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{g.n, g.n, g.n}, D: [3]float64{g.dx, g.dx, g.dx}, Origin: [3]float64{-g.half, -g.half, -g.half}})
	if err != nil {
		t.Fatal(err)
	}
	const e0, b0, l = 2.0, 3e-3, 0.01
	f.SetE(func(x, y, z float64) [3]float64 { return [3]float64{e0 * x / l, e0 * y / l, 0} })
	f.SetB(func(x, y, z float64) [3]float64 { return [3]float64{0, 0, b0} })
	var zero [3][]float64
	for c := 0; c < 3; c++ {
		zero[c] = make([]float64, len(f.E[c]))
	}
	jz, jzm, _, eb, _ := g.integrals(f.E, f.B, f.E, zero, zero)
	want := -Eps0 * e0 * b0 / l * 16 * math.Pow(g.b, 5) / 3
	// the midpoint rule on x² is low by dx²/12 per unit length: relative (dx/b)²/4
	if d := rel(jz, want); d > 0.3*sq(g.dx/g.b) {
		t.Errorf("J_z %.6e, closed form %.6e (rel %.1e)", jz, want, d)
	}
	if d := rel(jzm*Eps0, jz); d > 1e-12 { // D = E passed: Minkowski = Abraham/ε0
		t.Errorf("Minkowski path %.6e vs Abraham %.6e", jzm*Eps0, jz)
	}
	if math.Abs(eb) > 1e-15 {
		t.Errorf("∫E·B %.2e, want 0", eb)
	}
}

// Part 7b on the r/2 grid with short runs: the m = 1 mode turned 90° is
// orthogonal to itself (pure odd order); the travelling m = 1 carries
// J_z·ω/U = ±1 with the sign of its phase; the axisymmetric E + B carries
// no net J_z; its helicity changes sign with the relative phase; ledgers.
func TestSpinFieldCoarse(t *testing.T) {
	if testing.Short() {
		t.Skip("field-solver runs (about 2 minutes)")
	}
	sf := DefaultSpinField()
	sf.ProbeT, sf.FilterT, sf.ProbeLong, sf.FilterLong, sf.RunT = 8, 10, 30, 50, 4
	g := sf.grid()
	mE, err := sf.FindMode("E", "E", 1e9, 4e9)
	if err != nil {
		t.Fatal(err)
	}
	mB, err := sf.FindMode("B", "B", 1e9, 6e9)
	if err != nil {
		t.Fatal(err)
	}
	m1, err := sf.FindMode("m1", "m1", 1e9, 6e9)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(mE.Turn90-1) > 0.01 || math.Abs(m1.Turn90) > 0.01 || math.Abs(m1.Turn180+1) > 0.01 {
		t.Errorf("turn overlaps: E-ring 90° %.4f; m = 1 90° %.4f, 180° %.4f", mE.Turn90, m1.Turn90, m1.Turn180)
	}
	m1s := g.Rot90(m1)
	s2 := complex(math.Sqrt2/2, 0)
	i := complex(0, 1)
	tp, err := sf.RingDown("trav+", []Mode{m1, m1s}, []complex128{s2, i * s2})
	if err != nil {
		t.Fatal(err)
	}
	tm, err := sf.RingDown("trav−", []Mode{m1, m1s}, []complex128{s2, -i * s2})
	if err != nil {
		t.Fatal(err)
	}
	// The two senses are time reverses of each other at t = 0 (B flips):
	// J_z exactly opposite; net J_z clearly non-zero in the Minkowski form,
	// with the Abraham form of the same sign. (Its size is not a clean 1
	// on the staircased grid: fourfold symmetry mixes m with m ± 4.)
	if tp.JzMW0 < 0.5 || math.Abs(tp.JzMW0+tm.JzMW0) > 1e-9*tp.JzMW0 || tp.JzW0*tp.JzMW0 <= 0 || math.Abs(tp.JzW0+tm.JzW0) > 1e-9*math.Abs(tp.JzW0) {
		t.Errorf("travelling m = 1: Minkowski %+.4f / %+.4f, Abraham %+.4f / %+.4f", tp.JzMW0, tm.JzMW0, tp.JzW0, tm.JzW0)
	}
	ep, err := sf.RingDown("E+B +90", []Mode{mE, mB}, []complex128{s2, i * s2})
	if err != nil {
		t.Fatal(err)
	}
	em9, err := sf.RingDown("E+B −90", []Mode{mE, mB}, []complex128{s2, -i * s2})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ep.JzW0) > 0.02 || math.Abs(em9.JzW0) > 0.02 || math.Abs(ep.JzMW0) > 0.02 {
		t.Errorf("E + B (m = 0): J_z·ω/U %+.4f, %+.4f, want ≈ 0", ep.JzW0, em9.JzW0)
	}
	if ep.HW0*em9.HW0 >= 0 || math.Abs(ep.HW0+em9.HW0) > 0.05*math.Abs(ep.HW0) {
		t.Errorf("E + B helicity %+.4f and %+.4f: want opposite signs", ep.HW0, em9.HW0)
	}
	for _, r := range []SpinRun{tp, tm, ep, em9} {
		if math.Abs(r.Residual) > 1e-10 {
			t.Errorf("%s: field ledger %.1e", r.Name, r.Residual)
		}
	}
}
