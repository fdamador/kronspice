// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"testing"

	"github.com/fdamador/kronspice/analysis"
)

// Speed-controlled disks: the circuit's growth rate is (−R + M√(ω1ω2))/L.
func TestSpeedControlGrowth(t *testing.T) {
	for _, w := range [][2]float64{{2, 0.8}, {0.5, 1.5}, {2, -2}} {
		exact := -1.0
		if w[0]*w[1] > 0 {
			exact = -1 + math.Sqrt(w[0]*w[1])
		}
		c, r, err := speedRun(Rikitake(), constW(w[0]), constW(w[1]), 6/math.Max(math.Abs(exact), 0.05))
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(c-exact) > 1e-3 || math.Abs(r-exact) > 1e-3 {
			t.Errorf("ω = %v: circuit %g, RK4 %g, exact %g", w, c, r, exact)
		}
	}
}

// The feedback coil (G[0][0] = −K·M) and the modulated torque in the
// circuit agree with the RK4 equations.
func TestFeedbackAndModulationTwins(t *testing.T) {
	for _, p := range []Params{
		func() Params { p := Rikitake(); p.K = 0.5; return p }(),
		func() Params { p := Rikitake(); p.K = -0.25; return p }(),
		func() Params { p := Rikitake(); p.GMod, p.Wd = 0.3, 1.5; return p }(),
	} {
		tr, mis, err := p.RunCircuit(10, 0.01)
		if err != nil {
			t.Fatal(err)
		}
		ref := p.RunRK4(10, 0.01)
		for k := 1; k < len(tr.T); k++ {
			if d := math.Abs(tr.I1[k] - ref.I1[int(math.Round(tr.T[k]/0.01))]); d > 0.02 {
				t.Fatalf("K %g GMod %g: t = %g, |ΔI1| = %g", p.K, p.GMod, tr.T[k], d)
			}
		}
		if math.Abs(mis) > 1e-3 {
			t.Errorf("K %g: ledger mismatch %g J", p.K, mis)
		}
	}
}

// The exponential KS test accepts exponential samples and rejects
// Rikitake's intervals (which have a floor at the oscillation period).
func TestKSExponential(t *testing.T) {
	rng := rand.New(rand.NewPCG(1, 11)) // seed 7 draws an unlucky sample (D = 0.10); typical D ≈ 0.04
	x := make([]float64, 300)
	for i := range x {
		x[i] = 3 * rng.ExpFloat64()
	}
	if _, p := ksExp(x, 400, 3); p < 0.05 {
		t.Errorf("exponential sample rejected, p = %g", p)
	}
	r := Rikitake().StatRun(5000, 0.002, 100, 1e-2)
	if _, p := ksExp(r.Intervals, 400, 3); p > 0.01 {
		t.Errorf("Rikitake intervals accepted as exponential, p = %g", p)
	}
	if _, pp := ks2(x, x); pp < 0.99 {
		t.Errorf("two-sample KS of a sample with itself: p = %g", pp)
	}
}

// Two-height ionosphere: the fit puts all five modes within 3 % of the
// measured frequencies, with Q between 3 and 7.
func TestTwoHeightFit(t *testing.T) {
	th, _ := fitTwoHeight(5)
	for n := 1; n <= 5; n++ {
		f, q := th.Mode(n)
		if math.Abs(f/MeasuredSchumann[n-1]-1) > 0.03 || q < 3 || q > 7 {
			t.Errorf("n = %d: f %g Hz (measured %g), Q %g", n, f, MeasuredSchumann[n-1], q)
		}
	}
}

// Refined layer 2: storm books close, columns share one ionosphere, and
// the two-scale profile keeps the anchors.
func TestAtmosphereRefinements(t *testing.T) {
	a := StandardAtmosphere()
	s0, H := a.Calibrate()
	single := []Column{{"g", 1, func(z float64) float64 { return s0 * math.Exp(z/H) }}}
	st := Storms{Is: 2000, RUp: 100, RDown: 100, RCloud: 500, Distributed: true}
	g, err := a.SolveGlobal(single, st)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(g.PGen-g.PStorm-g.PFair) > 1e-5*g.PGen {
		t.Errorf("storm books %g ≠ %g + %g", g.PGen, g.PStorm, g.PFair)
	}
	sl, so := 1e-14, 2.2e-14
	cols := []Column{
		{"land", 0.29, func(z float64) float64 { return sl * math.Exp(z/H) }},
		{"ocean", 0.71, func(z float64) float64 { return so * math.Exp(z/H) }},
	}
	lo, err := a.SolveGlobal(cols, Storms{Is: a.GlobalCurrent()})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(lo.E0[0]/lo.E0[1]-1) > 1e-3 || math.Abs(lo.IFair[0]+lo.IFair[1]-a.GlobalCurrent()) > 1e-3*a.GlobalCurrent() {
		t.Errorf("land/ocean: E0 %v, I %v", lo.E0, lo.IFair)
	}
	two, _ := a.TwoScale(1500, 5000)
	ts, err := a.SolveGlobal([]Column{{"g", 1, two}}, Storms{Is: a.GlobalCurrent()})
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ts.VI/a.VI-1) > 2e-3 || math.Abs(ts.E0[0]/a.E0-1) > 5e-3 {
		t.Errorf("two-scale: V_I %g, E0 %g", ts.VI, ts.E0[0])
	}
}

// The unified network: one equipotential ionosphere at DC, Schumann-like
// peaks at ELF, from the same elements.
func TestUnifiedNetwork(t *testing.T) {
	prof, _ := StandardAtmosphere().TwoScale(1500, 5000)
	u := Unified{Bands: 31, Layers: 30, Top: 100e3, Sigma: prof, Hm: 90e3, RLoss: 0.15, FRef: 8}
	c, err := u.Circuit(1020, 1)
	if err != nil {
		t.Fatal(err)
	}
	if he := u.ElectricHeight(8); he < 35e3 || he > 60e3 {
		t.Errorf("electric height %g km", he/1e3)
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	v0, _ := op.Voltage("t0")
	v1, _ := op.Voltage(fmt.Sprintf("t%d", u.Bands-1))
	if math.Abs(v0/v1-1) > 1e-2 || math.Abs(v0/250e3-1) > 1e-2 {
		t.Errorf("DC ionosphere %g V at the pole, %g V at the other pole", v0, v1)
	}
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 196, FStart: 3, FStop: 42}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	vp, _ := ac.NodeVoltage("t0")
	var m []float64
	for k := range ac.Freq {
		m = append(m, cmplx.Abs(vp[k]))
	}
	pk := Peaks(ac.Freq, m)
	if len(pk) < 5 {
		t.Fatalf("peaks %v", pk)
	}
	for n := 0; n < 5; n++ {
		if math.Abs(pk[n]/MeasuredSchumann[n]-1) > 0.15 {
			t.Errorf("peak %d at %g Hz, measured %g", n+1, pk[n], MeasuredSchumann[n])
		}
	}
}
