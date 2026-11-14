// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// Magnetic sources are refused unless magnetic charge is enabled.
func TestMagneticNeedsHypotheticalMode(t *testing.T) {
	r, _ := New(Config{N: [3]int{8, 8, 8}, D: [3]float64{0.01, 0.01, 0.01}})
	one := func(float64) float64 { return 1 }
	if err := r.AddMagSource(MagFaceSource{Comp: 2, I: 3, J: 3, K: 3, Current: one}); !errors.Is(err, ErrMagneticDisabled) {
		t.Errorf("face source accepted without magnetic mode: %v", err)
	}
	if err := r.AddMagCellSource(MagCellSource{I: 3, J: 3, K: 3, Rate: one}); !errors.Is(err, ErrMagneticDisabled) {
		t.Errorf("cell source accepted without magnetic mode: %v", err)
	}
	if err := r.SetSm(func(x, y, z float64) float64 { return 0 }); !errors.Is(err, ErrMagneticDisabled) {
		t.Errorf("S_m set without MagAB: %v", err)
	}
}

// Gate 1, magnetic: with a conserved magnetic current and S_m = 0 at the
// start, MagAB is the dual-symmetric Maxwell of em/duality.go: S_m stays 0
// and E, B equal em.FDTD driven by the same current as a FaceSource
// (M = μ₀J_m), edge by edge; magnetic Gauss div B = μ₀ρ_m holds.
func TestMagneticReducesToDual(t *testing.T) {
	N, D := [3]int{14, 12, 10}, [3]float64{0.01, 0.012, 0.009}
	ref, err := em.NewFDTD(em.FDTDConfig{N: N, D: D, Courant: 0.99})
	if err != nil {
		t.Fatal(err)
	}
	dual, _ := New(Config{N: N, D: D, Courant: 0.99, Magnetic: MagDual})
	ab, _ := New(Config{N: N, D: D, Courant: 0.99, Magnetic: MagAB})
	g := gauss(40*ab.Dt, 10*ab.Dt)
	area := D[0] * D[1]
	cur := func(t float64) float64 { return 2e3 * g(t) * math.Sin(3e9*t) } // [A·m/s]
	if err := ref.AddMagneticSource(em.FaceSource{Comp: 2, I: 7, J: 6, K: 4,
		M: func(t float64) float64 { return em.Mu0 * cur(t) / area }}); err != nil {
		t.Fatal(err)
	}
	for _, r := range []*Region{dual, ab} {
		if err := r.AddMagSource(MagFaceSource{Comp: 2, I: 7, J: 6, K: 4, Current: cur}); err != nil {
			t.Fatal(err)
		}
	}
	var worstGauss float64
	for s := 0; s < 160; s++ {
		ref.Step()
		dual.Step()
		ab.Step()
		worstGauss = math.Max(worstGauss, dual.MagGaussDefect())
	}
	var dE, dB, dAB, mE, mB, mS float64
	for c := 0; c < 3; c++ {
		for n := range ref.E[c] {
			dE = math.Max(dE, math.Abs(ref.E[c][n]-dual.E[c][n]))
			dB = math.Max(dB, math.Abs(ref.B[c][n]-dual.B[c][n]))
			dAB = math.Max(dAB, math.Abs(ab.B[c][n]-dual.B[c][n]))
			mE = math.Max(mE, math.Abs(ref.E[c][n]))
			mB = math.Max(mB, math.Abs(ref.B[c][n]))
		}
	}
	for _, s := range ab.Sm {
		mS = math.Max(mS, math.Abs(s))
	}
	t.Logf("dual vs em.FDTD: E %.2g of %.3g V/m, B %.2g of %.3g T; AB vs dual B %.2g; max |S_m| %.2g T; magnetic Gauss defect %.2g",
		dE, mE, dB, mB, dAB, mS, worstGauss)
	if dE > 1e-12*mE || dB > 1e-12*mB {
		t.Errorf("MagDual differs from em.FDTD with a FaceSource: E %.3g, B %.3g", dE/mE, dB/mB)
	}
	if dAB > 1e-12*mB || mS > 1e-12*mB {
		t.Errorf("S_m grew without created magnetic charge: |S_m| %.3g, B diff %.3g", mS, dAB/mB)
	}
	if worstGauss > 1e-9 {
		t.Errorf("div B = μ₀ρ_m broken by a conserved magnetic current: %.3g", worstGauss)
	}
}

// Gate 3, both sectors: the energy theorem with S and S_m, electric and
// magnetic currents, created electric and magnetic charge, closes step by
// step, with every hypothetical column non-zero.
func TestEnergyTheoremBothScalars(t *testing.T) {
	r, err := New(Config{N: [3]int{17, 15, 13}, D: [3]float64{0.01, 0.01, 0.01}, Magnetic: MagAB})
	if err != nil {
		t.Fatal(err)
	}
	g := gauss(30*r.Dt, 8*r.Dt)
	_ = r.AddSource(EdgeSource{Comp: 0, I: 6, J: 7, K: 6, Current: func(t float64) float64 { return 3 * g(t) * math.Cos(2e9*t) }})
	_ = r.AddNodeSource(NodeSource{I: 9, J: 7, K: 6, Rate: func(t float64) float64 { return 1e-9 * g(t) / r.Dt }})
	if err := r.AddMagSource(MagFaceSource{Comp: 1, I: 5, J: 6, K: 7, Current: func(t float64) float64 { return 1e6 * g(t) * math.Sin(2.5e9*t) }}); err != nil {
		t.Fatal(err)
	}
	if err := r.AddMagCellSource(MagCellSource{I: 10, J: 5, K: 5, Rate: func(t float64) float64 { return 0.3 * g(t) / r.Dt }}); err != nil {
		t.Fatal(err)
	}
	r.SetS(func(x, y, z float64) float64 {
		return 1e-9 * math.Sin(math.Pi*x/0.16) * math.Sin(math.Pi*y/0.14) * math.Sin(math.Pi*z/0.12)
	})
	_ = r.SetSm(func(x, y, z float64) float64 {
		return 2e-9 * math.Sin(math.Pi*x/0.16) * math.Sin(2*math.Pi*y/0.14) * math.Sin(math.Pi*z/0.12)
	})
	worst := 0.0
	for s := 0; s < 400; s++ {
		r.Step()
		l := r.Ledger()
		scale := math.Max(math.Abs(l.W), math.Max(math.Abs(l.SourceWork)+math.Abs(l.MagSourceWork),
			math.Abs(l.ScalarWork)+math.Abs(l.MagScalarWork)))
		worst = math.Max(worst, math.Abs(l.Residual())/scale)
	}
	l := r.Ledger()
	t.Logf("W %.4g J (S %.3g, S_m %.3g); work: electric %.4g, magnetic %.4g, c²ρS %.4g, cρ_mS_m %.4g J; worst residual %.2g",
		l.W, l.WS, l.WSm, l.SourceWork, l.MagSourceWork, l.ScalarWork, l.MagScalarWork, worst)
	if worst > 1e-11 {
		t.Errorf("energy theorem residual %.3g", worst)
	}
	if l.MagScalarWork == 0 || l.MagSourceWork == 0 || l.ScalarWork == 0 {
		t.Error("a hypothetical column is zero with its source present")
	}
}

// Gate 2, magnetic, and □S_m = (μ₀/c)(∂ρ_m/∂t + div J_m): magnetic charge
// appearing in a cell at rate q̇_m(t) sends out the retarded
// S_m(r, t) = (μ₀/c) q̇_m(t − r/c)/(4πr), carried by B alone (E = 0, the
// dual of the electric S wave). The magnetic Gauss law in AB form holds;
// with S_m held at 0 (MagDual) the same source breaks it.
func TestRetardedMagneticScalarWave(t *testing.T) {
	const h = 1e-3
	r, err := New(Config{N: [3]int{42, 42, 42}, D: [3]float64{h, h, h}, Magnetic: MagAB})
	if err != nil {
		t.Fatal(err)
	}
	tau := 10 * r.Dt
	t0 := 5 * tau
	q0 := 300.0 // peak rate [A·m/s]
	rate := func(t float64) float64 { return q0 * gauss(t0, tau)(t) }
	if err := r.AddMagCellSource(MagCellSource{I: 20, J: 20, K: 20, Rate: rate}); err != nil {
		t.Fatal(err)
	}
	probes := [2][3]int{{30, 20, 20}, {27, 27, 20}}
	dist := [2]float64{10 * h, math.Hypot(7, 7) * h}
	var peaks, arrive [2]float64
	var maxB, maxE, worstGauss float64
	for s := 0; s < 260; s++ {
		r.Step()
		tS := r.Time()
		for p, q := range probes {
			if (tS-t0)*em.SpeedOfLight > 2*dist[p] { // direct wave only, as the electric test
				continue
			}
			if v := r.SmAt(q[0], q[1], q[2]); v > peaks[p] {
				peaks[p], arrive[p] = v, tS
			}
		}
		for c := 0; c < 3; c++ {
			for n := range r.B[c] {
				maxB = math.Max(maxB, math.Abs(r.B[c][n]))
				maxE = math.Max(maxE, math.Abs(r.E[c][n]))
			}
		}
		worstGauss = math.Max(worstGauss, r.MagGaussDefect())
	}
	c := em.SpeedOfLight
	for p := range probes {
		want := em.Mu0 / c * q0 / (4 * math.Pi * dist[p])
		lag := (arrive[p] - t0) * c / dist[p]
		t.Logf("r = %.2f mm: peak S_m %.4g T (retarded %.4g, %+.2f %%), arrival at %.3f r/c", dist[p]*1e3, peaks[p], want, 100*(peaks[p]/want-1), lag)
		if math.Abs(peaks[p]/want-1) > 0.05 {
			t.Errorf("probe %d: peak %.4g against retarded %.4g", p, peaks[p], want)
		}
		if math.Abs(lag-1) > 0.08 {
			t.Errorf("probe %d: arrives at %.3f r/c", p, lag)
		}
	}
	t.Logf("max |E| %.2g V/m against max c|B| %.2g V/m; magnetic Gauss defect (AB) %.2g", maxE, c*maxB, worstGauss)
	if maxE > 1e-12*c*maxB {
		t.Errorf("the S_m wave carries E: %.3g", maxE/(c*maxB))
	}
	if worstGauss > 1e-10 {
		t.Errorf("AB magnetic Gauss law defect %.3g", worstGauss)
	}
	dual, _ := New(Config{N: [3]int{21, 21, 21}, D: [3]float64{h, h, h}, Magnetic: MagDual})
	_ = dual.AddMagCellSource(MagCellSource{I: 10, J: 10, K: 10, Rate: rate})
	for s := 0; s < 120; s++ {
		dual.Step()
	}
	t.Logf("S_m held at 0, same source: magnetic Gauss defect %.3g", dual.MagGaussDefect())
	if dual.MagGaussDefect() < 0.5 {
		t.Errorf("MagDual should not hold div B = μ₀ρ_m for created magnetic charge: %.3g", dual.MagGaussDefect())
	}
}

// Duality between the two scalars: created electric charge at rate q̇ and
// created magnetic charge at rate c·q̇ send out the same scalar amplitude,
// S = S_m = μ₀q̇/(4πr), the electric one with E only and the magnetic one
// with B only.
func TestScalarDuality(t *testing.T) {
	const h = 1e-3
	el, _ := New(Config{N: [3]int{41, 41, 41}, D: [3]float64{h, h, h}})
	mg, _ := New(Config{N: [3]int{42, 42, 42}, D: [3]float64{h, h, h}, Magnetic: MagAB})
	tau := 10 * el.Dt
	t0 := 5 * tau
	q0 := 1e-6
	_ = el.AddNodeSource(NodeSource{I: 20, J: 20, K: 20, Rate: func(t float64) float64 { return q0 * gauss(t0, tau)(t) }})
	_ = mg.AddMagCellSource(MagCellSource{I: 20, J: 20, K: 20, Rate: func(t float64) float64 { return em.SpeedOfLight * q0 * gauss(t0, tau)(t) }})
	var pe, pm float64
	for s := 0; s < 90; s++ {
		el.Step()
		mg.Step()
		pe = math.Max(pe, el.SAt(30, 20, 20))
		pm = math.Max(pm, mg.SmAt(30, 20, 20))
	}
	t.Logf("peak at 10 cells: electric S %.5g T, magnetic S_m %.5g T (ratio %.4f)", pe, pm, pm/pe)
	if math.Abs(pm/pe-1) > 0.02 {
		t.Errorf("dual scalars differ: S %.4g, S_m %.4g", pe, pm)
	}
}
