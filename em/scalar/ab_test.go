// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

func gauss(t0, tau float64) func(float64) float64 {
	return func(t float64) float64 { x := (t - t0) / tau; return math.Exp(-x * x) }
}

// Gate 1 (§3k): with every current conserved and S = 0 at the start, AB
// electrodynamics is Maxwell's: S stays zero and E, B equal em.FDTD's edge
// by edge.
func TestReducesToYee(t *testing.T) {
	N, D := [3]int{14, 12, 10}, [3]float64{0.01, 0.012, 0.009}
	ref, err := em.NewFDTD(em.FDTDConfig{N: N, D: D, Courant: 0.99})
	if err != nil {
		t.Fatal(err)
	}
	ab, err := New(Config{N: N, D: D, Courant: 0.99, Mode: AB})
	if err != nil {
		t.Fatal(err)
	}
	g := gauss(40*ab.Dt, 10*ab.Dt)
	cur := func(t float64) float64 { return 2.5 * g(t) * math.Sin(3e9*t) }
	if err := ref.AddSource(em.EdgeSource{Comp: 2, I: 7, J: 6, K: 4, Current: cur}); err != nil {
		t.Fatal(err)
	}
	if err := ab.AddSource(EdgeSource{Comp: 2, I: 7, J: 6, K: 4, Current: cur}); err != nil {
		t.Fatal(err)
	}
	for s := 0; s < 160; s++ {
		ref.Step()
		ab.Step()
	}
	var dE, dB, mE, mB, mS float64
	for c := 0; c < 3; c++ {
		for n := range ref.E[c] {
			dE = math.Max(dE, math.Abs(ref.E[c][n]-ab.E[c][n]))
			dB = math.Max(dB, math.Abs(ref.B[c][n]-ab.B[c][n]))
			mE = math.Max(mE, math.Abs(ref.E[c][n]))
			mB = math.Max(mB, math.Abs(ref.B[c][n]))
		}
	}
	for _, s := range ab.S {
		mS = math.Max(mS, math.Abs(s))
	}
	t.Logf("E diff %.2g of %.3g V/m, B diff %.2g of %.3g T, max |S| %.2g T (c|S| / |E| = %.2g)",
		dE, mE, dB, mB, mS, em.SpeedOfLight*mS/mE)
	if dE > 1e-12*mE || dB > 1e-12*mB {
		t.Errorf("AB with S = 0 differs from em.FDTD: E %.3g, B %.3g", dE/mE, dB/mB)
	}
	if em.SpeedOfLight*mS > 1e-12*mE {
		t.Errorf("S grew without an extra current: c|S|/|E| = %.3g", em.SpeedOfLight*mS/mE)
	}
}

// Gate 3: the energy theorem closes step by step, with the hypothetical
// column c²ρS non-zero once S is sourced.
func TestEnergyTheorem(t *testing.T) {
	r, err := New(Config{N: [3]int{17, 15, 13}, D: [3]float64{0.01, 0.01, 0.01}})
	if err != nil {
		t.Fatal(err)
	}
	g := gauss(30*r.Dt, 8*r.Dt)
	_ = r.AddSource(EdgeSource{Comp: 0, I: 6, J: 7, K: 6, Current: func(t float64) float64 { return 3 * g(t) * math.Cos(2e9*t) }})
	_ = r.AddNodeSource(NodeSource{I: 9, J: 7, K: 6, Rate: func(t float64) float64 { return 1e-9 * g(t) / r.Dt }})
	r.SetS(func(x, y, z float64) float64 {
		return 1e-9 * math.Sin(math.Pi*x/0.16) * math.Sin(math.Pi*y/0.14) * math.Sin(math.Pi*z/0.12)
	})
	worst := 0.0
	for s := 0; s < 400; s++ {
		r.Step()
		l := r.Ledger()
		scale := math.Max(math.Abs(l.W), math.Max(math.Abs(l.SourceWork), math.Abs(l.ScalarWork)))
		worst = math.Max(worst, math.Abs(l.Residual())/scale)
	}
	l := r.Ledger()
	t.Logf("W %.4g J (S part %.3g), source work %.4g J, scalar work c²ρS %.4g J, worst residual %.2g",
		l.W, l.WS, l.SourceWork, l.ScalarWork, worst)
	if worst > 1e-11 {
		t.Errorf("energy theorem residual %.3g", worst)
	}
	if l.ScalarWork == 0 {
		t.Error("scalar column is zero with an extra current present")
	}
}

// Gate 2 and the source law □S = μ₀I: charge appearing at a node at rate
// q̇(t) sends out S(r, t) = μ₀ q̇(t − r/c)/(4πr), the retarded solution,
// with B = 0 (the wave is longitudinal). Gauss's law in AB form holds to
// round-off; in standard mode the same source breaks it.
func TestRetardedScalarWave(t *testing.T) {
	const h = 1e-3
	var peaks [2]float64
	var arrive [2]float64
	r, err := New(Config{N: [3]int{41, 41, 41}, D: [3]float64{h, h, h}})
	if err != nil {
		t.Fatal(err)
	}
	tau := 10 * r.Dt
	t0 := 5 * tau
	q0 := 1e-6 // peak rate [C/s]
	_ = r.AddNodeSource(NodeSource{I: 20, J: 20, K: 20, Rate: func(t float64) float64 { return q0 * gauss(t0, tau)(t) }})
	probes := [2][3]int{{30, 20, 20}, {27, 27, 20}}
	dist := [2]float64{10 * h, math.Hypot(7, 7) * h}
	var maxB, maxE float64
	for s := 0; s < 260; s++ {
		r.Step()
		tS := r.Time() + r.Dt/2
		for p, q := range probes {
			// direct wave only: the first wall echo (inverted) needs a path
			// of 30 cells, so stop the search at 2 r/c after the source peak.
			if (tS-t0)*em.SpeedOfLight > 2*dist[p] {
				continue
			}
			v := r.SAt(q[0], q[1], q[2])
			if v > peaks[p] {
				peaks[p], arrive[p] = v, tS
			}
		}
		for c := 0; c < 3; c++ {
			for n := range r.B[c] {
				maxB = math.Max(maxB, math.Abs(r.B[c][n]))
				maxE = math.Max(maxE, math.Abs(r.E[c][n]))
			}
		}
	}
	for p := range probes {
		want := em.Mu0 * q0 / (4 * math.Pi * dist[p])
		lag := (arrive[p] - t0) * em.SpeedOfLight / dist[p]
		t.Logf("r = %.2f mm: peak S %.4g T (retarded %.4g, %+.2f %%), arrival at %.3f r/c", dist[p]*1e3, peaks[p], want, 100*(peaks[p]/want-1), lag)
		if math.Abs(peaks[p]/want-1) > 0.05 {
			t.Errorf("probe %d: peak %.4g against retarded %.4g", p, peaks[p], want)
		}
		if math.Abs(lag-1) > 0.08 {
			t.Errorf("probe %d: arrives at %.3f r/c", p, lag)
		}
	}
	t.Logf("max |B| %.2g T against max |E|/c %.2g T; Gauss defect (AB) %.2g", maxB, maxE/em.SpeedOfLight, r.GaussDefect())
	if maxB > 1e-12*maxE/em.SpeedOfLight {
		t.Errorf("the S wave carries B: %.3g", maxB*em.SpeedOfLight/maxE)
	}
	if d := r.GaussDefect(); d > 1e-10 {
		t.Errorf("AB Gauss law defect %.3g", d)
	}
	std, _ := New(Config{N: [3]int{21, 21, 21}, D: [3]float64{h, h, h}, Mode: Standard})
	_ = std.AddNodeSource(NodeSource{I: 10, J: 10, K: 10, Rate: func(t float64) float64 { return q0 * gauss(t0, tau)(t) }})
	for s := 0; s < 120; s++ {
		std.Step()
	}
	t.Logf("standard mode, same source: Gauss defect %.3g (charge with no field to carry it)", std.GaussDefect())
	if std.GaussDefect() < 0.5 {
		t.Errorf("standard mode should not hold Gauss's law for created charge: %.3g", std.GaussDefect())
	}
}

// Free S modes in the conducting box: E longitudinal, B = 0, at the
// discrete Dirichlet frequency of the grid (the TM cutoff, Pinchuk et al.),
// energy constant, and S and E share it as a wave does.
func TestCavityScalarMode(t *testing.T) {
	N, D := [3]int{13, 11, 9}, [3]float64{0.01, 0.012, 0.015}
	r, err := New(Config{N: N, D: D})
	if err != nil {
		t.Fatal(err)
	}
	var L [3]float64
	lam := 0.0
	for a := 0; a < 3; a++ {
		L[a] = float64(N[a]-1) * D[a]
		s := math.Sin(math.Pi * D[a] / (2 * L[a]))
		lam += 4 / (D[a] * D[a]) * s * s
	}
	r.SetS(func(x, y, z float64) float64 {
		return 1e-6 * math.Sin(math.Pi*x/L[0]) * math.Sin(math.Pi*y/L[1]) * math.Sin(math.Pi*z/L[2])
	})
	c := em.SpeedOfLight
	wDisc := 2 / r.Dt * math.Asin(c*r.Dt*math.Sqrt(lam)/2)
	wCont := c * math.Pi * math.Sqrt(1/(L[0]*L[0])+1/(L[1]*L[1])+1/(L[2]*L[2]))
	ci, cj, ck := (N[0]-1)/2, (N[1]-1)/2, (N[2]-1)/2
	var crossings []float64
	prev, tPrev := r.SAt(ci, cj, ck), r.Time()+r.Dt/2
	var maxB, maxE, worst float64
	steps := int(40 * 2 * math.Pi / wDisc / r.Dt)
	for s := 0; s < steps; s++ {
		r.Step()
		v, tv := r.SAt(ci, cj, ck), r.Time()+r.Dt/2
		if prev < 0 && v >= 0 {
			crossings = append(crossings, tPrev+(tv-tPrev)*(-prev)/(v-prev))
		}
		prev, tPrev = v, tv
		for cc := 0; cc < 3; cc++ {
			for n := range r.B[cc] {
				maxB = math.Max(maxB, math.Abs(r.B[cc][n]))
				maxE = math.Max(maxE, math.Abs(r.E[cc][n]))
			}
		}
		l := r.Ledger()
		worst = math.Max(worst, math.Abs(l.Residual())/l.W)
	}
	if len(crossings) < 10 {
		t.Fatalf("only %d periods seen", len(crossings))
	}
	period := (crossings[len(crossings)-1] - crossings[0]) / float64(len(crossings)-1)
	wMeas := 2 * math.Pi / period
	t.Logf("ω measured %.7g, discrete Dirichlet %.7g (%.2g), continuum TM cutoff %.7g (%+.2f %%)",
		wMeas, wDisc, wMeas/wDisc-1, wCont, 100*(wMeas/wCont-1))
	t.Logf("max |B| %.2g T (|E|/c %.2g T); energy residual %.2g", maxB, maxE/c, worst)
	if math.Abs(wMeas/wDisc-1) > 1e-5 {
		t.Errorf("S mode at %.7g, discrete theory %.7g", wMeas, wDisc)
	}
	if math.Abs(wMeas/wCont-1) > 0.03 {
		t.Errorf("S mode %.3g %% from the TM cutoff", 100*(wMeas/wCont-1))
	}
	if maxB > 1e-12*maxE/c {
		t.Errorf("free S mode carries B: %.3g", maxB*c/maxE)
	}
	if worst > 1e-11 {
		t.Errorf("energy residual %.3g", worst)
	}
}
