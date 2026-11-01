// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"testing"
)

// The Hellwarth–Nouchi pulse is an exact solution: checked by finite
// differences of the formula itself (curl E = −∂B/∂t, curl B = ∂E/∂t/c²,
// div E = div B = 0).
func TestFlyingDoughnutMaxwell(t *testing.T) {
	d := FlyingDoughnut{Q1: 0.01, Q2: 0.1, F0: 1}
	c := SpeedOfLight
	h, ht := 1e-6, 1e-6/c
	pts := [][4]float64{{0.01, 0.02, 0.005, 0}, {-0.03, 0.01, -0.02, 0.05 / c}, {0.02, -0.02, 0.04, 0.02 / c}}
	for _, p := range pts {
		d1 := func(a int, field func(x, y, z, t float64) (E, B [3]float64)) (dE, dB [3]float64) {
			q := p
			q[a] += h
			if a == 3 {
				q[3] = p[3] + ht
			}
			Ep, Bp := field(q[0], q[1], q[2], q[3])
			q = p
			q[a] -= h
			if a == 3 {
				q[3] = p[3] - ht
			}
			Em, Bm := field(q[0], q[1], q[2], q[3])
			s := 2 * h
			if a == 3 {
				s = 2 * ht
			}
			for i := 0; i < 3; i++ {
				dE[i], dB[i] = (Ep[i]-Em[i])/s, (Bp[i]-Bm[i])/s
			}
			return
		}
		var gE, gB [4][3]float64
		for a := 0; a < 4; a++ {
			gE[a], gB[a] = d1(a, d.Fields)
		}
		curl := func(g [4][3]float64) [3]float64 {
			return [3]float64{g[1][2] - g[2][1], g[2][0] - g[0][2], g[0][1] - g[1][0]}
		}
		cE, cB := curl(gE), curl(gB)
		scale := 0.0
		for i := 0; i < 3; i++ {
			scale = math.Max(scale, math.Max(math.Abs(cE[i]), math.Abs(gB[3][i])))
		}
		for i := 0; i < 3; i++ {
			if r := math.Abs(cE[i]+gB[3][i]) / scale; r > 1e-6 {
				t.Errorf("Faraday residual %g at %v", r, p)
			}
			if r := math.Abs(cB[i]-gE[3][i]/(c*c)) / (scale / c); r > 1e-6 {
				t.Errorf("Ampère residual %g at %v", r, p)
			}
		}
		if r := math.Abs(gB[0][0]+gB[1][1]+gB[2][2]) / scale; r > 1e-6 {
			t.Errorf("div B %g at %v", r, p)
		}
	}
}

// A closed copper-walled box: the TM110 mode rings at the Yee scheme's
// own frequency (exactly) and the continuum frequency (to (k·dx)²), and
// the field energy stays constant to round-off.
func TestFDTDCavity(t *testing.T) {
	nx, ny, nz := 21, 17, 6
	dx := 0.01
	f, err := NewFDTD(FDTDConfig{N: [3]int{nx, ny, nz}, D: [3]float64{dx, dx, dx}, Ledger: true})
	if err != nil {
		t.Fatal(err)
	}
	a, b := float64(nx-1)*dx, float64(ny-1)*dx
	f.SetE(func(x, y, z float64) [3]float64 {
		return [3]float64{0, 0, math.Sin(math.Pi*x/a) * math.Sin(math.Pi*y/b)}
	})
	probe := f.index(nx/2, ny/2, 2)
	prev := f.E[2][probe]
	var cross []float64
	for f.Steps < 4000 {
		f.Step()
		v := f.E[2][probe]
		if (prev > 0) != (v > 0) {
			cross = append(cross, f.Time()-f.Dt*v/(v-prev))
		}
		prev = v
	}
	period := 2 * (cross[len(cross)-1] - cross[0]) / float64(len(cross)-1)
	w := 2 * math.Pi / period
	c := SpeedOfLight
	s := c * f.Dt * math.Sqrt(math.Pow(2/dx*math.Sin(math.Pi*dx/(2*a)), 2)+math.Pow(2/dx*math.Sin(math.Pi*dx/(2*b)), 2)) / 2
	wYee := 2 / f.Dt * math.Asin(s)
	wCont := c * math.Pi * math.Hypot(1/a, 1/b)
	if r := math.Abs(w-wYee) / wYee; r > 1e-5 {
		t.Errorf("ω %.8g, Yee %.8g (%.2g)", w, wYee, r)
	}
	if r := math.Abs(w-wCont) / wCont; r > 2e-3 {
		t.Errorf("ω %.8g, continuum %.8g (%.2g)", w, wCont, r)
	}
	l := f.Ledger()
	if r := math.Abs(l.Residual()) / l.W0; r > 1e-12 {
		t.Errorf("energy drift %.3g of %.3g J", l.Residual(), l.W0)
	}
	t.Logf("TM110: ω/ω_Yee − 1 = %.2g, ω/ω_cont − 1 = %.2g, energy residual %.2g", w/wYee-1, w/wCont-1, l.Residual()/l.W0)
}

// Sources in a closed box with a lossy slab: div B stays at round-off,
// Gauss's law holds at every node with the charge the sources carried,
// and the ledger (source work = field energy + heat) closes.
func TestFDTDGaussAndLedger(t *testing.T) {
	n := 18
	dx := 0.005
	f, err := NewFDTD(FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx}, Ledger: true,
		Material: func(x, y, z float64) (float64, float64) {
			if z > 0.06 {
				return 4, 0.05
			}
			return 1, 0
		}})
	if err != nil {
		t.Fatal(err)
	}
	pulse := func(t0, w float64) func(float64) float64 {
		return func(t float64) float64 { return math.Exp(-math.Pow((t-t0)/w, 2)) }
	}
	srcs := []EdgeSource{{2, 6, 7, 5, pulse(2e-10, 5e-11)}, {0, 9, 9, 9, pulse(3e-10, 8e-11)}, {1, 4, 10, 8, func(t float64) float64 { return math.Sin(3e10 * t) }}}
	for _, s := range srcs {
		if err := f.AddSource(s); err != nil {
			t.Fatal(err)
		}
	}
	if f.AddSource(EdgeSource{0, 3, 0, 4, pulse(0, 1)}) == nil {
		t.Error("source on a conducting wall accepted")
	}
	q := map[[3]int]float64{}
	for f.Steps < 600 {
		th := f.Time() + f.Dt/2
		for _, s := range srcs {
			a, b := [3]int{s.I, s.J, s.K}, [3]int{s.I, s.J, s.K}
			b[s.Comp]++
			q[a] -= s.Current(th) * f.Dt
			q[b] += s.Current(th) * f.Dt
		}
		f.Step()
	}
	worst, scale := f.DivB()
	if worst > 1e-12*scale {
		t.Errorf("div B %.3g against scale %.3g", worst, scale)
	}
	qmax, err2 := 0.0, 0.0
	for k := 1; k < 11; k++ { // nodes in vacuum (z < 0.05 m); the slab conducts charge itself
		for j := 1; j < n-1; j++ {
			for i := 1; i < n-1; i++ {
				qmax = math.Max(qmax, math.Abs(q[[3]int{i, j, k}]))
				err2 = math.Max(err2, math.Abs(f.NodeCharge(i, j, k)-q[[3]int{i, j, k}]))
			}
		}
	}
	if err2 > 1e-9*qmax {
		t.Errorf("Gauss's law off by %.3g C (charges up to %.3g C)", err2, qmax)
	}
	l := f.Ledger()
	if r := math.Abs(l.Residual()) / math.Abs(l.SourceWork); r > 1e-10 {
		t.Errorf("ledger residual %.3g of source work %.3g J", l.Residual(), l.SourceWork)
	}
	if l.MaterialLoss <= 0 {
		t.Error("no heat in the lossy slab")
	}
	t.Logf("div B %.2g of scale; Gauss %.2g of %.2g C; ledger: work %.4g J = field %.4g + heat %.4g (residual %.2g)",
		worst/scale, err2, qmax, l.SourceWork, l.W-l.W0, l.MaterialLoss, l.Residual())
}

// A short current element (Hertzian dipole) in an absorbing box radiates
// Larmor's power ω⁴p²/(12πε0c³), p = I0·dz/ω; the absorbed power is the
// radiated power and the ledger closes.
func TestFDTDDipoleRadiation(t *testing.T) {
	if testing.Short() {
		t.Skip("long")
	}
	dx := 0.005
	lambda := 40 * dx
	n := 2*24 + 2*14 + 1
	f, err := NewFDTD(FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx}, Absorber: 14, Ledger: true})
	if err != nil {
		t.Fatal(err)
	}
	w := 2 * math.Pi * SpeedOfLight / lambda
	I0 := 1.0
	T := 2 * math.Pi / w
	ramp := 3 * T
	if err := f.AddSource(EdgeSource{2, n / 2, n / 2, n / 2, func(t float64) float64 {
		env := 1.0
		if t < ramp {
			env = math.Pow(math.Sin(math.Pi/2*t/ramp), 2)
		}
		return I0 * env * math.Sin(w*t)
	}}); err != nil {
		t.Fatal(err)
	}
	f.Run(10 * T)
	a0, t0 := f.Ledger().Absorbed, f.Time()
	f.Run(14 * T)
	l := f.Ledger()
	P := (l.Absorbed - a0) / (f.Time() - t0)
	p := I0 * dx / w
	want := math.Pow(w, 4) * p * p / (12 * math.Pi * Eps0 * math.Pow(SpeedOfLight, 3)) // Larmor's time average for amplitude p
	if r := P/want - 1; math.Abs(r) > 0.05 {
		t.Errorf("radiated %.5g W, Larmor %.5g W (%+.3f)", P, want, r)
	}
	if r := math.Abs(l.Residual()) / l.SourceWork; r > 1e-10 {
		t.Errorf("ledger residual %.3g of %.3g J", l.Residual(), l.SourceWork)
	}
	t.Logf("dipole λ = 40 cells: radiated %.5g W, Larmor %.5g W (%+.2f %%); ledger residual %.2g of %.3g J", P, want, 100*(P/want-1), l.Residual(), l.SourceWork)
}

// The flying doughnut, started from the exact fields before its focus,
// propagates through the focus as the exact solution says. The error
// (3.8 % at 8 cells per q1) has two parts: grid dispersion, and the
// part of the converging pulse that starts outside the box and so never
// arrives; the second does not shrink with finer cells.
func TestFDTDFlyingDoughnut(t *testing.T) {
	if testing.Short() {
		t.Skip("long")
	}
	r, err := RunFlyingDoughnut(FlyingDoughnut{Q1: 0.01, Q2: 0.04, F0: 1}, 8, 8, 6)
	if err != nil {
		t.Fatal(err)
	}
	if r.RelError > 0.05 {
		t.Errorf("FDTD vs exact: relative L2 error %.3g", r.RelError)
	}
	if r.DivB > 1e-12 {
		t.Errorf("div B %.3g", r.DivB)
	}
	t.Logf("flying doughnut: L2 error %.3g after %.3g ns, div B %.2g, energy %.4g→%.4g J (absorbed %.2g)", r.RelError, r.Duration*1e9, r.DivB, r.EnergyStart, r.EnergyEnd, r.Absorbed)
}
