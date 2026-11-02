// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

const (
	mu0Test = 1.25663706212e-6
	c0Test  = 299792458.0
)

// Kron 1962 Fig. 1 in SI: the branch along each axis scales as the cell's
// cross-section over its length, times ε, σ or μ, with the scale factors
// of curvilinear coordinates; the inductor of the figure is a quarter of
// the face inductance.
func TestFieldCellValues(t *testing.T) {
	c := FieldCell{Dx: 2, Dy: 3, Dz: 5, Eps: 7, Mu: 11, Sigma: 13}
	near := func(what string, got, want float64) {
		t.Helper()
		if math.Abs(got-want) > 1e-12*math.Abs(want) {
			t.Errorf("%s = %g, want %g", what, got, want)
		}
	}
	near("C_x", c.Capacitance(0), 3*5/2.0*7)
	near("C_y", c.Capacitance(1), 2*5/3.0*7)
	near("G_z", c.Conductance(2), 2*3/5.0*13)
	near("L_x face", c.FaceInductance(0), 3*5/2.0*11)
	near("L_x quarter", c.QuarterInductance(0), 3*5/2.0*11/4)
	c.H1, c.H2, c.H3 = 2, 3, 5 // (h₂h₃/h₁) etc.
	near("C_x curvilinear", c.Capacitance(0), 3*5/2.0*3*5/2.0*7)
	near("L_y curvilinear", c.FaceInductance(1), 2*5/3.0*2*5/3.0*11)
}

// kronVoltages solves a lattice in Kron's node-pair frame at ω (no
// operating point, so lossless inductor loops and capacitor-only nodes
// are fine) and returns the voltage of any named element.
func kronVoltages(t *testing.T, c *circuit.Circuit, omega float64) func(elem string) complex128 {
	t.Helper()
	nw, err := Build(c)
	if err != nil {
		t.Fatal(err)
	}
	sol, err := nw.Solve(NodeForm, omega, false)
	if err != nil {
		t.Fatal(err)
	}
	return func(elem string) complex128 {
		k, ok := nw.ElementBranch(elem)
		if !ok {
			t.Fatalf("no element %s", elem)
		}
		return sol.E.At(k)
	}
}

// peak finds the frequency of the largest |V| across element elem of a
// circuit between f0 and f1: a scan, then a golden-section search (the
// driving-point response of a lossless lattice has a pole at each
// resonance).
func peak(t *testing.T, c *circuit.Circuit, elem string, f0, f1 float64) float64 {
	t.Helper()
	nw, err := Build(c)
	if err != nil {
		t.Fatal(err)
	}
	k, ok := nw.ElementBranch(elem)
	if !ok {
		t.Fatalf("no element %s", elem)
	}
	// A singular solve is the pole itself (a lossless lattice exactly at
	// resonance).
	resp := func(f float64) float64 {
		sol, err := nw.Solve(NodeForm, 2*math.Pi*f, false)
		if err != nil {
			return math.Inf(1)
		}
		return cmplx.Abs(sol.E.At(k))
	}
	const scan = 80
	best, bf := -1.0, f0
	for k := 0; k <= scan; k++ {
		f := f0 + (f1-f0)*(float64(k)+0.37)/scan
		if v := resp(f); v > best {
			best, bf = v, f
		}
	}
	lo, hi := bf-(f1-f0)/scan, bf+(f1-f0)/scan
	g := (math.Sqrt(5) - 1) / 2
	a, b := hi-g*(hi-lo), lo+g*(hi-lo)
	fa, fb := resp(a), resp(b)
	for k := 0; k < 60; k++ {
		if fa > fb {
			hi, b, fb = b, a, fa
			a = hi - g*(hi-lo)
			fa = resp(a)
		} else {
			lo, a, fa = a, b, fb
			b = lo + g*(hi-lo)
			fb = resp(b)
		}
	}
	return (lo + hi) / 2
}

// The TE lattice of a closed box resonates at the lattice's own
// dispersion value (the spatial part of Yee's, as em/fdtd.go's
// TestFDTDCavity and TestU1YeeCavity use) and at the continuum frequency
// to second order in Δ/λ; Kron's node-pair solve and MNA (with a 1e-9 Ω
// series leak so the operating point exists) agree to 1e-6.
func TestTELatticeResonance(t *testing.T) {
	eps0 := 1 / (mu0Test * c0Test * c0Test)
	l := TELattice{Nx: 12, Nz: 10, Cell: FieldCell{Dx: 0.01, Dy: 0.01, Dz: 0.01, Eps: eps0, Mu: mu0Test}, Walls: true}
	build := func(l TELattice) *circuit.Circuit {
		b := circuit.NewBuilder("TE11 box")
		if err := l.Build(b); err != nil {
			t.Fatal(err)
		}
		b.Add(circuit.Element{Name: "ID", Kind: circuit.CurrentSource, ACMag: 1e-3}, "0", l.Node(3, 4))
		c, err := b.Build()
		if err != nil {
			t.Fatal(err)
		}
		return c
	}
	c := build(l)
	wl, wc := l.Dispersion(1, 1)
	probe := "C_6_5"
	f := peak(t, c, probe, 0.9*wl/(2*math.Pi), 1.1*wl/(2*math.Pi))
	w := 2 * math.Pi * f
	if r := math.Abs(w-wl) / wl; r > 1e-6 {
		t.Errorf("TE11: ω %.8g, lattice %.8g (%.2g)", w, wl, r)
	}
	if r := math.Abs(w-wc) / wc; r > 5e-3 {
		t.Errorf("TE11: ω %.8g, continuum %.8g (%.2g)", w, wc, r)
	}
	t.Logf("TE11 on Kron's lattice: ω/ω_lattice − 1 = %.2g, ω/ω_continuum − 1 = %.2g", w/wl-1, w/wc-1)
	// Second order: halving the cell quarters the continuum error.
	l2 := l
	l2.Nx, l2.Nz = 24, 20
	l2.Cell.Dx, l2.Cell.Dy, l2.Cell.Dz = 0.005, 0.005, 0.005
	wl2, wc2 := l2.Dispersion(1, 1)
	if r := (wl2/wc2 - 1) / (wl/wc - 1); math.Abs(r-0.25) > 0.02 {
		t.Errorf("error ratio on halving the cell %.3g, want 0.25", r)
	}
	// MNA twin off resonance, with the series leak, and its ledger.
	lm := l
	lm.SeriesR = 1e-9
	cm := build(lm)
	fo := 0.8 * f
	r, err := analysis.AC(cm, analysis.ACSpec{Variation: "lin", Points: 1, FStart: fo, FStop: fo}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	kv := kronVoltages(t, c, 2*math.Pi*fo)
	worst := 0.0
	for i := 1; i < l.Nx; i++ {
		for j := 1; j < l.Nz; j++ {
			v, _ := r.NodeVoltage(l.Node(i, j))
			k := kv(fmt.Sprintf("C_%d_%d", i, j))
			worst = math.Max(worst, cmplx.Abs(v[0]-k)/cmplx.Abs(k))
		}
	}
	if worst > 1e-6 {
		t.Errorf("Kron node-pair solve vs MNA: %.2g", worst)
	}
	if lg := r.Ledger; lg.MaxRelativeResidual > 1e-9 {
		t.Errorf("AC ledger relative residual %.3g", lg.MaxRelativeResidual)
	}
	t.Logf("Kron vs MNA at 0.8 f₁₁: %.2g; AC ledger residual %.2g", worst, r.Ledger.MaxRelativeResidual)
}

// Kron's 1945 cavity (J. Appl. Phys. Fig. 14): the lattice solved exactly
// in Kron's frame against the junction potentials H_y Kron read on the
// network analyzer, whose unbalanced currents (Fig. 15) he put at
// 8 percent of the largest coil current.
func TestKron1945Cavity(t *testing.T) {
	c, err := Kron1945Cavity(1e9)
	if err != nil {
		t.Fatal(err)
	}
	kv := kronVoltages(t, c, 1)
	var got [5][5]float64
	vmax := 0.0
	// A lossless reactive lattice driven by in-phase currents has its
	// potentials in quadrature; the analyzer read their signed magnitudes.
	sign := 1.0
	if imag(kv("C_0_0")) < 0 {
		sign = -1
	}
	for i := 0; i < 5; i++ {
		for j := 0; j < 5; j++ {
			v := kv(fmt.Sprintf("C_%d_%d", i, j))
			if re := math.Abs(real(v)); re > 1e-6*math.Abs(imag(v))+1e-9 {
				t.Errorf("h_%d_%d not in quadrature with the drive: %v", i, j, v)
			}
			got[i][j] = sign * imag(v)
			vmax = math.Max(vmax, math.Abs(got[i][j]))
		}
	}
	// Antisymmetry about the centre column, exact on the lattice.
	for i := 0; i < 2; i++ {
		for j := 0; j < 5; j++ {
			if d := math.Abs(got[i][j] + got[4-i][j]); d > 1e-9*vmax {
				t.Errorf("antisymmetry h_%d_%d + h_%d_%d = %g", i, j, 4-i, j, d)
			}
		}
	}
	worst, sum, n := 0.0, 0.0, 0
	for j := 0; j < 5; j++ {
		var row string
		for i := 0; i < 5; i++ {
			k := Kron1945Potentials[j][i]
			d := math.Abs(got[i][j]-k) / vmax
			worst = math.Max(worst, d)
			sum += d
			n++
			row += fmt.Sprintf("  %8.2f (Kron %8.2f)", got[i][j], k)
		}
		t.Log(row)
	}
	t.Logf("H_y against Kron's analyzer values: worst %.3f, mean %.3f of the largest potential %.1f", worst, sum/float64(n), vmax)
	if worst > 0.12 {
		t.Errorf("worst deviation %.3f of the largest potential exceeds the analyzer's 8 %% unbalance by too much", worst)
	}
}

// Balasubramanian's waveguide (thesis Figs. 15–16): the error of the
// lattice's matching impedance against the theoretical TE₁₀ value grows
// with (Δz/λ)², about 1 % at Δz/λ = 5 % and 4 % at 10 %.
func TestWaveguideMatching(t *testing.T) {
	eps0 := 1 / (mu0Test * c0Test * c0Test)
	w := Waveguide{A: 0.11, Dy: 0.01, N: 5, Sections: 24, Freq: 2e9, Eps: eps0, Mu: mu0Test}
	if b := w.Beta(); math.Abs(b-0.731955) > 1e-5 {
		t.Errorf("βc/ω = %.6f, thesis 0.731955", b)
	}
	if z := w.Theory(); math.Abs(z/(2.3395e4*w.Dy)-1) > 2e-3 {
		t.Errorf("theoretical matching impedance %.5g, thesis 2.3395e4·Δy = %.5g", z, 2.3395e4*w.Dy)
	}
	ripple := func(r float64) float64 {
		c, _, err := w.Circuit(r)
		if err != nil {
			t.Fatal(err)
		}
		kv := kronVoltages(t, c, 2*math.Pi*w.Freq)
		lo, hi := math.Inf(1), 0.0
		for j := 3; j <= w.Sections-3; j++ {
			m := cmplx.Abs(kv(fmt.Sprintf("C_%d_%d", 2, j)))
			lo, hi = math.Min(lo, m), math.Max(hi, m)
		}
		return (hi - lo) / (hi + lo)
	}
	lambda := c0Test / w.Freq
	type pt struct{ dz, thesis float64 }
	for _, p := range []pt{{0.0075, 0.010}, {0.011, 0.016}, {0.0165, 0.034}, {0.022, 0.062}} {
		w.Dz = p.dz
		z0 := w.Theory()
		lo, hi := 0.7*z0, 1.3*z0
		g := (math.Sqrt(5) - 1) / 2
		a, b := hi-g*(hi-lo), lo+g*(hi-lo)
		fa, fb := ripple(a), ripple(b)
		for k := 0; k < 40; k++ {
			if fa < fb {
				hi, b, fb = b, a, fa
				a = hi - g*(hi-lo)
				fa = ripple(a)
			} else {
				lo, a, fa = a, b, fb
				b = lo + g*(hi-lo)
				fb = ripple(b)
			}
		}
		zm := (lo + hi) / 2
		zpi := w.LatticeImpedance()
		errPct := zm/z0 - 1
		rip := ripple(zm)
		// The transverse part of the lattice's error (five cells across the
		// width) is the same at every Δz; what varies with Δz is the
		// longitudinal part, which is what the thesis plotted.
		w0 := w
		w0.Dz = 1e-5
		transverse := w0.LatticeImpedance()/z0 - 1
		t.Logf("Δz/λ = %.1f %%: matched at %.5g (lattice Z_π %.5g), theory %.5g, error %+.2f %% = transverse %+.2f %% + longitudinal %+.2f %% (thesis ≈ %.1f %%), residual ripple %.1e",
			100*p.dz/lambda, zm, zpi, z0, 100*errPct, 100*transverse, 100*((1+errPct)/(1+transverse)-1), 100*p.thesis, rip)
		if math.Abs(zm/zpi-1) > 1e-4 {
			t.Errorf("Δz = %g: minimum-ripple termination %.5g, lattice Z_π %.5g", p.dz, zm, zpi)
		}
		if rip > 1e-4 {
			t.Errorf("Δz = %g: residual ripple %.2g at the match", p.dz, rip)
		}
		// The thesis's curve is the longitudinal (Δz/λ)² part (its points
		// read from Fig. 15 to about ±0.3 %); the lattice also carries the
		// fixed transverse discretization of the width into five cells.
		if long := (1+errPct)/(1+transverse) - 1; math.Abs(long-p.thesis) > 0.005 {
			t.Errorf("Δz = %g: longitudinal error %.4f, thesis %.3f", p.dz, long, p.thesis)
		}
	}
}
