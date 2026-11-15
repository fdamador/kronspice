// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"testing"
)

func maxDiff(a, b [][]float64) float64 {
	d := 0.0
	for i := range a {
		for j := range a[i] {
			d = math.Max(d, math.Abs(a[i][j]-b[i][j]))
		}
	}
	return d
}

// One conductor: the surge impedance is √(L/C) and the velocity 1/√(LC),
// as for TLine.
func TestLosslessModesSingle(t *testing.T) {
	tl := TLine{L: 250e-9, C: 100e-12, Length: 1}
	m, err := LosslessModes([][]float64{{tl.L}}, [][]float64{{tl.C}})
	if err != nil {
		t.Fatal(err)
	}
	if z := real(tl.Z0(0)); math.Abs(m.Zc[0][0]-z) > 1e-12*z {
		t.Errorf("Zc = %g, want %g", m.Zc[0][0], z)
	}
	if v := 1 / tl.Delay(); math.Abs(m.Velocity[0]-v) > 1e-12*v {
		t.Errorf("v = %g, want %g", m.Velocity[0], v)
	}
}

// Bewley Eq. 6‑15: conductors over ground in a homogeneous medium have
// L·K = 1/c², so every mode travels at c, while the surge impedance matrix
// still couples the conductors.
func TestLosslessModesHomogeneous(t *testing.T) {
	const c0 = 299792458.0
	k := [][]float64{{60e-12, -12e-12, -3e-12}, {-12e-12, 65e-12, -12e-12}, {-3e-12, -12e-12, 60e-12}}
	ki, err := invert(k)
	if err != nil {
		t.Fatal(err)
	}
	l := make([][]float64, 3)
	for i := range l {
		l[i] = make([]float64, 3)
		for j := range l[i] {
			l[i][j] = ki[i][j] / (c0 * c0)
		}
	}
	m, err := LosslessModes(l, k)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range m.Velocity {
		if math.Abs(v-c0) > 1e-9*c0 {
			t.Errorf("mode %d: v = %g, want c", i, v)
		}
	}
	if d := maxDiff(matMul(matMul(m.Zc, k), m.Zc), l); d > 1e-9*l[0][0] {
		t.Errorf("Zc·K·Zc − L = %g", d)
	}
}

// An inhomogeneous pair (a microstrip-like line): two velocities, and the
// identities L·K·T = T·diag(1/v²) and Zc·K·Zc = L hold.
func TestLosslessModesInhomogeneous(t *testing.T) {
	l := [][]float64{{311.7e-9, 71.0e-9}, {71.0e-9, 311.7e-9}}
	k := [][]float64{{144.9e-12, -18.6e-12}, {-18.6e-12, 144.9e-12}}
	m, err := LosslessModes(l, k)
	if err != nil {
		t.Fatal(err)
	}
	if m.Velocity[0] >= m.Velocity[1] {
		t.Errorf("velocities not two distinct ascending: %v", m.Velocity)
	}
	lk := matMul(l, k)
	lkt := matMul(lk, m.T)
	for i := 0; i < 2; i++ {
		for c := 0; c < 2; c++ {
			want := m.T[i][c] / (m.Velocity[c] * m.Velocity[c])
			if math.Abs(lkt[i][c]-want) > 1e-9*math.Abs(want)+1e-30 {
				t.Errorf("L·K·T[%d][%d] = %g, want %g", i, c, lkt[i][c], want)
			}
		}
	}
	if d := maxDiff(matMul(matMul(m.Zc, k), m.Zc), l); d > 1e-9*l[0][0] {
		t.Errorf("Zc·K·Zc − L = %g", d)
	}
	t.Logf("v = %.4g, %.4g m/s; Zc = %.4g", m.Velocity[0], m.Velocity[1], m.Zc)
}

// Bewley Fig. 6‑2 (pp. 79–82): two lines of surge impedance
// [[500, 125], [125, 500]] meet R₁ = 200, R₂ = 100 and two outgoing lines
// [[400, 100], [100, 400]]. Kron's elimination (Eq. 3‑18) of the 3×3
// network gives Z′ = [[580, 20], [20, 80]]; the reflected voltages are
// about (.10eₐ − .21e_b, −.72e_b) and the transmitted ones
// (.72eₐ − .11e_b, .28e_b). The book carries two digits (it prints −.20,
// 3.56 and 6.40 where full precision gives −.21, 3.57 and 6.34), so the
// check is to full precision against an independent computation and to
// ±0.02 against the print.
func TestBewleyFig62(t *testing.T) {
	z := [][]float64{{600, 100, -100}, {100, 400, -400}, {-100, -400, 500}}
	zr := [][]float64{
		{z[0][0] - z[0][2]*z[2][0]/z[2][2], z[0][1] - z[0][2]*z[2][1]/z[2][2]},
		{z[1][0] - z[1][2]*z[2][0]/z[2][2], z[1][1] - z[1][2]*z[2][1]/z[2][2]},
	}
	if d := maxDiff(zr, [][]float64{{580, 20}, {20, 80}}); d > 1e-12 {
		t.Fatalf("Z′ = %v", zr)
	}
	zc := [][]float64{{500, 125}, {125, 500}}
	r, err := Reflection(zc, zr)
	if err != nil {
		t.Fatal(err)
	}
	full := [][]float64{{0.10180, -0.20646}, {0, -0.72414}}
	printed := [][]float64{{0.10, -0.20}, {0, -0.72}}
	if d := maxDiff(r, full); d > 1e-4 {
		t.Errorf("R = %v, want %v", r, full)
	}
	if d := maxDiff(r, printed); d > 0.02 {
		t.Errorf("R = %v, Bewley prints %v", r, printed)
	}
	// Transmitted: branch currents from the total voltage (1 + R)·e, then
	// the outgoing lines' voltages.
	tot := [][]float64{{1 + r[0][0], r[0][1]}, {r[1][0], 1 + r[1][1]}}
	zi, err := invert(z)
	if err != nil {
		t.Fatal(err)
	}
	ib := matMul(zi, [][]float64{tot[0], tot[1], {0, 0}})
	out := [][]float64{ib[0], {ib[1][0] - ib[2][0], ib[1][1] - ib[2][1]}}
	tr := matMul([][]float64{{400, 100}, {100, 400}}, out)
	want := [][]float64{{0.71862, -0.11069}, {0, 0.27586}}
	if d := maxDiff(tr, want); d > 1e-4 {
		t.Errorf("transmitted %v, want %v", tr, want)
	}
	if d := maxDiff(tr, [][]float64{{0.72, -0.11}, {0, 0.28}}); d > 0.02 {
		t.Errorf("transmitted %v differs from Bewley's print", tr)
	}
}
