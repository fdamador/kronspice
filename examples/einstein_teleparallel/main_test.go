// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/kron"
)

// TestSameTorsion: Einstein's torsion (in frame components) is minus the
// frame's anholonomy, and Hoffmann's Ω is it lowered by L̄.
func TestSameTorsion(t *testing.T) {
	for _, f := range []Frame{{}, {S: 0.3}} {
		for _, th := range []float64{0.1, 0.7, 2.5} {
			A := f.Anholonomy(th)
			if d := maxDiff(A, f.ToFrame(f.EinsteinTorsion(th), th), -1); d > 1e-9 {
				t.Errorf("s %g θ %g: Einstein vs anholonomy %g", f.S, th, d)
			}
			om, err := kron.Torsion(Motor(), J, f.Moving(), th)
			if err != nil {
				t.Fatal(err)
			}
			lb := f.NewMetric(th)
			for l := 0; l < 4; l++ {
				for m := 0; m < 4; m++ {
					for n := 0; n < 4; n++ {
						v := 0.0
						for r := 0; r < 4; r++ {
							v += lb[n][r] * A[r][l][m]
						}
						if math.Abs(om[l][m][n]-v) > 1e-12 {
							t.Fatalf("Ω ≠ L̄Λ at %d%d%d", l, m, n)
						}
					}
				}
			}
			if d := TorsionOfW(f, th); d > 1e-8 {
				t.Errorf("Weitzenböck torsion off by %g", d)
			}
		}
	}
}

// TestPhiAndCurvature: Kron's frame has φ = 0, the scaled frame φ_θ = s;
// Einstein's connection is flat, the machine's is not.
func TestPhiAndCurvature(t *testing.T) {
	if p := Phi(Frame{}.EinsteinTorsion(0.7)); math.Abs(p[3]) > 1e-9 {
		t.Errorf("φ = %v", p)
	}
	if p := Phi(Frame{S: 0.3}.EinsteinTorsion(0.7)); math.Abs(p[3]-0.3) > 1e-8 {
		t.Errorf("scaled φ = %v", p)
	}
	ce, ck := Curvatures(Frame{}, 0.7)
	if ce > 1e-9 || ck < 1 {
		t.Errorf("curvatures %g, %g", ce, ck)
	}
}

// TestLinear: Einstein's eqs. 7–9 follow from eq. 5.
func TestLinear(t *testing.T) {
	for _, p := range [][4]complex128{{1, 0.3, -0.2, 0.5}, {0.4, 0.2 + 0.3i, -0.5, 0.1 - 0.7i}, {1, 1i, 0, 0}} {
		lc := CheckWave(p)
		if lc.Eq7 > 1e-12 || lc.Eq8 > 1e-12 || lc.Eq9 > 1e-12 || lc.Ricci > 1e-12 {
			t.Errorf("p %v: %+v", p, lc)
		}
	}
	if lc := CheckWave([4]complex128{1, 1i, 0, 0}); lc.MaxPhi < 0.1 || lc.Solutions != 10 {
		t.Errorf("light-like: %+v", lc)
	}
	if lc := CheckWave([4]complex128{1, 0.3, -0.2, 0.5}); lc.MaxPhi > 1e-12 || lc.Solutions != 4 {
		t.Errorf("non-null: %+v", lc)
	}
}

// TestMotor: the motor built from Einstein's torsion equals Kron's frame
// II(b); without the torsion it only coasts.
func TestMotor(t *testing.T) {
	_, s2, s3, s4, err := Specs()
	if err != nil {
		t.Fatal(err)
	}
	for i := range s2.G {
		for j := range s2.G[i] {
			if math.Abs(s2.G[i][j]-s3.G[i][j]) > 1e-12 || math.Abs(s2.L0[i][j]-s3.L0[i][j]) > 1e-12 {
				t.Fatalf("G or L̄ differ at %d,%d", i, j)
			}
		}
	}
	r4, err := Start("no torsion", s4, 40, 1)
	if err != nil {
		t.Fatal(err)
	}
	if want := 40 * math.Exp(-1e-4/J); math.Abs(r4.Final-want) > 1e-3*want {
		t.Errorf("no torsion: %g, want %g", r4.Final, want)
	}
}
