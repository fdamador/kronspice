// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math/rand"
	"testing"
)

func randSym(r *rand.Rand, m int) [][]complex128 {
	z := make([][]complex128, m)
	for i := range z {
		z[i] = make([]complex128, m)
	}
	for i := 0; i < m; i++ {
		for j := i; j < m; j++ {
			v := complex(r.Float64(), r.Float64())
			if i == j {
				v += 5
			}
			z[i][j], z[j][i] = v, v
		}
	}
	return z
}

// Bewley Eq. 5‑9 (built from pairwise leakages only) equals Kron's
// transformation of the full primitive impedance with the magnetizing
// current eliminated, C = [−n₁ … −n_n; I] (Eqs. 5‑13 to 5‑16).
func TestLeakageMatrixEq59(t *testing.T) {
	r := rand.New(rand.NewSource(5))
	n := []float64{1, 0.5, 2.3, 1.7}
	m := len(n)
	z := randSym(r, m)
	pair, err := PairLeakage(z, n)
	if err != nil {
		t.Fatal(err)
	}
	got, err := LeakageMatrix(pair, n)
	if err != nil {
		t.Fatal(err)
	}
	zp := NewMat[Primitive, Primitive](m, m)
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			zp.Set(i, j, z[i][j])
		}
	}
	c := NewMat[Primitive, Mesh](m, m-1)
	for k := 1; k < m; k++ {
		c.Set(0, k-1, complex(-n[k], 0))
		c.Set(k, k-1, 1)
	}
	want, err := Transform(c, zp)
	if err != nil {
		t.Fatal(err)
	}
	for j := 0; j < m-1; j++ {
		for k := 0; k < m-1; k++ {
			if !cnear(got[j][k], want.At(j, k), 1e-12) {
				t.Errorf("Z[%d][%d] = %v, Cₜ·z·C gives %v", j+1, k+1, got[j][k], want.At(j, k))
			}
		}
	}
}

// Eq. 5‑21: one primary N₀ and two secondaries N₁, N₂ in series give,
// referred to the load side, N₁/(N₁+N₂)·Z₀₋₁ + N₂/(N₁+N₂)·Z₀₋₂ −
// N₁N₂/(N₁+N₂)²·Z₁₋₂, which is n²·(Z₁₁ + 2Z₁₂ + Z₂₂) with n = N₀/(N₁+N₂).
func TestLeakageSeriesSecondaries(t *testing.T) {
	z01, z02, z12 := 0.1+2i, 0.12+2.5i, 0.05+0.9i
	N0, N1, N2 := 100.0, 40.0, 70.0
	n := []float64{1, N1 / N0, N2 / N0}
	pair := [][]complex128{{0, z01, z02}, {z01, 0, z12}, {z02, z12, 0}}
	zl, err := LeakageMatrix(pair, n)
	if err != nil {
		t.Fatal(err)
	}
	r := N0 / (N1 + N2)
	got := complex(r*r, 0) * (zl[0][0] + 2*zl[0][1] + zl[1][1])
	s := N1 + N2
	want := complex(N1/s, 0)*z01 + complex(N2/s, 0)*z02 - complex(N1*N2/(s*s), 0)*z12
	if !cnear(got, want, 1e-13) {
		t.Errorf("series leakage %v, Eq. 5‑21 gives %v", got, want)
	}
}

// Eq. 5‑42: the forked autotransformer of Fig. 5‑5 has the input
// impedance (Z + Z₂₂) − 2n·Z₁₂ + n²·Z₁₁ for C = (n, −1, 1) on the axes
// (1, 2, load).
func TestForkedAutotransformer(t *testing.T) {
	r := rand.New(rand.NewSource(9))
	zz := randSym(r, 2)
	zl := 0.3 + 2i
	nr := 1.7
	zp := NewMat[Primitive, Primitive](3, 3)
	zp.Set(0, 0, zz[0][0])
	zp.Set(0, 1, zz[0][1])
	zp.Set(1, 0, zz[1][0])
	zp.Set(1, 1, zz[1][1])
	zp.Set(2, 2, zl)
	c := NewMat[Primitive, Mesh](3, 1)
	c.Set(0, 0, complex(nr, 0))
	c.Set(1, 0, -1)
	c.Set(2, 0, 1)
	got, err := Transform(c, zp)
	if err != nil {
		t.Fatal(err)
	}
	n := complex(nr, 0)
	want := (zl + zz[1][1]) - 2*n*zz[0][1] + n*n*zz[0][0]
	if !cnear(got.At(0, 0), want, 1e-13) {
		t.Errorf("Z = %v, Eq. 5‑42 gives %v", got.At(0, 0), want)
	}
}
