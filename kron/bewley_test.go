// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math/rand"
	"testing"
)

// Worked examples from L. V. Bewley, Tensor Analysis of Electric Circuits
// and Machines, Ronald Press 1961 ("Bewley 1961"), reproduced with Kron's
// transformation law (Transform). Photographs of the pages are in
// reference/bewley_tensor_analysis/ (not tracked).

func matFrom[R, C Frame](rows [][]complex128) Mat[R, C] {
	m := NewMat[R, C](len(rows), len(rows[0]))
	for i, r := range rows {
		for j, v := range r {
			m.Set(i, j, v)
		}
	}
	return m
}

// Fig. 3‑1 (pp. 32–34): Z = [[4,2,0],[2,5,1],[0,1,6]], C = [[1,1],[1,0],
// [0,−1]] give Z′ = [[13,5],[5,10]] (det 105); with E = 105 on both new
// axes the new currents are (5, 8), the coil currents (13, 5, −8) and the
// coil voltages (62, 43, −43). All exact.
func TestBewleyFig31(t *testing.T) {
	z := matFrom[Primitive, Primitive]([][]complex128{{4, 2, 0}, {2, 5, 1}, {0, 1, 6}})
	c := matFrom[Primitive, Mesh]([][]complex128{{1, 1}, {1, 0}, {0, -1}})
	zp, err := Transform(c, z)
	if err != nil {
		t.Fatal(err)
	}
	want := [2][2]complex128{{13, 5}, {5, 10}}
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if zp.At(i, j) != want[i][j] {
				t.Errorf("Z′[%d][%d] = %v, want %v", i, j, zp.At(i, j), want[i][j])
			}
		}
	}
	det := zp.At(0, 0)*zp.At(1, 1) - zp.At(0, 1)*zp.At(1, 0)
	if det != 105 {
		t.Errorf("det Z′ = %v, want 105", det)
	}
	ip := NewVec[Mesh](2)
	ip.Set(0, (105*zp.At(1, 1)-105*zp.At(0, 1))/det)
	ip.Set(1, (105*zp.At(0, 0)-105*zp.At(1, 0))/det)
	i, err := MulVec(c, ip)
	if err != nil {
		t.Fatal(err)
	}
	e, err := MulVec(z, i)
	if err != nil {
		t.Fatal(err)
	}
	for k, w := range []complex128{5, 8} {
		if ip.At(k) != w {
			t.Errorf("new current %d = %v, want %v", k, ip.At(k), w)
		}
	}
	for k, w := range []complex128{13, 5, -8} {
		if i.At(k) != w {
			t.Errorf("coil current %d = %v, want %v", k, i.At(k), w)
		}
	}
	for k, w := range []complex128{62, 43, -43} {
		if e.At(k) != w {
			t.Errorf("coil voltage %d = %v, want %v", k, e.At(k), w)
		}
	}
}

// Fig. 3‑6 (pp. 40–41): the π-circuit's Z′ = [[12,5,1],[5,10,8],[1,8,11]]
// reduces by partitioned elimination (Eq. 3‑18) to (1/11)·[[131,47],
// [47,46]]. The driving-point impedance is exactly 3817/506 = 7.543 Ω and
// the transfer impedance −3817/517 = −7.383 Ω; the book prints 7.53 and
// −7.37 (a rounding slip in the last digit).
func TestBewleyFig36(t *testing.T) {
	z := matFrom[Primitive, Primitive]([][]complex128{
		{3, 1, 2, 2, 1}, {1, 2, 1, 1, 1}, {2, 1, 4, 3, 3}, {2, 1, 3, 5, 2}, {1, 1, 3, 2, 6}})
	c := matFrom[Primitive, Mesh]([][]complex128{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, -1}, {0, 1, 1}})
	zp, err := Transform(c, z)
	if err != nil {
		t.Fatal(err)
	}
	want := [3][3]complex128{{12, 5, 1}, {5, 10, 8}, {1, 8, 11}}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if zp.At(i, j) != want[i][j] {
				t.Errorf("Z′[%d][%d] = %v, want %v", i, j, zp.At(i, j), want[i][j])
			}
		}
	}
	var r [2][2]complex128
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			r[i][j] = zp.At(i, j) - zp.At(i, 2)*zp.At(2, j)/zp.At(2, 2)
		}
	}
	if !cnear(r[0][0], 131.0/11, 1e-14) || !cnear(r[0][1], 47.0/11, 1e-14) || !cnear(r[1][1], 46.0/11, 1e-14) {
		t.Errorf("reduced Z = %v, want (1/11)·[[131,47],[47,46]]", r)
	}
	det := r[0][0]*r[1][1] - r[0][1]*r[1][0]
	drive := det / r[1][1]     // 1/Y₁₁
	transfer := -det / r[0][1] // 1/Y₂₁
	if !cnear(drive, 3817.0/506, 1e-14) || !cnear(transfer, -3817.0/517, 1e-14) {
		t.Errorf("driving point %v, transfer %v; want 3817/506 and −3817/517", drive, transfer)
	}
}

// Fig. 7‑2 (pp. 85–86), the junction-pair (node) form: six coils, four
// junctions, junction-pair voltages (a′, b′, d′). Y′ = Cᵀ·Y·C gives the
// six printed entries, e.g. Y′a′a′ = Yaa − Yac − Yca + Ycc + Yee + Yff and
// Y′d′d′ = Ycc + Ydd + Yff.
func TestBewleyFig72(t *testing.T) {
	r := rand.New(rand.NewSource(72))
	y := NewMat[Primitive, Primitive](6, 6)
	for k := 0; k < 6; k++ {
		y.Set(k, k, complex(1+r.Float64(), r.Float64()))
	}
	const a, b, cc, d, e, f = 0, 1, 2, 3, 4, 5
	yab, yac := complex(r.Float64(), 0), complex(r.Float64(), 0)
	y.Set(a, b, yab)
	y.Set(b, a, yab)
	y.Set(a, cc, yac)
	y.Set(cc, a, yac)
	conn := matFrom[Primitive, Node]([][]complex128{{1, 0, 0}, {0, 1, 0}, {-1, -1, 1}, {0, 0, -1}, {1, 1, 0}, {1, 0, -1}})
	yp, err := Transform(conn, y)
	if err != nil {
		t.Fatal(err)
	}
	Y := func(i, j int) complex128 { return y.At(i, j) }
	want := map[[2]int]complex128{
		{0, 0}: Y(a, a) - Y(a, cc) - Y(cc, a) + Y(cc, cc) + Y(e, e) + Y(f, f),
		{0, 1}: Y(a, b) - Y(a, cc) + Y(cc, cc) + Y(e, e),
		{0, 2}: Y(a, cc) - Y(cc, cc) - Y(f, f),
		{1, 1}: Y(b, b) + Y(cc, cc) + Y(e, e),
		{1, 2}: -Y(cc, cc),
		{2, 2}: Y(cc, cc) + Y(d, d) + Y(f, f),
	}
	for ij, w := range want {
		if !cnear(yp.At(ij[0], ij[1]), w, 1e-14) {
			t.Errorf("Y′%v = %v, want %v", ij, yp.At(ij[0], ij[1]), w)
		}
	}
}

// Fig. 8‑3 (pp. 97–99), the complete ("all-mesh") network of five coils:
// the non-singular part of the connection tensor gives z′ with, among
// others, z′a′a′ = zaa + zbb + zff, z′b′b′ = zcc + zdd + zff and
// z′a′b′ = zfc + zff (coils c and f mutually coupled).
func TestBewleyFig83(t *testing.T) {
	r := rand.New(rand.NewSource(83))
	const a, b, c, d, f = 0, 1, 2, 3, 4
	z := NewMat[Primitive, Primitive](5, 5)
	for k := 0; k < 5; k++ {
		z.Set(k, k, complex(1+r.Float64(), r.Float64()))
	}
	zcf := complex(r.Float64(), r.Float64())
	z.Set(c, f, zcf)
	z.Set(f, c, zcf)
	// Columns a′, b′, p′, q′, r′ (Bewley p. 98): i_a = i_a′,
	// i_b = i_a′+i_p′+i_q′+i_r′, i_c = i_b′−i_q′, i_d = −i_b′,
	// i_f = i_a′+i_b′+i_p′.
	conn := matFrom[Primitive, Mesh]([][]complex128{
		{1, 0, 0, 0, 0}, {1, 0, 1, 1, 1}, {0, 1, 0, -1, 0}, {0, -1, 0, 0, 0}, {1, 1, 1, 0, 0}})
	zp, err := Transform(conn, z)
	if err != nil {
		t.Fatal(err)
	}
	Z := func(i, j int) complex128 { return z.At(i, j) }
	checks := map[[2]int]complex128{
		{0, 0}: Z(a, a) + Z(b, b) + Z(f, f),
		{1, 1}: Z(c, c) + Z(d, d) + Z(f, f) + 2*Z(c, f),
		{0, 1}: Z(f, c) + Z(f, f),
		{0, 2}: Z(b, b) + Z(f, f),
		{2, 2}: Z(b, b) + Z(f, f),
		{3, 3}: Z(b, b) + Z(c, c),
		{4, 4}: Z(b, b),
	}
	for ij, w := range checks {
		if !cnear(zp.At(ij[0], ij[1]), w, 1e-14) {
			t.Errorf("z′%v = %v, want %v", ij, zp.At(ij[0], ij[1]), w)
		}
	}
}
