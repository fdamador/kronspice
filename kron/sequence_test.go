// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math/cmplx"
	"testing"
)

func cnear(a, b complex128, tol float64) bool {
	return cmplx.Abs(a-b) <= tol*(1+cmplx.Abs(b))
}

// Bewley 1961 Eq. 4‑17 is unitary, and the sequence transform round-trips.
func TestFortescueUnitary(t *testing.T) {
	c := FortescueMatrix()
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			var s complex128
			for k := 0; k < 3; k++ {
				s += cmplx.Conj(c[k][i]) * c[k][j]
			}
			want := complex(0, 0)
			if i == j {
				want = 1
			}
			if !cnear(s, want, 1e-14) {
				t.Errorf("(C*ₜC)[%d][%d] = %v", i, j, s)
			}
		}
	}
	x := [3]complex128{1 + 2i, -0.3 + 0.1i, 4}
	y := FromSequence(ToSequence(x))
	for k := range x {
		if !cnear(y[k], x[k], 1e-14) {
			t.Errorf("round trip %d: %v, want %v", k, y[k], x[k])
		}
	}
	// Power invariance: e*·i is the same in both frames.
	e := [3]complex128{2 - 1i, 0.5i, -1 + 3i}
	var p, q complex128
	es, is := ToSequence(e), ToSequence(x)
	for k := 0; k < 3; k++ {
		p += cmplx.Conj(e[k]) * x[k]
		q += cmplx.Conj(es[k]) * is[k]
	}
	if !cnear(p, q, 1e-13) {
		t.Errorf("power %v in phase frame, %v in sequence frame", p, q)
	}
}

// Eq. 4‑24: a cyclic-symmetric impedance matrix is diagonal on sequence
// axes, with entries Z+Z′+Z″, Z+a²Z′+aZ″, Z+aZ′+a²Z″.
func TestSequenceCyclic(t *testing.T) {
	zs, zm1, zm2 := 2+3i, 0.5+1i, 0.2+0.4i
	z := [3][3]complex128{{zs, zm1, zm2}, {zm2, zs, zm1}, {zm1, zm2, zs}}
	a := SeqOperator
	want := [3]complex128{zs + zm1 + zm2, zs + a*a*zm1 + a*zm2, zs + a*zm1 + a*a*zm2}
	got := SequenceImpedance(z)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			w := complex(0, 0)
			if i == j {
				w = want[i]
			}
			if !cnear(got[i][j], w, 1e-13) {
				t.Errorf("Z′[%d][%d] = %v, want %v", i, j, got[i][j], w)
			}
		}
	}
}

// Eq. 4‑26: an unbalanced wye load diag(Zₐ, Z_b, Z_c) has first column
// (Zₐ+Z_b+Z_c)/3, (Zₐ+aZ_b+a²Z_c)/3, (Zₐ+a²Z_b+aZ_c)/3 on sequence axes.
func TestSequenceUnbalanced(t *testing.T) {
	za, zb, zc := 1+1i, 2+0.5i, complex(3, 0)
	got := SequenceImpedance([3][3]complex128{{za, 0, 0}, {0, zb, 0}, {0, 0, zc}})
	a := SeqOperator
	want := [3]complex128{(za + zb + zc) / 3, (za + a*zb + a*a*zc) / 3, (za + a*a*zb + a*zc) / 3}
	for i := 0; i < 3; i++ {
		if !cnear(got[i][0], want[i], 1e-13) {
			t.Errorf("Z′[%d][0] = %v, want %v", i, got[i][0], want[i])
		}
	}
}

// Bewley Fig. 4‑2, Eq. 4‑35: a grounded generator (Z₀, Z₁, Z₂) feeding a
// load Z per phase with neutral impedance Z_g reduces to the classic
// sequence network [[Z₀+Z₁+2Z+3Z_g, Z₀+Z+3Z_g], [Z₀+Z+3Z_g, Z₀+Z₂+2Z+3Z_g]].
func TestSequenceFig42(t *testing.T) {
	z0, z1, z2, zl, zg := 0.1+0.5i, 0.2+1i, 0.15+0.9i, 5+2i, 1+0.3i
	z := NewMat[Primitive, Primitive](7, 7)
	z.Set(0, 0, 3*z0)
	z.Set(1, 1, 3*z1)
	z.Set(2, 2, 3*z2)
	tri := [3][3]complex128{{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			z.Set(3+i, 3+j, zl*tri[i][j])
		}
	}
	z.Set(6, 6, 9*zg)
	for i := 0; i < 7; i++ {
		for j := 0; j < 7; j++ {
			z.Set(i, j, z.At(i, j)/3)
		}
	}
	c := NewMat[Primitive, Mesh](7, 2)
	for i, row := range [][2]complex128{{-1, -1}, {1, 0}, {0, 1}, {-1, -1}, {1, 0}, {0, 1}, {-1, -1}} {
		c.Set(i, 0, row[0])
		c.Set(i, 1, row[1])
	}
	got, err := Transform(c, z)
	if err != nil {
		t.Fatal(err)
	}
	want := [2][2]complex128{{z0 + z1 + 2*zl + 3*zg, z0 + zl + 3*zg}, {z0 + zl + 3*zg, z0 + z2 + 2*zl + 3*zg}}
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			if !cnear(got.At(i, j), want[i][j], 1e-13) {
				t.Errorf("Eq. 4‑35 [%d][%d] = %v, want %v", i, j, got.At(i, j), want[i][j])
			}
		}
	}
}

// Bewley Problem 4‑4, line-to-ground fault: phase a of a generator with
// sequence impedances (Z₀, Z₁, Z₂) shorted to ground, b and c open. Solved
// in the phase frame (Z = C·Z′·C*ₜ), the fault current is the classic
// 3E/(Z₀+Z₁+Z₂), and its sequence components are equal.
func TestSequenceLineToGround(t *testing.T) {
	z0, z1, z2 := 0.05+0.3i, 0.02+1.2i, 0.03+0.9i
	e := complex(1, 0)
	zp := PhaseImpedance([3][3]complex128{{z0, 0, 0}, {0, z1, 0}, {0, 0, z2}})
	// V = E − Z·I with I_b = I_c = 0 and V_a = 0 ⇒ I_a = E_a / Z_aa.
	ia := e / zp[0][0]
	want := 3 * e / (z0 + z1 + z2)
	if !cnear(ia, want, 1e-13) {
		t.Errorf("I_a = %v, want 3E/(Z₀+Z₁+Z₂) = %v", ia, want)
	}
	s := ToSequence([3]complex128{ia, 0, 0})
	for k := 1; k < 3; k++ {
		if !cnear(s[k], s[0], 1e-13) {
			t.Errorf("sequence currents %v not equal", s)
		}
	}
}
