// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
)

// Transmission lines from the field solution (CLAUDE.md §3c item 2, M9i).
// The telegrapher's equations ∂v/∂z = −(R + L∂ₜ)i, ∂i/∂z = −(G + C∂ₜ)v
// follow from Maxwell's equations on a uniform line in the TEM mode (Paul
// 2008, Ch. 3); their per-unit-length L and C come from a Section solve.
// In the circuit a line is a ladder of n sections — each a series R·Δz,
// L·Δz and a shunt C·Δz, G·Δz (a symmetric T: half the series element at
// each end of the section) — accurate while a section is short against the
// shortest wavelength (Δz < λ/20 keeps the phase error below about 1 %).

// TLine is a uniform line: per-unit-length R [Ω/m], L [H/m], G [S/m],
// C [F/m], and its length [m].
type TLine struct {
	R, L, G, C, Length float64
}

// FromSection makes a line from conductor k's self terms of a
// single-signal section result (R and G given separately).
func FromSection(s *SectionResult, k int, length, r, g float64) TLine {
	return TLine{R: r, L: s.L[k][k], G: g, C: s.C[k][k], Length: length}
}

// Z0 is the characteristic impedance √((R + jωL)/(G + jωC)) at ω (ω = 0
// with R = G = 0 gives the lossless √(L/C)).
func (t TLine) Z0(omega float64) complex128 {
	if omega == 0 && t.R == 0 && t.G == 0 {
		return complex(math.Sqrt(t.L/t.C), 0)
	}
	return cmplx.Sqrt(complex(t.R, omega*t.L) / complex(t.G, omega*t.C))
}

// Delay is the lossless one-way delay Length·√(LC) [s].
func (t TLine) Delay() float64 { return t.Length * math.Sqrt(t.L*t.C) }

// Ladder adds the line to a builder as n T-sections between nodes in and
// out, the return conductor being ref. Element names start with prefix.
func (t TLine) Ladder(b *circuit.Builder, prefix, in, out, ref string, n int) error {
	if n < 1 || !(t.Length > 0) || !(t.L > 0) || !(t.C > 0) {
		return fmt.Errorf("em: Ladder: %d sections of %+v", n, t)
	}
	dz := t.Length / float64(n)
	node := func(k int) string { return fmt.Sprintf("%s_n%d", prefix, k) }
	prev := in
	for k := 0; k < n; k++ {
		mid := fmt.Sprintf("%s_m%d", prefix, k)
		next := out
		if k < n-1 {
			next = node(k + 1)
		}
		half := func(tag, a, z string) {
			if t.R > 0 {
				r := fmt.Sprintf("%s_r%s%d", prefix, tag, k)
				b.Resistor("R"+r, a, r, t.R*dz/2)
				a = r
			}
			b.Inductor(fmt.Sprintf("L%s_%s%d", prefix, tag, k), a, z, t.L*dz/2)
		}
		half("a", prev, mid)
		b.Capacitor(fmt.Sprintf("C%s_%d", prefix, k), mid, ref, t.C*dz)
		if t.G > 0 {
			b.Resistor(fmt.Sprintf("G%s_%d", prefix, k), mid, ref, 1/(t.G*dz))
		}
		half("b", mid, next)
		prev = next
	}
	return nil
}

// CoupledLadder adds m coupled lossless lines (inductance matrix L and
// Maxwell capacitance matrix C per unit length, from a Section) as n
// sections: each line a T-section ladder like Ladder, the series inductors
// of the same half-section coupled by K = L_ij/√(L_ii·L_jj), and the
// mutual capacitance −C_ij between the lines' mid nodes, each line's
// capacitance to ref being its row sum. in[i] and out[i] are line i's
// ends.
func CoupledLadder(b *circuit.Builder, prefix string, in, out []string, ref string, L, C [][]float64, length float64, n int) error {
	m := len(in)
	if len(out) != m || len(L) != m || len(C) != m || n < 1 {
		return fmt.Errorf("em: CoupledLadder: %d lines, %d×%d matrices", m, len(L), len(C))
	}
	dz := length / float64(n)
	prev := append([]string(nil), in...)
	for k := 0; k < n; k++ {
		mids := make([]string, m)
		nexts := make([]string, m)
		for i := 0; i < m; i++ {
			mids[i] = fmt.Sprintf("%s%d_m%d", prefix, i, k)
			nexts[i] = out[i]
			if k < n-1 {
				nexts[i] = fmt.Sprintf("%s%d_n%d", prefix, i, k+1)
			}
		}
		for _, side := range []string{"a", "b"} {
			names := make([]string, m)
			for i := 0; i < m; i++ {
				names[i] = fmt.Sprintf("L%s%d_%s%d", prefix, i, side, k)
				a, z := prev[i], mids[i]
				if side == "b" {
					a, z = mids[i], nexts[i]
				}
				b.Inductor(names[i], a, z, L[i][i]*dz/2)
			}
			for i := 0; i < m; i++ {
				for j := i + 1; j < m; j++ {
					if kc := L[i][j] / math.Sqrt(L[i][i]*L[j][j]); kc != 0 {
						b.Mutual(fmt.Sprintf("K%s%d%d_%s%d", prefix, i, j, side, k), names[i], names[j], kc)
					}
				}
			}
		}
		for i := 0; i < m; i++ {
			row := 0.0
			for j := 0; j < m; j++ {
				row += C[i][j]
			}
			b.Capacitor(fmt.Sprintf("C%s%d_g%d", prefix, i, k), mids[i], ref, row*dz)
			for j := i + 1; j < m; j++ {
				if C[i][j] != 0 {
					b.Capacitor(fmt.Sprintf("C%s%d%d_%d", prefix, i, j, k), mids[i], mids[j], -C[i][j]*dz)
				}
			}
		}
		prev = nexts
	}
	return nil
}
