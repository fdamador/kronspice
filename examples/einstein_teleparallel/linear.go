// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
)

// Einstein's first approximation (§2): h_μα = δ_μα + k_μα with x⁴ imaginary,
// so all index sums are Euclidean. For a plane wave k_μα = ε_μα·e^{p·x}
// every ∂_μ becomes p_μ (complex p allows waves), and his field equation
// (5) becomes a 16×16 linear system M(p)·ε = 0:
//
//	p²ε_βα − p_μp_β ε_μα + p_μp_β ε_αμ − p_μp_α ε_βμ = 0          (5)
//
// with g_αβ − δ_αβ = ε_αβ + ε_βα (6) and 2φ_α = p_μ·ε_αμ − p_α·ε_μμ (2a).
// His claims, checked on every solution:
//
//	p·φ = 0                                                          (8)
//	p²·φ = 0                                                         (9)
//	½(−p²ḡ_αβ + p_μp_β ḡ_μα + p_μp_α ḡ_μβ − p_αp_β ḡ_μμ) = p_βφ_α + p_αφ_β   (7)
//
// and the left side of (7) is the linearized Ricci tensor, so with φ = 0
// it is R_αβ = 0, vacuum gravitation.

func idx(m, a int) int { return 4*m + a }

// fieldMatrix is M(p) acting on ε (16 components, ε[4μ+α] = ε_μα); row
// (β, α) is equation (5) for that index pair.
func fieldMatrix(p [4]complex128) [16][16]complex128 {
	var M [16][16]complex128
	p2 := complex(0, 0)
	for _, v := range p {
		p2 += v * v
	}
	for b := 0; b < 4; b++ {
		for a := 0; a < 4; a++ {
			row := idx(b, a)
			M[row][idx(b, a)] += p2
			for m := 0; m < 4; m++ {
				M[row][idx(m, a)] -= p[m] * p[b]
				M[row][idx(a, m)] += p[m] * p[b]
				M[row][idx(b, m)] -= p[m] * p[a]
			}
		}
	}
	return M
}

// nullSpace returns a basis of the null space of M (complex row
// reduction with partial pivoting; tol relative).
func nullSpace(M [16][16]complex128) [][16]complex128 {
	A := M
	piv := make([]int, 0, 16)
	r := 0
	scale := 0.0
	for i := range A {
		for j := range A[i] {
			scale = math.Max(scale, cmplx.Abs(A[i][j]))
		}
	}
	tol := 1e-10 * scale
	for c := 0; c < 16 && r < 16; c++ {
		p := r
		for i := r + 1; i < 16; i++ {
			if cmplx.Abs(A[i][c]) > cmplx.Abs(A[p][c]) {
				p = i
			}
		}
		if cmplx.Abs(A[p][c]) < tol {
			continue
		}
		A[r], A[p] = A[p], A[r]
		f := A[r][c]
		for k := range A[r] {
			A[r][k] /= f
		}
		for i := 0; i < 16; i++ {
			if i != r && A[i][c] != 0 {
				g := A[i][c]
				for k := range A[i] {
					A[i][k] -= g * A[r][k]
				}
			}
		}
		piv = append(piv, c)
		r++
	}
	isPiv := map[int]int{}
	for i, c := range piv {
		isPiv[c] = i
	}
	var out [][16]complex128
	for free := 0; free < 16; free++ {
		if _, ok := isPiv[free]; ok {
			continue
		}
		var v [16]complex128
		v[free] = 1
		for c, i := range isPiv {
			v[c] = -A[i][free]
		}
		out = append(out, v)
	}
	return out
}

// LinearCheck solves (5) for wave vector p and returns the number of
// independent solutions and the largest violations of (8), (9) and (7),
// and, among solutions with φ = 0, of R_αβ = 0.
type LinearCheck struct {
	Solutions                    int
	Eq8, Eq9, Eq7, Ricci, MaxPhi float64
}

// CheckWave runs the check for one wave vector.
func CheckWave(p [4]complex128) LinearCheck {
	sols := nullSpace(fieldMatrix(p))
	lc := LinearCheck{Solutions: len(sols)}
	p2 := complex(0, 0)
	for _, v := range p {
		p2 += v * v
	}
	for _, e := range sols {
		eps := func(m, a int) complex128 { return e[idx(m, a)] }
		var phi [4]complex128
		norm := 0.0
		for i := range e {
			norm = math.Max(norm, cmplx.Abs(e[i]))
		}
		for a := 0; a < 4; a++ { // (2a): 2φ_α = ∂k_αμ/∂x_μ − ∂k_μμ/∂x_α
			for m := 0; m < 4; m++ {
				phi[a] += 0.5 * (p[m]*eps(a, m) - p[a]*eps(m, m))
			}
		}
		g := func(a, b int) complex128 { return eps(a, b) + eps(b, a) }
		div := complex(0, 0)
		pn := 0.0
		for a := 0; a < 4; a++ {
			div += p[a] * phi[a]
			lc.Eq9 = math.Max(lc.Eq9, cmplx.Abs(p2*phi[a])/norm)
			pn = math.Max(pn, cmplx.Abs(phi[a])/norm)
		}
		lc.MaxPhi = math.Max(lc.MaxPhi, pn)
		lc.Eq8 = math.Max(lc.Eq8, cmplx.Abs(div)/norm)
		tr := complex(0, 0)
		for m := 0; m < 4; m++ {
			tr += g(m, m)
		}
		for a := 0; a < 4; a++ {
			for b := 0; b < 4; b++ {
				lhs := -p2 * g(a, b)
				for m := 0; m < 4; m++ {
					lhs += p[m]*p[b]*g(m, a) + p[m]*p[a]*g(m, b)
				}
				lhs = 0.5 * (lhs - p[a]*p[b]*tr)
				rhs := p[b]*phi[a] + p[a]*phi[b]
				lc.Eq7 = math.Max(lc.Eq7, cmplx.Abs(lhs-rhs)/norm)
				if pn < 1e-9 {
					lc.Ricci = math.Max(lc.Ricci, cmplx.Abs(lhs)/norm)
				}
			}
		}
	}
	return lc
}
