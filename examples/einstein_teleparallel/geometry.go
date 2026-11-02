// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Two connections on the same machine space (coordinates q1, q2, q3, θ;
// everything depends on θ only).
//
//   - Einstein's distant-parallelism (Weitzenböck) connection
//     W^α_μν = h_a^α ∂_ν h_μa: the frame legs are parallel everywhere, so
//     its curvature is zero identically and all the geometry is torsion.
//   - Kron's connection: the Christoffel symbols of the machine's metric
//     g(θ) (inductances and inertia), written in the moving frame
//     (Hoffmann eq. 17). A change of frame cannot change a curvature, so it
//     keeps the curvature of g; the torsion Ω only re-expresses it.

type conn [4][4][4]float64 // Γ^α_μν, index [α][μ][ν]

// weitzenbock returns W at θ.
func (f Frame) weitzenbock(th float64) conn {
	c, _ := f.C(th)
	h := 1e-5
	cp, _ := f.C(th + h)
	cm, _ := f.C(th - h)
	cip, cim := inv4(cp), inv4(cm)
	var W conn
	for al := 0; al < 4; al++ {
		for mu := 0; mu < 4; mu++ {
			// ∂_ν nonzero only for ν = θ (index 3)
			s := 0.0
			for a := 0; a < 4; a++ {
				s += c[al][a] * (cip[a][mu] - cim[a][mu]) / (2 * h)
			}
			W[al][mu][3] = s
		}
	}
	return W
}

// metric4 is the machine's metric g(θ) (inductances, inertia J).
func metric4(th float64) [4][4]float64 {
	s := Motor()
	l := make([][]float64, 3)
	for i := range l {
		l[i] = make([]float64, 3)
	}
	s.Inductance(th, l, nil, nil)
	var g [4][4]float64
	for i := 0; i < 3; i++ {
		copy(g[i][:3], l[i])
	}
	g[3][3] = J
	return g
}

// levi returns the Christoffel symbols Γ^α_μν of g at θ.
func levi(th float64) conn {
	h := 1e-5
	g := metric4(th)
	gp, gm := metric4(th+h), metric4(th-h)
	gi := inv4(g)
	dg := func(a, b, c int) float64 { // ∂_c g_ab
		if c != 3 {
			return 0
		}
		return (gp[a][b] - gm[a][b]) / (2 * h)
	}
	var G conn
	for al := 0; al < 4; al++ {
		for mu := 0; mu < 4; mu++ {
			for nu := 0; nu < 4; nu++ {
				s := 0.0
				for e := 0; e < 4; e++ {
					s += gi[al][e] * (dg(e, mu, nu) + dg(e, nu, mu) - dg(mu, nu, e))
				}
				G[al][mu][nu] = 0.5 * s
			}
		}
	}
	return G
}

// curvature returns the largest |R^α_βμν| of a connection given as a
// function of θ: R^α_βμν = ∂_μΓ^α_βν − ∂_νΓ^α_βμ + Γ^α_σμΓ^σ_βν − Γ^α_σνΓ^σ_βμ.
func curvature(G func(float64) conn, th float64) float64 {
	h := 1e-4
	g0, gp, gm := G(th), G(th+h), G(th-h)
	d := func(a, b, n, mu int) float64 { // ∂_μ Γ^a_bn
		if mu != 3 {
			return 0
		}
		return (gp[a][b][n] - gm[a][b][n]) / (2 * h)
	}
	worst := 0.0
	for a := 0; a < 4; a++ {
		for b := 0; b < 4; b++ {
			for mu := 0; mu < 4; mu++ {
				for nu := 0; nu < 4; nu++ {
					r := d(a, b, nu, mu) - d(a, b, mu, nu)
					for s := 0; s < 4; s++ {
						r += g0[a][s][mu]*g0[s][b][nu] - g0[a][s][nu]*g0[s][b][mu]
					}
					worst = math.Max(worst, math.Abs(r))
				}
			}
		}
	}
	return worst
}

// Curvatures returns the largest curvature component of Einstein's
// connection and of the machine's Riemannian connection at θ.
func Curvatures(f Frame, th float64) (einstein, kron float64) {
	return curvature(f.weitzenbock, th), curvature(levi, th)
}

// TorsionOfW checks that the torsion of Einstein's connection,
// W^α_μν − W^α_νμ, is twice his Λ^α_μν (sign as in his definition).
func TorsionOfW(f Frame, th float64) float64 {
	W := f.weitzenbock(th)
	L := f.EinsteinTorsion(th)
	worst := 0.0
	for a := 0; a < 4; a++ {
		for m := 0; m < 4; m++ {
			for n := 0; n < 4; n++ {
				worst = math.Max(worst, math.Abs((W[a][m][n]-W[a][n][m])-2*L[a][m][n]))
			}
		}
	}
	return worst
}
