// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
)

// Kron's single-phase motor (Hoffmann 1949, eq. 42; examples/hoffmann_motor)
// has four coordinates x^α = (q1, q2, q3, θ) — three winding charges and
// the rotor angle — so the change to stationary rotor coils (Hoffmann
// eq. 35) is a field of four frame vectors, a 4-bein, Einstein's case n = 4.
//
// Notation: C^α_λ(θ) (C[α][λ]) gives the old coordinates' differentials
// from the frame's, dx^α = C^α_λ dx̄^λ; c = C⁻¹. The frame vectors are
// e_λ = C^α_λ ∂_α; Einstein's legs (his h) are these, his h_{μa} the rows
// of c (A. Einstein, Sitzungsber. Preuss. Akad. Wiss. 1928, 217 and 224;
// A. Unzicker's translation).

// Motor is the frame-I machine (values of examples/hoffmann_motor, with a
// salient pole so the geometry is not trivial).
func Motor() *circuit.MachineSpec {
	s := circuit.NewMachineSpec(3)
	s.L0[0][0], s.L0[1][1], s.L0[2][2] = 0.25, 0.25, 0.25
	h1 := s.Harmonic(1)
	h1.Cos[0][1], h1.Cos[1][0] = 0.23, 0.23
	h1.Sin[0][2], h1.Sin[2][0] = -0.23, -0.23
	h2 := s.Harmonic(2)
	h2.Cos[1][1], h2.Cos[2][2] = -0.02, 0.02
	h2.Sin[1][2], h2.Sin[2][1] = 0.02, 0.02
	return s
}

// J is the rotor inertia [kg·m²].
const J = 1e-3

// Frame is Hoffmann's change to stationary rotor coils, optionally scaled
// by e^{sθ} on the rotor coils (a frame that is not a pure rotation).
type Frame struct {
	S float64 // scale rate (0 = Hoffmann's rotation)
}

// C returns C(θ) and dC/dθ, 4×4 (windings and θ).
func (f Frame) C(th float64) (c, d [4][4]float64) {
	r := kron.RotorFrame(3, [][2]int{{1, 2}}, 1)
	ce, de := r.C(th), r.D(th)
	k := math.Exp(f.S * th)
	for a := 0; a < 3; a++ {
		for b := 0; b < 3; b++ {
			sc, ds := 1.0, 0.0
			if b > 0 { // scale the rotor-coil legs
				sc, ds = k, f.S*k
			}
			c[a][b] = ce[a][b] * sc
			d[a][b] = de[a][b]*sc + ce[a][b]*ds
		}
	}
	c[3][3] = 1
	return c, d
}

// Moving returns the frame as a kron.MovingFrame (windings only).
func (f Frame) Moving() kron.MovingFrame {
	get := func(th float64, deriv bool) [][]float64 {
		c, d := f.C(th)
		m := c
		if deriv {
			m = d
		}
		out := make([][]float64, 3)
		for a := range out {
			out[a] = []float64{m[a][0], m[a][1], m[a][2]}
		}
		return out
	}
	return kron.MovingFrame{C: func(t float64) [][]float64 { return get(t, false) }, D: func(t float64) [][]float64 { return get(t, true) }}
}

func inv4(a [4][4]float64) [4][4]float64 {
	var m [4][8]float64
	for i := 0; i < 4; i++ {
		copy(m[i][:4], a[i][:])
		m[i][4+i] = 1
	}
	for c := 0; c < 4; c++ {
		p := c
		for r := c + 1; r < 4; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		m[c], m[p] = m[p], m[c]
		f := m[c][c]
		for k := range m[c] {
			m[c][k] /= f
		}
		for r := 0; r < 4; r++ {
			if r != c && m[r][c] != 0 {
				f := m[r][c]
				for k := range m[r] {
					m[r][k] -= f * m[c][k]
				}
			}
		}
	}
	var o [4][4]float64
	for i := 0; i < 4; i++ {
		copy(o[i][:], m[i][4:])
	}
	return o
}

// Anholonomy is the frame's torsion in frame components,
// Λ^ρ_λμ = ½ c^ρ_α (e_μ C^α_λ − e_λ C^α_μ); only θ is a coordinate the
// frame depends on, and e_μ(θ) = δ_μθ, so e_μ = δ_μθ ∂_θ.
func (f Frame) Anholonomy(th float64) [4][4][4]float64 {
	c, d := f.C(th)
	ci := inv4(c)
	var L [4][4][4]float64
	for r := 0; r < 4; r++ {
		for l := 0; l < 4; l++ {
			for m := 0; m < 4; m++ {
				s := 0.0
				for a := 0; a < 4; a++ {
					v := 0.0
					if m == 3 {
						v += d[a][l]
					}
					if l == 3 {
						v -= d[a][m]
					}
					s += ci[r][a] * v
				}
				L[r][l][m] = 0.5 * s
			}
		}
	}
	return L
}

// EinsteinTorsion is Einstein's Λ^α_μν = ½ h_a^α (∂h_μa/∂x^ν − ∂h_νa/∂x^μ)
// (his eqs. 9–10 of the June 7, 1928 note) in coordinate components, with
// h_a^α = C^α_a and h_μa = c^a_μ; only ∂/∂θ is non-zero.
func (f Frame) EinsteinTorsion(th float64) [4][4][4]float64 {
	c, _ := f.C(th)
	h := 1e-6
	cp, _ := f.C(th + h)
	cm, _ := f.C(th - h)
	ci, cip, cim := inv4(c), inv4(cp), inv4(cm)
	dci := func(a, m int) float64 { return (cip[a][m] - cim[a][m]) / (2 * h) } // ∂_θ c^a_μ
	_ = ci
	var L [4][4][4]float64
	for al := 0; al < 4; al++ {
		for mu := 0; mu < 4; mu++ {
			for nu := 0; nu < 4; nu++ {
				s := 0.0
				for a := 0; a < 4; a++ {
					v := 0.0
					if nu == 3 {
						v += dci(a, mu)
					}
					if mu == 3 {
						v -= dci(a, nu)
					}
					s += c[al][a] * v
				}
				L[al][mu][nu] = 0.5 * s
			}
		}
	}
	return L
}

// Phi is Einstein's φ_μ = Λ^α_μα (his eq. 2), the trace that he reads as
// the electromagnetic potential; zero means a "pure gravitational" field.
func Phi(L [4][4][4]float64) [4]float64 {
	var p [4]float64
	for mu := 0; mu < 4; mu++ {
		for a := 0; a < 4; a++ {
			p[mu] += L[a][mu][a]
		}
	}
	return p
}

// ToFrame transforms a coordinate tensor T^α_μν to frame components
// c^ρ_α T^α_μν C^μ_λ C^ν_κ.
func (f Frame) ToFrame(T [4][4][4]float64, th float64) [4][4][4]float64 {
	c, _ := f.C(th)
	ci := inv4(c)
	var o [4][4][4]float64
	for r := 0; r < 4; r++ {
		for l := 0; l < 4; l++ {
			for k := 0; k < 4; k++ {
				s := 0.0
				for a := 0; a < 4; a++ {
					for m := 0; m < 4; m++ {
						for n := 0; n < 4; n++ {
							s += ci[r][a] * T[a][m][n] * c[m][l] * c[n][k]
						}
					}
				}
				o[r][l][k] = s
			}
		}
	}
	return o
}

// NewMetric is L̄ = Cᵀ·L(θ)·C (Hoffmann eq. 15), 4×4 with J.
func (f Frame) NewMetric(th float64) [4][4]float64 {
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
	c, _ := f.C(th)
	var o [4][4]float64
	for a := 0; a < 4; a++ {
		for b := 0; b < 4; b++ {
			for i := 0; i < 4; i++ {
				for j := 0; j < 4; j++ {
					o[a][b] += c[i][a] * g[i][j] * c[j][b]
				}
			}
		}
	}
	return o
}
