// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"

	"github.com/fdamador/kronspice/circuit"
)

// Kron's non-Riemannian form of the machine equations, as set out by
// B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," Rev. Mod. Phys. 21,
// 535–540 (1949).
//
// Coordinates x^α = (q_1 … q_n, θ): the winding charges and the rotor angle
// (Hoffmann eq. 22). The metric is L_αβ(θ) with the inertia J in the θθ
// slot (eqs. 8, 24); the holonomic equations are Lagrange–Maxwell, eq. 5:
//
//	L_αβ·ẍ^β + [βγ, α]·ẋ^β·ẋ^γ + R_αβ·ẋ^β = F_α
//
// A change to current paths that do not turn with the rotor (Hoffmann's
// "stationary coils" 2̄, 3̄, Fig. 3) is a non-holonomic transformation
// dx^α = C_λ^α(θ)·dx̄^λ (eqs. 11, 35, 37–38). In the new frame the metric
// L̄_λμ = C_λ^α·C_μ^β·L_αβ (eq. 15) can be constant, but the connection is
// no longer the Christoffel symbols of L̄: it gains the torsion
//
//	Ω_λμν = ½·L_αγ·C_ν^γ·(∂C_λ^α/∂x̄^μ − ∂C_μ^α/∂x̄^λ)       (eq. 19)
//
// and the equations become (eq. 5′)
//
//	L̄_λν·ẍ^ν + ⟨[μν, λ]⟩·ẋ^μ·ẋ^ν + 2·Ω_μλν·ẋ^μ·ẋ^ν + R̄_λμ·ẋ^μ = F̄_λ
//
// "If we used these values of T̄ and R̄ in the original form of the
// Lagrange–Maxwell equations, we would obtain wrong results" (p. 540): the
// 2Ω term carries the speed voltages and the torque. In kronspice's
// MovingFrame, C_λ^α is C[α][λ] for the windings and 1 for θ.

// frameC returns Hoffmann's C_λ^α (index order [α][λ]) and its θ-derivative
// for the full coordinate set, windings plus θ.
func frameC(f MovingFrame, n int, theta float64) (c, dc [][]float64, err error) {
	ce, de := f.C(theta), f.D(theta)
	if len(ce) != n || len(de) != n {
		return nil, nil, fmt.Errorf("%w: frame of size %d for %d windings", ErrDimension, len(ce), n)
	}
	N := n + 1
	c, dc = square(N), square(N)
	for a := 0; a < n; a++ {
		copy(c[a][:n], ce[a])
		copy(dc[a][:n], de[a])
	}
	c[n][n] = 1
	return c, dc, nil
}

// metric returns L_αβ(θ) and ∂L_αβ/∂θ with the inertia J in the θθ slot.
func metric(s *circuit.MachineSpec, J, theta float64) (g, dg [][]float64) {
	n := s.N
	l, d1 := square(n), square(n)
	s.Inductance(theta, l, d1, nil)
	g, dg = square(n+1), square(n+1)
	for a := 0; a < n; a++ {
		copy(g[a][:n], l[a])
		copy(dg[a][:n], d1[a])
	}
	g[n][n] = J
	return g, dg
}

// Torsion returns Hoffmann's Ω_λμν (eq. 19) of the frame f at angle θ,
// index order [λ][μ][ν]. Only x̄^θ = θ appears in C, so only ∂/∂θ terms
// survive. Ω is antisymmetric in λ, μ.
func Torsion(s *circuit.MachineSpec, J float64, f MovingFrame, theta float64) ([][][]float64, error) {
	n, N := s.N, s.N+1
	c, dc, err := frameC(f, n, theta)
	if err != nil {
		return nil, err
	}
	g, _ := metric(s, J, theta)
	// h[α][ν] = L_αγ·C_ν^γ
	h := square(N)
	for a := 0; a < N; a++ {
		for nu := 0; nu < N; nu++ {
			for gm := 0; gm < N; gm++ {
				h[a][nu] += g[a][gm] * c[gm][nu]
			}
		}
	}
	out := make([][][]float64, N)
	for la := 0; la < N; la++ {
		out[la] = square(N)
		for mu := 0; mu < N; mu++ {
			for nu := 0; nu < N; nu++ {
				v := 0.0
				for a := 0; a < N; a++ {
					d := 0.0
					if mu == n {
						d += dc[a][la]
					}
					if la == n {
						d -= dc[a][mu]
					}
					v += h[a][nu] * d
				}
				out[la][mu][nu] = 0.5 * v
			}
		}
	}
	return out, nil
}

// NewMetric returns L̄_λμ = C_λ^α·C_μ^β·L_αβ (eq. 15) and its θ-derivative.
func NewMetric(s *circuit.MachineSpec, J float64, f MovingFrame, theta float64) (lb, dlb [][]float64, err error) {
	c, dc, err := frameC(f, s.N, theta)
	if err != nil {
		return nil, nil, err
	}
	g, dg := metric(s, J, theta)
	ct := transpose(c)
	lb = mul(ct, mul(g, c))
	dlb = add(add(mul(transpose(dc), mul(g, c)), mul(ct, mul(dg, c))), mul(ct, mul(g, dc)))
	return lb, dlb, nil
}

// NonRiemannian evaluates the left side of Hoffmann's eq. 5′ without
// resistance in frame f: the generalized forces F̄_λ (winding voltages in
// the new frame, then the shaft torque) for new-frame currents i, their
// derivatives di, speed ω and acceleration dω. With torsion false the
// 2Ω term is left out — the "wrong results" of ordinary Lagrange–Maxwell
// equations applied to the new metric.
func NonRiemannian(s *circuit.MachineSpec, J float64, f MovingFrame, theta float64, i, di []float64, omega, domega float64, torsion bool) ([]float64, error) {
	n, N := s.N, s.N+1
	if len(i) != n || len(di) != n {
		return nil, fmt.Errorf("%w: %d currents for %d windings", ErrDimension, len(i), n)
	}
	lb, dlb, err := NewMetric(s, J, f, theta)
	if err != nil {
		return nil, err
	}
	om, err := Torsion(s, J, f, theta)
	if err != nil {
		return nil, err
	}
	v := append(append([]float64(nil), i...), omega)
	a := append(append([]float64(nil), di...), domega)
	// Christoffel symbols of the new metric, ⟨[μν, λ]⟩: only θ-derivatives.
	dg := func(p, q, r int) float64 {
		if r != n {
			return 0
		}
		return dlb[p][q]
	}
	out := make([]float64, N)
	for la := 0; la < N; la++ {
		for nu := 0; nu < N; nu++ {
			out[la] += lb[la][nu] * a[nu]
		}
		for mu := 0; mu < N; mu++ {
			for nu := 0; nu < N; nu++ {
				chr := 0.5 * (dg(mu, la, nu) + dg(nu, la, mu) - dg(mu, nu, la))
				out[la] += chr * v[mu] * v[nu]
				if torsion {
					out[la] += 2 * om[mu][la][nu] * v[mu] * v[nu]
				}
			}
		}
	}
	return out, nil
}
