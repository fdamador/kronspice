// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math/cmplx"

// Hopfion is Rañada's knotted light in the Bateman form of Kedia et al.,
// PRL 111, 150404 (2013), copied from examples/knotted_light (main.go),
// where it is checked against Maxwell's equations and traced: units of L0
// and L0/c, F = E + icB = ∇α × ∇β, P = α∇β with F = ∇ × P.
type Hopfion struct{}

func (Hopfion) parts(x, y, z, t float64) (alpha complex128, gA, gB [3]complex128) {
	r2 := x*x + y*y + z*z
	tm := complex(t, -1)
	D := complex(r2, 0) - tm*tm
	N := complex(r2-t*t-1, 2*z)
	alpha = N / D
	M := complex(2*x, -2*y)
	beta := M / D
	_ = beta
	p := [3]float64{x, y, z}
	for j := 0; j < 3; j++ {
		dD := complex(2*p[j], 0)
		dN := dD
		if j == 2 {
			dN += 2i
		}
		var dM complex128
		switch j {
		case 0:
			dM = 2
		case 1:
			dM = -2i
		}
		gA[j] = dN/D - N*dD/(D*D)
		gB[j] = dM/D - M*dD/(D*D)
	}
	return
}

// F is ∇α × ∇β at (x, y, z, t).
func (h Hopfion) F(x, y, z, t float64) [3]complex128 {
	_, a, b := h.parts(x, y, z, t)
	return [3]complex128{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

// P is the potential α∇β, with F = ∇ × P.
func (h Hopfion) P(x, y, z, t float64) [3]complex128 {
	al, _, b := h.parts(x, y, z, t)
	return [3]complex128{al * b[0], al * b[1], al * b[2]}
}

// Sign is +1 when E = Re F, cB = Im F satisfy Maxwell's equations
// (i∂F/∂t = ∇ × F) and −1 when the conjugate does (then cB = −Im F). It is
// checked numerically at a few points, so the helicity convention of the
// formula cannot be taken on trust.
func (h Hopfion) Sign() float64 {
	const e = 1e-5
	var plus, minus float64
	for _, p := range [][4]float64{{0.3, 0.2, -0.1, 0.4}, {-0.5, 0.7, 0.3, -0.2}, {0.1, -0.4, 0.8, 0.9}} {
		x, y, z, t := p[0], p[1], p[2], p[3]
		dt := h.F(x, y, z, t+e)
		dm := h.F(x, y, z, t-e)
		// curl F by central differences
		d := func(c int, ax int) complex128 {
			q1, q2 := [3]float64{x, y, z}, [3]float64{x, y, z}
			q1[ax] += e
			q2[ax] -= e
			return (h.F(q1[0], q1[1], q1[2], t)[c] - h.F(q2[0], q2[1], q2[2], t)[c]) / complex(2*e, 0)
		}
		curl := [3]complex128{d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)}
		for c := 0; c < 3; c++ {
			ft := (dt[c] - dm[c]) / complex(2*e, 0)
			plus += cmplx.Abs(1i*ft - curl[c])
			minus += cmplx.Abs(-1i*ft - curl[c])
		}
	}
	if plus <= minus {
		return 1
	}
	return -1
}
