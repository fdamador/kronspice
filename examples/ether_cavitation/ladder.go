// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Ladder is the 1-D amplitude problem discretized as Kron's 1945 ladder
// network (Kron, J. Appl. Phys. 16, 172, 1945; the radial form as in
// examples/light_hydrogen): node i carries the amplitude a_i as a voltage,
// branch i joins nodes i and i+1 with conductance 2κ2·A_i/H, and node i has
// a nonlinear shunt drawing W_i·V′(a_i) to ground. The discrete energy
//
//	E = Σ_i κ2·A_i·(a_{i+1} − a_i)²/H + Σ_i W_i·(V(a_i) − Vref)
//
// has ∂E/∂a_i = −F_i, with F_i the node's KCL residual, so the network's
// operating points are exactly the stationary points of E (minima and
// saddles alike). Planar: A = 1, W = H. Radial (3-D, spherical): r_i = iH,
// A_i = (r_i + H/2)², W_i = r_i²·H, W_0 = H³/24, energies × 4π.

import (
	"errors"
	"math"
)

// Ladder holds the geometry; FixL/FixR hold the end amplitudes.
type Ladder struct {
	H          float64
	A          []float64 // len N−1
	W          []float64 // len N
	FixL, FixR bool
	Scale      float64 // 1 planar (energy per area), 4π radial
}

// Planar returns a planar ladder of n nodes, both ends fixed.
func Planar(n int, h float64) Ladder {
	l := Ladder{H: h, A: make([]float64, n-1), W: make([]float64, n), FixL: true, FixR: true, Scale: 1}
	for i := range l.A {
		l.A[i] = 1
	}
	for i := range l.W {
		l.W[i] = h
	}
	l.W[0], l.W[n-1] = h/2, h/2
	return l
}

// Radial returns a spherical ladder of n nodes on [0, (n−1)h]: centre free,
// outer end fixed.
func Radial(n int, h float64) Ladder {
	l := Ladder{H: h, A: make([]float64, n-1), W: make([]float64, n), FixR: true, Scale: 4 * math.Pi}
	for i := range l.A {
		r := (float64(i) + 0.5) * h
		l.A[i] = r * r
	}
	l.W[0] = h * h * h / 24
	for i := 1; i < n; i++ {
		r := float64(i) * h
		l.W[i] = r * r * h
	}
	return l
}

// Residual returns F (zero at fixed nodes).
func (l Ladder) Residual(e EOS, a []float64) []float64 {
	n := len(a)
	f := make([]float64, n)
	for i := 0; i < n; i++ {
		if (i == 0 && l.FixL) || (i == n-1 && l.FixR) {
			continue
		}
		s := -l.W[i] * e.D1(a[i])
		if i > 0 {
			s -= 2 * kappa2 * l.A[i-1] * (a[i] - a[i-1]) / l.H
		}
		if i < n-1 {
			s += 2 * kappa2 * l.A[i] * (a[i+1] - a[i]) / l.H
		}
		f[i] = s
	}
	return f
}

// Energy is the discrete energy above, relative to the density Vref.
func (l Ladder) Energy(e EOS, a []float64, vref float64) float64 {
	var s float64
	for i := 0; i+1 < len(a); i++ {
		d := a[i+1] - a[i]
		s += kappa2 * l.A[i] * d * d / l.H
	}
	for i, v := range a {
		s += l.W[i] * (e.V(v) - vref)
	}
	return l.Scale * s
}

var errNoConverge = errors.New("ether_cavitation: Newton did not converge")

// Solve runs damped Newton on F = 0 from a (modified in place); it finds
// minima and saddles alike, whichever is nearest. Returns the iterations.
func (l Ladder) Solve(e EOS, a []float64, tol float64, maxit int) (int, error) {
	n := len(a)
	i0, i1 := 0, n-1
	if l.FixL {
		i0 = 1
	}
	if l.FixR {
		i1 = n - 2
	}
	m := i1 - i0 + 1
	lo, di, up, rhs := make([]float64, m), make([]float64, m), make([]float64, m), make([]float64, m)
	for it := 1; it <= maxit; it++ {
		f := l.Residual(e, a)
		var norm float64
		for i := i0; i <= i1; i++ {
			norm = math.Max(norm, math.Abs(f[i])/math.Max(l.W[i], 1e-300))
		}
		if norm < tol {
			return it, nil
		}
		// Jacobian ∂F/∂a (tridiagonal).
		for k := 0; k < m; k++ {
			i := i0 + k
			d := -l.W[i] * e.D2(a[i])
			if i > 0 {
				d -= 2 * kappa2 * l.A[i-1] / l.H
				lo[k] = 2 * kappa2 * l.A[i-1] / l.H
			}
			if i < n-1 {
				d -= 2 * kappa2 * l.A[i] / l.H
				up[k] = 2 * kappa2 * l.A[i] / l.H
			}
			di[k] = d
			rhs[k] = -f[i]
		}
		dx := thomas(lo, di, up, rhs)
		var big float64
		for _, v := range dx {
			big = math.Max(big, math.Abs(v))
		}
		s := 1.0
		if big > 0.2 {
			s = 0.2 / big
		}
		for k := 0; k < m; k++ {
			a[i0+k] += s * dx[k]
		}
	}
	return maxit, errNoConverge
}

// thomas solves a tridiagonal system (lo[0], up[m−1] unused).
func thomas(lo, di, up, rhs []float64) []float64 {
	m := len(di)
	c, d := make([]float64, m), make([]float64, m)
	c[0], d[0] = up[0]/di[0], rhs[0]/di[0]
	for k := 1; k < m; k++ {
		den := di[k] - lo[k]*c[k-1]
		if k < m-1 {
			c[k] = up[k] / den
		}
		d[k] = (rhs[k] - lo[k]*d[k-1]) / den
	}
	x := make([]float64, m)
	x[m-1] = d[m-1]
	for k := m - 2; k >= 0; k-- {
		x[k] = d[k] - c[k]*x[k+1]
	}
	return x
}

// negatives counts the negative eigenvalues of the symmetric tridiagonal
// matrix (di, off) by Sturm's sequence (LDLᵀ inertia, Sylvester).
func negatives(di, off []float64) int {
	cnt := 0
	q := di[0]
	if q < 0 {
		cnt++
	}
	for k := 1; k < len(di); k++ {
		if q == 0 {
			q = 1e-300
		}
		q = di[k] - off[k-1]*off[k-1]/q
		if q < 0 {
			cnt++
		}
	}
	return cnt
}
