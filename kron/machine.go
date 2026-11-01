// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// Kron's rotating electric machinery (G. Kron, Equivalent Circuits of
// Electric Machinery, Wiley 1951, Epilogue "The Electrodynamics of
// Equivalent Circuits", pp. 253–259; G. Kron, "Non-Riemannian Dynamics of
// Rotating Electrical Machinery," J. Math. Phys. 13, 103–194, 1934).
//
// A machine's windings and its rotor angle are the coordinates
// x = (q_1 … q_n, θ) of one dynamical system. Its "metric tensor" g holds
// every self- and mutual inductance L(θ) and the rotor inertia J; the
// holonomic frame I (axes fixed to the iron) gives Kron's eq. 3
//
//	e_α = R_αβ·ẋ^β + g_αβ·ẍ^β + [βγ, α]·ẋ^β·ẋ^γ
//	[βγ, α] = ½·(∂g_βα/∂x^γ + ∂g_γα/∂x^β − ∂g_βγ/∂x^α)     (eq. 4)
//
// the motion of a point in an n-dimensional Riemannian space: the winding
// rows are e = Ri + p(Li) (speed voltages from the Christoffel symbols),
// the θ row is the torque equation. Axes that do not follow the iron (the
// two-reaction, d-q, cross-field frames) are non-holonomic: the connection
// tensor C(θ) depends on the coordinate itself, and transforming adds the
// speed term G·i·pθ with G from ∂(L·C)/∂θ (eqs. 6–9). With such a C the
// metric becomes constant — the trigonometric terms disappear — and the
// geometry carries the motion instead.
//
// Christoffel and RotatingFrame implement both views for kronspice's N
// element (circuit.Machine), whose equations are the frame-I or
// frame-II(b) form depending on its spec.

// Christoffel returns the Christoffel symbols of the first kind
// [βγ, α] (eq. 4) of a machine's metric at angle θ, for the coordinates
// (q_1 … q_n, θ): only θ-derivatives of the inductances are non-zero (the
// inertia is constant). Index order: out[β][γ][α].
func Christoffel(s *circuit.MachineSpec, theta float64) [][][]float64 {
	n := s.N
	d1 := make([][]float64, n)
	for i := range d1 {
		d1[i] = make([]float64, n)
	}
	s.Inductance(theta, nil, d1, nil)
	N := n + 1
	dg := func(a, b, c int) float64 { // ∂g_ab/∂x^c
		if c != n || a == n || b == n {
			return 0
		}
		return d1[a][b]
	}
	out := make([][][]float64, N)
	for b := 0; b < N; b++ {
		out[b] = make([][]float64, N)
		for c := 0; c < N; c++ {
			out[b][c] = make([]float64, N)
			for a := 0; a < N; a++ {
				out[b][c][a] = 0.5 * (dg(b, a, c) + dg(c, a, b) - dg(b, c, a))
			}
		}
	}
	return out
}

// Geodesic evaluates the right-hand side of Kron's eq. 3 without
// resistance: e_α = g_αβ·ẍ^β + [βγ, α]·ẋ^β·ẋ^γ, for winding currents i,
// their derivatives di, rotor speed ω and acceleration dω, with inertia J.
// The first n entries are the winding voltages (beyond resistance), the
// last is the torque that must be applied to the shaft.
func Geodesic(s *circuit.MachineSpec, J, theta float64, i, di []float64, omega, domega float64) []float64 {
	n := s.N
	l := make([][]float64, n)
	for k := range l {
		l[k] = make([]float64, n)
	}
	s.Inductance(theta, l, nil, nil)
	v := append(append([]float64(nil), i...), omega)
	a := append(append([]float64(nil), di...), domega)
	gam := Christoffel(s, theta)
	N := n + 1
	g := func(p, q int) float64 {
		switch {
		case p < n && q < n:
			return l[p][q]
		case p == n && q == n:
			return J
		}
		return 0
	}
	e := make([]float64, N)
	for al := 0; al < N; al++ {
		for b := 0; b < N; b++ {
			e[al] += g(al, b) * a[b]
			for c := 0; c < N; c++ {
				e[al] += gam[b][c][al] * v[b] * v[c]
			}
		}
	}
	return e
}

// MovingFrame is a change of winding axes that depends on the rotor angle:
// winding currents i = C(θ)·i′, voltages e′ = Cᵀ·e (power invariant).
// D is dC/dθ.
type MovingFrame struct {
	C, D func(theta float64) [][]float64
}

// RotorFrame returns the frame that refers the rotor winding pairs (each
// a two-phase winding a, b displaced by 90 electrical degrees, at
// electrical angle p·θ) to axes fixed to the stator: the pair's block is
// the rotation R(−p·θ), the identity elsewhere (Park's transformation in
// two-phase form; Kron 1951 Ch. 4).
func RotorFrame(n int, pairs [][2]int, p float64) MovingFrame {
	build := func(theta float64, deriv bool) [][]float64 {
		c := make([][]float64, n)
		for i := range c {
			c[i] = make([]float64, n)
			if !deriv {
				c[i][i] = 1
			}
		}
		a := -p * theta
		co, si := math.Cos(a), math.Sin(a)
		for _, pr := range pairs {
			x, y := pr[0], pr[1]
			if deriv {
				// d/dθ of [[cos a, −sin a], [sin a, cos a]] with a = −pθ.
				c[x][x], c[x][y], c[y][x], c[y][y] = p*si, p*co, -p*co, p*si
			} else {
				c[x][x], c[x][y], c[y][x], c[y][y] = co, -si, si, co
			}
		}
		return c
	}
	return MovingFrame{C: func(t float64) [][]float64 { return build(t, false) }, D: func(t float64) [][]float64 { return build(t, true) }}
}

// RotatingFrame transforms a machine to frame f (Kron 1951, Epilogue eqs.
// 6–9): the new inductances are Cᵀ·L(θ)·C, which must be constant over a
// turn, and the new speed matrix is G′ = Cᵀ·(∂L/∂θ·C + L·∂C/∂θ) + Cᵀ·G·C.
// The result is a frame-II(b) machine (constant L, speed voltages G′·i′·ω,
// torque i′ᵀ·G′·i′) that, connected to the same circuit through e′ = Cᵀ·e,
// behaves exactly as the original. C must be orthogonal (power
// invariant).
func RotatingFrame(s *circuit.MachineSpec, f MovingFrame) (*circuit.MachineSpec, error) {
	n := s.N
	kmax := 0
	for _, h := range s.Harmonics {
		if h.K > kmax {
			kmax = h.K
		}
	}
	steps := 4*kmax + 16
	l := square(n)
	d1 := square(n)
	out := circuit.NewMachineSpec(n)
	out.Theta0 = s.Theta0
	scale := 0.0
	for a := 0; a < steps; a++ {
		th := 2 * math.Pi * float64(a) / float64(steps)
		c, dc := f.C(th), f.D(th)
		if len(c) != n || len(dc) != n {
			return nil, fmt.Errorf("%w: frame of size %d for %d windings", ErrDimension, len(c), n)
		}
		ctc := mul(transpose(c), c)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				want := 0.0
				if i == j {
					want = 1
				}
				if math.Abs(ctc[i][j]-want) > 1e-12 {
					return nil, fmt.Errorf("%w: the frame is not orthogonal at θ = %.4g (CᵀC ≠ I)", ErrUnsupported, th)
				}
			}
		}
		s.Inductance(th, l, d1, nil)
		lp := mul(transpose(c), mul(l, c))
		g := add(mul(transpose(c), add(mul(d1, c), mul(l, dc))), mul(transpose(c), mul(s.G, c)))
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				scale = math.Max(scale, math.Abs(lp[i][j]))
			}
		}
		if a == 0 {
			out.L0, out.G = lp, g
			continue
		}
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if math.Abs(lp[i][j]-out.L0[i][j]) > 1e-9*scale || math.Abs(g[i][j]-out.G[i][j]) > 1e-9*scale*float64(kmax+1) {
					return nil, fmt.Errorf("%w: in this frame the inductances still depend on θ (entry %d,%d)", ErrUnsupported, i+1, j+1)
				}
			}
		}
	}
	// Symmetrize and drop round-off.
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if math.Abs(out.L0[i][j]) < 1e-15*scale {
				out.L0[i][j] = 0
			}
			if math.Abs(out.G[i][j]) < 1e-15*scale {
				out.G[i][j] = 0
			}
		}
		for j := i + 1; j < n; j++ {
			v := 0.5 * (out.L0[i][j] + out.L0[j][i])
			out.L0[i][j], out.L0[j][i] = v, v
		}
	}
	if err := out.CheckPositive(); err != nil {
		return nil, err
	}
	return out, nil
}

func square(n int) [][]float64 {
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
	}
	return a
}

func transpose(a [][]float64) [][]float64 {
	t := square(len(a))
	for i := range a {
		for j := range a[i] {
			t[j][i] = a[i][j]
		}
	}
	return t
}

func mul(a, b [][]float64) [][]float64 {
	n := len(a)
	c := square(n)
	for i := 0; i < n; i++ {
		for k := 0; k < n; k++ {
			if a[i][k] == 0 {
				continue
			}
			for j := 0; j < n; j++ {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}

func add(a, b [][]float64) [][]float64 {
	c := square(len(a))
	for i := range a {
		for j := range a[i] {
			c[i][j] = a[i][j] + b[i][j]
		}
	}
	return c
}
