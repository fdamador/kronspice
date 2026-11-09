// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
)

// Does a plasmoid of "resonant size" stop radiating (Nachamkin §IV, eq. 35)?
//
// Inside, the plasma (permittivity ε = 1 − ωp²/ω²) carries the l = 1
// fields with wavenumber k_in = √(ω² − ωp²)/c; outside is vacuum, k = ω/c.
// Without surface currents (a plasma boundary carries none) the tangential
// E and H must be continuous, and with nothing feeding it from outside the
// exterior field can only be outgoing (spherical Hankel h1). That makes
// the l = 1 modes leaky resonances with complex frequency ω = ω′ − iγ
// (time factor e^{−iωt}); the stored energy decays as e^{−2γt}, Q = ω′/2γ.
// Mode conditions (Mie theory; J. A. Stratton, Electromagnetic Theory,
// §9.22):
//
//	TE (Nachamkin's M part):  x·j1′(x)/j1(x) = y·h1′(y)/h1(y)
//	TM (his N part):          [1 + x·j1′(x)/j1(x)]/ε = 1 + y·h1′(y)/h1(y)
//
// with x = k_in·R, y = k·R.

func cj1(z complex128) complex128 { return cmplx.Sin(z)/(z*z) - cmplx.Cos(z)/z }
func cj1p(z complex128) complex128 {
	return cmplx.Sin(z)/z - 2*cj1(z)/z
}
func ch1(z complex128) complex128 { // spherical Hankel h1^(1)
	return -cmplx.Exp(1i*z) * (z + 1i) / (z * z)
}
func ch1p(z complex128) complex128 {
	h0 := -1i * cmplx.Exp(1i*z) / z
	return h0 - 2*ch1(z)/z
}

// Sphere is a uniform cold-plasma sphere.
type Sphere struct {
	R, Wp float64
}

// residual of the TE (te = true) or TM mode condition at complex ω, in a
// form without poles (multiplied through by j1·h1).
func (s Sphere) residual(w complex128, te bool) complex128 {
	kin := cmplx.Sqrt(w*w-complex(s.Wp*s.Wp, 0)) / C0
	x := kin * complex(s.R, 0)
	y := w / C0 * complex(s.R, 0)
	jx, jpx := cj1(x), cj1p(x)
	hy, hpy := ch1(y), ch1p(y)
	if te {
		return x*jpx*hy - y*hpy*jx
	}
	eps := 1 - complex(s.Wp*s.Wp, 0)/(w*w)
	return (jx+x*jpx)*hy - eps*(hy+y*hpy)*jx
}

// scaled residual, normalized by |h1(y)| so its size does not depend on
// the exterior argument.
func (s Sphere) scaled(w complex128, te bool) complex128 {
	y := w / C0 * complex(s.R, 0)
	return s.residual(w, te) / complex(cmplx.Abs(ch1(y)), 0)
}

// Mode is a complex resonance.
type Mode struct {
	TE       bool
	Omega    complex128 // ω′ − iγ
	Q        float64
	Lifetime float64 // energy e-folding time 1/(2γ) [s]
	Found    bool
}

// TMZero is the first zero of d(x·j1(x))/dx, where the l = 1 TM mode of
// an ε → 0 sphere sits.
func TMZero() float64 {
	f := func(x float64) float64 { return j1(x) + x*j1p(x) }
	lo, hi := 1.0, J1Zero
	for i := 0; i < 100; i++ {
		m := 0.5 * (lo + hi)
		if f(lo)*f(m) <= 0 {
			hi = m
		} else {
			lo = m
		}
	}
	return 0.5 * (lo + hi)
}

// FindMode finds the lowest l = 1 TE or TM resonance of a sphere whose
// plasma frequency is just below the wave frequency (the plasmoid regime).
// Start: the interior argument x0 at the first TE zero (j1 = 0, Nachamkin's
// αR) or the first TM zero (d(xj1)/dx = 0), shifted by the outgoing-wave
// loss −i·x0/(kR) for TE; ω from ω² = ωp² + (c·x/R)²; then Newton.
func (s Sphere) FindMode(te bool) Mode {
	x := J1Zero
	if !te {
		x = TMZero()
	}
	return s.FindModeFrom(te, x)
}

// FindModeFrom finds the l = 1 TE or TM resonance whose interior argument
// starts at xStart (a zero of j1 for TE, of d(x·j1)/dx for TM): the
// overtones of the sphere.
func (s Sphere) FindModeFrom(te bool, xStart float64) Mode {
	x0 := complex(xStart, 0)
	y0 := s.Wp * s.R / C0
	if te {
		x0 -= 1i * x0 / complex(y0, 0)
	}
	k := x0 / complex(s.R, 0)
	w := cmplx.Sqrt(complex(s.Wp*s.Wp, 0) + C0*C0*k*k)
	w0 := s.Wp
	m := Mode{TE: te}
	for it := 0; it < 200; it++ {
		f := s.scaled(w, te)
		h := complex(1e-7*w0, 0)
		d := (s.scaled(w+h, te) - s.scaled(w-h, te)) / (2 * h)
		step := f / d
		// damp large steps
		if cmplx.Abs(step) > 0.02*w0 {
			step *= complex(0.02*w0/cmplx.Abs(step), 0)
		}
		w -= step
		if cmplx.Abs(step) < 1e-13*w0 {
			m.Found = true
			break
		}
	}
	m.Omega = w
	g := -imag(w)
	m.Q = real(w) / (2 * g)
	m.Lifetime = 1 / (2 * g)
	return m
}

func cmplxAbs(z complex128) float64 { return cmplx.Abs(z) }

// Large-kR lifetimes. For y = kR ≫ 1 the outgoing-wave factor is
// 1 + y·h1′(y)/h1(y) → i·y, and the mode conditions above close:
//
//	TE: x·j1′/j1 = i·y  →  δx = −i·x1/y  (x1 = 4.4934)
//	TM: 1 + x·j1′/j1 = ε·i·y, ε = x²/y²  →  δx = i·x_t²/(y·g′(x_t))
//
// with g(x) = 1 + x·j1′(x)/j1(x) and x_t its first zero (2.7437). Since
// ω² = ωp² + c²x²/R², δω = c·x·δx/(R·y), so the energy lifetimes 1/(2γ)
// are
//
//	τ_TE = R·y²/(2c·x1²),   τ_TM = R·y²·|g′(x_t)|/(2c·x_t³).
//
// The complex root finder above checks these where it converges.

// TauTE is the large-kR TE energy lifetime [s].
func TauTE(R, y float64) float64 { return R * y * y / (2 * C0 * J1Zero * J1Zero) }

// TauTM is the large-kR TM energy lifetime [s].
func TauTM(R, y float64) float64 {
	xt := TMZero()
	g := func(x float64) float64 { return 1 + x*j1p(x)/j1(x) }
	h := 1e-6
	gp := math.Abs((g(xt+h) - g(xt-h)) / (2 * h))
	return R * y * y * gp / (2 * C0 * xt * xt * xt)
}
