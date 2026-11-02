// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// The l = 1 force-free (Beltrami, Chandrasekhar–Kendall) field that
// Nachamkin uses (his eq. 8): ∇×F = αF inside a sphere whose radius is the
// first zero of j1(αr), αR = 4.4934. The same field is the spheromak
// (Rosenbluth & Bussac, Nucl. Fusion 19, 489, 1979), in spherical
// components (r, θ, φ):
//
//	F_r = 2·j1(αr)/(αr)·cos θ
//	F_θ = −[j1(αr)/(αr) + j1′(αr)]·sin θ
//	F_φ = j1(αr)·sin θ
//
// Its field lines lie on nested tori around the ring where F_r = F_θ = 0
// (Nachamkin Figs. 2–4): the poloidal flux ψ = r·j1(αr)·sin²θ is constant
// along every line.

// J1Zero is the first zero of j1 (αR for a sphere of radius R).
const J1Zero = 4.493409457909064

func j1(x float64) float64 {
	if math.Abs(x) < 1e-4 {
		return x / 3
	}
	return math.Sin(x)/(x*x) - math.Cos(x)/x
}

func j1p(x float64) float64 { // dj1/dx
	if math.Abs(x) < 1e-4 {
		return 1.0 / 3
	}
	return math.Sin(x)/x - 2*j1(x)/x
}

// Beltrami returns (F_r, F_θ, F_φ) at (r, θ) for wavenumber α.
func Beltrami(alpha, r, th float64) (fr, ft, fp float64) {
	x := alpha * r
	var jx float64
	if x < 1e-4 {
		jx = 1.0 / 3
	} else {
		jx = j1(x) / x
	}
	return 2 * jx * math.Cos(th), -(jx + j1p(x)) * math.Sin(th), j1(x) * math.Sin(th)
}

// cart returns the field in Cartesian components at (x, y, z).
func cart(alpha float64, p [3]float64) [3]float64 {
	r := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
	rho := math.Hypot(p[0], p[1])
	th := math.Atan2(rho, p[2])
	ph := math.Atan2(p[1], p[0])
	fr, ft, fp := Beltrami(alpha, r, th)
	st, ct, sp, cp := math.Sin(th), math.Cos(th), math.Sin(ph), math.Cos(ph)
	return [3]float64{
		fr*st*cp + ft*ct*cp - fp*sp,
		fr*st*sp + ft*ct*sp + fp*cp,
		fr*ct - ft*st,
	}
}

// CurlCheck returns the largest |∇×F − αF| ÷ |F| at sample points inside
// the sphere (central differences) and the largest |∇·F|.
func CurlCheck(alpha float64) (curl, div float64) {
	R := J1Zero / alpha
	h := 1e-5 * R
	pts := [][3]float64{{0.3, 0.2, 0.1}, {-0.5, 0.1, 0.4}, {0.2, -0.6, -0.3}, {0.05, 0.7, 0.2}, {-0.3, -0.3, 0.6}}
	for _, q := range pts {
		p := [3]float64{q[0] * R, q[1] * R, q[2] * R}
		d := func(i, j int) float64 { // ∂F_i/∂x_j
			pp, pm := p, p
			pp[j] += h
			pm[j] -= h
			return (cart(alpha, pp)[i] - cart(alpha, pm)[i]) / (2 * h)
		}
		f := cart(alpha, p)
		c := [3]float64{d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)}
		n := math.Sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2])
		for i := 0; i < 3; i++ {
			curl = math.Max(curl, math.Abs(c[i]-alpha*f[i])/(alpha*n))
		}
		div = math.Max(div, math.Abs(d(0, 0)+d(1, 1)+d(2, 2))/(alpha*n))
	}
	return curl, div
}

// Psi is the poloidal flux function r·j1(αr)·sin²θ.
func Psi(alpha float64, p [3]float64) float64 {
	r := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
	rho := math.Hypot(p[0], p[1])
	st := rho / r
	return r * j1(alpha*r) * st * st
}

// Trace follows a field line from p0 (RK4, arc-length steps) and returns
// the points and the largest relative change of ψ along it.
func Trace(alpha float64, p0 [3]float64, steps int, ds float64) (pts [][3]float64, drift float64) {
	p := p0
	psi0 := Psi(alpha, p)
	dir := func(q [3]float64) [3]float64 {
		f := cart(alpha, q)
		n := math.Sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2])
		return [3]float64{f[0] / n, f[1] / n, f[2] / n}
	}
	add := func(a [3]float64, s float64, b [3]float64) [3]float64 {
		return [3]float64{a[0] + s*b[0], a[1] + s*b[1], a[2] + s*b[2]}
	}
	for i := 0; i < steps; i++ {
		k1 := dir(p)
		k2 := dir(add(p, ds/2, k1))
		k3 := dir(add(p, ds/2, k2))
		k4 := dir(add(p, ds, k3))
		for j := 0; j < 3; j++ {
			p[j] += ds / 6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
		}
		if i%4 == 0 {
			pts = append(pts, p)
		}
		drift = math.Max(drift, math.Abs(Psi(alpha, p)-psi0)/math.Abs(psi0))
	}
	return pts, drift
}

// MagneticAxis is the radius of the ring (in the equatorial plane) where
// the poloidal field vanishes: d(r·j1(αr))/dr = 0, i.e. the ring of the
// "electron ring" carrying the toroidal current.
func MagneticAxis(alpha float64) float64 {
	// solve j1(x) + x·j1′(x) = 0 for x in (1, 4.49) by bisection
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
	return 0.5 * (lo + hi) / alpha
}
