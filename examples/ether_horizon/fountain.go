// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
)

// Part 8a: the fountain. Instead of destroying the ether (Part 5), the
// core lets it circulate: in round the equator, swirling, out along the
// axis at the poles, and back round the outside. In an ideal ether
// (inviscid, incompressible away from the void) a steady axisymmetric flow
// with swirl obeys the Bragg–Hawthorne equation for the Stokes stream
// function ψ in spherical (r, θ):
//
//	Δ*ψ = s²·H'(ψ) − C(ψ)·C'(ψ),   Δ* = ∂²/∂r² + (sin θ/r²)·∂/∂θ(1/sin θ·∂/∂θ),
//
// s = r sin θ, H the Bernoulli head, C = s·u_φ the swirl (Bragg & Hawthorne,
// J. Aeronaut. Sci. 17, 243, 1950). With H constant and
//
//	C(ψ) = √(C0² + k²ψ²)
//
// the right side is −k²ψ, so the equation is linear: Δ*ψ + k²ψ = 0, solved
// by ψ = A·(kr)·j_l(kr)·Θ_l(θ) with Θ_1 = sin²θ, Θ_2 = sin²θ·cos θ (the
// Chandrasekhar–Kendall modes, Astrophys. J. 126, 457, 1957). The
// sphere r = a is a stream surface when j_l(ka) = 0. C0 is the circulation
// of a line vortex on the axis: its core is the void of Part 1, and the
// swirl is the same above and below the equator because C depends on ψ².
//
//	l = 1: one cell, in at the south pole, out at the north (Hicks, Moffatt).
//	l = 2: two cells, in round the equator, out at both poles — the
//	       question's flow.
//
// Units: a = 1, A = 1, ρ = 1.

// Fountain is one of the two flows.
type Fountain struct {
	L  int     // 1 or 2
	K  float64 // k·a: the first zero of j_l
	C0 float64 // axis circulation (line vortex), units of k·ψ
}

// Zeros of j_1 and j_2 (Abramowitz & Stegun, Table 10.6).
const (
	J1Zero = 4.493409457909064
	J2Zero = 5.763459196894550
)

// NewFountain makes the l = 1 or l = 2 flow with axis circulation c0.
func NewFountain(l int, c0 float64) Fountain {
	k := J1Zero
	if l == 2 {
		k = J2Zero
	}
	return Fountain{L: l, K: k, C0: c0}
}

func sphj(l int, x float64) float64 {
	if x < 1e-4 {
		if l == 1 {
			return x / 3
		}
		return x * x / 15
	}
	s, c := math.Sin(x), math.Cos(x)
	if l == 1 {
		return s/(x*x) - c/x
	}
	return (3/(x*x)-1)*s/x - 3*c/(x*x)
}

// Radial and angular parts and their derivatives.
func (f Fountain) rad(r float64) (R, dR float64) {
	h := 1e-6
	g := func(r float64) float64 { return f.K * r * sphj(f.L, f.K*r) }
	return g(r), (g(r+h) - g(r-h)) / (2 * h)
}

func (f Fountain) ang(th float64) (T, dT float64) {
	s, c := math.Sin(th), math.Cos(th)
	if f.L == 1 {
		return s * s, 2 * s * c
	}
	return s * s * c, 2*s*c*c - s*s*s
}

// Psi is the stream function ψ(r, θ).
func (f Fountain) Psi(r, th float64) float64 {
	R, _ := f.rad(r)
	T, _ := f.ang(th)
	return R * T
}

// Swirl is C(ψ) = √(C0² + k²ψ²).
func (f Fountain) Swirl(psi float64) float64 { return math.Hypot(f.C0, f.K*psi) }

// Velocity returns (u_r, u_θ, u_φ) from u_r = ψ_θ/(r² sin θ), u_θ = −ψ_r/(r sin θ), u_φ = C/(r sin θ).
func (f Fountain) Velocity(r, th float64) (ur, ut, up float64) {
	R, dR := f.rad(r)
	T, dT := f.ang(th)
	s := math.Sin(th)
	return R * dT / (r * r * s), -dR * T / (r * s), f.Swirl(R*T) / (r * s)
}

// Residual checks the Bragg–Hawthorne equation by finite differences:
// |Δ*ψ + k²ψ| relative to |k²ψ| at (r, θ).
func (f Fountain) Residual(r, th float64) float64 {
	h := 1e-4
	p := f.Psi(r, th)
	prr := (f.Psi(r+h, th) - 2*p + f.Psi(r-h, th)) / (h * h)
	g := func(t float64) float64 { return (f.Psi(r, t+h/2) - f.Psi(r, t-h/2)) / h / math.Sin(t) }
	pth := math.Sin(th) / (r * r) * (g(th+h/2) - g(th-h/2)) / h
	return math.Abs(prr+pth+f.K*f.K*p) / math.Abs(f.K*f.K*p)
}

// SphereFlux is the net outward flux through the sphere of radius r,
// ∫u_r dA = 2π[ψ(r, π) − ψ(r, 0)]·… evaluated by quadrature (it must be 0:
// nothing is destroyed).
func (f Fountain) SphereFlux(r float64, n int) float64 {
	s := 0.0
	for i := 0; i < n; i++ {
		th := (float64(i) + 0.5) * math.Pi / float64(n)
		ur, _, _ := f.Velocity(r, th)
		s += ur * 2 * math.Pi * r * r * math.Sin(th) * math.Pi / float64(n)
	}
	return s
}

// CellFlux is the volume flux circulating in one cell, 2π·max|ψ|.
func (f Fountain) CellFlux() float64 {
	ro, to := f.OPoint()
	return 2 * math.Pi * math.Abs(f.Psi(ro, to))
}

// OPoint is the centre of the (upper) cell, where ∇ψ = 0: θ from Θ' = 0,
// r from (kr·j_l)' = 0.
func (f Fountain) OPoint() (r, th float64) {
	th = math.Pi / 2
	if f.L == 2 {
		th = math.Acos(1 / math.Sqrt(3))
	}
	lo, hi := 0.05, 0.999
	for i := 0; i < 80; i++ {
		m := 0.5 * (lo + hi)
		if _, d := f.rad(m); d > 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi), th
}

// Surface is one closed streamline of the poloidal flow, traced once round.
type Surface struct {
	Start   float64 // starting radius on the ray from the O-point (fraction of the way to the edge)
	Q       float64 // toroidal turns per poloidal turn
	Major   float64 // distance of the O-point from the axis (in a)
	Minor   float64 // mean distance of the streamline from the O-point (in a)
	Period  float64 // poloidal transit time
	Closure float64 // distance between start and end point of the traced loop (in a)
	S, Z    []float64
}

// Trace follows the streamline that starts at fraction x of the way from
// the O-point to the sphere along the cell's own ray (θ = θ_O), once round
// the poloidal loop, with RK4 in (s, z) and dφ/dt = C/s².
func (f Fountain) Trace(x float64, dt float64) Surface {
	ro, to := f.OPoint()
	r0 := ro + x*(1-ro)
	if f.L == 1 {
		r0 = ro + x*(1-ro) // on the equator, outward
	}
	s, z := r0*math.Sin(to), r0*math.Cos(to)
	so, zo := ro*math.Sin(to), ro*math.Cos(to)
	vel := func(s, z float64) (float64, float64, float64) {
		r := math.Hypot(s, z)
		th := math.Atan2(s, z)
		ur, ut, up := f.Velocity(r, th)
		return ur*math.Sin(th) + ut*math.Cos(th), ur*math.Cos(th) - ut*math.Sin(th), up / s
	}
	out := Surface{Start: x, Major: so}
	ang0 := math.Atan2(z-zo, s-so)
	phi, t, wound := 0.0, 0.0, 0.0
	prev := ang0
	sumR, n := 0.0, 0
	for step := 0; step < 2000000; step++ {
		k1s, k1z, k1p := vel(s, z)
		k2s, k2z, k2p := vel(s+dt/2*k1s, z+dt/2*k1z)
		k3s, k3z, k3p := vel(s+dt/2*k2s, z+dt/2*k2z)
		k4s, k4z, k4p := vel(s+dt*k3s, z+dt*k3z)
		s += dt / 6 * (k1s + 2*k2s + 2*k3s + k4s)
		z += dt / 6 * (k1z + 2*k2z + 2*k3z + k4z)
		phi += dt / 6 * (k1p + 2*k2p + 2*k3p + k4p)
		t += dt
		out.S = append(out.S, s)
		out.Z = append(out.Z, z)
		sumR += math.Hypot(s-so, z-zo)
		n++
		a := math.Atan2(z-zo, s-so)
		d := a - prev
		if d > math.Pi {
			d -= 2 * math.Pi
		}
		if d < -math.Pi {
			d += 2 * math.Pi
		}
		wound += d
		prev = a
		if math.Abs(wound) >= 2*math.Pi {
			// Interpolate back to exactly one turn.
			over := (math.Abs(wound) - 2*math.Pi) / math.Abs(d)
			t -= over * dt
			phi -= over * dt * k4p
			break
		}
	}
	out.Q = math.Abs(phi) / (2 * math.Pi)
	out.Period = t
	out.Minor = sumR / float64(n)
	out.Closure = math.Hypot(s-r0*math.Sin(to), z-r0*math.Cos(to))
	return out
}

// QProfile traces surfaces at the given fractions.
func (f Fountain) QProfile(xs []float64, dt float64) []Surface {
	out := make([]Surface, len(xs))
	for i, x := range xs {
		out[i] = f.Trace(x, dt)
	}
	return out
}

// SurfaceForQ finds by bisection the surface whose winding is q (q grows
// outward from the O-point to the separatrix, where it diverges).
func (f Fountain) SurfaceForQ(q, dt float64) (Surface, bool) {
	lo, hi := 0.02, 0.98
	slo, shi := f.Trace(lo, dt), f.Trace(hi, dt)
	if (slo.Q-q)*(shi.Q-q) > 0 {
		return Surface{}, false
	}
	var m Surface
	for i := 0; i < 40; i++ {
		x := 0.5 * (lo + hi)
		m = f.Trace(x, dt)
		if (m.Q-q)*(slo.Q-q) > 0 {
			lo, slo = x, m
		} else {
			hi = x
		}
	}
	return m, true
}

// AngularMomentum is ρ∫C dV inside the sphere (ρ = 1), by quadrature;
// the line vortex's C0 counts in full.
func (f Fountain) AngularMomentum(n int) float64 {
	s := 0.0
	for i := 0; i < n; i++ {
		r := (float64(i) + 0.5) / float64(n)
		for j := 0; j < n; j++ {
			th := (float64(j) + 0.5) * math.Pi / float64(n)
			dv := 2 * math.Pi * r * r * math.Sin(th) / float64(n) * math.Pi / float64(n)
			s += f.Swirl(f.Psi(r, th)) * dv
		}
	}
	return s
}

// Energy is the kinetic energy ½∫|u|² dV inside the sphere, outside a
// void of radius xi round the axis (the swirl C0/s diverges as ln(1/xi)).
func (f Fountain) Energy(xi float64, n int) float64 {
	e := 0.0
	for i := 0; i < n; i++ {
		r := (float64(i) + 0.5) / float64(n)
		for j := 0; j < 2*n; j++ {
			th := (float64(j) + 0.5) * math.Pi / float64(2*n)
			if r*math.Sin(th) < xi {
				continue
			}
			ur, ut, up := f.Velocity(r, th)
			dv := 2 * math.Pi * r * r * math.Sin(th) / float64(n) * math.Pi / float64(2*n)
			e += 0.5 * (ur*ur + ut*ut + up*up) * dv
		}
	}
	return e
}

// PsiO is |ψ| at the cell's centre: the natural unit for C0 (C0 = β·k·ψ_O).
func (f Fountain) PsiO() float64 {
	ro, to := f.OPoint()
	return math.Abs(f.Psi(ro, to))
}

// FitTorus finds the axis circulation β (C0 = β·k·ψ_O) for which the
// streamline surface winding q has tube/ring ratio aspect, by bisection
// (more swirl moves the q surface inward). It returns β and that surface.
func FitTorus(l int, q, aspect, dt float64) (float64, Surface, bool) {
	unit := NewFountain(l, 0)
	unit.C0 = 0
	ko := unit.K * unit.PsiO()
	ratio := func(beta float64) (float64, Surface, bool) {
		s, ok := NewFountain(l, beta*ko).SurfaceForQ(q, dt)
		return s.Minor / s.Major, s, ok
	}
	lo, hi := 0.3, 6.0
	var best Surface
	for i := 0; i < 40; i++ {
		m := 0.5 * (lo + hi)
		r, s, ok := ratio(m)
		if !ok {
			// q not reached inside the cell: too little swirl near the edge, or too much at the centre.
			if NewFountain(l, m*ko).Trace(0.02, dt).Q > q {
				hi = m
			} else {
				lo = m
			}
			continue
		}
		best = s
		if r > aspect {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi), best, best.Q > 0
}
