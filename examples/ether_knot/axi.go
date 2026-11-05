// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The Q = 1 knot is axially symmetric (Battye & Sutcliffe 1998: "an
// axially symmetric lump"), with the field turning once round the axis:
// n(ρ, φ, z) = R_z(φ) m(ρ, z), m a unit vector on the half-plane ρ > 0.
// Rañada's Hopf map has this form exactly (n_x + i n_y = 2αβ̄ ∝ e^{iφ}).
// With ∂_φ n = R_z(φ)(ẑ × m), the energy densities become
//
//	e2 = |m_ρ|² + |m_z|² + |ẑ × m|²/ρ²
//	e4 = |m_ρ × m_z|² + (|m_ρ × (ẑ × m)|² + |m_z × (ẑ × m)|²)/ρ²
//
// and E = 2π ∫∫ (κ2 e2 + κ4 e4) ρ dρ dz. A unit vector m (rather than two
// angles) has no coordinate singularity at the knot's core or on the axis.

import (
	"math"
)

// Axi is m on nodes ρ_i = (i + ½)h, i = 0..NR−1, and z_j = −Z + jh,
// j = 0..NZ−1. The outer ρ row and both z rows are held at the vacuum
// (0, 0, 1); next to the axis, ρ_0 = h/2, the 1/ρ² terms force m → ẑ.
type Axi struct {
	NR, NZ int
	H      float64
	K2, K4 float64
	m      [3][]float64
}

func (a *Axi) at(i, j int) int   { return i + a.NR*j }
func (a *Axi) rho(i int) float64 { return (float64(i) + 0.5) * a.H }
func (a *Axi) z(j int) float64   { return -float64(a.NZ-1)/2*a.H + float64(j)*a.H }
func (a *Axi) vec(n int) [3]float64 {
	return [3]float64{a.m[0][n], a.m[1][n], a.m[2][n]}
}

// NewAxi fills the half-plane from the compact Hopf map at scale s (the
// map at φ = 0 is m), reaching the vacuum at 0.9·min(R, Z).
func NewAxi(R, Z, h, s, k2, k4 float64) *Axi {
	a := &Axi{NR: int(math.Round(R / h)), NZ: 2*int(math.Round(Z/h)) + 1, H: h, K2: k2, K4: k4}
	for c := 0; c < 3; c++ {
		a.m[c] = make([]float64, a.NR*a.NZ)
	}
	for j := 0; j < a.NZ; j++ {
		for i := 0; i < a.NR; i++ {
			v := [3]float64{0, 0, 1}
			if i < a.NR-1 && j > 0 && j < a.NZ-1 {
				v = CompactHopf(a.rho(i), 0, a.z(j), s, 0.9*math.Min(R, Z))
			}
			for c := 0; c < 3; c++ {
				a.m[c][a.at(i, j)] = v[c]
			}
		}
	}
	return a
}

func zc(m [3]float64) [3]float64 { return [3]float64{-m[1], m[0], 0} } // ẑ × m
func sub(p, q [3]float64) [3]float64 {
	return [3]float64{p[0] - q[0], p[1] - q[1], p[2] - q[2]}
}
func scale(p [3]float64, s float64) [3]float64 { return [3]float64{p[0] * s, p[1] * s, p[2] * s} }

// Energy returns E2, E4 (2π ∫∫ … ρ dρ dz, couplings excluded) and, with
// grad, adds ∂(κ2E2 + κ4E4)/∂m into it. Discretization: the |m_ρ|², |m_z|²
// terms on links (forward differences, link-midpoint ρ), |ẑ×m|²/ρ² at
// nodes, and e4 at interior nodes with central differences (as the 3-D
// grid). Axis side: node i = 0 has a mirror neighbour i = −1 with
// m(−1) = R_z(π) m(0) (the same physical point set across the axis), so
// m_ρ at the axis is (m(1) − R_z(π)m(0))/2h for e4 and the link 0 – (−1)
// carries e2.
func (a *Axi) Energy(grad *[3][]float64) (e2, e4 float64) {
	h := a.H
	w := 2 * math.Pi * h * h
	flip := func(p [3]float64) [3]float64 { return [3]float64{-p[0], -p[1], p[2]} }
	add := func(n int, g [3]float64, k float64) {
		if grad == nil {
			return
		}
		for c := 0; c < 3; c++ {
			(*grad)[c][n] += k * g[c]
		}
	}
	free := func(i, j int) bool { return i < a.NR-1 && j > 0 && j < a.NZ-1 }
	for j := 0; j < a.NZ; j++ {
		for i := 0; i < a.NR; i++ {
			n := a.at(i, j)
			p := a.vec(n)
			r := a.rho(i)
			// ρ-link to i+1 at ρ + h/2.
			if i+1 < a.NR {
				q := a.vec(a.at(i+1, j))
				d := sub(q, p)
				wl := w * (r + h/2) / (h * h)
				e2 += wl * dot(d, d)
				if free(i, j) {
					add(n, d, -2*wl*a.K2)
				}
				if free(i+1, j) {
					add(a.at(i+1, j), d, 2*wl*a.K2)
				}
			}
			// mirror link across the axis (i = 0 to its image), weight ρ = 0
			// at the axis: zero; the 1/ρ² term below does the work there.
			// z-link to j+1.
			if j+1 < a.NZ {
				q := a.vec(a.at(i, j+1))
				d := sub(q, p)
				wl := w * r / (h * h)
				e2 += wl * dot(d, d)
				if free(i, j) {
					add(n, d, -2*wl*a.K2)
				}
				if free(i, j+1) {
					add(a.at(i, j+1), d, 2*wl*a.K2)
				}
			}
			// |ẑ × m|²/ρ².
			c := zc(p)
			wn := w * r / (r * r)
			e2 += wn * dot(c, c)
			if free(i, j) {
				add(n, [3]float64{2 * p[0], 2 * p[1], 0}, wn*a.K2)
			}
			// e4 at interior nodes.
			if j == 0 || j == a.NZ-1 || i == a.NR-1 {
				continue
			}
			var pm [3]float64
			im := -1
			if i == 0 {
				pm = flip(p)
			} else {
				im = a.at(i-1, j)
				pm = a.vec(im)
			}
			ip := a.at(i+1, j)
			jp, jm := a.at(i, j+1), a.at(i, j-1)
			mr := scale(sub(a.vec(ip), pm), 1/(2*h))
			mz := scale(sub(a.vec(jp), a.vec(jm)), 1/(2*h))
			x1 := cross(mr, mz)
			x2 := cross(mr, c)
			x3 := cross(mz, c)
			wv := w * r
			e4 += wv * (dot(x1, x1) + (dot(x2, x2)+dot(x3, x3))/(r*r))
			if grad == nil {
				continue
			}
			// ∂/∂mr, ∂/∂mz, ∂/∂c of the node's e4·wv·κ4.
			k := wv * a.K4
			gr := cross(mz, cross(mr, mz))
			gz := cross(mr, cross(mz, mr))
			t := cross(c, cross(mr, c))
			gr = [3]float64{2 * (gr[0] + t[0]/(r*r)), 2 * (gr[1] + t[1]/(r*r)), 2 * (gr[2] + t[2]/(r*r))}
			t = cross(c, cross(mz, c))
			gz = [3]float64{2 * (gz[0] + t[0]/(r*r)), 2 * (gz[1] + t[1]/(r*r)), 2 * (gz[2] + t[2]/(r*r))}
			u2, u3 := cross(mr, cross(c, mr)), cross(mz, cross(c, mz))
			gc := [3]float64{2 * (u2[0] + u3[0]) / (r * r), 2 * (u2[1] + u3[1]) / (r * r), 2 * (u2[2] + u3[2]) / (r * r)}
			gm := cross(gc, [3]float64{0, 0, 1}) // ∂E/∂m from c = ẑ × m
			add(n, gm, k)
			if free(i+1, j) {
				add(ip, gr, k/(2*h))
			}
			if i == 0 {
				// pm = flip(p): ∂/∂p of −gr/(2h)·flip(p)
				g := scale(gr, -1/(2*h))
				add(n, [3]float64{-g[0], -g[1], g[2]}, k)
			} else if free(i-1, j) {
				add(im, gr, -k/(2*h))
			}
			if free(i, j+1) {
				add(jp, gz, k/(2*h))
			}
			if free(i, j-1) {
				add(jm, gz, -k/(2*h))
			}
		}
	}
	return e2, e4
}

// Total is κ2E2 + κ4E4.
func (a *Axi) Total() float64 {
	e2, e4 := a.Energy(nil)
	return a.K2*e2 + a.K4*e4
}

// Size returns the core radius (where m_z is least on z = 0, the circle
// the knot's preimage of the south pole runs round) and the rms radius of
// the energy density.
func (a *Axi) Size() (core, rms float64) {
	j0 := (a.NZ - 1) / 2
	best := 2.0
	for i := 0; i < a.NR-1; i++ {
		if v := a.m[2][a.at(i, j0)]; v < best {
			best, core = v, a.rho(i)
		}
	}
	// energy density per node (e2 + e4 with node-local differences)
	var s, sr float64
	h := a.H
	for j := 1; j < a.NZ-1; j++ {
		for i := 1; i < a.NR-1; i++ {
			p := a.vec(a.at(i, j))
			r := a.rho(i)
			mr := scale(sub(a.vec(a.at(i+1, j)), a.vec(a.at(i-1, j))), 1/(2*h))
			mz := scale(sub(a.vec(a.at(i, j+1)), a.vec(a.at(i, j-1))), 1/(2*h))
			c := zc(p)
			x1, x2, x3 := cross(mr, mz), cross(mr, c), cross(mz, c)
			e := a.K2*(dot(mr, mr)+dot(mz, mz)+dot(c, c)/(r*r)) + a.K4*(dot(x1, x1)+(dot(x2, x2)+dot(x3, x3))/(r*r))
			zz := a.z(j)
			s += e * r
			sr += e * r * (r*r + zz*zz)
		}
	}
	return core, math.Sqrt(sr / s)
}

// Relax: arrested Newton flow (Battye & Sutcliffe): m̈ = −∂E/∂m per unit
// "mass" 2πρh², projected on the sphere; velocity zeroed when E rises.
// Stops after steps, or when E changes by less than tol (relative) over
// 1000 steps.
func (a *Axi) Relax(steps int, dt, tol float64, every int, report func(step int, e2, e4 float64)) (e2, e4 float64, done int) {
	N := a.NR * a.NZ
	var v, g [3][]float64
	for c := 0; c < 3; c++ {
		v[c] = make([]float64, N)
		g[c] = make([]float64, N)
	}
	tau := dt * a.H
	prev, ref := math.Inf(1), math.Inf(1)
	for s := 0; ; s++ {
		for c := 0; c < 3; c++ {
			for n := range g[c] {
				g[c][n] = 0
			}
		}
		e2, e4 = a.Energy(&g)
		e := a.K2*e2 + a.K4*e4
		if report != nil && every > 0 && s%every == 0 {
			report(s, e2, e4)
		}
		if s%1000 == 0 {
			if math.Abs(ref-e) < tol*e {
				return e2, e4, s
			}
			ref = e
		}
		if s >= steps {
			return e2, e4, s
		}
		if e > prev {
			for c := 0; c < 3; c++ {
				for n := range v[c] {
					v[c][n] = 0
				}
			}
		}
		prev = e
		for j := 1; j < a.NZ-1; j++ {
			for i := 0; i < a.NR-1; i++ {
				n := a.at(i, j)
				mass := 2 * math.Pi * a.rho(i) * a.H * a.H
				p := a.vec(n)
				f := [3]float64{-g[0][n] / mass, -g[1][n] / mass, -g[2][n] / mass}
				fp := dot(f, p)
				u := [3]float64{v[0][n], v[1][n], v[2][n]}
				for c := 0; c < 3; c++ {
					u[c] += tau * (f[c] - fp*p[c])
				}
				up := dot(u, p)
				for c := 0; c < 3; c++ {
					u[c] -= up * p[c]
					p[c] += tau * u[c]
				}
				l := math.Sqrt(dot(p, p))
				for c := 0; c < 3; c++ {
					a.m[c][n] = p[c] / l
					v[c][n] = u[c]
				}
			}
		}
	}
}

// To3D samples the axisymmetric field on a 3-D grid (nearest-neighbour
// bilinear interpolation in (ρ, z), rotated by φ).
func (a *Axi) To3D(N int, half float64) *Grid {
	g := NewGrid(N, half, 1, a.K2, a.K4)
	for k := 1; k < N-1; k++ {
		for j := 1; j < N-1; j++ {
			for i := 1; i < N-1; i++ {
				x, y, z := g.pos(i), g.pos(j), g.pos(k)
				r := math.Hypot(x, y)
				m := a.sample(r, z)
				cs, sn := 1.0, 0.0
				if r > 0 {
					cs, sn = x/r, y/r
				}
				v := [3]float64{cs*m[0] - sn*m[1], sn*m[0] + cs*m[1], m[2]}
				l := math.Sqrt(dot(v, v))
				n := g.idx(i, j, k)
				for c := 0; c < 3; c++ {
					g.n[c][n] = v[c] / l
				}
			}
		}
	}
	return g
}

func (a *Axi) sample(r, z float64) [3]float64 {
	fi := r/a.H - 0.5
	fj := (z - a.z(0)) / a.H
	if fi < 0 { // between the axis image and node 0: m(−1) = R_z(π)m(0)
		u := fi + 1
		p := a.vec(a.at(0, clampi(int(math.Floor(fj)), 0, a.NZ-2)))
		q := a.vec(a.at(0, clampi(int(math.Floor(fj))+1, 0, a.NZ-1)))
		tj := fj - math.Floor(fj)
		m := [3]float64{p[0]*(1-tj) + q[0]*tj, p[1]*(1-tj) + q[1]*tj, p[2]*(1-tj) + q[2]*tj}
		return [3]float64{m[0] * (2*u - 1), m[1] * (2*u - 1), m[2]}
	}
	i0, j0 := int(math.Floor(fi)), int(math.Floor(fj))
	if i0 >= a.NR-1 || j0 < 0 || j0 >= a.NZ-1 {
		return [3]float64{0, 0, 1}
	}
	ti, tj := fi-float64(i0), fj-float64(j0)
	var out [3]float64
	for _, c := range [][3]float64{{0, 0, (1 - ti) * (1 - tj)}, {1, 0, ti * (1 - tj)}, {0, 1, (1 - ti) * tj}, {1, 1, ti * tj}} {
		p := a.vec(a.at(i0+int(c[0]), j0+int(c[1])))
		for k := 0; k < 3; k++ {
			out[k] += c[2] * p[k]
		}
	}
	return out
}

func clampi(v, lo, hi int) int {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}
