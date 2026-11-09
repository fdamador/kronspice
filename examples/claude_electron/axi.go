// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Axisymmetric Faddeev knots, copied from examples/ether_knot/axi.go and
// generalized to Battye & Sutcliffe's A_{M,K} family (Proc. R. Soc. A 455,
// 4305, 1999): the field turns M times round the axis,
//
//	n(ρ, φ, z) = R_z(Mφ) m(ρ, z),
//
// and its phase winds K times round the core circle, so the Hopf charge
// is Q = M·K. The preimage curves (the field lines of B_n) lie on nested
// tori and wind K times the long way round for every M times through the
// hole. Williamson & van der Mark's double loop (twice round the ring,
// once round the tube) is A_{1,2}: Q = 2. Rañada's knot is A_{1,1}.
// With ∂_φ n = R_z(Mφ)(M ẑ × m),
//
//	e2 = |m_ρ|² + |m_z|² + M²|ẑ × m|²/ρ²
//	e4 = |m_ρ × m_z|² + M²(|m_ρ × (ẑ × m)|² + |m_z × (ẑ × m)|²)/ρ²
//
// and E = 2π ∫∫ (κ2 e2 + κ4 e4) ρ dρ dz.

import "math"

// Axi is m on nodes ρ_i = (i + ½)h and z_j; the outer ρ row and both z
// rows are held at the vacuum (0, 0, 1).
type Axi struct {
	NR, NZ int
	H      float64
	K2, K4 float64
	M, K   int // azimuthal and poloidal windings
	m      [3][]float64
}

func (a *Axi) at(i, j int) int   { return i + a.NR*j }
func (a *Axi) rho(i int) float64 { return (float64(i) + 0.5) * a.H }
func (a *Axi) z(j int) float64   { return -float64(a.NZ-1)/2*a.H + float64(j)*a.H }
func (a *Axi) vec(n int) [3]float64 {
	return [3]float64{a.m[0][n], a.m[1][n], a.m[2][n]}
}
func (a *Axi) free(i, j int) bool { return i < a.NR-1 && j > 0 && j < a.NZ-1 }

// RationalMap is the A_{M,K} start: with (Z1, Z0) the point of the
// 3-sphere over x (inverse stereographic projection, scale s),
// Z1 = 2(x + iy)/(r² + 1), Z0 = (2z + i(r² − 1))/(r² + 1), the value
// W = Z1^M / Z0^K and n = (2 Re W, 2 Im W, 1 − |W|²)/(1 + |W|²). It tends
// to the vacuum (0, 0, 1) at infinity, and its Hopf charge is M·K.
func RationalMap(x, y, z, s float64, M, K int) [3]float64 {
	x, y, z = x/s, y/s, z/s
	r2 := x*x + y*y + z*z
	d := r2 + 1
	z1 := complex(2*x/d, 2*y/d)
	z0 := complex(2*z/d, (r2-1)/d)
	w := cpow(z1, M) / cpow(z0, K)
	w2 := real(w)*real(w) + imag(w)*imag(w)
	return [3]float64{2 * real(w) / (1 + w2), 2 * imag(w) / (1 + w2), (1 - w2) / (1 + w2)}
}

func cpow(c complex128, n int) complex128 {
	p := complex(1, 0)
	for k := 0; k < n; k++ {
		p *= c
	}
	return p
}

// compact pulls space in so the map reaches the vacuum smoothly at rc
// (as ether_knot's CompactHopf; the map and its charge are unchanged).
func compact(x, y, z, s, rc float64, M, K int) [3]float64 {
	r2 := (x*x + y*y + z*z) / (rc * rc)
	if r2 >= 1 {
		return [3]float64{0, 0, 1}
	}
	f := 1 / (1 - r2)
	return RationalMap(x*f, y*f, z*f, s, M, K)
}

// NewAxi fills the half-plane with the A_{M,K} map at φ = 0 (scale s),
// reaching the vacuum at 0.9·min(R, Z).
func NewAxi(R, Z, h, s, k2, k4 float64, M, K int) *Axi {
	a := &Axi{NR: int(math.Round(R / h)), NZ: 2*int(math.Round(Z/h)) + 1, H: h, K2: k2, K4: k4, M: M, K: K}
	for c := 0; c < 3; c++ {
		a.m[c] = make([]float64, a.NR*a.NZ)
	}
	for j := 0; j < a.NZ; j++ {
		for i := 0; i < a.NR; i++ {
			v := [3]float64{0, 0, 1}
			if a.free(i, j) {
				v = compact(a.rho(i), 0, a.z(j), s, 0.9*math.Min(R, Z), M, K)
			}
			for c := 0; c < 3; c++ {
				a.m[c][a.at(i, j)] = v[c]
			}
		}
	}
	return a
}

func cross(a, b [3]float64) [3]float64 {
	return [3]float64{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}
func dot(a, b [3]float64) float64 { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
func sub(p, q [3]float64) [3]float64 {
	return [3]float64{p[0] - q[0], p[1] - q[1], p[2] - q[2]}
}
func scale(p [3]float64, s float64) [3]float64 { return [3]float64{p[0] * s, p[1] * s, p[2] * s} }

// zc is M ẑ × m.
func (a *Axi) zc(m [3]float64) [3]float64 {
	M := float64(a.M)
	return [3]float64{-M * m[1], M * m[0], 0}
}

// mirror is the axis image m(−1) = R_z(Mπ) m(0).
func (a *Axi) mirror(p [3]float64) [3]float64 {
	if a.M%2 == 0 {
		return p
	}
	return [3]float64{-p[0], -p[1], p[2]}
}

// Energy is ether_knot's discretization with the winding M: E2, E4
// (couplings excluded) and, with grad, ∂(κ2E2 + κ4E4)/∂m added into it.
func (a *Axi) Energy(grad *[3][]float64) (e2, e4 float64) {
	h := a.H
	w := 2 * math.Pi * h * h
	M2 := float64(a.M * a.M)
	add := func(n int, g [3]float64, k float64) {
		if grad == nil {
			return
		}
		for c := 0; c < 3; c++ {
			(*grad)[c][n] += k * g[c]
		}
	}
	for j := 0; j < a.NZ; j++ {
		for i := 0; i < a.NR; i++ {
			n := a.at(i, j)
			p := a.vec(n)
			r := a.rho(i)
			if i+1 < a.NR {
				q := a.vec(a.at(i+1, j))
				d := sub(q, p)
				wl := w * (r + h/2) / (h * h)
				e2 += wl * dot(d, d)
				if a.free(i, j) {
					add(n, d, -2*wl*a.K2)
				}
				if a.free(i+1, j) {
					add(a.at(i+1, j), d, 2*wl*a.K2)
				}
			}
			if j+1 < a.NZ {
				q := a.vec(a.at(i, j+1))
				d := sub(q, p)
				wl := w * r / (h * h)
				e2 += wl * dot(d, d)
				if a.free(i, j) {
					add(n, d, -2*wl*a.K2)
				}
				if a.free(i, j+1) {
					add(a.at(i, j+1), d, 2*wl*a.K2)
				}
			}
			c := a.zc(p)
			wn := w / r
			e2 += wn * dot(c, c)
			if a.free(i, j) {
				add(n, [3]float64{2 * M2 * p[0], 2 * M2 * p[1], 0}, wn*a.K2)
			}
			if j == 0 || j == a.NZ-1 || i == a.NR-1 {
				continue
			}
			var pm [3]float64
			im := -1
			if i == 0 {
				pm = a.mirror(p)
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
			k := wv * a.K4
			gr := cross(mz, cross(mr, mz))
			gz := cross(mr, cross(mz, mr))
			t := cross(c, cross(mr, c))
			gr = [3]float64{2 * (gr[0] + t[0]/(r*r)), 2 * (gr[1] + t[1]/(r*r)), 2 * (gr[2] + t[2]/(r*r))}
			t = cross(c, cross(mz, c))
			gz = [3]float64{2 * (gz[0] + t[0]/(r*r)), 2 * (gz[1] + t[1]/(r*r)), 2 * (gz[2] + t[2]/(r*r))}
			u2, u3 := cross(mr, cross(c, mr)), cross(mz, cross(c, mz))
			gc := [3]float64{2 * (u2[0] + u3[0]) / (r * r), 2 * (u2[1] + u3[1]) / (r * r), 2 * (u2[2] + u3[2]) / (r * r)}
			// c = M ẑ × m: ∂/∂m_x = M ∂/∂c_y, ∂/∂m_y = −M ∂/∂c_x
			gm := [3]float64{float64(a.M) * gc[1], -float64(a.M) * gc[0], 0}
			add(n, gm, k)
			if a.free(i+1, j) {
				add(ip, gr, k/(2*h))
			}
			if i == 0 {
				g := scale(gr, -1/(2*h))
				add(n, a.mirror(g), k)
			} else if a.free(i-1, j) {
				add(im, gr, -k/(2*h))
			}
			if a.free(i, j+1) {
				add(jp, gz, k/(2*h))
			}
			if a.free(i, j-1) {
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

// Relax is ether_knot's arrested Newton flow (Battye & Sutcliffe):
// velocity zeroed whenever the energy rises; stops after steps or when E
// changes by less than tol (relative) over 1000 steps.
func (a *Axi) Relax(steps int, dt, tol float64) (e2, e4 float64, done int) {
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

// derivs returns m, m_ρ, m_z and M ẑ×m at interior node (i, j).
func (a *Axi) derivs(i, j int) (p, mr, mz, c [3]float64) {
	h := a.H
	p = a.vec(a.at(i, j))
	pm := a.mirror(p)
	if i > 0 {
		pm = a.vec(a.at(i-1, j))
	}
	mr = scale(sub(a.vec(a.at(i+1, j)), pm), 1/(2*h))
	mz = scale(sub(a.vec(a.at(i, j+1)), a.vec(a.at(i, j-1))), 1/(2*h))
	c = a.zc(p)
	return
}

// Core is the core radius: where m_z is least on z = 0 (the circle the
// preimage of the south pole runs round).
func (a *Axi) Core() float64 {
	j0 := (a.NZ - 1) / 2
	best, core := 2.0, 0.0
	for i := 0; i < a.NR-1; i++ {
		if v := a.m[2][a.at(i, j0)]; v < best {
			best, core = v, a.rho(i)
		}
	}
	return core
}

// Winding is the number of turns of m's phase (m_x + i m_y) round a
// circle of radius r about the core point (core, 0) in the (ρ, z)
// plane: the poloidal winding K, which with M gives the Hopf charge.
func (a *Axi) Winding(r float64) float64 {
	core := a.Core()
	const n = 720
	tot, prev := 0.0, 0.0
	for k := 0; k <= n; k++ {
		th := 2 * math.Pi * float64(k) / n
		m := a.sample(core+r*math.Cos(th), r*math.Sin(th))
		ph := math.Atan2(m[1], m[0])
		if k > 0 {
			d := ph - prev
			for d > math.Pi {
				d -= 2 * math.Pi
			}
			for d < -math.Pi {
				d += 2 * math.Pi
			}
			tot += d
		}
		prev = ph
	}
	return tot / (2 * math.Pi)
}

// Degree is the degree of m on the (ρ, z) half-plane, whose axis, box
// edge and infinity all map to the vacuum: (1/4π)∫∫ m·(m_ρ × m_z) dρ dz
// (sign made positive for the A_{M,K} starts). The Hopf charge is
// Q = M·|Degree|; it is an integer in the continuum and changes on a
// lattice only when the knot unties through the grid.
func (a *Axi) Degree() float64 {
	h := a.H
	s := 0.0
	for j := 1; j < a.NZ-1; j++ {
		for i := 0; i < a.NR-1; i++ {
			p, mr, mz, _ := a.derivs(i, j)
			s += dot(p, cross(mr, mz))
		}
	}
	return -s * h * h / (4 * math.Pi)
}

// sample is m at (ρ, z) by bilinear interpolation (ρ ≥ h/2).
func (a *Axi) sample(r, z float64) [3]float64 {
	fi := r/a.H - 0.5
	fj := (z - a.z(0)) / a.H
	i0, j0 := int(math.Floor(fi)), int(math.Floor(fj))
	if i0 < 0 || i0 >= a.NR-1 || j0 < 0 || j0 >= a.NZ-1 {
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
