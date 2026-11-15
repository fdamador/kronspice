// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
)

// Part 13: the full ether simulation of the donut void, with no thin-ring
// formulas.
//
// Two fields on an axisymmetric (ρ, z) grid, in healing-length units:
//
//	ether  ψ(ρ, z): Gross–Pitaevskii, E_e = ∫ ½|∇ψ|² + ¼(|ψ|² − 1)² dV;
//	       a vortex ring is a phase winding round a point (R, 0) of the
//	       meridional plane;
//	light  χ(ρ, z)·e^{iφ/2}: the double-loop mode (m = ½, antiperiodic
//	       round the axis), lowest eigenmode of
//	       H = −∇² + m²/ρ² + V0·|ψ|²,  E_L = ε·√λ  (ε plays ħc).
//
// The ether is opaque to light where it is dense (V0|ψ|², like a plasma
// frequency), so light lives in the void. The light pushes the ether out
// through ∂E_L/∂ψ* = (ε V0 / 2√λ)·χ²·ψ (Hellmann–Feynman; χ normalized).
// dV = 2πρ dρ dz. ρ is cell-centred, ρ_i = (i + ½)h; the outer walls are
// Neumann for ψ and Dirichlet for χ; m²/ρ² keeps χ off the axis.

// EtherGrid is the (ρ, z) grid.
type EtherGrid struct {
	NR, NZ int
	H, LZ  float64 // spacing; z runs over (−LZ, LZ)
	Rho, Z []float64
}

// NewEtherGrid makes ρ ∈ (0, lr), z ∈ (−lz, lz) with spacing h.
func NewEtherGrid(h, lr, lz float64) *EtherGrid {
	g := &EtherGrid{NR: int(lr / h), NZ: int(2 * lz / h), H: h, LZ: lz}
	for i := 0; i < g.NR; i++ {
		g.Rho = append(g.Rho, (float64(i)+0.5)*h)
	}
	for j := 0; j < g.NZ; j++ {
		g.Z = append(g.Z, (float64(j)+0.5)*h-lz)
	}
	return g
}

func (g *EtherGrid) at(i, j int) int { return i*g.NZ + j }

// DV is the cell volume 2πρh².
func (g *EtherGrid) DV(i int) float64 { return 2 * math.Pi * g.Rho[i] * g.H * g.H }

// LapC is the axisymmetric Laplacian of a complex field, Neumann walls.
func (g *EtherGrid) LapC(f, out []complex128) {
	h2 := g.H * g.H
	for i := 0; i < g.NR; i++ {
		rp, rm, r := g.Rho[i]+g.H/2, g.Rho[i]-g.H/2, g.Rho[i]
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			c := f[k]
			ip, im := c, c
			if i+1 < g.NR {
				ip = f[g.at(i+1, j)]
			}
			if i > 0 {
				im = f[g.at(i-1, j)]
			}
			jp, jm := c, c
			if j+1 < g.NZ {
				jp = f[k+1]
			}
			if j > 0 {
				jm = f[k-1]
			}
			out[k] = (complex(rp, 0)*(ip-c)-complex(rm, 0)*(c-im))/complex(r*h2, 0) + (jp-2*c+jm)/complex(h2, 0)
		}
	}
}

// LightOp applies the ρ-weighted light operator A x = ρ·(H − σ)x (symmetric):
// stiffness with Dirichlet walls, ρ(m²/ρ² + V0·dens − σ) on the diagonal.
func (g *EtherGrid) LightOp(x, out, dens []float64, v0, sigma float64) {
	h2 := g.H * g.H
	for i := 0; i < g.NR; i++ {
		rp, rm, r := g.Rho[i]+g.H/2, g.Rho[i]-g.H/2, g.Rho[i]
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			c := x[k]
			s := (rp + rm + 2*r) / h2 * c
			if i+1 < g.NR {
				s -= rp / h2 * x[g.at(i+1, j)]
			}
			if i > 0 {
				s -= rm / h2 * x[g.at(i-1, j)]
			}
			if j+1 < g.NZ {
				s -= r / h2 * x[k+1]
			}
			if j > 0 {
				s -= r / h2 * x[k-1]
			}
			s += r * (0.25/(r*r) + v0*dens[k] - sigma) * c
			out[k] = s
		}
	}
}

// Normalize scales χ to ∫χ² dV = 1.
func (g *EtherGrid) Normalize(chi []float64) {
	s := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			c := chi[g.at(i, j)]
			s += c * c * dv
		}
	}
	s = math.Sqrt(s)
	for k := range chi {
		chi[k] /= s
	}
}

// Rayleigh is λ = ⟨χ|H|χ⟩ for normalized χ.
func (g *EtherGrid) Rayleigh(chi, dens []float64, v0 float64) float64 {
	out := make([]float64, len(chi))
	g.LightOp(chi, out, dens, v0, 0)
	s := 0.0
	for k := range chi {
		s += chi[k] * out[k]
	}
	return s * 2 * math.Pi * g.H * g.H // the operator carries ρ; dV carries 2πρh²
}

// LightMode refines the lowest mode by inverse iteration with shift
// σ < λ: each step solves (A − σB)y = Bχ by conjugate gradients (A − σB is
// positive definite there). Returns λ; chi is updated in place.
func (g *EtherGrid) LightMode(chi, dens []float64, v0, sigma float64, iters int) float64 {
	n := len(chi)
	b := make([]float64, n)
	y := make([]float64, n)
	r := make([]float64, n)
	p := make([]float64, n)
	ap := make([]float64, n)
	for it := 0; it < iters; it++ {
		for i := 0; i < g.NR; i++ {
			for j := 0; j < g.NZ; j++ {
				k := g.at(i, j)
				b[k] = g.Rho[i] * chi[k]
			}
		}
		copy(y, chi)
		g.LightOp(y, ap, dens, v0, sigma)
		rr := 0.0
		for k := range r {
			r[k] = b[k] - ap[k]
			p[k] = r[k]
			rr += r[k] * r[k]
		}
		bb := 0.0
		for k := range b {
			bb += b[k] * b[k]
		}
		for cg := 0; cg < 2000 && rr > 1e-22*bb; cg++ {
			g.LightOp(p, ap, dens, v0, sigma)
			pap := 0.0
			for k := range p {
				pap += p[k] * ap[k]
			}
			a := rr / pap
			nr := 0.0
			for k := range y {
				y[k] += a * p[k]
				r[k] -= a * ap[k]
				nr += r[k] * r[k]
			}
			bet := nr / rr
			rr = nr
			for k := range p {
				p[k] = r[k] + bet*p[k]
			}
		}
		copy(chi, y)
		g.Normalize(chi)
	}
	return g.Rayleigh(chi, dens, v0)
}

// EtherState is the two fields with their parameters.
type EtherState struct {
	G        *EtherGrid
	Psi      []complex128
	Chi      []float64
	Eps, V0  float64
	Lam      float64
	dens, lp []float64
	lap      []complex128
	k        [][]complex128
	inVoid   []bool
}

// NewRingState starts a vortex ring of radius r0 with an emptied core of
// radius a0 (r0 = a0 is the horn torus) and the light's mode inside it.
func NewRingState(g *EtherGrid, eps, v0, r0, a0 float64) *EtherState {
	n := g.NR * g.NZ
	s := &EtherState{G: g, Eps: eps, V0: v0, Psi: make([]complex128, n), Chi: make([]float64, n),
		dens: make([]float64, n), lap: make([]complex128, n)}
	for i, r := range g.Rho {
		for j, z := range g.Z {
			d1 := math.Hypot(r-r0, z)
			d2 := math.Hypot(r+r0, z)
			ph := complex(r-r0, z) / complex(math.Max(d1, 1e-12), 0) * complex(r+r0, -z) / complex(math.Max(d2, 1e-12), 0)
			amp := math.Tanh(math.Max(d1-a0, 0) + 1e-3)
			k := g.at(i, j)
			s.Psi[k] = ph * complex(amp, 0)
			if d1 < a0 {
				s.Chi[k] = math.Sqrt(r) * math.Cos(0.5*math.Pi*d1/a0)
			}
			s.Chi[k] += 1e-9
		}
	}
	g.Normalize(s.Chi)
	s.updateDens()
	s.Lam = g.Rayleigh(s.Chi, s.dens, v0)
	s.Lam = g.LightMode(s.Chi, s.dens, v0, 0, 6)
	return s
}

func (s *EtherState) updateDens() {
	for k, p := range s.Psi {
		s.dens[k] = real(p)*real(p) + imag(p)*imag(p)
	}
}

// Force is δE/δψ* (ether plus light).
func (s *EtherState) Force(psi, out []complex128, chi []float64, lam float64) {
	g := s.G
	g.LapC(psi, s.lap)
	c := s.Eps * s.V0 / (2 * math.Sqrt(lam))
	for k, p := range psi {
		d := real(p)*real(p) + imag(p)*imag(p)
		out[k] = -0.5*s.lap[k] + complex(0.5*(d-1)+c*chi[k]*chi[k], 0)*p
	}
}

// Relax runs dissipative (imaginary-time) flow: energy only goes down.
// The light is re-solved every `every` steps.
func (s *EtherState) Relax(steps, every int) {
	dt := 0.15 * s.G.H * s.G.H
	f := make([]complex128, len(s.Psi))
	for it := 0; it < steps; it++ {
		s.Force(s.Psi, f, s.Chi, s.Lam)
		for k := range s.Psi {
			s.Psi[k] -= complex(dt, 0) * f[k]
		}
		if it%every == every-1 {
			s.updateDens()
			s.Lam = s.G.LightMode(s.Chi, s.dens, s.V0, 0.9*s.Lam, 2)
		}
	}
	s.updateDens()
	s.Lam = s.G.LightMode(s.Chi, s.dens, s.V0, 0.9*s.Lam, 3)
}

// Step advances real time by dt with no friction: classical RK4 with the
// light frozen over the step and re-solved after it (Born–Oppenheimer: in
// the electron's scaling the light is ~10⁸ times faster than the ether's
// sound).
func (s *EtherState) Step(dt float64) {
	n := len(s.Psi)
	if s.k == nil {
		s.k = make([][]complex128, 5)
		for i := range s.k {
			s.k[i] = make([]complex128, n)
		}
	}
	k1, k2, k3, k4, tmp := s.k[0], s.k[1], s.k[2], s.k[3], s.k[4]
	eval := func(p, out []complex128) {
		s.Force(p, out, s.Chi, s.Lam)
		for i := range out {
			out[i] *= -1i
		}
	}
	eval(s.Psi, k1)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(0.5*dt, 0)*k1[i]
	}
	eval(tmp, k2)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(0.5*dt, 0)*k2[i]
	}
	eval(tmp, k3)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(dt, 0)*k3[i]
	}
	eval(tmp, k4)
	for i := range s.Psi {
		s.Psi[i] += complex(dt/6, 0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
	}
	s.updateDens()
	s.Lam = s.G.LightMode(s.Chi, s.dens, s.V0, 0.95*s.Lam, 1)
}

// Energies returns the ether energy, the light energy and the impulse ∫Im(ψ*∂zψ)dV.
func (s *EtherState) Energies() (ee, el, pz float64) {
	g := s.G
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			p := s.Psi[k]
			var gr, gz complex128
			if i+1 < g.NR {
				gr = (s.Psi[g.at(i+1, j)] - p) / complex(g.H, 0)
			}
			if j+1 < g.NZ {
				gz = (s.Psi[k+1] - p) / complex(g.H, 0)
			}
			d := real(p)*real(p) + imag(p)*imag(p)
			ee += (0.5*(sq(gr)+sq(gz)) + 0.25*(d-1)*(d-1)) * dv
			if j > 0 && j+1 < g.NZ {
				dzp := (s.Psi[k+1] - s.Psi[k-1]) / complex(2*g.H, 0)
				pz += imag(cmplx.Conj(p)*dzp) * dv
			}
		}
	}
	return ee, s.Eps * math.Sqrt(s.Lam), pz
}

func sq(c complex128) float64 { return real(c)*real(c) + imag(c)*imag(c) }

// VoidShape describes the region |ψ|² < ½.
type VoidShape struct {
	Any                bool
	Inner, Outer       float64 // ρ extent
	ZC, Height, Volume float64
	OnAxis             bool    // the void reaches the axis: the hole is closed
	Winding            float64 // phase winding round the whole void (1: the vortex survives)
	LightIn            float64 // share of the light inside the void
}

// Shape measures the void that holds the light: the connected region of
// |ψ|² < ½ containing the light's brightest cell (sound waves can dip
// below ½ elsewhere; they are not the void).
func (s *EtherState) Shape() VoidShape {
	g := s.G
	v := VoidShape{Inner: math.Inf(1)}
	best, kb := -1.0, -1
	for k, c := range s.Chi {
		if c*c > best && s.dens[k] < 0.5 {
			best, kb = c*c, k
		}
	}
	if kb < 0 {
		return v
	}
	in := make([]bool, len(s.dens))
	stack := []int{kb}
	in[kb] = true
	for len(stack) > 0 {
		k := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		i, j := k/g.NZ, k%g.NZ
		for _, nb := range [][2]int{{i + 1, j}, {i - 1, j}, {i, j + 1}, {i, j - 1}} {
			if nb[0] < 0 || nb[1] < 0 || nb[0] >= g.NR || nb[1] >= g.NZ {
				continue
			}
			q := g.at(nb[0], nb[1])
			if !in[q] && s.dens[q] < 0.5 {
				in[q] = true
				stack = append(stack, q)
			}
		}
	}
	s.inVoid = in
	zmin, zmax := math.Inf(1), math.Inf(-1)
	zs := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			if !in[k] {
				continue
			}
			v.Any = true
			v.Inner = math.Min(v.Inner, g.Rho[i])
			v.Outer = math.Max(v.Outer, g.Rho[i])
			zmin, zmax = math.Min(zmin, g.Z[j]), math.Max(zmax, g.Z[j])
			v.Volume += dv
			zs += g.Z[j] * dv
			v.LightIn += s.Chi[k] * s.Chi[k] * dv
			if i == 0 {
				v.OnAxis = true
			}
		}
	}
	v.ZC = zs / v.Volume
	v.Height = zmax - zmin + g.H
	v.Winding = s.windingAround(math.Min(v.Outer+4, g.Rho[g.NR-1]-g.H), v.ZC, v.Height/2+4)
	return v
}

// windingAround counts the phase winding along the boundary of the
// meridional box ρ ∈ (0, b), z ∈ zc ± c.
func (s *EtherState) windingAround(b, zc, c float64) float64 {
	g := s.G
	i1 := int(b / g.H)
	if i1 > g.NR-1 {
		i1 = g.NR - 1
	}
	j0 := int((zc - c + g.LZ) / g.H)
	j1 := int((zc + c + g.LZ) / g.H)
	if j0 < 1 {
		j0 = 1
	}
	if j1 > g.NZ-2 {
		j1 = g.NZ - 2
	}
	var path []complex128
	for j := j1; j > j0; j-- {
		path = append(path, s.Psi[g.at(0, j)])
	}
	for i := 0; i < i1; i++ {
		path = append(path, s.Psi[g.at(i, j0)])
	}
	for j := j0; j < j1; j++ {
		path = append(path, s.Psi[g.at(i1, j)])
	}
	for i := i1; i > 0; i-- {
		path = append(path, s.Psi[g.at(i, j1)])
	}
	tot := 0.0
	for k := range path {
		d := cmplx.Phase(path[(k+1)%len(path)]) - cmplx.Phase(path[k])
		d = math.Mod(d+3*math.Pi, 2*math.Pi) - math.Pi
		tot += d
	}
	return math.Round(100*tot/(2*math.Pi)) / 100
}

// WallInvDGrid is ⟨1/d⟩ for charge spread evenly over the void's wall
// (cells with |ψ|² < ½ next to one with |ψ|² ≥ ½), each wall cell a ring
// round the axis weighted by its ρ; rings interact through 4K(m)/√(...)
// averaged over 2π, the same-cell term from a ring offset by h/4.
func (s *EtherState) WallInvDGrid() float64 {
	g := s.G
	type pt struct{ r, z, w float64 }
	var pts []pt
	for i := 0; i < g.NR; i++ {
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			if s.inVoid == nil || !s.inVoid[k] {
				continue
			}
			edge := false
			for _, nb := range [][2]int{{i + 1, j}, {i - 1, j}, {i, j + 1}, {i, j - 1}} {
				if nb[0] < 0 || nb[1] < 0 || nb[0] >= g.NR || nb[1] >= g.NZ {
					continue
				}
				if s.dens[g.at(nb[0], nb[1])] >= 0.5 {
					edge = true
				}
			}
			if edge {
				pts = append(pts, pt{g.Rho[i], g.Z[j], g.Rho[i]})
			}
		}
	}
	wsum := 0.0
	for _, p := range pts {
		wsum += p.w
	}
	tot := 0.0
	for a, p := range pts {
		for b, q := range pts {
			r2, dz := q.r, q.z-p.z
			if a == b {
				dz = g.H / 4
			}
			den := (p.r+r2)*(p.r+r2) + dz*dz
			m := 4 * p.r * r2 / den
			if m >= 1 {
				m = 1 - 1e-15
			}
			tot += p.w * q.w * (2 / math.Pi) * ellipK(m) / math.Sqrt(den)
		}
	}
	return tot / (wsum * wsum)
}
