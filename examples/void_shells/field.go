// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"sort"
)

// The two fields, on an axisymmetric (ρ, z) grid in healing-length units
// (the solver of examples/ether_horizon Part 13, copied so the examples
// stay independent, with several light modes):
//
//	ether ψ: E_e = ∫ ½|∇ψ|² + ¼(|ψ|² − 1)² dV  (Gross–Pitaevskii);
//	light: double-loop modes χ_q·e^{iφ/2} (m = ½), eigenmodes of
//	       H = −∇² + m²/ρ² + V0·|ψ|²,  each occupied mode holds one photon
//	       of energy ε√λ_q (ε plays ħc); spin ½ photons, so one per mode.
//
// Light pushes ether: ∂E_L/∂ψ* = Σ_q (εV0/2√λ_q)·χ_q²·ψ (Hellmann–Feynman).
// dV = 2πρ dρ dz, ρ_i = (i + ½)h; outer walls Neumann for ψ, Dirichlet for χ.

// Grid is the (ρ, z) grid.
type Grid struct {
	NR, NZ int
	H, LZ  float64
	Rho, Z []float64
}

// NewGrid makes ρ ∈ (0, lr), z ∈ (−lz, lz).
func NewGrid(h, lr, lz float64) *Grid {
	g := &Grid{NR: int(lr / h), NZ: int(2 * lz / h), H: h, LZ: lz}
	for i := 0; i < g.NR; i++ {
		g.Rho = append(g.Rho, (float64(i)+0.5)*h)
	}
	for j := 0; j < g.NZ; j++ {
		g.Z = append(g.Z, (float64(j)+0.5)*h-lz)
	}
	return g
}

func (g *Grid) at(i, j int) int { return i*g.NZ + j }

// DV is the cell volume 2πρh².
func (g *Grid) DV(i int) float64 { return 2 * math.Pi * g.Rho[i] * g.H * g.H }

// lapC is the axisymmetric Laplacian, Neumann walls.
func (g *Grid) lapC(f, out []complex128) {
	h2 := g.H * g.H
	for i := 0; i < g.NR; i++ {
		rp, rm, r := g.Rho[i]+g.H/2, g.Rho[i]-g.H/2, g.Rho[i]
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			c := f[k]
			ip, im, jp, jm := c, c, c, c
			if i+1 < g.NR {
				ip = f[k+g.NZ]
			}
			if i > 0 {
				im = f[k-g.NZ]
			}
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

// lightOp applies A = ρ·H (symmetric, positive definite): Dirichlet
// stiffness plus ρ(m²/ρ² + V0·dens).
func (g *Grid) lightOp(x, out, dens []float64, v0 float64) {
	h2 := g.H * g.H
	for i := 0; i < g.NR; i++ {
		rp, rm, r := g.Rho[i]+g.H/2, g.Rho[i]-g.H/2, g.Rho[i]
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			c := x[k]
			s := (rp + rm + 2*r) / h2 * c
			if i+1 < g.NR {
				s -= rp / h2 * x[k+g.NZ]
			}
			if i > 0 {
				s -= rm / h2 * x[k-g.NZ]
			}
			if j+1 < g.NZ {
				s -= r / h2 * x[k+1]
			}
			if j > 0 {
				s -= r / h2 * x[k-1]
			}
			out[k] = s + r*(0.25/(r*r)+v0*dens[k])*c
		}
	}
}

// bdot is ∫ x·y dV.
func (g *Grid) bdot(x, y []float64) float64 {
	s := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			s += x[k] * y[k] * dv
		}
	}
	return s
}

// solveCG solves A y = ρ·x by conjugate gradients from the guess y.
func (g *Grid) solveCG(x, y, dens []float64, v0 float64, ws [4][]float64) {
	b, r, p, ap := ws[0], ws[1], ws[2], ws[3]
	for i := 0; i < g.NR; i++ {
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			b[k] = g.Rho[i] * x[k]
		}
	}
	g.lightOp(y, ap, dens, v0)
	rr, bb := 0.0, 0.0
	for k := range r {
		r[k] = b[k] - ap[k]
		p[k] = r[k]
		rr += r[k] * r[k]
		bb += b[k] * b[k]
	}
	for it := 0; it < 3000 && rr > 1e-22*bb; it++ {
		g.lightOp(p, ap, dens, v0)
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
}

// LightModes refines the lowest len(chi) modes in place by subspace
// (inverse) iteration with Rayleigh–Ritz, and returns their eigenvalues in
// ascending order. Modes are normalized, ∫χ² dV = 1.
func (g *Grid) LightModes(chi [][]float64, dens []float64, v0 float64, iters int) []float64 {
	n := len(chi[0])
	nm := len(chi)
	var ws [4][]float64
	for i := range ws {
		ws[i] = make([]float64, n)
	}
	for it := 0; it < iters; it++ {
		for q := range chi {
			y := append([]float64(nil), chi[q]...)
			g.solveCG(chi[q], y, dens, v0, ws)
			chi[q] = y
		}
		g.orthonormalize(chi)
	}
	// Rayleigh–Ritz in the subspace
	hc := make([][]float64, nm)
	for q := range chi {
		hc[q] = make([]float64, n)
		g.lightOp(chi[q], hc[q], dens, v0)
	}
	a := make([][]float64, nm)
	for p := range a {
		a[p] = make([]float64, nm)
		for q := range a[p] {
			s := 0.0
			for k := range hc[q] {
				s += chi[p][k] * hc[q][k]
			}
			a[p][q] = s * 2 * math.Pi * g.H * g.H
		}
	}
	vals, vecs := jacobiEigen(a)
	rot := make([][]float64, nm)
	for q := 0; q < nm; q++ {
		rot[q] = make([]float64, n)
		for p := 0; p < nm; p++ {
			c := vecs[p][q]
			for k := range rot[q] {
				rot[q][k] += c * chi[p][k]
			}
		}
	}
	for q := range chi {
		chi[q] = rot[q]
	}
	g.orthonormalize(chi)
	return vals
}

func (g *Grid) orthonormalize(chi [][]float64) {
	for q := range chi {
		for p := 0; p < q; p++ {
			c := g.bdot(chi[p], chi[q])
			for k := range chi[q] {
				chi[q][k] -= c * chi[p][k]
			}
		}
		nn := math.Sqrt(g.bdot(chi[q], chi[q]))
		for k := range chi[q] {
			chi[q][k] /= nn
		}
	}
}

// jacobiEigen diagonalizes a small symmetric matrix: eigenvalues ascending,
// eigenvectors as the columns of vecs.
func jacobiEigen(a0 [][]float64) ([]float64, [][]float64) {
	n := len(a0)
	a := make([][]float64, n)
	v := make([][]float64, n)
	for i := range a {
		a[i] = append([]float64(nil), a0[i]...)
		v[i] = make([]float64, n)
		v[i][i] = 1
	}
	for sweep := 0; sweep < 100; sweep++ {
		off := 0.0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(a[p][q]) < 1e-300 {
					continue
				}
				th := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := math.Copysign(1, th) / (math.Abs(th) + math.Sqrt(th*th+1))
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p], a[k][q] = c*akp-s*akq, s*akp+c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k], a[q][k] = c*apk-s*aqk, s*apk+c*aqk
				}
				for k := 0; k < n; k++ {
					vkp, vkq := v[k][p], v[k][q]
					v[k][p], v[k][q] = c*vkp-s*vkq, s*vkp+c*vkq
				}
			}
		}
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(x, y int) bool { return a[idx[x]][idx[x]] < a[idx[y]][idx[y]] })
	vals := make([]float64, n)
	vecs := make([][]float64, n)
	for i := range vecs {
		vecs[i] = make([]float64, n)
	}
	for c, o := range idx {
		vals[c] = a[o][o]
		for r := 0; r < n; r++ {
			vecs[r][c] = v[r][o]
		}
	}
	return vals, vecs
}

// State is the ether, the light modes and which modes hold a photon.
type State struct {
	G       *Grid
	Psi     []complex128
	Chi     [][]float64 // the lowest NModes modes, ascending
	Lam     []float64
	Occ     []int // occupied mode indices (one photon each)
	Eps, V0 float64
	dens    []float64
	lap     []complex128
	k       [5][]complex128
}

// NewState makes a state from an ether amplitude profile amp(ρ, z) (the
// phase is 1: no vortex) with photons in modes occ.
func NewState(g *Grid, eps, v0 float64, occ []int, amp func(r, z float64) float64) *State {
	n := g.NR * g.NZ
	s := &State{G: g, Eps: eps, V0: v0, Occ: occ, Psi: make([]complex128, n), dens: make([]float64, n), lap: make([]complex128, n)}
	for i, r := range g.Rho {
		for j, z := range g.Z {
			s.Psi[g.at(i, j)] = complex(amp(r, z), 0)
		}
	}
	nm := 0
	for _, q := range occ {
		if q+1 > nm {
			nm = q + 1
		}
	}
	nm++ // one spare mode keeps the subspace iteration honest at the top
	s.Chi = make([][]float64, nm)
	for q := range s.Chi {
		s.Chi[q] = make([]float64, n)
		for i, r := range g.Rho {
			for j, z := range g.Z {
				// smooth, distinct starting vectors: √ρ times low polynomials in z and ρ
				s.Chi[q][g.at(i, j)] = math.Sqrt(r) * math.Exp(-(r*r+z*z)/60) * math.Pow(z/5, float64(q%3)) * math.Pow(r/5, float64(q/3))
			}
		}
	}
	g.orthonormalize(s.Chi)
	s.updateDens()
	s.Lam = g.LightModes(s.Chi, s.dens, v0, 12)
	return s
}

// Ball is a void of radius r0 (ether density rising over one healing length).
func Ball(r0 float64) func(r, z float64) float64 {
	return func(r, z float64) float64 { return math.Tanh(math.Max(math.Hypot(r, z)-r0, 0) + 1e-3) }
}

// ShellAmp is an ether core of radius rc inside a void shell out to ro.
func ShellAmp(rc, ro float64) func(r, z float64) float64 {
	return func(r, z float64) float64 {
		d := math.Hypot(r, z)
		if d < rc {
			return math.Tanh(rc - d + 1e-3)
		}
		return math.Tanh(math.Max(d-ro, 0) + 1e-3)
	}
}

// OnionAmp is a void core r < r1, an ether layer r1 … r2, a void shell r2 … r3.
func OnionAmp(r1, r2, r3 float64) func(r, z float64) float64 {
	return func(r, z float64) float64 {
		d := math.Hypot(r, z)
		switch {
		case d < r1:
			return 0
		case d < r2:
			return math.Min(math.Tanh(d-r1+1e-3), math.Tanh(r2-d+1e-3))
		default:
			return math.Tanh(math.Max(d-r3, 0) + 1e-3)
		}
	}
}

func (s *State) updateDens() {
	for k, p := range s.Psi {
		s.dens[k] = real(p)*real(p) + imag(p)*imag(p)
	}
}

// force is δE/δψ*.
func (s *State) force(psi, out []complex128) {
	s.G.lapC(psi, s.lap)
	coup := make([]float64, len(psi))
	for _, q := range s.Occ {
		c := s.Eps * s.V0 / (2 * math.Sqrt(s.Lam[q]))
		for k, x := range s.Chi[q] {
			coup[k] += c * x * x
		}
	}
	for k, p := range psi {
		d := real(p)*real(p) + imag(p)*imag(p)
		out[k] = -0.5*s.lap[k] + complex(0.5*(d-1)+coup[k], 0)*p
	}
}

func (s *State) refreshLight(iters int) {
	s.updateDens()
	s.Lam = s.G.LightModes(s.Chi, s.dens, s.V0, iters)
}

// Relax runs dissipative flow (energy only falls), light refreshed every 25 steps.
func (s *State) Relax(steps int) {
	f := make([]complex128, len(s.Psi))
	dt := s.relaxDT()
	for it := 0; it < steps; it++ {
		if it%25 == 0 {
			dt = s.relaxDT()
		}
		s.force(s.Psi, f)
		for k := range s.Psi {
			s.Psi[k] -= complex(dt, 0) * f[k]
		}
		if it%25 == 24 {
			s.refreshLight(1)
		}
	}
	s.refreshLight(4)
}

// relaxDT is the stable step of the dissipative flow: the stencil limits it
// to ~0.15h², the light's push (εV0/2√λ)·χ² to below 1/push.
func (s *State) relaxDT() float64 {
	dt := 0.15 * s.G.H * s.G.H
	pmax := 0.0
	for k := range s.Psi {
		p := 0.0
		for _, q := range s.Occ {
			p += s.Eps * s.V0 / (2 * math.Sqrt(s.Lam[q])) * s.Chi[q][k] * s.Chi[q][k]
		}
		pmax = math.Max(pmax, p)
	}
	if pmax*dt > 0.5 {
		dt = 0.5 / pmax
	}
	return dt
}

// Step advances real time (no friction) by RK4 with the light frozen over
// the step and refreshed after it.
func (s *State) Step(dt float64) {
	n := len(s.Psi)
	for i := range s.k {
		if s.k[i] == nil {
			s.k[i] = make([]complex128, n)
		}
	}
	k1, k2, k3, k4, tmp := s.k[0], s.k[1], s.k[2], s.k[3], s.k[4]
	ev := func(p, out []complex128) {
		s.force(p, out)
		for i := range out {
			out[i] *= -1i
		}
	}
	ev(s.Psi, k1)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(0.5*dt, 0)*k1[i]
	}
	ev(tmp, k2)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(0.5*dt, 0)*k2[i]
	}
	ev(tmp, k3)
	for i := range tmp {
		tmp[i] = s.Psi[i] + complex(dt, 0)*k3[i]
	}
	ev(tmp, k4)
	for i := range s.Psi {
		s.Psi[i] += complex(dt/6, 0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
	}
	s.refreshLight(1)
}

// Energies: ether, light (occupied modes).
func (s *State) Energies() (ee, el float64) {
	g := s.G
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			p := s.Psi[k]
			var gr, gz complex128
			if i+1 < g.NR {
				gr = (s.Psi[k+g.NZ] - p) / complex(g.H, 0)
			}
			if j+1 < g.NZ {
				gz = (s.Psi[k+1] - p) / complex(g.H, 0)
			}
			d := real(p)*real(p) + imag(p)*imag(p)
			ee += (0.5*(sq(gr)+sq(gz)) + 0.25*(d-1)*(d-1)) * dv
		}
	}
	for _, q := range s.Occ {
		el += s.Eps * math.Sqrt(s.Lam[q])
	}
	return
}

func sq(c complex128) float64 { return real(c)*real(c) + imag(c)*imag(c) }

// Layout describes the voids (|ψ|² < ½).
type Layout struct {
	Voids      int     // separate void regions (4-connected in the cut)
	CoreDens   float64 // ether density at the centre
	REq, RAxis float64 // void's reach along ρ at z = 0 and along the axis
	LightIn    float64 // share of the occupied light inside voids
	Volume     float64
}

// Measure counts the voids and their extent.
func (s *State) Measure() Layout {
	g := s.G
	var L Layout
	lab := make([]int, len(s.dens))
	for k := range s.dens {
		if s.dens[k] >= 0.5 || lab[k] != 0 {
			continue
		}
		L.Voids++
		stack := []int{k}
		lab[k] = L.Voids
		for len(stack) > 0 {
			c := stack[len(stack)-1]
			stack = stack[:len(stack)-1]
			i, j := c/g.NZ, c%g.NZ
			for _, nb := range [][2]int{{i + 1, j}, {i - 1, j}, {i, j + 1}, {i, j - 1}} {
				if nb[0] < 0 || nb[1] < 0 || nb[0] >= g.NR || nb[1] >= g.NZ {
					continue
				}
				q := g.at(nb[0], nb[1])
				if lab[q] == 0 && s.dens[q] < 0.5 {
					lab[q] = L.Voids
					stack = append(stack, q)
				}
			}
		}
	}
	jc := g.NZ / 2
	L.CoreDens = 0.25 * (s.dens[g.at(0, jc)] + s.dens[g.at(0, jc-1)] + s.dens[g.at(1, jc)] + s.dens[g.at(1, jc-1)])
	tot := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			w := 0.0
			for _, q := range s.Occ {
				w += s.Chi[q][k] * s.Chi[q][k]
			}
			tot += w * dv
			if s.dens[k] < 0.5 {
				L.LightIn += w * dv
				L.Volume += dv
				if j == jc || j == jc-1 {
					L.REq = math.Max(L.REq, g.Rho[i]+g.H/2)
				}
				if i == 0 {
					L.RAxis = math.Max(L.RAxis, math.Abs(g.Z[j])+g.H/2)
				}
			}
		}
	}
	L.LightIn /= tot
	return L
}

// G0 is g = Mc²/E_L for the whole object (ether energy plus light).
func (s *State) G0() float64 {
	ee, el := s.Energies()
	return (ee + el) / el
}
