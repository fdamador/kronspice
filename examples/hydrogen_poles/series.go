// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Series is the exact field of an axisymmetric azimuthal current density
// J_φ(r, θ), from the Legendre expansion of the Green's function
// (Jackson §5.5, eq. 5.33 with the addition theorem):
//
//	A_φ(r,θ) = (μ0/2) Σ_l P_l¹(cosθ)/(l(l+1)) [r^−(l+1) In_l(r) + r^l Out_l(r)],
//	In_l(r)  = ∫_0^r j_l(r') r'^(l+2) dr',  Out_l(r) = ∫_r^∞ j_l(r') r'^(1−l) dr',
//	j_l(r')  = ∫_0^π J_φ(r', θ') P_l¹(cosθ') sinθ' dθ'.
//
// B_r = −(μ0/2) Σ_l P_l(cosθ) G_l(r)/r with G_l = r^−(l+1) In_l + r^l Out_l.
// Outside all sources the potential is Σ_l (μ0/4π) m_l (l+1) P_l/r^(l+2)
// in B_r, with the 2^l-pole moments m_l = −2π In_l(∞)/(l+1) (m_1 the
// dipole moment, m_3 the octupole). Fields are returned in B0 = μ0μB/(4πa0³).
type Series struct {
	Ls    []int
	Dr    float64
	Edges []float64   // radial cell edges
	In    [][]float64 // In[l index][edge]
	Out   [][]float64
	Jl    [][]float64 // j_l at cell centres
}

// NewSeries projects J onto P_l¹ for each l in ls, on n radial cells of
// width dr (midpoint rule) and nt Gauss–Legendre points in θ.
func NewSeries(J func(r, th float64) float64, ls []int, dr float64, n, nt int) *Series {
	s := &Series{Ls: ls, Dr: dr}
	x, w := GaussLegendre(nt)
	jl := make([][]float64, len(ls))
	for k := range jl {
		jl[k] = make([]float64, n)
	}
	parallel(n, func(i int) {
		r := (float64(i) + 0.5) * dr
		for q := range x {
			th := math.Acos(x[q])
			j := J(r, th) * w[q] // sinθ dθ = dx
			for k, l := range ls {
				jl[k][i] += j * Plm1(l, x[q])
			}
		}
	})
	s.Jl = jl
	s.Edges = make([]float64, n+1)
	for i := range s.Edges {
		s.Edges[i] = float64(i) * dr
	}
	s.In, s.Out = make([][]float64, len(ls)), make([][]float64, len(ls))
	for k, l := range ls {
		in, out := make([]float64, n+1), make([]float64, n+1)
		for i := 0; i < n; i++ {
			r := (float64(i) + 0.5) * dr
			in[i+1] = in[i] + jl[k][i]*math.Pow(r, float64(l+2))*dr
		}
		for i := n - 1; i >= 0; i-- {
			r := (float64(i) + 0.5) * dr
			out[i] = out[i+1] + jl[k][i]*math.Pow(r, float64(1-l))*dr
		}
		s.In[k], s.Out[k] = in, out
	}
	return s
}

func (s *Series) at(k int, r float64) (in, out float64) {
	f := r / s.Dr
	i := int(f)
	if i >= len(s.Edges)-1 {
		n := len(s.Edges) - 1
		return s.In[k][n], 0
	}
	t := f - float64(i)
	return s.In[k][i]*(1-t) + s.In[k][i+1]*t, s.Out[k][i]*(1-t) + s.Out[k][i+1]*t
}

// Moment is the 2^l-pole moment m_l = −2π In_l(∞)/(l+1) [μB a0^(l−1)].
func (s *Series) Moment(l int) float64 {
	for k, ll := range s.Ls {
		if ll == l {
			return -2 * math.Pi * s.In[k][len(s.Edges)-1] / float64(l+1)
		}
	}
	return 0
}

// Br is B_r at (r, θ) in B0.
func (s *Series) Br(r, th float64) float64 {
	x := math.Cos(th)
	var b float64
	for k, l := range s.Ls {
		in, out := s.at(k, r)
		g := in/math.Pow(r, float64(l+1)) + out*math.Pow(r, float64(l))
		b += Pl(l, x) * g
	}
	// −(μ0/2) Σ … /r in units of μ0μB/(4π a0³): × 4π
	return -2 * math.Pi * b / r
}

// Flux is the flux function ψ = r sinθ A_φ (in μ0μB/(4π a0)); field lines are
// its contours.
func (s *Series) Flux(r, th float64) float64 {
	x := math.Cos(th)
	var a float64
	for k, l := range s.Ls {
		in, out := s.at(k, r)
		g := in/math.Pow(r, float64(l+1)) + out*math.Pow(r, float64(l))
		a += Plm1(l, x) / float64(l*(l+1)) * g
	}
	return 2 * math.Pi * a * r * math.Sin(th)
}

// Pl is the Legendre polynomial P_l(x).
func Pl(l int, x float64) float64 {
	p0, p1 := 1.0, x
	if l == 0 {
		return p0
	}
	for n := 1; n < l; n++ {
		p0, p1 = p1, (float64(2*n+1)*x*p1-float64(n)*p0)/float64(n+1)
	}
	return p1
}

// Plm1 is the associated Legendre function P_l¹(x) = −(1−x²)^½ P_l'(x)
// (Condon–Shortley phase).
func Plm1(l int, x float64) float64 {
	if l == 0 {
		return 0
	}
	// P_l'(x) = l (x P_l − P_{l−1})/(x² − 1)
	s := math.Sqrt(math.Max(0, 1-x*x))
	if s < 1e-12 {
		return 0
	}
	d := float64(l) * (Pl(l-1, x) - x*Pl(l, x)) / (1 - x*x)
	return -s * d
}

// GaussLegendre returns n nodes and weights on [−1, 1] (Newton on P_n).
func GaussLegendre(n int) (x, w []float64) {
	x, w = make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		for it := 0; it < 100; it++ {
			p := Pl(n, z)
			dp := float64(n) * (z*p - Pl(n-1, z)) / (z*z - 1)
			dz := p / dp
			z -= dz
			if math.Abs(dz) < 1e-15 {
				break
			}
		}
		dp := float64(n) * (z*Pl(n, z) - Pl(n-1, z)) / (z*z - 1)
		x[i], w[i] = z, 2/((1-z*z)*dp*dp)
	}
	return
}

// SignChanges counts sign changes of v.
func SignChanges(v []float64) int {
	n := 0
	for i := 1; i < len(v); i++ {
		if v[i-1]*v[i] < 0 {
			n++
		}
	}
	return n
}

// Bt is B_θ at (r, θ) in B0:
// B_θ = −(μ0/2) Σ_l P_l¹/(l(l+1)) (1/r) d(r G_l)/dr,
// d(r G_l)/dr = −l In_l/r^(l+1) + (l+1) r^l Out_l.
func (s *Series) Bt(r, th float64) float64 {
	x := math.Cos(th)
	var b float64
	for k, l := range s.Ls {
		in, out := s.at(k, r)
		d := -float64(l)*in/math.Pow(r, float64(l+1)) + float64(l+1)*out*math.Pow(r, float64(l))
		b += Plm1(l, x) / float64(l*(l+1)) * d
	}
	return -2 * math.Pi * b / r
}
