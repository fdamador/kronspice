// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// Plain vector calculus on the same grid (Gibbs–Heaviside form), written
// independently of the quaternion code so the triple-twin tests compare
// two implementations, not one. Same fourth-order stencil.

// VField is a 3-vector field.
type VField struct {
	G Grid
	V [][3]float64
}

// SampleV samples a vector field.
func SampleV(g Grid, f func(x, y, z float64) [3]float64) VField {
	out := VField{G: g, V: make([][3]float64, g.Len())}
	for n := range out.V {
		i, j, k := g.IJK(n)
		p := g.Point(i, j, k)
		out.V[n] = f(p[0], p[1], p[2])
	}
	return out
}

// ddx is ∂(component c of v)/∂x_a at (i, j, k), from explicit neighbours.
func (v VField) ddx(c, a, i, j, k int) float64 {
	g := v.G
	at := func(o int) float64 {
		ii, jj, kk := i, j, k
		switch a {
		case 0:
			ii += o
		case 1:
			jj += o
		default:
			kk += o
		}
		return v.V[g.Index(ii, jj, kk)][c]
	}
	return (at(-2) - 8*at(-1) + 8*at(1) - at(2)) / (12 * g.H[a])
}

// Div is ∇·v (NaN at the edges).
func Div(v VField) []float64 {
	g := v.G
	out := make([]float64, g.Len())
	for n := range out {
		if !g.Interior(n) {
			out[n] = math.NaN()
			continue
		}
		i, j, k := g.IJK(n)
		out[n] = v.ddx(0, 0, i, j, k) + v.ddx(1, 1, i, j, k) + v.ddx(2, 2, i, j, k)
	}
	return out
}

// Curl is ∇×v (NaN at the edges).
func Curl(v VField) VField {
	g := v.G
	out := VField{G: g, V: make([][3]float64, g.Len())}
	nan := math.NaN()
	for n := range out.V {
		if !g.Interior(n) {
			out.V[n] = [3]float64{nan, nan, nan}
			continue
		}
		i, j, k := g.IJK(n)
		out.V[n] = [3]float64{
			v.ddx(2, 1, i, j, k) - v.ddx(1, 2, i, j, k),
			v.ddx(0, 2, i, j, k) - v.ddx(2, 0, i, j, k),
			v.ddx(1, 0, i, j, k) - v.ddx(0, 1, i, j, k),
		}
	}
	return out
}

// Grad is ∇ψ of a scalar field (NaN at the edges).
func Grad(g Grid, psi []float64) VField {
	s := VField{G: g, V: make([][3]float64, g.Len())}
	for n, p := range psi {
		s.V[n] = [3]float64{p, 0, 0}
	}
	out := VField{G: g, V: make([][3]float64, g.Len())}
	nan := math.NaN()
	for n := range out.V {
		if !g.Interior(n) {
			out.V[n] = [3]float64{nan, nan, nan}
			continue
		}
		i, j, k := g.IJK(n)
		out.V[n] = [3]float64{s.ddx(0, 0, i, j, k), s.ddx(0, 1, i, j, k), s.ddx(0, 2, i, j, k)}
	}
	return out
}

// VectorResidual holds the residuals of Maxwell's four equations in vector
// form at each point: Gauss div E − ρ/ε₀ [V/m²], Faraday curl E + ∂B/∂t
// [T/s], div B [T/m], Ampère–Maxwell curl B − μ₀J − ∂E/∂t/c² [T/m].
type VectorResidual struct {
	Gauss   []float64
	Faraday [][3]float64
	DivB    []float64
	Ampere  [][3]float64
}

// MaxwellVector evaluates the residuals of Maxwell's equations in vector
// form for fields E, B, their time derivatives, charge density rho and
// current density J (rho or J may be nil for 0).
func MaxwellVector(E, B, dEdt, dBdt VField, rho []float64, J *VField) VectorResidual {
	g := E.G
	dE, cE := Div(E), Curl(E)
	dB, cB := Div(B), Curl(B)
	r := VectorResidual{Gauss: make([]float64, g.Len()), Faraday: make([][3]float64, g.Len()),
		DivB: dB, Ampere: make([][3]float64, g.Len())}
	c2 := SpeedOfLight * SpeedOfLight
	for n := range r.Gauss {
		q := 0.0
		if rho != nil {
			q = rho[n]
		}
		r.Gauss[n] = dE[n] - q/Eps0
		for a := 0; a < 3; a++ {
			j := 0.0
			if J != nil {
				j = J.V[n][a]
			}
			r.Faraday[n][a] = cE.V[n][a] + dBdt.V[n][a]
			r.Ampere[n][a] = cB.V[n][a] - Mu0*j - dEdt.V[n][a]/c2
		}
	}
	return r
}
