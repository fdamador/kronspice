// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"sort"
)

// Dirac's equation, solved numerically. Natural units ħ = c = m = 1
// (energies in m c², lengths in the reduced Compton wavelength ƛ = ħ/mc).
//
//	H = α·p + β,  α_i = [[0, σ_i], [σ_i, 0]],  β = diag(1, 1, −1, −1)
//
// (P. A. M. Dirac, Proc. R. Soc. A 117, 610 (1928); Bjorken & Drell §1–3).

type cmat [][]complex128

func zeros(n int) cmat {
	m := make(cmat, n)
	for i := range m {
		m[i] = make([]complex128, n)
	}
	return m
}

var pauli = [3][2][2]complex128{
	{{0, 1}, {1, 0}},
	{{0, -1i}, {1i, 0}},
	{{1, 0}, {0, -1}},
}

// FreeDirac is the 4×4 Hamiltonian at momentum p.
func FreeDirac(p [3]float64) cmat {
	h := zeros(4)
	for k := 0; k < 3; k++ {
		for a := 0; a < 2; a++ {
			for b := 0; b < 2; b++ {
				h[a][2+b] += complex(p[k], 0) * pauli[k][a][b]
				h[2+a][b] += complex(p[k], 0) * pauli[k][a][b]
			}
		}
	}
	h[0][0], h[1][1], h[2][2], h[3][3] = 1, 1, -1, -1
	return h
}

// SpinOp is Σ_k/2 = ½ diag(σ_k, σ_k).
func SpinOp(k int) cmat {
	s := zeros(4)
	for a := 0; a < 2; a++ {
		for b := 0; b < 2; b++ {
			s[a][b] = 0.5 * pauli[k][a][b]
			s[2+a][2+b] = 0.5 * pauli[k][a][b]
		}
	}
	return s
}

// HermEigen diagonalizes a Hermitian matrix through its real symmetric
// embedding [[A, −B], [B, A]] (each eigenvalue appears twice there) and
// returns the eigenvalues ascending with one complex eigenvector each.
func HermEigen(h cmat) ([]float64, [][]complex128) {
	n := len(h)
	r := make([][]float64, 2*n)
	for i := range r {
		r[i] = make([]float64, 2*n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			a, b := real(h[i][j]), imag(h[i][j])
			r[i][j], r[i][n+j], r[n+i][j], r[n+i][n+j] = a, -b, b, a
		}
	}
	vals, vecs := jacobiEigen(r)
	// take every second eigenpair (the pairs are degenerate copies), orthonormalized
	var outV []float64
	var outX [][]complex128
	for c := 0; c < 2*n; c++ {
		v := make([]complex128, n)
		for i := 0; i < n; i++ {
			v[i] = complex(vecs[i][c], vecs[n+i][c])
		}
		for _, u := range outX {
			ov := cdot(u, v)
			for i := range v {
				v[i] -= ov * u[i]
			}
		}
		nn := math.Sqrt(real(cdot(v, v)))
		if nn < 1e-6 || len(outX) == n {
			continue
		}
		for i := range v {
			v[i] /= complex(nn, 0)
		}
		outV = append(outV, vals[c])
		outX = append(outX, v)
	}
	return outV, outX
}

func cdot(a, b []complex128) complex128 {
	var s complex128
	for i := range a {
		s += cmplx.Conj(a[i]) * b[i]
	}
	return s
}

func apply(m cmat, v []complex128) []complex128 {
	out := make([]complex128, len(v))
	for i := range m {
		for j := range v {
			out[i] += m[i][j] * v[j]
		}
	}
	return out
}

// HelicityOnPositive returns the eigenvalues of the helicity S·p̂ in the
// two positive-energy states at momentum p, and the energies found.
func HelicityOnPositive(p [3]float64) (energies, helicities []float64) {
	vals, vecs := HermEigen(FreeDirac(p))
	pn := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
	hel := zeros(4)
	for k := 0; k < 3; k++ {
		s := SpinOp(k)
		for i := range hel {
			for j := range hel {
				hel[i][j] += complex(p[k]/pn, 0) * s[i][j]
			}
		}
	}
	var pos [][]complex128
	for i, e := range vals {
		energies = append(energies, e)
		if e > 0 {
			pos = append(pos, vecs[i])
		}
	}
	m := zeros(2)
	for a := 0; a < 2; a++ {
		for b := 0; b < 2; b++ {
			m[a][b] = cdot(pos[a], apply(hel, pos[b]))
		}
	}
	hv, _ := HermEigen(m)
	return energies, hv
}

// SpinSquared is ⟨S²⟩ = s(s + 1) on any state: Σ_k (Σ_k/2)² = ¾·1.
func SpinSquared() float64 {
	tot := 0.0
	for k := 0; k < 3; k++ {
		s := SpinOp(k)
		for i := 0; i < 4; i++ {
			for j := 0; j < 4; j++ {
				tot += real(s[i][j] * s[j][i])
			}
		}
	}
	return tot / 4 // trace over 4 states, so the per-state value
}

// ---------------------------------------------------------------- Landau levels and g

// LandauLevel is one eigenstate of the Dirac electron in a field B along z
// (p_z = 0), built in the oscillator basis of the cyclotron motion.
type LandauLevel struct {
	E  float64 // energy [m c²]
	N  int     // dominant oscillator level of the upper (large) components
	Up bool    // dominant spin of the upper components
}

// LandauDirac diagonalizes H = [[1, σ·π], [σ·π, −1]] with
// σ·π = √(2b)·[[0, a], [a†, 0]], b = eB ħ/(m²c²) (the field over the
// critical field), truncated at nmax oscillator levels. Basis order:
// upper ↑ (n), upper ↓ (n), lower ↑ (n), lower ↓ (n).
func LandauDirac(b float64, nmax int) []LandauLevel {
	n := nmax
	dim := 4 * n
	h := make([][]float64, dim)
	for i := range h {
		h[i] = make([]float64, dim)
	}
	uU, uD, lU, lD := 0, n, 2*n, 3*n
	s := math.Sqrt(2 * b)
	for k := 0; k < n; k++ {
		h[uU+k][uU+k], h[uD+k][uD+k] = 1, 1
		h[lU+k][lU+k], h[lD+k][lD+k] = -1, -1
		// (σ·π)_{↑k, ↓k+1} = √(2b)·√(k+1) (a lowers), and its transpose
		if k+1 < n {
			v := s * math.Sqrt(float64(k+1))
			h[uU+k][lD+k+1], h[lD+k+1][uU+k] = v, v
			h[lU+k][uD+k+1], h[uD+k+1][lU+k] = v, v
		}
	}
	vals, vecs := jacobiEigen(h)
	var out []LandauLevel
	for c := range vals {
		if vals[c] < 0 {
			continue
		}
		best, bi := 0.0, 0
		for i := 0; i < 2*n; i++ {
			if w := vecs[i][c] * vecs[i][c]; w > best {
				best, bi = w, i
			}
		}
		if bi%n >= n-2 { // states at the basis edge are truncation artifacts
			continue
		}
		out = append(out, LandauLevel{E: vals[c], N: bi % n, Up: bi < n})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].E < out[j].E })
	return out
}

// GFromLandau reads g from the levels: the spin splitting of the lowest
// cyclotron level over the cyclotron spacing gives g/2.
func GFromLandau(lv []LandauLevel) (g, spin, cyc float64) {
	find := func(n int, up bool) float64 {
		for _, l := range lv {
			if l.N == n && l.Up == up {
				return l.E
			}
		}
		return math.NaN()
	}
	e0d, e0u, e1d := find(0, false), find(0, true), find(1, false)
	spin, cyc = e0u-e0d, e1d-e0d
	return 2 * spin / cyc, spin, cyc
}

// ---------------------------------------------------------------- hydrogen

// DiracExact is the exact Coulomb level (infinite nuclear mass) in m c²:
// E = [1 + (Zα/(n − |κ| + √(κ² − Z²α²)))²]^(−½).
func DiracExact(n, kappa int, za float64) float64 {
	k := math.Abs(float64(kappa))
	g := math.Sqrt(k*k - za*za)
	d := float64(n) - k + g
	return 1 / math.Sqrt(1+za*za/(d*d))
}

// radialEnd integrates the radial Dirac equations
//
//	G′ = −(κ/r)G + (E + 1 − V)F,  F′ = (κ/r)F − (E − 1 − V)G,  V = −Zα/r,
//
// outward in x = ln r by RK4 from r0, starting on the regular solution
// G = r^γ, F = r^γ(γ + κ)/(Zα), and returns G at rmax (its sign flips at
// each eigenvalue).
func radialEnd(e float64, kappa int, za, r0, rmax float64, steps int) float64 {
	k := float64(kappa)
	gam := math.Sqrt(k*k - za*za)
	x0, x1 := math.Log(r0), math.Log(rmax)
	h := (x1 - x0) / float64(steps)
	f := func(x, G, F float64) (float64, float64) {
		r := math.Exp(x)
		v := -za / r
		return r * (-(k/r)*G + (e+1-v)*F), r * ((k/r)*F - (e-1-v)*G)
	}
	G := math.Pow(r0, gam)
	F := G * (gam + k) / za
	x := x0
	for i := 0; i < steps; i++ {
		a1, b1 := f(x, G, F)
		a2, b2 := f(x+h/2, G+h/2*a1, F+h/2*b1)
		a3, b3 := f(x+h/2, G+h/2*a2, F+h/2*b2)
		a4, b4 := f(x+h, G+h*a3, F+h*b3)
		G += h / 6 * (a1 + 2*a2 + 2*a3 + a4)
		F += h / 6 * (b1 + 2*b2 + 2*b3 + b4)
		x += h
		if s := math.Abs(G) + math.Abs(F); s > 1e100 { // keep the growing solution finite
			G /= 1e100
			F /= 1e100
		}
	}
	return G
}

// HydrogenLevel finds the level with principal number n and κ by
// bisection on the sign of the outward solution, bracketing the binding
// energy between 0.5 and 1.5 times the Bohr value (Zα)²/(2n²) — a
// bracket that holds no other level of the same κ.
func HydrogenLevel(n, kappa int, za float64) float64 { return hydrogenLevelSteps(n, kappa, za, 60000) }

func hydrogenLevelSteps(n, kappa int, za float64, steps int) float64 {
	bohr := za * za / (2 * float64(n*n))
	lo, hi := 1-1.5*bohr, 1-0.5*bohr
	rmax := 40 * float64(n*n) / za
	end := func(e float64) float64 { return radialEnd(e, kappa, za, 1e-3, rmax, steps) }
	flo := end(lo)
	for it := 0; it < 80; it++ {
		mid := 0.5 * (lo + hi)
		fm := end(mid)
		if math.Signbit(fm) == math.Signbit(flo) {
			lo, flo = mid, fm
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// ---------------------------------------------------------------- small symmetric eigen

// jacobiEigen diagonalizes a real symmetric matrix (cyclic Jacobi):
// eigenvalues (unsorted) and eigenvectors as the columns of vecs.
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
	vals := make([]float64, n)
	for i := range vals {
		vals[i] = a[i][i]
	}
	return vals, v
}
