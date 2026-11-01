// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/hyper/ga"
)

// Oscillating point dipole: the coupling between a lumped (much smaller
// than a wavelength) charge distribution and light. It is how a vibrating
// molecule absorbs and releases infrared light (examples/dna_vibration):
// the molecule is ~1 nm, the light 3–30 µm, so it sits deep in the lumped
// limit of §3c and couples only through its dipole moment.
//
// Fields of p(t) = Re(p·e^{−iωt}) at r, k = ω/c (J. D. Jackson, Classical
// Electrodynamics, 3rd ed., eqs. 9.18, SI):
//
//	H = (c k²/4π) (n×p) (e^{ikr}/r) (1 − 1/(ikr))
//	E = (1/4πε₀) { k² (n×p)×n e^{ikr}/r + [3n(n·p) − p] (1/r³ − ik/r²) e^{ikr} }
//
// and the time-averaged radiated power, Larmor's formula for a harmonic
// dipole (Jackson eq. 9.24), P = c² Z₀ k⁴ |p|²/(12π) = ω⁴|p|²/(12π ε₀ c³).

// DipoleFields returns the phasors E [V/m] and B [T] of a dipole p [C·m]
// at the origin, observed at r [m], angular frequency omega [rad/s].
func DipoleFields(p [3]complex128, r [3]float64, omega float64) (E, B [3]complex128) {
	c := SpeedOfLight
	k := omega / c
	rr := math.Sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
	n := [3]float64{r[0] / rr, r[1] / rr, r[2] / rr}
	ph := cmplx.Exp(complex(0, k*rr))
	cr := func(a [3]float64, b [3]complex128) [3]complex128 {
		ac := [3]complex128{complex(a[0], 0), complex(a[1], 0), complex(a[2], 0)}
		return [3]complex128{ac[1]*b[2] - ac[2]*b[1], ac[2]*b[0] - ac[0]*b[2], ac[0]*b[1] - ac[1]*b[0]}
	}
	np := cr(n, p)                                                        // n×p
	npn := [3]complex128{np[1]*complex(n[2], 0) - np[2]*complex(n[1], 0), // (n×p)×n
		np[2]*complex(n[0], 0) - np[0]*complex(n[2], 0),
		np[0]*complex(n[1], 0) - np[1]*complex(n[0], 0)}
	ndp := complex(n[0], 0)*p[0] + complex(n[1], 0)*p[1] + complex(n[2], 0)*p[2]
	near := complex(1/(rr*rr*rr), -k/(rr*rr))
	hf := complex(c*k*k/(4*math.Pi), 0) * ph / complex(rr, 0) * (1 - 1/complex(0, k*rr))
	ef := complex(1/(4*math.Pi*Eps0), 0)
	for d := 0; d < 3; d++ {
		H := hf * np[d]
		B[d] = complex(Mu0, 0) * H
		E[d] = ef * (complex(k*k, 0)*npn[d]*ph/complex(rr, 0) + (3*complex(n[d], 0)*ndp-p[d])*near*ph)
	}
	return E, B
}

// LarmorPower is the time-averaged power [W] radiated by a dipole of
// amplitude |p| [C·m] oscillating at omega [rad/s].
func LarmorPower(p [3]complex128, omega float64) float64 {
	pp := 0.0
	for _, v := range p {
		pp += real(v * cmplx.Conj(v))
	}
	c := SpeedOfLight
	return math.Pow(omega, 4) * pp / (12 * math.Pi * Eps0 * c * c * c)
}

// PoyntingGA is the time-averaged Poynting vector ½·Re(E×H*) [W/m²] from
// field phasors, computed as the grade-1 part of ½·F·F† for the real and
// imaginary parts of the field multivector F = E + I·c·B (§3f; with the
// reverse F† = E − I·c·B, ½FF† = ½(E² + c²B²) + c·E×B).
func PoyntingGA(E, B [3]complex128) [3]float64 {
	var s [3]float64
	for _, part := range []func(complex128) float64{func(z complex128) float64 { return real(z) }, func(z complex128) float64 { return imag(z) }} {
		f := FieldMV([3]float64{part(E[0]), part(E[1]), part(E[2])}, [3]float64{part(B[0]), part(B[1]), part(B[2])})
		m := f.Mul(f.Reverse()).Scale(0.5)
		g := m.Grade(1)
		s[0] += g[ga.G3E1]
		s[1] += g[ga.G3E2]
		s[2] += g[ga.G3E3]
	}
	// grade 1 of ½FF† is c·E×B; ⟨S⟩ = ½·Re(E×B*)/μ₀.
	k := 0.5 / (SpeedOfLight * Mu0)
	return [3]float64{k * s[0], k * s[1], k * s[2]}
}

// PoyntingVector is the same time average from the vector cross product,
// ½·Re(E×B*)/μ₀ (the plain-vector twin).
func PoyntingVector(E, B [3]complex128) [3]float64 {
	bc := [3]complex128{cmplx.Conj(B[0]), cmplx.Conj(B[1]), cmplx.Conj(B[2])}
	x := [3]complex128{E[1]*bc[2] - E[2]*bc[1], E[2]*bc[0] - E[0]*bc[2], E[0]*bc[1] - E[1]*bc[0]}
	return [3]float64{0.5 * real(x[0]) / Mu0, 0.5 * real(x[1]) / Mu0, 0.5 * real(x[2]) / Mu0}
}

// RadiatedPower integrates the time-averaged Poynting flux of a dipole
// through a sphere of radius R [m] (Gauss–Legendre in cos θ with nTheta
// points, nPhi equally spaced azimuths), using the GA form when useGA is
// true and the vector form otherwise.
func RadiatedPower(p [3]complex128, omega, R float64, nTheta, nPhi int, useGA bool) float64 {
	x, w := gaussLegendre(nTheta)
	total := 0.0
	for i := range x {
		st := math.Sqrt(1 - x[i]*x[i])
		for j := 0; j < nPhi; j++ {
			phi := 2 * math.Pi * (float64(j) + 0.5) / float64(nPhi)
			n := [3]float64{st * math.Cos(phi), st * math.Sin(phi), x[i]}
			r := [3]float64{R * n[0], R * n[1], R * n[2]}
			E, B := DipoleFields(p, r, omega)
			var S [3]float64
			if useGA {
				S = PoyntingGA(E, B)
			} else {
				S = PoyntingVector(E, B)
			}
			total += (S[0]*n[0] + S[1]*n[1] + S[2]*n[2]) * w[i] * (2 * math.Pi / float64(nPhi)) * R * R
		}
	}
	return total
}

// gaussLegendre returns n nodes and weights on [−1, 1] (Newton on the
// Legendre recurrence; Golub & Welsch's values for small n).
func gaussLegendre(n int) ([]float64, []float64) {
	x := make([]float64, n)
	w := make([]float64, n)
	for i := 0; i < n; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		for it := 0; it < 100; it++ {
			p0, p1 := 1.0, z
			for k := 2; k <= n; k++ {
				p0, p1 = p1, ((2*float64(k)-1)*z*p1-(float64(k)-1)*p0)/float64(k)
			}
			dp := float64(n) * (z*p1 - p0) / (z*z - 1)
			dz := p1 / dp
			z -= dz
			if math.Abs(dz) < 1e-16 {
				p0, p1 = 1.0, z
				for k := 2; k <= n; k++ {
					p0, p1 = p1, ((2*float64(k)-1)*z*p1-(float64(k)-1)*p0)/float64(k)
				}
				dp = float64(n) * (z*p1 - p0) / (z*z - 1)
				x[i] = z
				w[i] = 2 / ((1 - z*z) * dp * dp)
				break
			}
		}
	}
	return x, w
}
