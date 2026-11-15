// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

// The spin connection of ECE theory (M. W. Evans; H. Eckardt), hypothetical
// mode (H48). Primary source: H. Eckardt, *ECE Unified Field Theory* Vol. 1
// (2022), Ch. 4–6 (equation numbers below are the book's); the slide decks
// "The ECE Engineering Model" and its linearized companion give the same
// equations without the factor c on ω₀. Copies in examples/AIAS/papers/other/.
//
// ECE makes the potential proportional to Cartan's tetrad, A^a = A⁽⁰⁾q^a,
// and the field to its torsion, so the field gets the spin connection
// through the first structure equation (4.172). With one polarization
// (4.207–4.212), ω₀ the scalar and ω the vector spin connection, both [1/m]:
//
//	E = −∇φ − ∂A/∂t − cω₀A + ωφ          B = ∇×A − ω×A          (ECE)
//
// Summing the polarization index instead of dropping it gives a factor −2
// (6.189–6.190):
//
//	E = −∇φ − ∂A/∂t + 2(cω₀A − ωφ)       B = ∇×A + 2ω×A         (ECE2)
//
// The field equations for E and B are Maxwell's in form (4.72–4.75), with the
// homogeneous ("magnetic") currents set to zero in practice. Written in the
// potentials they become constraints on the spin connection (4.213–4.214):
//
//	Gauss    ∇·(ω×A) = 0
//	Faraday  −c∇×(ω₀A) + ∇×(ωφ) − ∂(ω×A)/∂t = 0
//
// and ECE replaces the gauge by "antisymmetry" laws (5.9, 5.14, 5.18–5.20):
//
//	electric  ∇φ − ∂A/∂t − cω₀A − ωφ = 0
//	magnetic  ∂_jA_i + ∂_iA_j + ω_jA_i + ω_iA_j = 0   (i ≠ j)
//
// With the electric law E = −2(∂A/∂t + cω₀A) (5.24); with the further split
// cω₀A = ωφ = ½(∇φ − ∂A/∂t) of the slides (slide 17) the two spin-connection
// terms cancel and E is Maxwell's −∇φ − ∂A/∂t. Where the Gauss and Faraday
// constraints fail, ∇·B ≠ 0 and ∇×E + ∂B/∂t ≠ 0: magnetic charge and
// current, the Magnetic modes of magnetic.go (H2).
//
// ECE2 writes the field equations with "wave numbers" (6.164–6.169)
//
//	κ₀ = φ/(cW⁽⁰⁾) − ω₀        κ = A/W⁽⁰⁾ − ω        [1/m]
//	∇·E = −2κ·E = ρ/ε₀        ∇×B − (1/c²)∂E/∂t = 2(κ₀E/c + κ×B) = μ₀J
//
// so that charge and current are geometry: a conductivity σ = 2κ₀/(μ₀c) and a
// Hall-like current (2/μ₀)κ×B. (The slides, 34, print ∇·E = +κ·E with κ
// carrying the 2.) In free space the book sets κ₀ = κ = 0 (6.178–6.179), which
// is Maxwell with W = W⁽⁰⁾ω as the potential. Region evolves the vacuum with a
// declared κ₀, κ (spin.go), the current's work in its own column, SpinWork;
// the book gives no energy theorem with the spin connection (§5.3.1 states
// Poynting's per polarization index).
//
// This file evaluates the relations at a point from declared potentials and
// spin connections, by central differences.

import (
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// Potentials are declared fields: φ [V], A [V·s/m], ω₀ [1/m], ω [1/m],
// each a function of position x [m] and time t [s]. Nil spin connections
// are zero.
type Potentials struct {
	Phi    func(x [3]float64, t float64) float64
	A      func(x [3]float64, t float64) [3]float64
	Omega0 func(x [3]float64, t float64) float64
	Omega  func(x [3]float64, t float64) [3]float64
}

// ECEPoint holds the fields and constraint residuals at one point.
type ECEPoint struct {
	EStd, BStd   [3]float64 // Maxwell: −∇Φ − ∂A/∂t, ∇×A
	EECE, BECE   [3]float64 // 4.211–4.212: − cω₀A + ωφ, − ω×A
	EECE2, BECE2 [3]float64 // 6.189–6.190: + 2(cω₀A − ωφ), + 2ω×A

	// Magnetic charge and current that the ECE fields carry: ρ_m = ∇·B/μ₀,
	// μ₀J_m = −(∇×E + ∂B/∂t) [the units of em/duality.go]. Zero is ECE's
	// Gauss and Faraday constraint (4.213–4.214).
	DivBECE     float64    // ∇·B_ECE = −∇·(ω×A)
	FaradayECE  [3]float64 // ∇×E_ECE + ∂B_ECE/∂t
	ElecAnti    [3]float64 // ∇φ − ∂A/∂t − cω₀A − ωφ (5.14)
	MagAnti     [3]float64 // pairs (2,3), (1,3), (1,2) of 5.18–5.20
	SpinFromPot [3]float64 // ω = (∇φ − ∂A/∂t)/(2φ) (slide 17)
	Omega0FromP float64    // ω₀ = (∇φ − ∂A/∂t)·A/(2cA²) (slide 17, with the book's c)
}

// Kappa is ECE2's (κ₀, κ) [1/m] (6.164–6.165), for a declared W⁽⁰⁾ [Wb].
func (p Potentials) Kappa(x [3]float64, t, w0 float64) (k0 float64, k [3]float64) {
	o0, o := p.spin(x, t)
	a := p.A(x, t)
	k0 = p.Phi(x, t)/(em.SpeedOfLight*w0) - o0
	for c := 0; c < 3; c++ {
		k[c] = a[c]/w0 - o[c]
	}
	return
}

func (p Potentials) spin(x [3]float64, t float64) (o0 float64, o [3]float64) {
	if p.Omega0 != nil {
		o0 = p.Omega0(x, t)
	}
	if p.Omega != nil {
		o = p.Omega(x, t)
	}
	return
}

func cross(a, b [3]float64) [3]float64 {
	return [3]float64{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

// fields is (E, B) in one convention: s = 0 Maxwell, 1 ECE, 2 ECE2.
func (p Potentials) fields(x [3]float64, t, h, ht float64, s int) (e, b [3]float64) {
	d := func(f func(y [3]float64) [3]float64, c int) [3]float64 {
		xp, xm := x, x
		xp[c] += h
		xm[c] -= h
		fp, fm := f(xp), f(xm)
		return [3]float64{(fp[0] - fm[0]) / (2 * h), (fp[1] - fm[1]) / (2 * h), (fp[2] - fm[2]) / (2 * h)}
	}
	av := func(y [3]float64) [3]float64 { return p.A(y, t) }
	var dA [3][3]float64 // dA[c][i] = ∂A_i/∂x_c
	for c := 0; c < 3; c++ {
		dA[c] = d(av, c)
	}
	ap, am := p.A(x, t+ht), p.A(x, t-ht)
	a, phi := p.A(x, t), p.Phi(x, t)
	o0, o := p.spin(x, t)
	oa := cross(o, a)
	cl := em.SpeedOfLight
	for c := 0; c < 3; c++ {
		xp, xm := x, x
		xp[c] += h
		xm[c] -= h
		gphi := (p.Phi(xp, t) - p.Phi(xm, t)) / (2 * h)
		e[c] = -gphi - (ap[c]-am[c])/(2*ht)
		switch s {
		case 1:
			e[c] += -cl*o0*a[c] + o[c]*phi
		case 2:
			e[c] += 2 * (cl*o0*a[c] - phi*o[c])
		}
	}
	b = [3]float64{dA[1][2] - dA[2][1], dA[2][0] - dA[0][2], dA[0][1] - dA[1][0]}
	for c := 0; c < 3; c++ {
		switch s {
		case 1:
			b[c] -= oa[c]
		case 2:
			b[c] += 2 * oa[c]
		}
	}
	return
}

// At evaluates the three field conventions and ECE's constraints at x, t,
// with spatial step h [m] and time step ht [s] for the central differences
// (second-order accurate; exact for quadratic fields).
func (p Potentials) At(x [3]float64, t, h, ht float64) ECEPoint {
	var r ECEPoint
	r.EStd, r.BStd = p.fields(x, t, h, ht, 0)
	r.EECE, r.BECE = p.fields(x, t, h, ht, 1)
	r.EECE2, r.BECE2 = p.fields(x, t, h, ht, 2)

	// ∇·B and ∇×E + ∂B/∂t of the ECE fields, by a second, outer difference.
	H := 8 * h
	var curlE [3]float64
	for c := 0; c < 3; c++ {
		xp, xm := x, x
		xp[c] += H
		xm[c] -= H
		ep, bp := p.fields(xp, t, h, ht, 1)
		em_, bm := p.fields(xm, t, h, ht, 1)
		r.DivBECE += (bp[c] - bm[c]) / (2 * H)
		// contributions of ∂/∂x_c to the curl
		c1, c2 := (c+1)%3, (c+2)%3
		curlE[c2] += (ep[c1] - em_[c1]) / (2 * H)
		curlE[c1] -= (ep[c2] - em_[c2]) / (2 * H)
	}
	HT := 8 * ht
	_, bt1 := p.fields(x, t+HT, h, ht, 1)
	_, bt0 := p.fields(x, t-HT, h, ht, 1)
	for c := 0; c < 3; c++ {
		r.FaradayECE[c] = curlE[c] + (bt1[c]-bt0[c])/(2*HT)
	}

	// Antisymmetry (slide 16) and the spin connection it implies (slide 17).
	a, phi := p.A(x, t), p.Phi(x, t)
	o0, o := p.spin(x, t)
	ap, am := p.A(x, t+ht), p.A(x, t-ht)
	var g [3]float64 // ∇Φ − ∂A/∂t
	var dA [3][3]float64
	for c := 0; c < 3; c++ {
		xp, xm := x, x
		xp[c] += h
		xm[c] -= h
		g[c] = (p.Phi(xp, t)-p.Phi(xm, t))/(2*h) - (ap[c]-am[c])/(2*ht)
		fp, fm := p.A(xp, t), p.A(xm, t)
		for i := 0; i < 3; i++ {
			dA[c][i] = (fp[i] - fm[i]) / (2 * h)
		}
	}
	a2 := a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
	var ga float64
	for c := 0; c < 3; c++ {
		r.ElecAnti[c] = g[c] - em.SpeedOfLight*o0*a[c] - o[c]*phi
		ga += g[c] * a[c]
		if phi != 0 {
			r.SpinFromPot[c] = g[c] / (2 * phi)
		}
	}
	if a2 != 0 {
		r.Omega0FromP = ga / (2 * em.SpeedOfLight * a2)
	}
	pairs := [3][2]int{{1, 2}, {0, 2}, {0, 1}}
	for m, ij := range pairs {
		i, j := ij[0], ij[1]
		r.MagAnti[m] = dA[j][i] + dA[i][j] + o[j]*a[i] + o[i]*a[j]
	}
	return r
}

// UFT 456 (Eckardt & Goel, "The homogeneous and inhomogeneous ECE current,
// Part III", 9 Mar 2024, rev. 3 Jun 2025) writes both ECE2 current pairs as
// fields (eqs. 28–31; UFT 454 eqs. 24–27):
//
//	ρ_p = −2κ_Λ·B                 [V·s/m³]   ∇·B = ρ_p            (20)
//	V   = 2(cκ_Λ0 B − κ_Λ×E)      [V/m²]     ∂B/∂t + ∇×E = V      (21)
//	ρ_e = −2ε₀κ·E                 [A·s/m³]   ∇·D = ρ_e            (22)
//	J   = (2/μ₀)(κ₀E/c + κ×B)     [A/m²]     −∂D/∂t + ∇×H = J     (23)
//
// (20)–(23) are Maxwell's equations with magnetic sources in the weber
// convention, V = −J_m (Jackson 3rd ed. §6.11). Checked
// (examples/AIAS/uft456/check.py, sympy):
//   - Gauss's law (20)+(28) is carried forward by Faraday's law (21)+(29):
//     ∂ρ_p/∂t − ∇·V = −2cκ_Λ0(∇·B − ρ_p).
//   - Coulomb's law (22)+(30) is not carried forward by Ampère's law
//     (23)+(31): ∂ρ_e/∂t + ∇·J = −4κ·(∇×H) + 2cκ₀(∇·D − ρ_e), so it breaks at
//     the rate −4κ·∇×H (CoulombDefectRate; Region.CoulombDefect measures it
//     on the grid). With the sign of κ×B in (31) reversed it would hold: the
//     two Lorentz-type terms are not dual images.
//   - The Poynting theorem is ∂u/∂t + ∇·(E×H) = −J·E + V·H (eq. 57 with the
//     conductivity parts): V·H is a source density on the right, not a term
//     of the Poynting vector.

// ECE2Sources are the field-built charges and currents of UFT 456 eqs. 28–31.
type ECE2Sources struct {
	RhoP float64    // ρ_p = −2κ_Λ·B [V·s/m³] (28)
	V    [3]float64 // 2(cκ_Λ0 B − κ_Λ×E) [V/m²] (29)
	RhoE float64    // −2ε₀κ·E [A·s/m³] (30)
	J    [3]float64 // (2/μ₀)(κ₀E/c + κ×B) [A/m²] (31)
}

// ECE2SourcesAt evaluates eqs. 28–31 of UFT 456 for fields E [V/m], B [T]
// and wave numbers κ₀, κ (inhomogeneous) and κ_Λ0, κ_Λ (homogeneous) [1/m].
func ECE2SourcesAt(k0 float64, k [3]float64, kL0 float64, kL [3]float64, e, b [3]float64) ECE2Sources {
	cl := em.SpeedOfLight
	kxB, kLxE := cross(k, b), cross(kL, e)
	var s ECE2Sources
	for i := 0; i < 3; i++ {
		s.RhoP -= 2 * kL[i] * b[i]
		s.RhoE -= 2 * em.Eps0 * k[i] * e[i]
		s.V[i] = 2 * (cl*kL0*b[i] - kLxE[i])
		s.J[i] = 2 / em.Mu0 * (k0*e[i]/cl + kxB[i])
	}
	return s
}

// ECE2Conductivities returns UFT 456's conductivities and rates:
// σ₀ = 2κ₀/(μ₀c) [S/m] (33), η₀ = 2cκ_Λ0μ₀ [Ω/m] (39, a magnetic
// conductivity, not the impedance of free space), ω_e = σ₀/ε₀ = 2cκ₀ (36),
// ω_p = η₀/μ₀ = 2cκ_Λ0 (42) [1/s]. In the usual lossy-medium form
// ∇×E = −μ₀∂H/∂t − σ*H the paper's V = η₀H is σ* = −η₀: η₀ > 0 is gain.
func ECE2Conductivities(k0, kL0 float64) (sigma0, eta0, omegaE, omegaP float64) {
	cl := em.SpeedOfLight
	sigma0 = 2 * k0 / (em.Mu0 * cl)
	eta0 = 2 * cl * kL0 * em.Mu0
	return sigma0, eta0, sigma0 / em.Eps0, eta0 / em.Mu0
}

// ECE2Dispersion returns the two frequencies ω [1/s] of a plane wave
// exp(i(kx − ωt)) in UFT 456's medium, eq. 73:
//
//	ω² + iω(ω_e − ω_p) + ω_eω_p = c²k²
//
// Im ω > 0 grows. The principal part of eq. 73 is (1/c²)∂²/∂t² − ∇² for
// every ω_e, ω_p, so the signal (front) speed is c: the rates change
// amplitudes and phase velocities, never the front (Sommerfeld, Brillouin
// 1914). The "equilibrium" ω_e = ω_p = ω₀ (eq. 109) gives ω² = c²k² − ω₀²,
// growing for k < ω₀/c; ω_p = −ω_e is Heaviside's distortionless medium,
// ω = ±ck − iω_e.
func ECE2Dispersion(k, omegaE, omegaP float64) (w1, w2 complex128) {
	cl := em.SpeedOfLight
	b := complex(0, omegaE-omegaP)
	cc := complex(omegaE*omegaP-cl*cl*k*k, 0)
	d := cmplx.Sqrt(b*b - 4*cc)
	return (-b + d) / 2, (-b - d) / 2
}

// CoulombDefectRate is ∂ρ_e/∂t + ∇·J on fields that satisfy Coulomb's law
// (22)+(30) at the instant: −4κ·(∇×H) [A/m³], with H = B/μ₀. Zero only where
// κ ⊥ ∇×H.
func CoulombDefectRate(k, curlH [3]float64) float64 {
	return -4 * (k[0]*curlH[0] + k[1]*curlH[1] + k[2]*curlH[2])
}
