// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Harmonic and size analysis of Nachamkin's plasmoid, as the paper states
// it (PL-TR-92-3044, 1992).
//
// A plasmoid of radius R must satisfy his eq. 35,
//
//	j1(αR) = 0  and  j1(kR) = 0,
//
// with α the spatial wavenumber of the force-free field and k = ω/c, and
// his eq. 48, k > α (above the critical frequency ω0 = αc, eq. 15). He
// takes αR at the first zero, 4.4934 (the minimum radius, §II), so kR
// must be a higher zero x_n of j1: the allowed frequencies of a plasmoid
// of radius R are a comb, f_n = x_n·c/(2πR), n = 2, 3, … ("harmonics").
// The free-electron density follows from his eq. 49 (ω² = ωp² + α²c²).
//
// Appendix B picks the harmonic from the stored energy and his
// small-oscillation factor C (eqs. 50, 58): the field amplitude E0 holds
// half the energy, and ω² = |e/m|·E0/(C·R). The electrons then oscillate
// with amplitude ξ = C·R at speed ω·ξ, so
//
//	v/c = C·kR = C·x_n,
//
// and his non-relativistic analysis (§V–VII) needs C·x_n ≪ 1.

import "math"

// Loschmidt is the number density of an ideal gas at 0 °C, 1 atm [m⁻³]
// (Nachamkin's 1.126×10²³ molecules in a 10 cm sphere).
const Loschmidt = 2.686780111e25

// AirBreakdown is the DC breakdown field of air at STP, order of
// magnitude [V/m].
const AirBreakdown = 3e6

// J1Zeros returns the zeros of j1 in (0, xmax], in order (the n-th entry
// is x_n; x_1 = 4.4934).
func J1Zeros(xmax float64) []float64 {
	var z []float64
	// zeros of j1 lie in ((n+½)π − π/2, (n+½)π + π/2); tan x = x
	for n := 1; ; n++ {
		lo := float64(n)*math.Pi + 1e-9
		hi := (float64(n) + 0.5) * math.Pi
		if lo > xmax {
			break
		}
		for i := 0; i < 200; i++ {
			m := 0.5 * (lo + hi)
			if j1(lo)*j1(m) <= 0 {
				hi = m
			} else {
				lo = m
			}
		}
		x := 0.5 * (lo + hi)
		if x <= xmax {
			z = append(z, x)
		}
	}
	return z
}

// jl is the spherical Bessel function j_l(x) by upward recurrence (x > 0).
func jl(l int, x float64) float64 {
	j0 := math.Sin(x) / x
	if l == 0 {
		return j0
	}
	j1v := math.Sin(x)/(x*x) - math.Cos(x)/x
	for k := 1; k < l; k++ {
		j0, j1v = j1v, float64(2*k+1)/x*j1v-j0
	}
	return j1v
}

// FirstZero is the first zero of j_l: αR for an order-l plasmoid of
// minimum radius (Nachamkin eq. 8 with general l).
func FirstZero(l int) float64 {
	step := 0.01
	for x := float64(l) + 0.5; x < 50; x += step {
		if jl(l, x)*jl(l, x+step) <= 0 {
			lo, hi := x, x+step
			for i := 0; i < 100; i++ {
				m := 0.5 * (lo + hi)
				if jl(l, lo)*jl(l, m) <= 0 {
					hi = m
				} else {
					lo = m
				}
			}
			return 0.5 * (lo + hi)
		}
	}
	return math.NaN()
}

// Design is one plasmoid: radius R, harmonic x_n = kR, and what follows.
type Design struct {
	R, W, C           float64 // radius [m], total energy [J] (field + fluid), small-oscillation factor
	Alpha, X          float64 // α [1/m], kR (a zero of j1)
	N                 int     // index of the zero (x_N)
	Omega, F, F0, Lam float64 // ω [rad/s], f [Hz], critical frequency αc/2π [Hz], free-space wavelength [m]
	E0                float64 // peak field holding half of W [V/m]
	Ne, Electrons     float64 // free-electron density [m⁻³] and number
	IonFrac           float64 // Ne ÷ Loschmidt
	Beta              float64 // electron speed ÷ c = C·X
	EOverBreakdown    float64 // peak field ÷ air breakdown
	ModeTE, ModeTM    Mode    // leaky l = 1 resonances (root finder, kR ≤ 1000)
	TauTE, TauTM      float64 // energy lifetimes [s]
}

// IM is ∫|M|² dV over the sphere for the l = 1 field of his eq. 8 (the
// normalization E0 multiplies): his eq. B3 gives 35.87/α³.
func IM(alpha float64) float64 {
	R := J1Zero / alpha
	vol := 4.0 / 3 * math.Pi * R * R * R
	return meanSquare() * MMax2 * vol
}

// MMax2 is max|M|² of eq. 8 (l = 1), reached at the centre: |M(0)| = 2/3.
const MMax2 = 4.0 / 9

// fieldE0 is the amplitude E0 for a time-averaged field energy wf, his
// eq. B2: (ε0·E0²/4)·∫|M|² dV = wf (the magnetic energy is smaller by
// (α/k)², which he neglects).
func fieldE0(R, wf float64) float64 {
	return math.Sqrt(4 * wf / (Eps0 * IM(J1Zero/R)))
}

// finish fills the quantities that follow from R, X and C.
func (d *Design) finish() {
	d.Alpha = J1Zero / d.R
	d.Omega = d.X * C0 / d.R
	d.F = d.Omega / (2 * math.Pi)
	d.F0 = d.Alpha * C0 / (2 * math.Pi)
	d.Lam = 2 * math.Pi * d.R / d.X
	d.Ne = Eps0 * Me * (d.Omega*d.Omega - d.Alpha*d.Alpha*C0*C0) / (Qe * Qe) // his eq. 49
	d.Electrons = d.Ne * 4.0 / 3 * math.Pi * d.R * d.R * d.R
	d.IonFrac = d.Ne / Loschmidt
	d.Beta = d.C * d.X
	d.EOverBreakdown = d.E0 * math.Sqrt(MMax2) / AirBreakdown // peak field |E0·M| at the centre
	// energy lifetimes: the complex root finder up to kR = 1000, the
	// large-kR closed forms beyond (radiation.go; they agree at 1000)
	d.TauTE, d.TauTM = TauTE(d.R, d.X), TauTM(d.R, d.X)
	if d.X <= 1000 {
		s := Sphere{R: d.R, Wp: math.Sqrt(d.Ne * Qe * Qe / (Me * Eps0))}
		d.ModeTE, d.ModeTM = s.FindMode(true), s.FindMode(false)
		if d.ModeTE.Found {
			d.TauTE = d.ModeTE.Lifetime
		}
		if d.ModeTM.Found {
			d.TauTM = d.ModeTM.Lifetime
		}
	}
}

// DesignFromC follows Appendix B: energy W split equally between field
// and fluid, field amplitude from W, ω from C (eq. 58), then kR moved to
// the nearest zero of j1 above αR (eq. 35); C is recomputed for that k at
// the same field.
func DesignFromC(R, W, C float64) Design {
	d := Design{R: R, W: W}
	d.E0 = fieldE0(R, W/2)
	w := math.Sqrt(EoM * d.E0 / (C * R))
	d.X, d.N = nearestZero(w * R / C0)
	om := d.X * C0 / R
	d.C = EoM * d.E0 / (R * om * om)
	d.finish()
	return d
}

// DesignFromBeta picks the lowest harmonic at which the electrons move at
// no more than beta·c, for energy W: with E0 fixed by W,
// v/c = |e/m|·E0·R/(X·c²), so X ≥ |e/m|·E0·R/(beta·c²).
func DesignFromBeta(R, W, beta float64) Design {
	d := Design{R: R, W: W}
	d.E0 = fieldE0(R, W/2)
	xmin := math.Max(EoM*d.E0*R/(beta*C0*C0), J1Zero*1.0001)
	d.X, d.N = zeroAbove(xmin)
	om := d.X * C0 / R
	d.C = EoM * d.E0 / (R * om * om)
	d.finish()
	return d
}

// DesignAt is a plasmoid of radius R at the harmonic X = kR, with energy W.
func DesignAt(R, W, X float64) Design {
	d := Design{R: R, W: W, X: X}
	d.E0 = fieldE0(R, W/2)
	om := X * C0 / R
	d.C = EoM * d.E0 / (R * om * om)
	d.finish()
	return d
}

// nearestZero is the zero of j1 nearest to x, above the first.
func nearestZero(x float64) (float64, int) {
	z := J1Zeros(math.Max(x, J1Zero) + 2*math.Pi)
	best, n := z[1], 2
	for i := 1; i < len(z); i++ {
		if math.Abs(z[i]-x) < math.Abs(best-x) {
			best, n = z[i], i+1
		}
	}
	return best, n
}

// zeroAbove is the first zero of j1 at or above x (and above x_1).
func zeroAbove(x float64) (float64, int) {
	z := J1Zeros(x + 2*math.Pi)
	for i := 1; i < len(z); i++ {
		if z[i] >= x {
			return z[i], i + 1
		}
	}
	return z[len(z)-1], len(z)
}

// AirElectronsPerMolecule is the number of electrons in an average air
// molecule (78 % N₂ with 14, 21 % O₂ with 16, 1 % Ar with 18): the
// free-electron density of fully stripped sea-level air is this many
// times Loschmidt's number.
const AirElectronsPerMolecule = 0.78*14 + 0.21*16 + 0.01*18

// ForLifetime is the plasmoid of radius R whose shorter (TE) radiative
// lifetime is tau: from τ_TE = R·y²/(2c·x1²), y = kR = √(2c·x1²·τ/R), and
// with eq. 49 the free-electron density. Neither depends on the energy.
func ForLifetime(R, tau float64) (y, f, ne float64) {
	y = math.Sqrt(2 * C0 * J1Zero * J1Zero * tau / R)
	om := y * C0 / R
	f = om / (2 * math.Pi)
	ne = Eps0 * Me * (om*om - J1Zero*J1Zero*C0*C0/(R*R)) / (Qe * Qe)
	return
}
