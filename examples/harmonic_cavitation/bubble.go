// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Stage 5: one cavity's collapse, and many collapses gathering on the axis.
//
// The cavity is taken as a sphere of the same vapour volume, collapsing
// from rest under the driving pressure Δp = p∞ − p_v, with a little
// non-condensable gas (adiabatic, γ = 1.4), surface tension and viscosity:
// the Rayleigh–Plesset equation (Brennen 1995, eq. 2.12),
//   ρ(R R̈ + 3/2 Ṙ²) = p_g(R) + p_v − p∞ − 2S/R − 4μṘ/R.
// Its ledger: the work p∞·V and the gas, surface and kinetic terms,
//   Δ[2πρR³Ṙ² + 4πSR² + (p∞ − p_v)·4πR³/3] − Δ(gas energy) = −viscous heat,
// is checked along the run. With no gas, S = μ = 0 the collapse time is
// Rayleigh's 0.915·R0·√(ρ/Δp) (tested). The pressure radiated to distance d
// is p = ρ(R²R̈ + 2RṘ²)/d (incompressible near field, Brennen eq. 2.37).

// Collapse is the result of one run.
type Collapse struct {
	R0, Rmin, T       float64 // start radius, minimum radius, collapse time
	PGasMax           float64 // peak gas pressure [Pa]
	EPot              float64 // (p∞ − p_v)·4πR0³/3 [J]
	KMax              float64 // peak kinetic energy of the liquid [J]
	PRadPeakAt1m      float64 // peak radiated pressure × distance [Pa·m]
	LedgerRelResidual float64
}

// RunCollapse integrates from rest at R0 to the first minimum.
func RunCollapse(r0, dp, gasFrac, surface, mu float64) Collapse {
	return RunCollapseTrace(r0, dp, gasFrac, surface, mu, nil)
}

// RunCollapseTrace is RunCollapse, calling rec(t, R, Ṙ, ρ(R²R̈ + 2RṘ²)) at
// every step when rec is not nil.
func RunCollapseTrace(r0, dp, gasFrac, surface, mu float64, rec func(t, r, v, prad float64)) Collapse {
	const gamma = 1.4
	pg0 := gasFrac * dp // gas partial pressure at R0 (declared fraction)
	pInf := dp + pVap
	c := Collapse{R0: r0, EPot: dp * 4 * math.Pi * r0 * r0 * r0 / 3}
	acc := func(r, v float64) float64 {
		pg := pg0 * math.Pow(r0/r, 3*gamma)
		return ((pg+pVap-pInf-2*surface/r-4*mu*v/r)/rho - 1.5*v*v) / r
	}
	energy := func(r, v float64) float64 {
		gas := pg0 * 4 * math.Pi * r0 * r0 * r0 / 3 / (gamma - 1) * math.Pow(r0/r, 3*(gamma-1))
		return 2*math.Pi*rho*r*r*r*v*v + 4*math.Pi*surface*r*r + dp*4*math.Pi*r*r*r/3 + gas
	}
	r, v, t := r0, 0.0, 0.0
	e0 := energy(r, v)
	heat := 0.0
	tr := 0.915 * r0 * math.Sqrt(rho/dp)
	h := tr * 1e-4
	for k := 0; k < 5_000_000; k++ {
		// RK4 with step scaled to the local time scale
		hs := math.Min(h, 0.01*r/math.Max(math.Abs(v), 1e-9))
		k1r, k1v := v, acc(r, v)
		k2r, k2v := v+0.5*hs*k1v, acc(r+0.5*hs*k1r, v+0.5*hs*k1v)
		k3r, k3v := v+0.5*hs*k2v, acc(r+0.5*hs*k2r, v+0.5*hs*k2v)
		k4r, k4v := v+hs*k3v, acc(r+hs*k3r, v+hs*k3v)
		nr := r + hs/6*(k1r+2*k2r+2*k3r+k4r)
		nv := v + hs/6*(k1v+2*k2v+2*k3v+k4v)
		heat += hs * 16 * math.Pi * mu * r * 0.5 * (v*v + nv*nv) // viscous: 4μṘ/R · Ṙ · 4πR²
		if nv >= 0 && v < 0 {                                    // minimum reached
			c.Rmin, c.T = r, t
			break
		}
		r, v, t = nr, nv, t+hs
		if r <= 1e-4*r0 || math.IsNaN(r) { // empty cavity: stop at 10⁻⁴ R0
			c.Rmin, c.T = r, t
			break
		}
		kin := 2 * math.Pi * rho * r * r * r * v * v
		c.KMax = math.Max(c.KMax, kin)
		pg := pg0 * math.Pow(r0/r, 3*gamma)
		c.PGasMax = math.Max(c.PGasMax, pg)
		rdd := acc(r, v)
		c.PRadPeakAt1m = math.Max(c.PRadPeakAt1m, rho*(r*r*rdd+2*r*v*v))
		if rec != nil {
			rec(t, r, v, rho*(r*r*rdd+2*r*v*v))
		}
	}
	c.LedgerRelResidual = math.Abs(energy(r, v)+heat-e0) / e0
	return c
}

// Focus gives the peak pressure on the axis when k collapses on a ring of
// radius d happen together: equidistant, their pulses arrive together and
// add, k·p1(d) (linear superposition, an upper bound). The energy at the
// focus is at most the k collapse energies, all paid by the shaft.
func Focus(c Collapse, k int, d float64) (peak float64) {
	return float64(k) * c.PRadPeakAt1m / d
}

// Simultaneous returns how many cavities collapse together on one ring:
// rotor and housing cavities line up gcd(Nr, Ns) at a time.
func Simultaneous(nr, ns int) int {
	if ns == 0 {
		return 1
	}
	for ns != 0 {
		nr, ns = ns, nr%ns
	}
	return nr
}
