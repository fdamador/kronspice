// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Eq3Check builds Kaluza's metric from smooth sample fields (a weak 4-D
// metric near −δ, a potential q_κ and corner g, all quadratic in x) and
// compares his three-index symbols (eq. 1) with his closed forms (eq. 3):
// Γ_0κλ = αF_κλ, Γ_κλ0 = −αΣ_κλ, Γ_00κ = −Γ_0κ0 = g,κ, Γ_000 = 0, and
// Γ_κλμ unchanged by the electromagnetic framing. Returns the largest
// difference.
func Eq3Check(alpha float64) float64 {
	c := [4][4]float64{{0.3, -0.2, 0.1, 0.4}, {0.05, 0.2, -0.3, 0.1}, {-0.1, 0.15, 0.25, -0.2}, {0.2, 0.1, -0.05, 0.3}}
	q := func(x []float64) [4]float64 {
		var o [4]float64
		for k := 0; k < 4; k++ {
			for l := 0; l < 4; l++ {
				o[k] += c[k][l]*x[l+1] + 0.1*c[l][k]*x[l+1]*x[k+1]
			}
		}
		return o
	}
	g4 := func(x []float64) [4][4]float64 {
		var m [4][4]float64
		for k := 0; k < 4; k++ {
			m[k][k] = -1
			for l := 0; l < 4; l++ {
				m[k][l] += 0.01 * (c[k][l] + c[l][k]) * x[1+(k+l)%4]
			}
		}
		return m
	}
	g := func(x []float64) float64 { return -0.5 + 0.02*x[1]*x[2] + 0.03*x[3] - 0.01*x[4]*x[4] }
	metric := KaluzaMetric(g4, q, g, alpha)
	plain := KaluzaMetric(g4, func([]float64) [4]float64 { return [4]float64{} }, g, alpha)
	x := []float64{7, 0.4, -0.3, 0.2, 0.6} // x⁰ arbitrary: nothing depends on it
	h := 1e-4
	G := KaluzaGamma(metric, x, h)
	G0 := KaluzaGamma(plain, x, h)
	dq := func(k, l int) float64 { // q_κ,λ
		xp := append([]float64(nil), x...)
		xm := append([]float64(nil), x...)
		xp[l+1] += h
		xm[l+1] -= h
		return (q(xp)[k] - q(xm)[k]) / (2 * h)
	}
	dg := func(k int) float64 {
		xp := append([]float64(nil), x...)
		xm := append([]float64(nil), x...)
		xp[k+1] += h
		xm[k+1] -= h
		return (g(xp) - g(xm)) / (2 * h)
	}
	worst := math.Abs(G[0][0][0])
	for k := 0; k < 4; k++ {
		worst = math.Max(worst, math.Abs(G[0][0][k+1]-dg(k)))
		worst = math.Max(worst, math.Abs(G[0][k+1][0]+dg(k)))
		for l := 0; l < 4; l++ {
			F := dq(k, l) - dq(l, k)
			S := dq(k, l) + dq(l, k)
			worst = math.Max(worst, math.Abs(G[0][k+1][l+1]-alpha*F))
			worst = math.Max(worst, math.Abs(G[k+1][l+1][0]+alpha*S))
			for mu := 0; mu < 4; mu++ {
				worst = math.Max(worst, math.Abs(G[k+1][l+1][mu+1]-G0[k+1][l+1][mu+1]))
			}
		}
	}
	return worst
}

// Constants (CGS, CODATA 2018) for Kaluza's numbers.
const (
	Gcgs      = 6.67430e-8      // cm³ g⁻¹ s⁻²
	Ccgs      = 2.99792458e10   // cm/s
	Ecgs      = 4.803204712e-10 // esu
	MeCgs     = 9.1093837015e-28
	MpCgs     = 1.67262192369e-24
	HbarC     = 1.054571817e-27 // erg·s
	AlphaFS   = 7.2973525693e-3
	LPlanckM  = 1.616255e-35 // m
	GeVPerErg = 624.150907
)

// KaluzaAlpha is his α = √(κ/2) with κ = 8πG/c² (eq. before 8a).
func KaluzaAlpha() float64 { return math.Sqrt(8 * math.Pi * Gcgs / (Ccgs * Ccgs) / 2) }

// ChargeToMass returns q/(m√G) (dimensionless, Gaussian) for charge q and mass m.
func ChargeToMass(q, m float64) float64 { return q / (m * math.Sqrt(Gcgs)) }

// KleinRadius is the radius of Klein's circle (Z. Phys. 37, 895, 1926)
// that makes p_w = n·ħ/R carry charge n·e: q = k·P_w/c with k = 2√G
// gives R = 2√G·ħ/(c·e) = 2·l_P/√α [m].
func KleinRadius() float64 {
	return 2 * math.Sqrt(Gcgs) * HbarC / (Ccgs * Ecgs) / 100
}

// KKMassGeV is the mass of the first charged Kaluza–Klein mode, ħ/(R·c),
// in GeV/c².
func KKMassGeV() float64 {
	return HbarC * Ccgs / (KleinRadius() * 100) * GeVPerErg
}
