// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// QED's correction to g: a_e = (g − 2)/2 = Σ_n A^(2n) (α/π)^n + hadronic + weak.

// SchwingerF2 is the one-loop vertex's magnetic form factor in Feynman
// parameters (Peskin & Schroeder eqs. 6.47, 6.56), as the coefficient of α/π:
//
//	F2(q²) = (α/2π) ∫dx dy dz δ(x+y+z−1) · 2m²z(1−z)/Δ,  Δ = (1−z)²m² + xy·Q²,
//
// with Q² = −q² ≥ 0 (spacelike) in units of m². A midpoint double sum over
// z and x (y = 1 − x − z). At Q² = 0 it is ½: Schwinger's α/2π (1948).
func SchwingerF2(Q2 float64, n int) float64 {
	s := 0.0
	hz := 1 / float64(n)
	for i := 0; i < n; i++ {
		z := (float64(i) + 0.5) * hz
		w := 1 - z
		hx := w / float64(n)
		for j := 0; j < n; j++ {
			x := (float64(j) + 0.5) * hx
			y := w - x
			d := w*w + x*y*Q2
			s += 2 * z * (1 - z) / d * hx * hz
		}
	}
	return s / 2
}

// Zeta3 is ζ(3) by its alternating series ζ(3) = (4/3)·Σ (−1)^{k+1}/k³,
// summed with averaging of partial sums.
func Zeta3() float64 {
	// direct sum with Euler–Maclaurin tail: Σ_{k≤N} 1/k³ + 1/(2N²) − 1/(2N³) + 1/(4N⁴)
	N := 2000
	s := 0.0
	for k := N; k >= 1; k-- {
		fk := float64(k)
		s += 1 / (fk * fk * fk)
	}
	fN := float64(N)
	return s + 1/(2*fN*fN) - 1/(2*fN*fN*fN) + 1/(4*fN*fN*fN*fN)
}

// A4Petermann is the two-loop coefficient in closed form (Petermann 1957,
// Sommerfield 1957): 197/144 + π²/12 − (π²/2)·ln 2 + (3/4)·ζ(3).
func A4Petermann() float64 {
	return 197.0/144 + math.Pi*math.Pi/12 - math.Pi*math.Pi/2*math.Ln2 + 0.75*Zeta3()
}

// Literature terms (not computed here).
const (
	A6     = 1.181241456587 // Laporta & Remiddi, PLB 379, 283 (1996)
	A8     = -1.9122457649  // Laporta, PLB 772, 232 (2017)
	A10    = 6.737          // Aoyama, Kinoshita, Nio, Atoms 7, 28 (2019); ±0.159, disputed by Volkov (arXiv:2404.00649)
	A10Err = 0.159
	// mass-dependent (muon, tau loops): Aoyama, Hayakawa, Kinoshita, Nio, PRL 109, 111807 (2012)
	A4Mu, A4Tau = 5.19738667e-7, 1.83798e-9
	A6Mu, A6Tau = -7.37394155e-6, -6.5830e-8
	A8Mu, A8Tau = 9.222e-4, 8.24e-6
	A8MuTau     = 7.465e-7
	A10Mu       = -0.00382
	// hadronic (LO + NLO vacuum polarization + light-by-light) and weak, same paper
	AHadLO, AHadNLO = 1.875e-12, -0.225e-12
	AHadLbL         = 0.035e-12 // light-by-light (Prades, de Rafael & Vainshtein 2009), quoted by AHKN
	AWeak           = 0.0297e-12
)

// Alpha inverse from atom recoil.
const (
	AlphaInvRb = 137.035999206 // Morel et al., Nature 588, 61 (2020), ±0.000000011
	AlphaInvCs = 137.035999046 // Parker et al., Science 360, 191 (2018), ±0.000000027
	// AeMeasured is (g − 2)/2 from Fan et al., PRL 130, 071801 (2023), ±0.13e-12.
	AeMeasured    = 1.15965218059e-3
	AeMeasuredErr = 0.13e-12
)

// QEDTerm is one order of the series at a given α.
type QEDTerm struct {
	Name  string
	Coef  float64 // coefficient of (α/π)^n (mass-dependent parts included)
	Order int
	Value float64
}

// AeSeries evaluates a_e(α) term by term.
func AeSeries(alphaInv float64) ([]QEDTerm, float64) {
	x := 1 / alphaInv / math.Pi
	terms := []QEDTerm{
		{"one loop (Schwinger), computed", SchwingerF2(0, 400), 1, 0},
		{"two loops (Petermann closed form), computed", A4Petermann() + A4Mu + A4Tau, 2, 0},
		{"three loops (Laporta & Remiddi)", A6 + A6Mu + A6Tau, 3, 0},
		{"four loops (Laporta; AHKN)", A8 + A8Mu + A8Tau + A8MuTau, 4, 0},
		{"five loops (AHKN, disputed)", A10 + A10Mu, 5, 0},
	}
	tot := 0.0
	for i := range terms {
		terms[i].Value = terms[i].Coef * math.Pow(x, float64(terms[i].Order))
		tot += terms[i].Value
	}
	had := AHadLO + AHadNLO + AHadLbL
	terms = append(terms, QEDTerm{"hadronic", 0, 0, had}, QEDTerm{"weak", 0, 0, AWeak})
	tot += had + AWeak
	return terms, tot
}

// GFromAe is g = 2(1 + a_e).
func GFromAe(a float64) float64 { return 2 * (1 + a) }
