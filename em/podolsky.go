// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// Bopp–Podolsky generalized electrodynamics — HYPOTHETICAL (register H19).
//
// F. Bopp, Ann. Phys. 430, 345 (1940); B. Podolsky, "A generalized
// electrodynamics," Phys. Rev. 62, 68 (1942). The Lagrangian adds a term
// with higher derivatives of the field, −¼F_μνF^μν + (a²/2)∂_μF^μσ∂^νF_νσ,
// with one new constant: a length a (the "Podolsky mass" is m = ħ/(ac)).
// In Lorenz gauge the field equation becomes fourth order,
//
//	(1 + a²□) □ A^μ = μ0 J^μ        (□ = ∂²/c²∂t² − ∇²),
//
// so the second, third and fourth time derivatives of the potentials appear
// in the field equations themselves, not only through the material. The
// free field splits into the ordinary Maxwell wave (□A = 0, speed c) and a
// massive mode ((□ + 1/a²)A = 0) (form as given by R. Bufalo, B. M.
// Pimentel, D. E. Soto, Phys. Rev. D 90, 085012, 2014, arXiv:1407.1476,
// Sec. 2.2). As a → 0 it is Maxwell exactly. The same paper bounds the
// Podolsky mass from measured Bhabha (e⁺e⁻) scattering; that value is to be
// entered here from the paper's results section, not assumed.
//
// Rules (CLAUDE.md §1, §3f): off by default, every result labeled
// hypothetical, the standard (Maxwell, a = 0) prediction returned beside
// it. Point-charge results below are the static solutions.

// PodolskyPotential is the potential of a point charge q at distance r
// [V]: q/(4πε0 r)·(1 − e^(−r/a)); finite at r = 0, equal to q/(4πε0 a).
// a = 0 gives Coulomb.
func PodolskyPotential(q, r, a float64) float64 {
	k := q / (4 * math.Pi * Eps0)
	if a == 0 {
		return k / r
	}
	if r == 0 {
		return k / a
	}
	return k / r * -math.Expm1(-r/a)
}

// PodolskyField is the radial electric field of a point charge [V/m]:
// q/(4πε0 r²)·[1 − (1 + r/a)e^(−r/a)], which tends to the finite value
// q/(8πε0 a²) near the centre instead of diverging.
func PodolskyField(q, r, a float64) float64 {
	k := q / (4 * math.Pi * Eps0)
	if a == 0 {
		return k / (r * r)
	}
	x := r / a
	if x < 1e-4 { // series: 1 − (1+x)e^(−x) = x²/2 − x³/3 + …
		return k / (a * a) * (0.5 - x/3)
	}
	return k / (r * r) * (1 - (1+x)*math.Exp(-x))
}

// PodolskyRelativeDeviation is how far the potential at r falls below
// Coulomb's: e^(−r/a). An experiment that finds Coulomb's law to a
// fractional precision ε at distance r bounds a < r / ln(1/ε).
func PodolskyRelativeDeviation(r, a float64) float64 {
	if a == 0 {
		return 0
	}
	return math.Exp(-r / a)
}

// PodolskyBound is the largest a compatible with Coulomb's law holding to
// fractional precision eps at distance r.
func PodolskyBound(r, eps float64) float64 { return r / math.Log(1/eps) }

// PodolskySelfEnergy is the field energy of a point charge [J],
// q²/(8πε0 a): finite, where Maxwell's is infinite.
func PodolskySelfEnergy(q, a float64) float64 { return q * q / (8 * math.Pi * Eps0 * a) }

// PodolskyModes returns the two angular frequencies at wavenumber k: the
// Maxwell wave ck and the massive mode c√(k² + 1/a²).
func PodolskyModes(k, a float64) (maxwell, massive float64) {
	maxwell = SpeedOfLight * k
	if a == 0 {
		return maxwell, math.Inf(1)
	}
	return maxwell, SpeedOfLight * math.Sqrt(k*k+1/(a*a))
}

// PodolskyLength converts a Podolsky mass m (energy, in eV) to the length
// a = ħc/m [m].
func PodolskyLength(massEV float64) float64 {
	const hbarC = 197.3269804e-9 // eV·m (CODATA 2018)
	return hbarC / massEV
}
