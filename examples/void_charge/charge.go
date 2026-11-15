// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Measured constants (CODATA 2018 unless noted).
const (
	Alpha     = 1 / 137.035999084
	LbarM     = 3.8615926796e-13 // reduced Compton length ƛ = ħ/mc, m
	MC2Hz     = 1.2355899e20     // mc²/h, Hz
	MC2keV    = 510.99895        // mc², keV
	ECharge   = 4.80320471e-10   // e, statC (Gaussian)
	CLightCGS = 2.99792458e10    // c, cm/s
	HbarCGeV  = 0.1973269804     // ħc, GeV·fm
	// LambBound is how well the 2S Lamb shift pins any extra contact term:
	// measured 1057.845(9) MHz (Lundeen & Pipkin 1981; ±9 kHz).
	LambBound = 9e3
	// LEP2: e+e− → μ+μ−, ττ and Bhabha at √s ≈ 200 GeV agree with the standard
	// model to about 1–2 % (LEP EW Working Group 2013), so the electron's
	// form factor obeys |F|² ≳ 0.97, F ≳ 0.985.
	LEPsqrtS = 200.0
	LEPFmin  = 0.985
)

// DarwinR2 is the mean square charge radius that Dirac's Darwin term acts
// like, in units of ƛ²: the zitterbewegung smearing ¾ƛ² (rms 0.866 ƛ).
const DarwinR2 = 0.75

// ContactShiftHz is the shift of an nS level of hydrogen from a charge
// spread with mean square radius r2 (units of ƛ²), first order:
// ΔE = (2π/3)Zα ħc |ψ(0)|² ⟨r²⟩, |ψ_nS(0)|² = 1/(πn³a₀³), a₀ = ƛ/α, so
// ΔE = (2/3) α⁴ mc² (⟨r²⟩/ƛ²) / n³. With r2 = ¾ it is the Darwin term.
func ContactShiftHz(r2 float64, n int) float64 {
	return 2.0 / 3 * math.Pow(Alpha, 4) * MC2Hz * r2 / float64(n*n*n)
}

// LEPq is the momentum transfer at LEP2 in units of 1/ƛ.
func LEPq() float64 { return LEPsqrtS / HbarCGeV * (LbarM * 1e15) }

// LEPrmsBound is the largest rms charge radius (units of ƛ) LEP2 allows:
// F ≈ 1 − q²⟨r²⟩/6 ≥ F_min.
func LEPrmsBound() float64 { return math.Sqrt(6*(1-LEPFmin)) / LEPq() }

// UehlingRatio is QED's effective charge inside radius r (r in ƛ), relative
// to e, for r ≪ ƛ (Uehling 1935; Berestetskii, Lifshitz & Pitaevskii §114):
// α(r)/α = 1 + (2α/3π)(ln(1/r) − γ − 5/6). For r ≳ ƛ it tends to 1.
func UehlingRatio(r float64) float64 {
	if r >= 1 {
		return 1
	}
	v := 1 + 2*Alpha/(3*math.Pi)*(math.Log(1/r)-0.5772156649-5.0/6)
	return math.Max(v, 1)
}

// Tewari's charged void (Universal Principles of Space and Matter, 2017;
// §2.3 eq. 2.4, §2.6 eq. 2.14, §3.1 eqs. 3.1.3–3.1.7, §8.0–8.1).
type Tewari struct {
	ReCm      float64 // void radius from q = π² r_e² c (eq. 2.8), cm
	ReLbar    float64 // r_e / ƛ
	EnergyMC2 float64 // energy released per particle in annihilation, eq. 2.14: (4/5) mc²
	WaveOverL float64 // his annihilation light's wavelength r_e over the measured λ = h/mc
	LoopWork  float64 // work around a closed loop in his field, in units of e²/r_e (zero for any electrostatic field)
	Ly2pSplit float64 // spread of hydrogen 2p by orientation to the spin axis, eV
}

// TewariVoid computes the numbers that test his charged void.
func TewariVoid() Tewari {
	var t Tewari
	t.ReCm = math.Sqrt(ECharge / (math.Pi * math.Pi * CLightCGS))
	t.ReLbar = t.ReCm * 1e-2 / LbarM
	t.EnergyMC2 = 4.0 / 5
	t.WaveOverL = t.ReCm * 1e-2 / (2 * math.Pi * LbarM)
	// His field E_r = −A sin²θ / r², E_θ = 0 (eq. 3.1.3): a loop out along
	// θ = π/2 from r_e to 2r_e, along an arc to θ = 0 (E_θ = 0), back along
	// the axis (E = 0) and the arc home. Work ∫E·dl = −A(1/r_e − 1/2r_e) per
	// turn; with A normalised so the angle average is Coulomb's e/r², A = 3e/2:
	// |W| = (3/2)(1/2) e²/r_e.
	t.LoopWork = 0.75
	// Hydrogen 2p: his field is zero along the spin axis and 3/2 Coulomb's in
	// the equator (normalised to the same average). Along radial lines, as he
	// computes his potential (eq. 3.1.6), V = V_C·(3/2)sin²θ = V_C(1 − P₂).
	// For l = 1, ⟨P₂⟩ = 2/5 (m = 0) and −1/5 (|m| = 1) about the spin axis;
	// ⟨e²/r⟩_2p = 27.211/4 eV, so the m levels spread by (3/5)·6.80 eV.
	t.Ly2pSplit = 3.0 / 5 * 27.211386 / 4
	return t
}
