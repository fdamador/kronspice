// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Steps 0, 3 and 4: the electron test bench numbers, the charge options,
// and the size windows that fix κ2 and κ4.

import "math"

// CODATA 2018 unless stated.
const (
	hbarC    = 197.3269804     // MeV·fm
	meC2     = 0.51099895000   // MeV
	mpOverMe = 1836.15267343   //
	alpha    = 7.2973525693e-3 //
	lbarC    = hbarC / meC2    // reduced Compton wavelength, fm (386.159)
	muBh     = 13.996244936    // μ_B/h, GHz/T
	gQED     = 2.00231930436   // measured g (Fan et al., PRL 130, 071801, 2023)
	aMeas    = 0.00115965218059
	lymanA   = 121.56701 // nm, hydrogen Lyman-α (NIST)

	// Knot tone from examples/ether_knot (squashing mode, extrapolated
	// h → 0): ω = 1.36 c/ℓ, Q ≈ 6–7. Not recomputed here.
	toneW = 1.36

	// Bounds used (declared, cited):
	eStarMinMeV = 5.6e6 // excited electron m_e* > 5.6 TeV, CMS, JHEP 05 (2020) 052, Λ = m_e*
	sizeMaxFm   = 1e-3  // electron size < 10⁻¹⁸ m (register, H21 step 3)
	h1s2sTolMHz = 0.010 // hydrogen 1S–2S: theory–experiment agreement taken as 10 kHz (declared)
)

// Step 0: the standard twin. Dirac's electron read by Hestenes (Found.
// Phys. 20, 1213, 1990): a point charge circulating at c with radius
// ħ/2mc at angular frequency 2mc²/ħ. g = 2 comes from the Dirac current
// (Gordon decomposition: the spin magnetization), not from classical
// circulation: a charge e moving at c on a circle of radius r has
// μ = ecr/2, so with S = ħ/2 its g is 2μ/μ_B = 2r/λ̄.
type Twin struct {
	ZitterRadiusFm, ZitterOmega float64 // fm, rad/s
	GDirac, GQED                float64
	GClassicalAtZitter          float64 // 2r/λ̄ at r = ħ/2mc
	RadiusForG2Fm               float64 // circulation radius that gives g = 2 at S = ħ/2
}

func StandardTwin() Twin {
	r := lbarC / 2
	return Twin{
		ZitterRadiusFm:     r,
		ZitterOmega:        2 * meC2 / hbarC * 2.99792458e23, // 2mc²/ħ in rad/s (c in fm/s)
		GDirac:             2,
		GQED:               gQED,
		GClassicalAtZitter: GCirculation(r),
		RadiusForG2Fm:      lbarC,
	}
}

// GCirculation is the g of a charge e circulating at c on radius r (fm)
// with spin ħ/2: g = 2μ/μ_B, μ = ecr/2.
func GCirculation(rFm float64) float64 { return 2 * rFm / lbarC }

// WvdM is the Williamson–van der Mark ring: loop radius λ/4π = λ̄/2 and
// charge √(3ε0ħc)/2π (their eq. 5) in units of e.
func WvdMRadiusFm() float64 { return lbarC / 2 }
func WvdMCharge() float64   { return math.Sqrt(3/(4*math.Pi*alpha)) / (2 * math.Pi) }

// Charge options of step 3, each as q/e.
type ChargeOption struct {
	Name          string
	Q             float64 // q/e
	Anomaly       float64 // (q/e)²·α/2π (H23 eq. 12 with α′ = (q/e)²α)
	LymanNm       float64 // hydrogen Lyman-α with the electron charge q (levels ∝ (q/e)²)
	NetAtomCharge float64 // (1 − q/e) e for an atom with a proton of charge e
	Derived       bool
}

func ChargeOptions() []ChargeOption {
	mk := func(name string, q float64, derived bool) ChargeOption {
		return ChargeOption{name, q, q * q * alpha / (2 * math.Pi), lymanA / (q * q), 1 - q, derived}
	}
	return []ChargeOption{
		mk("H23 paper, √(3ε0ħc)/2π", WvdMCharge(), true),
		mk("Rañada topological unit √(ε0ħc)", 1/math.Sqrt(4*math.Pi*alpha), true),
		mk("declared e", 1, false),
	}
}

// Knot is what step 1–2 give for one knot, in E0 and ℓ units.
type Knot struct {
	Name       string
	Q          int
	Fermion    bool    // Finkelstein–Rubinstein: odd Q only
	E          float64 // static energy / E0
	Core       float64 // core radius / ℓ
	RmsRho     float64 // √⟨ρ²⟩ / ℓ (energy-weighted)
	RmsR       float64 // √⟨r²⟩ / ℓ
	IPar, IPer float64 // / (E0 ℓ²/c²)
}

// SpinC is c_r for the J = ½ state (only meaningful when Fermion).
func (k Knot) SpinC() float64 { return 1/(4*k.IPer) + 1/(8*k.IPar) }

// MinSizeSpin is the least ℓ (fm) for which the spin-½ knot can weigh
// m_e: m_e c² = E·E0 + c_r(ħc/ℓ)²/E0 has a real E0 only for
// ℓ ≥ 2√(E·c_r)·λ̄. At that ℓ half the mass is rotation.
func (k Knot) MinSizeSpin() float64 { return 2 * math.Sqrt(k.E*k.SpinC()) * lbarC }

// E0At solves for E0 (MeV) at ℓ (fm) with the spin energy included; ok is
// false when no real solution exists. The smaller-rotation root is taken.
func (k Knot) E0At(lFm float64) (e0 float64, ok bool) {
	b := meC2
	c := k.SpinC() * (hbarC / lFm) * (hbarC / lFm)
	disc := b*b - 4*k.E*c
	if disc < -1e-9*b*b {
		return 0, false
	}
	return (b + math.Sqrt(math.Max(disc, 0))) / (2 * k.E), true
}

// GRigid is the g of the J = ½ knot with its charge e spread like its
// energy and turning rigidly with it: μ = (e/2)ω⟨ρ²⟩, ω = ħ/(2I∥), so
// g = m⟨ρ²⟩/I∥ = E·⟨ρ²⟩/I∥ (static mass). Equal mass and charge
// distributions with a mechanical inertia would give exactly 1.
func (k Knot) GRigid() float64 { return k.E * k.RmsRho * k.RmsRho / k.IPar }

// Reading is one choice of ℓ with everything it predicts.
type Reading struct {
	Name          string
	LFm           float64
	RingFm, RmsFm float64
	E0MeV         float64 // static-only E0 = m_e c²/E
	ERotOverM     float64 // spin energy / m_e c² at that E0
	E0SpinOK      bool    // a real E0 exists with the spin included
	ToneMeV       float64 // ħω of the squashing tone
	CoulombOverM  float64 // αħc/(2·ring)/m_e c² (charge e on a shell at the ring)
	RimSpeedC     float64 // ω⊥·ring/c in the J = ½ state
	H1S2SMHz      float64 // finite-size shift of 1S–2S with ⟨r²⟩ = RmsFm²
	G, ESRGHzT    float64
}

func (k Knot) At(name string, lFm float64) Reading {
	e0 := meC2 / k.E
	r := Reading{Name: name, LFm: lFm, RingFm: k.Core * lFm, RmsFm: k.RmsR * lFm, E0MeV: e0}
	r.ERotOverM = k.SpinC() * (hbarC / lFm) * (hbarC / lFm) / e0 / meC2
	_, r.E0SpinOK = k.E0At(lFm)
	r.ToneMeV = toneW * hbarC / lFm
	r.CoulombOverM = alpha * hbarC / (2 * r.RingFm) / meC2
	// ω⊥ = ħ√(J(J+1))/I⊥, I⊥ = IPer·E0ℓ²/c²  →  ω⊥·ring/c = √¾·ħc·Core/(IPer·E0·ℓ)
	r.RimSpeedC = math.Sqrt(0.75) * hbarC * k.Core / (k.IPer * e0 * lFm)
	r.H1S2SMHz = H1S2SShiftMHz(r.RmsFm)
	r.G = k.GRigid()
	r.ESRGHzT = r.G * muBh
	return r
}

// H1S2SShiftMHz is the first-order finite-size shift of hydrogen's 1S–2S
// interval for a charge distribution of rms radius rFm (the electron's,
// adding to the proton's): ΔE_nS = (2/3n³)α⁴ m_r c² ⟨r²⟩/λ̄_r², and
// 1S–2S takes 7/8 of the 1S shift. For r = 0.8409 fm (the proton) it
// gives the familiar ≈ 1 MHz.
func H1S2SShiftMHz(rFm float64) float64 {
	mr := meC2 / (1 + 1/mpOverMe)
	lr := hbarC / mr
	d1s := 2.0 / 3 * math.Pow(alpha, 4) * mr * (rFm / lr) * (rFm / lr) // MeV
	return 7.0 / 8 * d1s * 1e6 / 4.135667696e-9                        // MeV → eV → MHz (h = 4.1357e-15 eV·s)
}

// Window is the range of ℓ (fm) each measured constraint allows.
type Window struct {
	Name   string
	Lo, Hi float64 // ℓ ∈ [Lo, Hi]; 0 / +Inf when open
}

func (k Knot) Windows() []Window {
	inf := math.Inf(1)
	// 1S–2S: shift ∝ ℓ²; find ℓ where it equals the tolerance.
	l1 := 1.0 // fm
	s1 := H1S2SShiftMHz(k.RmsR * l1)
	hMax := l1 * math.Sqrt(h1s2sTolMHz/s1)
	w := []Window{
		{"scattering: rms size < 10⁻¹⁸ m", 0, sizeMaxFm / k.RmsR},
		{"no excited electron below 5.6 TeV (CMS 2020)", 0, toneW * hbarC / eStarMinMeV},
		{"hydrogen 1S–2S finite size < 10 kHz", 0, hMax},
		{"Coulomb energy of e at the ring < m_e c²", alpha * lbarC / (2 * k.Core), inf},
	}
	if k.Fermion {
		w = append(w, Window{"spin ½ rotation fits in m_e c²", k.MinSizeSpin(), inf})
	}
	return w
}

// Overlap intersects the windows.
func Overlap(ws []Window) (lo, hi float64) {
	lo, hi = 0, math.Inf(1)
	for _, w := range ws {
		lo = math.Max(lo, w.Lo)
		hi = math.Min(hi, w.Hi)
	}
	return
}
