// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
)

// The Euler–Heisenberg vacuum: the nonlinear electrodynamics of the QED
// vacuum at fields well below the critical field, as a constitutive law.
//
// Source: W. Heisenberg, H. Euler, "Folgerungen aus der Diracschen
// Theorie des Positrons," Z. Phys. 98, 714–732 (1936); in SI (e.g.
// J. S. Toll, PhD thesis, Princeton 1952; R. Battesti, C. Rizzo, Rep.
// Prog. Phys. 76, 016401, 2013, eqs. 9–12):
//
//	L = ½ε₀(E² − c²B²) + ξ·[(E² − c²B²)² + 7c²(E·B)²],
//	ξ = 2α²ε₀²ħ³ / (45 m_e⁴ c⁵),
//
// the one-loop effective Lagrangian of a constant field, to fourth order;
// the next terms are smaller by (E/E_s)². It depends only on the two
// Lorentz invariants F = E² − c²B² and G = c·E·B, so it is the same in
// every inertial frame. The constitutive relations are its derivatives,
//
//	D = ∂L/∂E = ε₀E + ξ·[4F·E + 14c²G′·B]       (G′ = E·B)
//	H = −∂L/∂B = B/μ₀ + ξ·[4c²F·B − 14c²G′·E],
//
// and Maxwell's equations keep their form (∇·D = ρ, ∇×H = J + ∂D/∂t,
// ∇×E = −∂B/∂t, ∇·B = 0): the vacuum becomes a nonlinear medium. That is
// why it fits Kron's network directly: the topology (Faraday, Ampère) is
// unchanged and only the C and L values move (kron.MaxwellNetwork).
//
// What is measured. Light-by-light scattering, the process this
// Lagrangian describes at low energy, was observed in ultra-peripheral
// Pb–Pb collisions (ATLAS, Nature Phys. 13, 852, 2017; CMS, Phys. Lett.
// B 797, 134826, 2019). Vacuum magnetic birefringence, its laboratory
// signature Δn = 3a ≈ 4.0×10⁻²⁴ (B/1 T)², is not yet resolved: PVLAS
// reached Δn/B² = (1.9 ± 2.7)×10⁻²⁴ T⁻² (Ejlli et al., Phys. Rep. 871, 1,
// 2020).

// CODATA 2018 constants for the vacuum's QED scale.
const (
	ElementaryCharge = 1.602176634e-19  // e [C] (exact)
	ReducedPlanck    = 1.054571817e-34  // ħ [J·s] (exact)
	ElectronMass     = 9.1093837015e-31 // m_e [kg]
)

// FineStructure is α = e²/(4πε₀ħc).
var FineStructure = ElementaryCharge * ElementaryCharge / (4 * math.Pi * Eps0 * ReducedPlanck * SpeedOfLight)

// SchwingerField is the critical field E_s = m_e²c³/(eħ) ≈ 1.32×10¹⁸ V/m;
// CriticalB = E_s/c ≈ 4.41×10⁹ T.
var (
	SchwingerField = ElectronMass * ElectronMass * SpeedOfLight * SpeedOfLight * SpeedOfLight / (ElementaryCharge * ReducedPlanck)
	CriticalB      = SchwingerField / SpeedOfLight
)

// EHXi is ξ = 2α²ε₀²ħ³/(45 m_e⁴c⁵) [(J/m³)/(V/m)⁴], the coefficient of
// the quartic terms (≈ 1.3×10⁻⁵² in SI).
var EHXi = 2 * FineStructure * FineStructure * Eps0 * Eps0 * math.Pow(ReducedPlanck, 3) /
	(45 * math.Pow(ElectronMass, 4) * math.Pow(SpeedOfLight, 5))

// ErrEHStrongField reports a field outside the weak-field expansion.
var ErrEHStrongField = errors.New("em: field outside the Euler–Heisenberg weak-field range")

// EHWeakLimit is the largest max(|E|, c|B|)/E_s accepted: the neglected
// sixth-order terms are then below about 1 % of the quartic ones.
const EHWeakLimit = 0.1

// EHStrength is max(|E|, c|B|)/E_s.
func EHStrength(e, b [3]float64) float64 {
	return math.Max(math.Sqrt(dot(e, e)), SpeedOfLight*math.Sqrt(dot(b, b))) / SchwingerField
}

// EHLagrangian is the Lagrangian density L(E, B) [J/m³].
func EHLagrangian(e, b [3]float64) float64 {
	c2 := SpeedOfLight * SpeedOfLight
	f := dot(e, e) - c2*dot(b, b)
	g := dot(e, b)
	return 0.5*Eps0*f + EHXi*(f*f+7*c2*g*g)
}

// EHConstitutive returns D [C/m²] and H [A/m] for E [V/m] and B [T].
func EHConstitutive(e, b [3]float64) (d, h [3]float64) {
	c2 := SpeedOfLight * SpeedOfLight
	f := dot(e, e) - c2*dot(b, b)
	g := dot(e, b)
	for i := 0; i < 3; i++ {
		d[i] = Eps0*e[i] + EHXi*(4*f*e[i]+14*c2*g*b[i])
		h[i] = b[i]/Mu0 + EHXi*(4*c2*f*b[i]-14*c2*g*e[i])
	}
	return d, h
}

// EHEnergyDensity is the field energy u = E·D − L [J/m³] (for a static
// field, the work ∫E·dD + ∫H·dB done to build it).
func EHEnergyDensity(e, b [3]float64) float64 {
	d, _ := EHConstitutive(e, b)
	return dot(e, d) - EHLagrangian(e, b)
}

// EHSmallSignal is the Euler–Heisenberg vacuum linearized about a static
// background (E0, B0): for small fields e, b on top of it,
//
//	δD = Eps·e + Xeb·b,   δH = Xhe·e + Nu·b,
//
// Eps = ∂D/∂E the permittivity tensor [F/m], Nu = ∂H/∂B the inverse
// permeability tensor [m/H], and Xeb = ∂D/∂B, Xhe = ∂H/∂E the
// magnetoelectric terms (both zero when E0 or B0 is zero).
type EHSmallSignal struct {
	Eps, Nu, Xeb, Xhe [3][3]float64
}

// EHLinearize returns the small-signal tensors about (E0, B0), or
// ErrEHStrongField beyond EHWeakLimit.
func EHLinearize(e0, b0 [3]float64) (EHSmallSignal, error) {
	var s EHSmallSignal
	if r := EHStrength(e0, b0); !(r <= EHWeakLimit) {
		return s, fmt.Errorf("%w: max(|E|, c|B|)/E_s = %.3g > %g", ErrEHStrongField, r, EHWeakLimit)
	}
	c2 := SpeedOfLight * SpeedOfLight
	f := dot(e0, e0) - c2*dot(b0, b0)
	g := dot(e0, b0)
	x := EHXi
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			dij := 0.0
			if i == j {
				dij = 1
			}
			// D_i = ε₀E_i + ξ(4F·E_i + 14c²G·B_i), F = E² − c²B², G = E·B
			s.Eps[i][j] = Eps0*dij + x*(4*f*dij+8*e0[i]*e0[j]+14*c2*b0[i]*b0[j])
			s.Xeb[i][j] = x * (-8*c2*e0[i]*b0[j] + 14*c2*(e0[j]*b0[i]+g*dij))
			// H_i = B_i/μ₀ + ξ(4c²F·B_i − 14c²G·E_i)
			s.Nu[i][j] = dij/Mu0 + x*(4*c2*f*dij-8*c2*c2*b0[i]*b0[j]-14*c2*e0[i]*e0[j])
			s.Xhe[i][j] = x * (8*c2*b0[i]*e0[j] - 14*c2*(b0[j]*e0[i]+g*dij))
		}
	}
	return s, nil
}

// EHBirefringenceParameter is a = (α/90π)(B/B_c)² = ξc²B²/ε₀ for a
// magnetic background B [T].
func EHBirefringenceParameter(b float64) float64 {
	return FineStructure / (90 * math.Pi) * (b / CriticalB) * (b / CriticalB)
}

// EHIndices returns the refractive indices of light crossing a magnetic
// background B [T] at right angles, to first order in a: n∥ = 1 + 7a with
// the light's electric field along B, n⊥ = 1 + 4a across it (Toll 1952;
// S. L. Adler, Ann. Phys. 67, 599, 1971).
func EHIndices(b float64) (nPar, nPerp float64) {
	a := EHBirefringenceParameter(b)
	return 1 + 7*a, 1 + 4*a
}
