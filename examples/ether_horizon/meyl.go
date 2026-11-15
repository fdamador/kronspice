// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Part 10: Konstantin Meyl's vortex picture beside H49.
//
// Meyl (Potential Vortex, vol. 2: Theory of Objectivity, 2nd ed.) adds a
// vortex of the electric field, the potential vortex, which contracts,
// against the eddy current's vortex, which expands; particles are
// spherical vortices in that balance. Two of his laws are testable here:
//
//   - the field sets the length scale: E, H ∝ 1/L² (his eq. 2.20, "the
//     field determines the dimensions!", §2.5);
//   - a particle's measured radius depends on the probe, because the
//     fields of probe and target add, (H_p + H_e) ∝ 1/R² (§3.5).
//
// His chapters on the structure and computation of the particles were
// moved to a third volume and are not quoted here; nothing below uses
// numbers from them.

// ---------------------------------------------------------------- 10a, 10b

// Contraction is a squeezing energy S(R) added to the light's ħc·x/R.
type Contraction struct {
	Name string
	S    func(R float64) float64 // energy [units of ħc·x/R0 at R = R0 = 1]
}

// Contractions are the laws scanned: the ether's pressure (R³) and
// tension (R²), a linear and a logarithmic squeeze, and two that weaken as
// the void grows, −1/R (Meyl's self-similar field: E ∝ 1/L² gives a field
// energy ∫E²dV ∝ L⁻⁴·L³ = 1/L) and −1/R².
func Contractions() []Contraction {
	return []Contraction{
		{"pressure, +p·R³", func(R float64) float64 { return R * R * R / 3 }},
		{"tension, +σ·R²", func(R float64) float64 { return R * R / 2 }},
		{"linear, +k·R", func(R float64) float64 { return R }},
		{"logarithmic, +k·ln R", func(R float64) float64 { return math.Log(R) }},
		{"Meyl self-similar, −k/R (k = ½)", func(R float64) float64 { return -0.5 / R }},
		{"Meyl self-similar, −k/R (k = 1)", func(R float64) float64 { return -1 / R }},
		{"steeper, −k/R² (k = ½)", func(R float64) float64 { return -0.5 / (R * R) }},
	}
}

// Balance classifies E(R) = 1/R + S(R) on R ∈ [1e-3, 1e3]: it returns the
// stationary points with the sign of E” at each (+1 stable, −1 unstable),
// and "none" if E' has no zero (E' ≡ 0 is reported as neutral).
func Balance(c Contraction) (roots []float64, kinds []int, neutral bool) {
	E := func(R float64) float64 { return 1/R + c.S(R) }
	d := func(R float64) float64 { h := 1e-6 * R; return (E(R+h) - E(R-h)) / (2 * h) }
	dd := func(R float64) float64 { h := 1e-4 * R; return (E(R+h) - 2*E(R) + E(R-h)) / (h * h) }
	n := 6000
	allZero := true
	lo0 := math.Log(1.37e-3) // off-centre grid, so no node sits on a root
	prev := lo0
	for i := 1; i <= n; i++ {
		lr := lo0 + float64(i)*math.Log(1e6)/float64(n)
		a, b := math.Exp(prev), math.Exp(lr)
		da, db := d(a), d(b)
		if math.Abs(da) > 1e-9*math.Abs(1/(a*a)) {
			allZero = false
		}
		if da*db < 0 {
			lo, hi := a, b
			for k := 0; k < 100; k++ {
				m := math.Sqrt(lo * hi)
				if d(m)*d(lo) > 0 {
					lo = m
				} else {
					hi = m
				}
			}
			r := math.Sqrt(lo * hi)
			roots = append(roots, r)
			if dd(r) > 0 {
				kinds = append(kinds, 1)
			} else {
				kinds = append(kinds, -1)
			}
		}
		prev = lr
	}
	return roots, kinds, allZero
}

// ---------------------------------------------------------------- 10c

// Measured hadron data.
const (
	Mp       = 1.67262192369e-27 // kg
	Mn       = 1.67492749804e-27 // kg
	RpCODATA = 0.84075e-15       // m, CODATA 2022 (PDG MOHR 2025)
	RpMuonic = 0.84087e-15       // m, μp Lamb shift (Antognini 2013)
	RpH2026  = 0.8406e-15        // m, hydrogen 2s–6p (2026)
	RpH2026U = 0.0015e-15
	RpMuU    = 0.00039e-15
	RpOld    = 0.8775e-15  // m, CODATA 2010 (e–p data), the value Meyl compares with
	RnMag    = 0.864e-15   // m, neutron magnetic rms radius (PDG)
	MuP      = 2.792847344 // μ_N
	MuN      = -1.91304276 // μ_N
)

// VoidAtMass is the void radius (TM₁₁ or TE₁₁; pressure- or
// tension-held) for a particle of mass m: the 9a balance scales with the
// reduced Compton wavelength ħ/(mc).
func VoidAtMass(m, x float64, tension bool) float64 {
	lc := Hbar / (m * C0)
	if tension {
		return 1.5 * x * lc
	}
	return 4.0 / 3 * x * lc
}

// UniformRMS is the rms radius of a uniform ball of radius R, √(3/5)·R —
// what a charge radius would read if the void were uniformly filled.
func UniformRMS(R float64) float64 { return math.Sqrt(0.6) * R }

// BagMomentRatio is the quark-model (SU(6)) ratio μ_n/μ_p = −2/3, the bag
// model's moments at the same radius.
const BagMomentRatio = -2.0 / 3

// ProbeBound is the 2σ bound on a probe-dependent proton radius from the
// electronic (hydrogen 2026) and muonic measurements: |Δr/r| ≤ bound.
func ProbeBound() (diff, bound float64) {
	d := RpH2026 - RpMuonic
	s := math.Hypot(RpH2026U, RpMuU)
	return d / RpMuonic, (math.Abs(d) + 2*s) / RpMuonic
}

// ---------------------------------------------------------------- 10d

// Onion is a chain of nested shells (void, light, void, light, …) as
// coupled resonators: shell k holds a cavity mode at ω_k = ω0·q^(−k) (a
// self-similar onion, radius ratio q < 1 inward), each an LC tank with L
// fixed and C_k = 1/(L·ω_k²), neighbours coupled through their shared wall,
// a capacitor of K·min(C_k, C_{k+1}) (so the smaller shell is detuned
// by only ~K).
type Onion struct {
	N     int
	Q     float64 // radius ratio between neighbouring shells (inner/outer)
	K     float64 // wall coupling
	L, F0 float64 // inductance [H] and outer shell's frequency [Hz]
}

// Freqs are the uncoupled shell frequencies [Hz].
func (o Onion) Freqs() []float64 {
	f := make([]float64, o.N)
	for k := range f {
		f[k] = o.F0 * math.Pow(o.Q, -float64(k))
	}
	return f
}

// Modes runs the chain in kronspice AC, driven by a 1 A source into the
// outer shell through a 10 MΩ coupling, and returns the peak frequencies
// of |V| summed over the shells (the chain's notes).
func (o Onion) Modes() ([]float64, error) {
	b := circuit.NewBuilder("onion of shells")
	fs := o.Freqs()
	cs := make([]float64, o.N)
	for k := 0; k < o.N; k++ {
		ks := strconv.Itoa(k)
		w := 2 * math.Pi * fs[k]
		cs[k] = 1 / (o.L * w * w)
		b.Inductor("L"+ks, "s"+ks, "0", o.L)
		b.Capacitor("C"+ks, "s"+ks, "0", cs[k])
		b.Resistor("R"+ks, "s"+ks, "0", 1e3*math.Sqrt(o.L/cs[k])) // Q ≈ 1000
	}
	for k := 0; k+1 < o.N; k++ {
		ks := strconv.Itoa(k)
		b.Capacitor("CW"+ks, "s"+ks, "s"+strconv.Itoa(k+1), o.K*math.Min(cs[k], cs[k+1]))
	}
	b.CurrentSource("IIN", "0", "s0", 0)
	c, err := b.Build()
	if err != nil {
		return nil, err
	}
	c.Elements[c.Index("IIN")].ACMag = 1
	res, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 4000, FStart: fs[0] / 3, FStop: fs[o.N-1] * 3}, analysis.Options{})
	if err != nil {
		return nil, err
	}
	amp := make([]float64, len(res.Freq))
	for k := 0; k < o.N; k++ {
		v, _ := res.NodeVoltage("s" + strconv.Itoa(k))
		for i := range amp {
			x := real(v[i])*real(v[i]) + imag(v[i])*imag(v[i])
			amp[i] += x * math.Pow(fs[k]/fs[0], 2) // weight inner shells up: their impedance is lower
		}
	}
	var peaks []float64
	for i := 1; i+1 < len(amp); i++ {
		if amp[i] > amp[i-1] && amp[i] > amp[i+1] {
			peaks = append(peaks, res.Freq[i])
		}
	}
	return peaks, nil
}

// Lepton mass ratios (PDG): m_μ/m_e and m_τ/m_μ.
const (
	MuOverE  = 206.7682830
	TauOverM = 1776.93 / 105.6583755
)
