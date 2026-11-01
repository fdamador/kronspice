// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// Unit conversions.
const (
	avogadro = 6.02214076e23
	amu      = 1.66053906660e-27 // kg
	cCm      = 2.99792458e10     // speed of light, cm/s
	// kcalMolA2ToMdynA converts kcal/mol/Å² to mdyn/Å (= 100 N/m).
	kcalMolA2ToMdynA = 4184 / avogadro / 1e-20 / 100
)

// wavenumber converts √(stiffness[mdyn/Å]/mass[u]) to cm⁻¹ (as in
// examples/water_vibration).
var wavenumber = math.Sqrt(100/amu) / (2 * math.Pi * cCm)

// Mode is one normal mode.
type Mode struct {
	Wavenumber float64   // cm⁻¹ (negative: imaginary, a saddle direction)
	Vector     []float64 // Cartesian displacement per unit normal coordinate, L_k/√m (Å·u^-½)
	Intensity  float64   // IR intensity with the fixed AMBER charges [km/mol]
	DipoleDer  [3]float64
	Inter      float64 // fraction of the mode's potential energy in the pairing (H-bond) terms
	Stack      float64 // fraction in the forces between stacked pairs
}

// State is a minimized molecule with its stiffness.
type State struct {
	Mol   *Molecule
	X     []float64   // minimized coordinates, Å
	K     [][]float64 // Hessian, mdyn/Å
	KPart map[int][][]float64
	RMS   float64 // RMS gradient after minimization, kcal/mol/Å
	Modes []Mode  // the 3N−6 vibrations, ascending
}

// Prepare minimizes m and computes its Hessian and normal modes.
func Prepare(m *Molecule) (*State, error) {
	x := m.Coords()
	rms, err := m.Minimize(x)
	if err != nil {
		return nil, err
	}
	if rms > 1e-6 {
		return nil, fmt.Errorf("%s: minimization ended at rms gradient %.3g kcal/mol/Å", m.Title, rms)
	}
	s := &State{Mol: m, X: x, RMS: rms, KPart: map[int][][]float64{}}
	conv := func(h [][]float64) [][]float64 {
		for i := range h {
			for j := range h[i] {
				h[i][j] *= kcalMolA2ToMdynA
			}
		}
		return h
	}
	s.K = conv(m.Hessian(x, PartAll, 1e-4))
	for _, p := range []int{PartBond, PartAngle, PartTorsion, PartIntra, PartInter, PartStack} {
		s.KPart[p] = conv(m.Hessian(x, p, 1e-4))
	}
	s.Modes = NormalModes(m, s.K, s.KPart[PartInter])
	for k := range s.Modes {
		s.Modes[k].Stack = share(s.Modes[k].Vector, s.KPart[PartStack], s.K)
	}
	return s, nil
}

// NormalModes diagonalizes the mass-weighted Hessian K (mdyn/Å) and drops
// the six rigid-body modes (the six smallest |λ|). kInter, when not nil,
// gives each mode's share of potential energy in the inter-base terms.
func NormalModes(m *Molecule, K [][]float64, kInter [][]float64) []Mode {
	n := len(K)
	mass := make([]float64, n)
	for i := range mass {
		mass[i] = m.Atoms[i/3].Mass
	}
	h := make([][]float64, n)
	for i := range h {
		h[i] = make([]float64, n)
		for j := range h[i] {
			h[i][j] = K[i][j] / math.Sqrt(mass[i]*mass[j])
		}
	}
	vals, vecs := Eigen(h)
	// Rigid-body modes: the six eigenvalues closest to zero.
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	for i := 1; i < n; i++ {
		for j := i; j > 0 && math.Abs(vals[idx[j]]) < math.Abs(vals[idx[j-1]]); j-- {
			idx[j], idx[j-1] = idx[j-1], idx[j]
		}
	}
	rigid := map[int]bool{}
	for _, k := range idx[:6] {
		rigid[k] = true
	}
	var out []Mode
	for k := 0; k < n; k++ {
		if rigid[k] {
			continue
		}
		md := Mode{Vector: make([]float64, n)}
		w := math.Sqrt(math.Abs(vals[k])) * wavenumber
		if vals[k] < 0 {
			w = -w
		}
		md.Wavenumber = w
		for i := 0; i < n; i++ {
			md.Vector[i] = vecs[i][k] / math.Sqrt(mass[i])
		}
		// Dipole derivative with fixed charges: dμ/dQ = Σ q_a L_a/√m_a [e/√u].
		for a := range m.Atoms {
			for d := 0; d < 3; d++ {
				md.DipoleDer[d] += m.Atoms[a].Charge * md.Vector[3*a+d]
			}
		}
		md.Intensity = IRIntensity(md.DipoleDer)
		if kInter != nil {
			num, den := 0.0, 0.0
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					num += md.Vector[i] * kInter[i][j] * md.Vector[j]
					den += md.Vector[i] * K[i][j] * md.Vector[j]
				}
			}
			md.Inter = num / den
		}
		out = append(out, md)
	}
	return out
}

// share is v·P·v / v·K·v: the fraction of a mode's potential energy in part P.
func share(v []float64, P, K [][]float64) float64 {
	num, den := 0.0, 0.0
	for i := range v {
		for j := range v {
			num += v[i] * P[i][j] * v[j]
			den += v[i] * K[i][j] * v[j]
		}
	}
	return num / den
}

// IRIntensity returns the integrated absorption intensity [km/mol] of a
// fundamental with dipole derivative dμ/dQ [e/√u]:
// A = N_A π |dμ/dQ|² / (3 c² · 4πε₀) (Wilson, Decius & Cross §7-3; in
// these units 974.86 km/mol per (e²/u), the factor used by quantum
// chemistry codes for intensities from atomic polar tensors).
func IRIntensity(d [3]float64) float64 {
	const e = 1.602176634e-19
	const eps0 = 8.8541878128e-12
	const c = 2.99792458e8
	s := (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) * e * e / amu         // C²/kg
	a := avogadro * math.Pi * s / (3 * c * c * 4 * math.Pi * eps0) // m/mol
	return a / 1000                                                // km/mol
}

// EinsteinA returns the spontaneous emission rate [1/s] of the v=1→0
// fundamental of a mode at ν̃ [cm⁻¹] with intensity S [km/mol]:
// A = 8π c ν̃² ∫σ dν̃, ∫σ dν̃ = S·10⁵ cm/mol ÷ N_A (Hilborn, Am. J. Phys.
// 50, 982, 1982, eq. 8, non-degenerate levels).
func EinsteinA(nu, S float64) float64 {
	sigma := S * 1e5 / avogadro // cm per molecule
	return 8 * math.Pi * cCm * nu * nu * sigma
}
