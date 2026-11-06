// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// The hydrogen molecule in the Heitler–London picture (Heitler & London,
// Z. Phys. 44, 455, 1927): two 1s orbitals a, b on protons at z = ∓R/2,
// spatial part a(1)b(2) ± b(1)a(2). The + (singlet, electron spins
// opposite) binds; the − (triplet, spins parallel) repels. Atomic units
// (hartree, a0) throughout.

// Hartree in eV (CODATA 2018).
const Hartree = 27.211386245988

// Overlap and one- and two-electron integrals over 1s orbitals at
// separation R (Slater, Quantum Theory of Molecules and Solids I, 1963,
// App. 6; Sugiura, Z. Phys. 45, 484, 1927, for the exchange integral).
func hlIntegrals(R float64) (S, jab, kab, J2, K2 float64) {
	e := math.Exp(-R)
	S = e * (1 + R + R*R/3)
	jab = 1/R - math.Exp(-2*R)*(1+1/R) // ⟨a|1/r_b|a⟩
	kab = e * (1 + R)                  // ⟨a|1/r_a|b⟩
	J2 = 1/R - math.Exp(-2*R)*(1/R+11.0/8+3*R/4+R*R/6)
	Sp := math.Exp(R) * (1 - R + R*R/3)
	const gamma = 0.57721566490153286
	K2 = (-math.Exp(-2*R)*(-25.0/8+23*R/4+3*R*R+R*R*R/3) +
		6/R*(S*S*(gamma+math.Log(R))+Sp*Sp*Ei(-4*R)-2*S*Sp*Ei(-2*R))) / 5
	return
}

// HLEnergy is the Heitler–London energy [hartree] of the singlet (+1) or
// triplet (−1) at separation R, including the proton repulsion:
// E = (H_aa ± H_ab)/(1 ± S²),
// H_aa = 2E_H − 2⟨a|1/r_b|a⟩ + J′ + 1/R,
// H_ab = 2E_H S² − 2S⟨a|1/r_a|b⟩ + K′ + S²/R.
func HLEnergy(R float64, sign float64) float64 {
	S, jab, kab, J2, K2 := hlIntegrals(R)
	const eh = -0.5
	haa := 2*eh - 2*jab + J2 + 1/R
	hab := 2*eh*S*S - 2*S*kab + K2 + S*S/R
	return (haa + sign*hab) / (1 + sign*S*S)
}

// Ei is the exponential integral for x < 0, Ei(x) = −E1(−x) (series for
// small |x|, continued fraction otherwise; Numerical Recipes §6.3).
func Ei(x float64) float64 {
	return -E1(-x)
}

// E1 is the exponential integral E1(x), x > 0.
func E1(x float64) float64 {
	const gamma = 0.57721566490153286
	if x <= 1 {
		sum, term := 0.0, 1.0
		for k := 1; k < 60; k++ {
			term *= -x / float64(k)
			sum -= term / float64(k)
		}
		return -gamma - math.Log(x) + sum
	}
	b := x + 1
	c := 1e300
	d := 1 / b
	h := d
	for i := 1; i < 200; i++ {
		an := -float64(i * i)
		b += 2
		d = 1 / (an*d + b)
		c = b + an/c
		del := c * d
		h *= del
		if math.Abs(del-1) < 1e-16 {
			break
		}
	}
	return h * math.Exp(-x)
}

// Minimum finds the bond length and well depth of the singlet curve
// (golden section on R ∈ [1, 3]).
func HLMinimum() (re, de float64) {
	f := func(r float64) float64 { return HLEnergy(r, 1) }
	a, b := 1.0, 3.0
	g := (math.Sqrt(5) - 1) / 2
	for i := 0; i < 200; i++ {
		c, d := b-g*(b-a), a+g*(b-a)
		if f(c) < f(d) {
			b = d
		} else {
			a = c
		}
	}
	re = (a + b) / 2
	return re, -1 - f(re)
}

// Molecule is H₂ at separation R in a chosen spin state, spins quantized
// along the bond (z) axis so the field is axisymmetric.
type Molecule struct {
	Name    string
	R       float64
	Triplet bool // electrons M_S = +1 (both up) — else singlet, spin density 0
	Ortho   bool // protons M_I = +1 (both up) — else para (I = 0), moment 0
	Bound   bool
}

// Molecules are the four spin arrangements of the four small magnets.
func Molecules(R float64) []Molecule {
	return []Molecule{
		{"singlet electrons + para protons", R, false, false, true},
		{"singlet electrons + ortho protons", R, false, true, true},
		{"triplet electrons + para protons", R, true, false, false},
		{"triplet electrons + ortho protons (all four aligned)", R, true, true, false},
	}
}

// orbitals: 1s on each proton, a at z = −R/2, b at z = +R/2.
func (m Molecule) orb(rho, z float64) (a, b, ra, rb float64) {
	ra = math.Hypot(rho, z+m.R/2)
	rb = math.Hypot(rho, z-m.R/2)
	return math.Exp(-ra) / math.Sqrt(math.Pi), math.Exp(-rb) / math.Sqrt(math.Pi), ra, rb
}

func (m Molecule) overlap() float64 {
	R := m.R
	return math.Exp(-R) * (1 + R + R*R/3)
}

// Density is the electron number density of the HL state (2 electrons):
// (a² + b² ± 2S ab)/(1 ± S²), + singlet, − triplet.
func (m Molecule) Density(rho, z float64) float64 {
	a, b, _, _ := m.orb(rho, z)
	S := m.overlap()
	if m.Triplet {
		return (a*a + b*b - 2*S*a*b) / (1 - S*S)
	}
	return (a*a + b*b + 2*S*a*b) / (1 + S*S)
}

// SpinDensity is ⟨Σ σ_z δ(r − r_i)⟩: the number density for M_S = +1
// (both spins up), and exactly zero for the singlet, whose spin part
// (↑↓ − ↓↑)/√2 gives ⟨σ_z⁽¹⁾⟩ = ⟨σ_z⁽²⁾⟩ = 0 at every point.
func (m Molecule) SpinDensity(rho, z float64) float64 {
	if !m.Triplet {
		return 0
	}
	return m.Density(rho, z)
}

// J is the magnetization current of the electron spins, M = −(g_s/2)μB s_z ẑ:
// J_φ = (∇×M)_φ = −∂M_z/∂ρ = (g_s/2) ∂s_z/∂ρ (analytic: ∂a/∂ρ = −a ρ/r_a).
func (m Molecule) J(rho, z float64) float64 {
	if !m.Triplet {
		return 0
	}
	a, b, ra, rb := m.orb(rho, z)
	S := m.overlap()
	da, db := -a*rho/ra, -b*rho/rb
	ds := (2*a*da + 2*b*db - 2*S*(da*b+a*db)) / (1 - S*S)
	return GS / 2 * ds
}

// ProtonMoments are the two proton dipoles (μB) at z = ∓R/2: +μp each for
// ortho M_I = +1; zero for para, whose nuclear singlet has ⟨I_z⟩ = 0 on each
// proton.
func (m Molecule) ProtonMoments() (za, ma, zb, mb float64) {
	if m.Ortho {
		return -m.R / 2, MuPr, m.R / 2, MuPr
	}
	return -m.R / 2, 0, m.R / 2, 0
}

// MolField is the molecule's field: the electron spin current as coaxial
// loops, plus the two proton dipoles (analytic). Inside the cloud the loop
// field is evaluated on the corners of the source cells, where the four
// neighbouring loops' near fields cancel by symmetry, and interpolated
// bilinearly; outside it the loops are summed directly.
type MolField struct {
	M     Molecule
	Loops Loops
	G     *Grid
}

// Grid is a loop field sampled on cell corners.
type Grid struct {
	H, Rmax, Zmax float64
	Nr, Nz        int
	Br, Bz, Psi   []float64
}

var gridCache = map[string]*Grid{}
var loopCache = map[string]Loops{}

// NewMolField builds the loop grid (cell h, extent ext) and the corner grid
// (extent gext). The electron part depends only on the electron spin state
// and R, so it is shared between the para and ortho cases.
func NewMolField(m Molecule, h, ext, gext float64) MolField {
	f := MolField{M: m}
	if !m.Triplet {
		return f
	}
	key := fmt.Sprintf("%g/%g/%g/%g", m.R, h, ext, gext)
	if l, ok := loopCache[key]; ok {
		f.Loops, f.G = l, gridCache[key]
		return f
	}
	f.Loops = NewLoops(m.J, h, ext, -ext, ext)
	g := &Grid{H: h, Rmax: gext, Zmax: gext, Nr: int(math.Round(gext/h)) + 1, Nz: int(math.Round(2*gext/h)) + 1}
	g.Br, g.Bz, g.Psi = make([]float64, g.Nr*g.Nz), make([]float64, g.Nr*g.Nz), make([]float64, g.Nr*g.Nz)
	parallel(g.Nr*g.Nz, func(k int) {
		i, j := k/g.Nz, k%g.Nz
		g.Br[k], g.Bz[k], g.Psi[k] = f.Loops.Field(float64(i)*h, -gext+float64(j)*h)
	})
	f.G = g
	loopCache[key], gridCache[key] = f.Loops, g
	return f
}

func (g *Grid) at(rho, z float64) (br, bz, psi float64, ok bool) {
	x, y := rho/g.H, (z+g.Zmax)/g.H
	i, j := int(x), int(y)
	if i < 0 || j < 0 || i >= g.Nr-1 || j >= g.Nz-1 {
		return 0, 0, 0, false
	}
	tx, ty := x-float64(i), y-float64(j)
	w := [4]float64{(1 - tx) * (1 - ty), tx * (1 - ty), (1 - tx) * ty, tx * ty}
	idx := [4]int{i*g.Nz + j, (i+1)*g.Nz + j, i*g.Nz + j + 1, (i+1)*g.Nz + j + 1}
	for q := range w {
		br += w[q] * g.Br[idx[q]]
		bz += w[q] * g.Bz[idx[q]]
		psi += w[q] * g.Psi[idx[q]]
	}
	return br, bz, psi, true
}

// Field returns B_ρ, B_z [B0] and the flux function at (ρ, z).
func (f MolField) Field(rho, z float64) (br, bz, psi float64) {
	if f.M.Triplet {
		var ok bool
		if br, bz, psi, ok = f.G.at(rho, z); !ok {
			br, bz, psi = f.Loops.Field(rho, z)
		}
	}
	za, ma, zb, mb := f.M.ProtonMoments()
	for _, p := range [][2]float64{{za, ma}, {zb, mb}} {
		if p[1] == 0 {
			continue
		}
		a, b, c := Dipole(p[1], p[0], rho, z)
		br, bz, psi = br+a, bz+b, psi+c
	}
	return
}

// Moment is the total magnetic moment along z [μB]: −(g_s/2)·2 for the
// triplet electrons plus the protons.
func (m Molecule) Moment() float64 {
	var mu float64
	if m.Triplet {
		mu -= GS
	}
	if m.Ortho {
		mu += 2 * MuPr
	}
	return mu
}
