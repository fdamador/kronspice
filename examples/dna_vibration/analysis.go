// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/em"
)

// Physical constants for the light coupling (CODATA 2018).
const (
	hPlanck   = 6.62607015e-34
	kBoltz    = 1.380649e-23
	eCharge   = 1.602176634e-19
	TempK     = 300.0 // room temperature for the thermal release spectrum
	FWHMcm    = 6.0   // display linewidth of every line [cm⁻¹]
	SweepFrom = 0.0
	SweepTo   = 3700.0
	SweepPts  = 1851 // 2 cm⁻¹ steps
)

// Release is how a mode gives its energy back as infrared light.
type Release struct {
	A        float64 // spontaneous emission rate of v=1→0 [1/s] (from the IR intensity)
	ALarmor  float64 // the same from Maxwell: Larmor power of the classical dipole / its energy
	Lifetime float64 // radiative lifetime 1/A [s]
	Thermal  float64 // power radiated at TempK in thermal equilibrium [W per molecule]
}

// ReleaseOf computes a mode's release by both routes.
func ReleaseOf(md Mode) Release {
	nu := md.Wavenumber
	A := EinsteinA(nu, md.Intensity)
	// Maxwell route: a classical oscillation with one quantum of energy E = hcν̃
	// has amplitude Q0 = √(2E)/ω in the mass-weighted coordinate; its dipole
	// amplitude is p0 = (dμ/dQ)·Q0, it radiates Larmor's power P, and P/E is
	// the rate at which the mode loses energy to light.
	omega := 2 * math.Pi * cCm * nu
	E := hPlanck * cCm * nu
	q0 := math.Sqrt(2*E) / omega // √kg·m
	conv := eCharge / math.Sqrt(amu)
	var p [3]complex128
	for d := 0; d < 3; d++ {
		p[d] = complex(md.DipoleDer[d]*conv*q0, 0)
	}
	aL := em.LarmorPower(p, omega) / E
	x := hPlanck * cCm * nu / (kBoltz * TempK)
	nbar := 1 / math.Expm1(x)
	return Release{A: A, ALarmor: aL, Lifetime: 1 / A, Thermal: A * E * nbar}
}

// Label names a mode by its potential-energy distribution over the
// circuit's springs (Wilson, Decius & Cross §5-5): the share of the mode's
// potential energy in each bond or angle spring, in the out-of-plane and
// nonbonded residual, and in the base-pairing couplings.
func (s *State) Label(md Mode, springs []Spring, res []Entry) (string, []float64) {
	type part struct {
		name string
		e    float64
	}
	var parts []part
	total := 0.0
	for i := range md.Vector {
		for j := range md.Vector {
			total += md.Vector[i] * s.K[i][j] * md.Vector[j]
		}
	}
	for _, sp := range springs {
		q := 0.0
		for i, b := range sp.B {
			q += b * md.Vector[i]
		}
		parts = append(parts, part{sp.Name, sp.K * q * q})
	}
	rIntra, rInter, rStack := 0.0, 0.0, 0.0
	for _, e := range res {
		v := md.Vector[e.I] * e.V * md.Vector[e.J]
		switch e.Group {
		case "HB":
			rInter += v
		case "ST":
			rStack += v
		default:
			rIntra += v
		}
	}
	parts = append(parts, part{"torsion/out-of-plane/nonbonded", rIntra}, part{"base pairing (H-bonds)", rInter}, part{"stacking", rStack})
	sort.Slice(parts, func(i, j int) bool { return math.Abs(parts[i].e) > math.Abs(parts[j].e) })
	var names []string
	for _, p := range parts[:3] {
		sh := p.e / total
		if sh < 0.15 {
			break
		}
		names = append(names, fmt.Sprintf("%s %.0f%%", pretty(p.name), 100*sh))
	}
	if len(names) == 0 {
		names = append(names, "delocalized")
	}
	return strings.Join(names, ", "), nil
}

// pretty turns a spring name into chemistry: b_a_n6_a_h61 → "A N6–H61
// stretch", a_t_c2_t_n3_t_c4 → "T C2–N3–C4 bend".
func pretty(n string) string {
	f := strings.Split(n, "_")
	if len(f) < 3 || (f[0] != "b" && f[0] != "a") {
		return n
	}
	var atoms []string
	base := ""
	for i := 1; i+1 < len(f); i += 2 {
		if base == "" {
			base = strings.ToUpper(f[i])
		}
		atoms = append(atoms, strings.ToUpper(f[i+1]))
	}
	if f[0] == "b" {
		return base + " " + strings.Join(atoms, "–") + " stretch"
	}
	return base + " " + strings.Join(atoms, "–") + " bend"
}

// HBond is a hydrogen bond between the two bases of a pair.
type HBond struct {
	Donor, H, Acceptor string
	HA, DA             float64 // H···A and D···A distances [Å]
	Angle              float64 // D–H···A [deg]
}

// HBonds finds the inter-base hydrogen bonds: H on N or O within 2.5 Å of
// an N or O of the other base.
func (s *State) HBonds() []HBond {
	m := s.Mol
	var out []HBond
	base := func(i int) string { return strings.SplitN(m.Atoms[i].Name, ":", 2)[0] }
	for h := range m.Atoms {
		if m.Atoms[h].Element != "H" {
			continue
		}
		donor := -1
		for _, b := range m.Bonds {
			if b.I == h {
				donor = b.J
			} else if b.J == h {
				donor = b.I
			}
		}
		if donor < 0 || (m.Atoms[donor].Element != "N" && m.Atoms[donor].Element != "O") {
			continue
		}
		for a := range m.Atoms {
			if base(a) == base(h) || (m.Atoms[a].Element != "N" && m.Atoms[a].Element != "O") {
				continue
			}
			ha := norm(sub(at(s.X, h), at(s.X, a)))
			if ha > 2.5 {
				continue
			}
			u := sub(at(s.X, donor), at(s.X, h))
			v := sub(at(s.X, a), at(s.X, h))
			ang := math.Acos(dot(u, v)/(norm(u)*norm(v))) * 180 / math.Pi
			out = append(out, HBond{m.Atoms[donor].Name, m.Atoms[h].Name, m.Atoms[a].Name, ha, norm(sub(at(s.X, donor), at(s.X, a))), ang})
		}
	}
	return out
}

// Size is the largest distance between two atoms [Å].
func (s *State) Size() float64 {
	d := 0.0
	for i := range s.Mol.Atoms {
		for j := range s.Mol.Atoms {
			d = math.Max(d, norm(sub(at(s.X, i), at(s.X, j))))
		}
	}
	return d
}

// Emission turns the circuit's absorption curve into the release spectra:
// one quantum in every mode radiates in proportion to A·hν ∝ ν³·(absorbed
// power line), and a molecule at TempK radiates that times the thermal
// occupation 1/(e^{hcν/kT} − 1) (Kirchhoff's law). Both are normalized to
// a peak of 1.
func Emission(nu, abs []float64) (quantum, thermal []float64) {
	quantum = make([]float64, len(nu))
	thermal = make([]float64, len(nu))
	for k, f := range nu {
		if f <= 0 {
			continue
		}
		quantum[k] = abs[k] * f * f * f
		thermal[k] = quantum[k] / math.Expm1(hPlanck*cCm*f/(kBoltz*TempK))
	}
	normalize(quantum)
	normalize(thermal)
	return
}

func normalize(v []float64) {
	m := 0.0
	for _, x := range v {
		m = math.Max(m, x)
	}
	if m > 0 {
		for i := range v {
			v[i] /= m
		}
	}
}

// Band is a named group vibration and the model's frequency for it: the
// mode with the largest potential-energy share in the named bond
// stretches (for an NH₂ group, the two modes with the largest share in its
// two N–H stretches; the higher is antisymmetric).
type Band struct {
	Molecule, Name string
	Wavenumber     float64
	Share          float64
}

// share returns, per mode, the potential-energy fraction in the springs.
func (s *State) share(names ...string) []float64 {
	sp := map[string]Spring{}
	for _, x := range s.Springs() {
		sp[x.Name] = x
	}
	out := make([]float64, len(s.Modes))
	for k, md := range s.Modes {
		tot := 0.0
		for i := range md.Vector {
			for j := range md.Vector {
				tot += md.Vector[i] * s.K[i][j] * md.Vector[j]
			}
		}
		for _, n := range names {
			x, ok := sp[n]
			if !ok {
				continue
			}
			q := 0.0
			for i, b := range x.B {
				q += b * md.Vector[i]
			}
			out[k] += x.K * q * q / tot
		}
	}
	return out
}

// Stretch is the mode with the largest share in one bond.
func (s *State) Stretch(bond string) (float64, float64) {
	sh := s.share(bond)
	best := 0
	for k := range sh {
		if sh[k] > sh[best] {
			best = k
		}
	}
	return s.Modes[best].Wavenumber, sh[best]
}

// NH2 returns the antisymmetric and symmetric stretch of an amino group.
func (s *State) NH2(b1, b2 string) (asym, sym, share float64) {
	sh := s.share(b1, b2)
	i1, i2 := -1, -1
	for k := range sh {
		if i1 < 0 || sh[k] > sh[i1] {
			i1, i2 = k, i1
		} else if i2 < 0 || sh[k] > sh[i2] {
			i2 = k
		}
	}
	f1, f2 := s.Modes[i1].Wavenumber, s.Modes[i2].Wavenumber
	return math.Max(f1, f2), math.Min(f1, f2), math.Min(sh[i1], sh[i2])
}

// Bands lists the group frequencies compared with measurement.
func Bands(st map[string]*State) []Band {
	var out []Band
	add := func(mol, name string, f, sh float64) { out = append(out, Band{mol, name, f, sh}) }
	nh2 := func(mol, pre, n, h1, h2, label string) {
		a, s, sh := st[mol].NH2("b_"+pre+"_"+n+"_"+pre+"_"+h1, "b_"+pre+"_"+n+"_"+pre+"_"+h2)
		add(mol, label+" NH₂ antisymmetric stretch", a, sh)
		add(mol, label+" NH₂ symmetric stretch", s, sh)
	}
	one := func(mol, bond, label string) {
		f, sh := st[mol].Stretch(bond)
		add(mol, label, f, sh)
	}
	for _, mol := range []string{"A", "AT"} {
		nh2(mol, "a", "n6", "h61", "h62", "A")
	}
	for _, mol := range []string{"T", "AT"} {
		one(mol, "b_t_n3_t_h3", "T N3–H stretch")
		one(mol, "b_t_c4_t_o4", "T C4=O stretch")
		one(mol, "b_t_c2_t_o2", "T C2=O stretch")
	}
	for _, mol := range []string{"G", "GC"} {
		nh2(mol, "g", "n2", "h21", "h22", "G")
		one(mol, "b_g_n1_g_h1", "G N1–H stretch")
		one(mol, "b_g_c6_g_o6", "G C6=O stretch")
	}
	for _, mol := range []string{"C", "GC"} {
		nh2(mol, "c", "n4", "h41", "h42", "C")
		one(mol, "b_c_c2_c_o2", "C C2=O stretch")
	}
	return out
}

// PairMotion names how a pair mode moves the hydrogen bonds: the change of
// each donor···acceptor distance per unit of the mode's displacement.
// All lengthening together is the H-bond stretch; opposite signs is a
// shear (opening) in the plane; no change is an out-of-plane motion
// (buckle or propeller twist).
func (s *State) PairMotion(md Mode) string {
	hb := s.HBonds()
	if len(hb) == 0 {
		return ""
	}
	idx := map[string]int{}
	for i, a := range s.Mol.Atoms {
		idx[a.Name] = i
	}
	tot := 0.0
	for _, v := range md.Vector {
		tot += v * v
	}
	tot = math.Sqrt(tot)
	var d []float64
	pos, neg, big := 0, 0, 0.0
	for _, h := range hb {
		a, b := idx[h.Donor], idx[h.Acceptor]
		e := sub(at(s.X, b), at(s.X, a))
		e = scale(e, 1/norm(e))
		dv := sub(at(md.Vector, b), at(md.Vector, a))
		x := dot(dv, e) / tot
		d = append(d, x)
		if x > 0.05 {
			pos++
		} else if x < -0.05 {
			neg++
		}
		big = math.Max(big, math.Abs(x))
	}
	switch {
	case big < 0.05:
		return "out of plane (buckle / propeller)"
	case pos == len(d) || neg == len(d):
		return "H-bond stretch (all bonds together)"
	default:
		return "in-plane shear / opening"
	}
}
