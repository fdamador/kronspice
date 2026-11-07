// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"strings"
)

// Circuit form of a molecule (the mobility analogy of examples/water_vibration):
//
//	velocity of an atom coordinate  ↔ node voltage   (node a_n9_x …)
//	mass                            ↔ capacitor C = m [u] to ground
//	bond or angle spring            ↔ inductor L = 1/k [Å/mdyn]
//	force                           ↔ current
//
// Every AMBER bond and angle is a spring branch: a chain of VCVS forms the
// rate of change of its internal coordinate, q′ = Σ B_j v_j (Wilson's B
// row), an inductor carries the spring force, and CCCS hand the force back
// to each atom, B_j·f (Kron's connection tensor between the bond frame and
// the atom frame). Everything else in the force field — torsions,
// out-of-plane impropers, the electrostatic and van der Waals forces inside
// a base, the pre-stress (first-derivative) part of bonds and angles at the
// minimum, and, in a pair, every force between the two bases (the hydrogen
// bonds and the rest of their electrostatics) — is the residual stiffness
// R = K − Σ k·BᵀB. It acts through displacement integrators (a VCCS feeds
// each node voltage into a 1 u capacitor, so that capacitor's voltage is
// the coordinate's displacement; a 1e14 Ω leak fixes the DC operating
// point and changes the response by less than 1e−9 above 1 cm⁻¹) and one
// VCCS per non-zero entry. The inter-base block of R
// is written as its own group, so the base pairing is a separate,
// named part of the circuit (G_HB_…). The circuit's stiffness is K exactly, up to
// entries below 1e−9 of the largest (finite-difference noise), which are
// dropped and counted.
//
// Frequencies: circuit frequency f_c [Hz] × 2π·(wavenumber) = cm⁻¹, as in
// the water example.

// Spring is one bond or angle spring branch.
type Spring struct {
	Name string
	B    map[int]float64 // coordinate → B coefficient
	K    float64         // mdyn/Å (bond) or mdyn·Å/rad² (angle)
}

// Springs returns the bond and angle springs of s at its minimum.
func (s *State) Springs() []Spring {
	m := s.Mol
	var out []Spring
	for _, b := range m.Bonds {
		d := sub(at(s.X, b.I), at(s.X, b.J))
		r := norm(d)
		B := map[int]float64{}
		for k := 0; k < 3; k++ {
			B[3*b.I+k] += d[k] / r
			B[3*b.J+k] -= d[k] / r
		}
		out = append(out, Spring{Name: "b_" + short(m.Atoms[b.I].Name) + "_" + short(m.Atoms[b.J].Name), B: B, K: b.K * kcalMolA2ToMdynA})
	}
	for _, a := range m.RawAngles {
		u := sub(at(s.X, a.I), at(s.X, a.J))
		v := sub(at(s.X, a.K), at(s.X, a.J))
		nu, nv := norm(u), norm(v)
		c := dot(u, v) / (nu * nv)
		sn := math.Sqrt(1 - c*c)
		B := map[int]float64{}
		for k := 0; k < 3; k++ {
			gi := (c*u[k]/nu - v[k]/nv) / (nu * sn)
			gk := (c*v[k]/nv - u[k]/nu) / (nv * sn)
			B[3*a.I+k] += gi
			B[3*a.K+k] += gk
			B[3*a.J+k] -= gi + gk
		}
		out = append(out, Spring{Name: "a_" + short(m.Atoms[a.I].Name) + "_" + short(m.Atoms[a.J].Name) + "_" + short(m.Atoms[a.K].Name), B: B, K: a.Kf * kcalMolA2ToMdynA})
	}
	return out
}

// Residual returns R = K − Σ k·BᵀB split into its intra-base part, the
// pairing part (between the two bases of a pair) and the stacking part
// (between bases of different pairs; zero for one base or one pair).
func (s *State) Residual(springs []Spring) (intra, inter, stack [][]float64) {
	n := len(s.K)
	intra = make([][]float64, n)
	inter = make([][]float64, n)
	stack = make([][]float64, n)
	for i := range intra {
		intra[i] = make([]float64, n)
		inter[i] = append([]float64(nil), s.KPart[PartInter][i]...)
		stack[i] = append([]float64(nil), s.KPart[PartStack][i]...)
		for j := range intra[i] {
			intra[i][j] = s.K[i][j] - inter[i][j] - stack[i][j]
		}
	}
	for _, sp := range springs {
		for i, bi := range sp.B {
			for j, bj := range sp.B {
				intra[i][j] -= sp.K * bi * bj
			}
		}
	}
	return intra, inter, stack
}

// short turns "A:N9" into "a_n9" (a SPICE-safe name).
func short(name string) string {
	r := strings.NewReplacer(":", "_", "'", "p", "*", "s")
	return strings.ToLower(r.Replace(name))
}

func (s *State) node(i int) string { return short(s.Mol.Atoms[i/3].Name) + "_" + string("xyz"[i%3]) }

// NetlistOptions control the circuit.
type NetlistOptions struct {
	// FWHM is the display linewidth [cm⁻¹] given to every line by a damping
	// resistor proportional to mass on each node (mass-proportional damping
	// keeps the normal modes uncoupled). It is not a computed linewidth.
	FWHM float64
	// Field selects the drive: "x", "y" or "z" puts a current q_a (the
	// atom's charge, a unit field along that axis) into every atom's node —
	// the force of an infrared light field on the fixed AMBER charges; ""
	// drives nothing.
	Field string
	// Analysis is the control card(s) appended before .end.
	Analysis string
}

// Netlist writes the molecule's circuit.
func (s *State) Netlist(o NetlistOptions) (string, Stats) {
	m := s.Mol
	n := len(s.K)
	var st Stats
	var b strings.Builder
	fmt.Fprintf(&b, "%s vibration circuit\n", m.Title)
	fmt.Fprintf(&b, "* examples/rna_vibration. Units: C in u, L in Angstrom/mdyn; f_c [Hz] x %.6f = wavenumber [1/cm].\n", 2*math.Pi*wavenumber)
	fmt.Fprintf(&b, "* Force field: %s\n", m.Note)
	gamma := 2 * math.Pi * o.FWHM / (2 * math.Pi * wavenumber) // damping rate in circuit rad/s
	for i := 0; i < n; i++ {
		mass := m.Atoms[i/3].Mass
		fmt.Fprintf(&b, "CM_%s %s 0 %.12g\n", s.node(i), s.node(i), mass)
		if o.FWHM > 0 {
			fmt.Fprintf(&b, "RD_%s %s 0 %.12g\n", s.node(i), s.node(i), 1/(mass*gamma))
		}
	}
	b.WriteString("* bond and angle springs\n")
	springs := s.Springs()
	for _, sp := range springs {
		prev := "0"
		for t, i := range sortedKeys(sp.B) {
			next := fmt.Sprintf("%s_s%d", sp.Name, t)
			fmt.Fprintf(&b, "E_%s_%d %s %s %s 0 %.15g\n", sp.Name, t, next, prev, s.node(i), sp.B[i])
			prev = next
		}
		// RL: the bond and angle coordinates are redundant (ring closures), so at
		// DC their inductor loops would leave a current undetermined; SpringR
		// in series fixes the operating point. It changes a spring's force by
		// R/(ωL) ≲ 1e−6 at the lowest modes; the direct twin includes it.
		fmt.Fprintf(&b, "VS_%s %s %s_l 0\nL_%s %s_l %s_r %.15g\nRL_%s %s_r 0 %g\n", sp.Name, prev, sp.Name, sp.Name, sp.Name, sp.Name, 1/sp.K, sp.Name, sp.Name, SpringR)
		for t, i := range sortedKeys(sp.B) {
			fmt.Fprintf(&b, "F_%s_%d %s 0 VS_%s %.15g\n", sp.Name, t, s.node(i), sp.Name, sp.B[i])
		}
		st.Springs++
	}
	res, rst := s.ResidualEntries(springs)
	st.Intra, st.Inter, st.Stack, st.Dropped, st.DroppedMax = rst.Intra, rst.Inter, rst.Stack, rst.Dropped, rst.DroppedMax
	b.WriteString("* displacement integrators: V(x_<node>) = displacement\n")
	for i := 0; i < n; i++ {
		fmt.Fprintf(&b, "G_X_%s 0 x_%s %s 0 1\nCX_%s x_%s 0 1\nRX_%s x_%s 0 %g\n", s.node(i), s.node(i), s.node(i), s.node(i), s.node(i), s.node(i), s.node(i), LeakR)
	}
	b.WriteString("* residual stiffness inside each base (G_R): torsions, impropers, nonbonded, pre-stress;\n")
	b.WriteString("* base pairing (G_HB): every force between the two bases (hydrogen bonds and the rest of their electrostatics)\n")
	b.WriteString("* stacking (G_ST): every force between bases of different stacked pairs\n")
	for _, e := range res {
		fmt.Fprintf(&b, "G_%s_%d_%d %s 0 x_%s 0 %.15g\n", e.Group, e.I, e.J, s.node(e.I), s.node(e.J), e.V)
	}
	if o.Field != "" {
		d := strings.Index("xyz", o.Field)
		fmt.Fprintf(&b, "* infrared drive: unit field along %s, force q_a on each atom\n", o.Field)
		for a := range m.Atoms {
			if q := m.Atoms[a].Charge; q != 0 {
				fmt.Fprintf(&b, "I_E_%s 0 %s AC %.12g\n", short(m.Atoms[a].Name), s.node(3*a+d), q)
			}
		}
	}
	if o.Analysis != "" {
		b.WriteString(o.Analysis + "\n")
	}
	b.WriteString(".end\n")
	return b.String(), st
}

// SpringR and LeakR fix the DC operating point (see Netlist); both are in
// the direct twin (DirectResponse) too.
const (
	SpringR = 1e-8
	LeakR   = 1e8
)

// Entry is one coupling of the residual stiffness.
type Entry struct {
	Group string // "R" inside a base, "HB" between the bases of a pair, "ST" between stacked pairs
	I, J  int
	V     float64 // mdyn/Å
}

// ResidualEntries returns the non-zero residual couplings the circuit
// carries, dropping entries below 1e−9 of the largest stiffness
// (finite-difference noise of the Hessian).
func (s *State) ResidualEntries(springs []Spring) ([]Entry, Stats) {
	var st Stats
	intra, inter, stack := s.Residual(springs)
	kmax := 0.0
	for i := range s.K {
		for j := range s.K[i] {
			kmax = math.Max(kmax, math.Abs(s.K[i][j]))
		}
	}
	tol := 1e-9 * kmax
	var out []Entry
	for _, g := range []struct {
		name string
		r    [][]float64
	}{{"R", intra}, {"HB", inter}, {"ST", stack}} {
		for i := range g.r {
			for j := range g.r[i] {
				v := g.r[i][j]
				if v == 0 {
					continue
				}
				if math.Abs(v) < tol {
					st.Dropped++
					st.DroppedMax = math.Max(st.DroppedMax, math.Abs(v))
					continue
				}
				out = append(out, Entry{g.name, i, j, v})
				switch g.name {
				case "HB":
					st.Inter++
				case "ST":
					st.Stack++
				default:
					st.Intra++
				}
			}
		}
	}
	return out, st
}

// Stats counts the circuit's parts.
type Stats struct {
	Springs, Intra, Inter, Stack, Dropped int
	DroppedMax                            float64
}

func sortedKeys(m map[int]float64) []int {
	out := make([]int, 0, len(m))
	for k := range m {
		out = append(out, k)
	}
	for i := 1; i < len(out); i++ {
		for j := i; j > 0 && out[j] < out[j-1]; j-- {
			out[j], out[j-1] = out[j-1], out[j]
		}
	}
	return out
}
