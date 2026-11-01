// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	_ "embed"
	"encoding/json"
	"fmt"
	"math"
)

//go:embed molecules.json
var moleculesJSON []byte

// Atom of a molecule (molecules.json, written by prep/prep.py).
type Atom struct {
	Name    string     `json:"name"`
	Element string     `json:"element"`
	Mass    float64    `json:"mass"`   // u
	Charge  float64    `json:"charge"` // e
	Type    string     `json:"type"`
	Sigma   float64    `json:"sigma"`   // Å
	Epsilon float64    `json:"epsilon"` // kcal/mol
	XYZ     [3]float64 `json:"xyz"`     // Å
}

// Torsion term k·(1 + cos(n·φ − phase)) [kcal/mol].
type TorsionTerm struct {
	N     int     `json:"n"`
	K     float64 `json:"k"`
	Phase float64 `json:"phase"`
}

// Molecule with its AMBER terms. Bonds and angles are ½k(x − x0)² (OpenMM's
// form, k in kcal/mol/Å² or kcal/mol/rad²).
type Molecule struct {
	Title string `json:"title"`
	Note  string `json:"note"`
	Atoms []Atom `json:"atoms"`
	Bonds []struct {
		I, J  int
		K, R0 float64
	} `json:"bonds"`
	RawAngles []struct {
		I      int     `json:"i"`
		J      int     `json:"j"`
		K      int     `json:"k"`
		Kf     float64 `json:"kf"`
		Theta0 float64 `json:"theta0"`
	} `json:"angles"`
	Propers []struct {
		Atoms [4]int        `json:"atoms"`
		Terms []TorsionTerm `json:"terms"`
	} `json:"propers"`
	Impropers []struct {
		Atoms [4]int        `json:"atoms"`
		Terms []TorsionTerm `json:"terms"`
	} `json:"impropers"`
	Pairs []struct {
		I       int  `json:"i"`
		J       int  `json:"j"`
		Scale14 bool `json:"scale14"`
		Inter   bool `json:"inter"` // between the two bases of one pair
		Stack   bool `json:"stack"` // between bases of different pairs of a stack
	} `json:"pairs"`
	Coulomb14 float64 `json:"coulomb14"`
	LJ14      float64 `json:"lj14"`
}

// Molecules returns the six molecules: "A", "T", "G", "C" (the methylated
// bases) and "AT", "GC" (the Watson–Crick pairs).
func Molecules() (map[string]*Molecule, error) {
	var raw struct {
		Source    string               `json:"source"`
		Molecules map[string]*Molecule `json:"molecules"`
	}
	if err := json.Unmarshal(moleculesJSON, &raw); err != nil {
		return nil, err
	}
	for k, m := range raw.Molecules {
		if len(m.Atoms) == 0 {
			return nil, fmt.Errorf("molecule %s: no atoms", k)
		}
	}
	return raw.Molecules, nil
}

// Coulomb's constant in kcal·Å/(mol·e²): OpenMM's ONE_4PI_EPS0 =
// 138.935456 kJ·nm/(mol·e²), so the energies match the force field as
// parameterized.
const coulombK = 138.935456 / 4.184 * 10

// Parts of the energy, for the circuit's split into springs and couplings.
const (
	PartBond = 1 << iota
	PartAngle
	PartTorsion
	PartIntra // nonbonded within a base
	PartInter // nonbonded between the two bases of a pair
	PartStack // nonbonded between bases of different stacked pairs
	PartAll   = PartBond | PartAngle | PartTorsion | PartIntra | PartInter | PartStack
)

func sub(a, b [3]float64) [3]float64 { return [3]float64{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func dot(a, b [3]float64) float64    { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
func cross(a, b [3]float64) [3]float64 {
	return [3]float64{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}
func scale(a [3]float64, s float64) [3]float64 { return [3]float64{a[0] * s, a[1] * s, a[2] * s} }
func norm(a [3]float64) float64                { return math.Sqrt(dot(a, a)) }

func at(x []float64, i int) [3]float64 { return [3]float64{x[3*i], x[3*i+1], x[3*i+2]} }
func add(g []float64, i int, v [3]float64, s float64) {
	g[3*i] += s * v[0]
	g[3*i+1] += s * v[1]
	g[3*i+2] += s * v[2]
}

// Dihedral returns the torsion angle i–j–k–l and its gradient with
// respect to the four positions (Blondel & Karplus, J. Comput. Chem. 17,
// 1132, 1996, eqs. 27–29).
func Dihedral(p1, p2, p3, p4 [3]float64) (phi float64, g [4][3]float64) {
	f := sub(p1, p2)
	gv := sub(p2, p3)
	h := sub(p4, p3)
	a := cross(f, gv)
	b := cross(h, gv)
	gn := norm(gv)
	aa, bb := dot(a, a), dot(b, b)
	phi = math.Atan2(dot(cross(b, a), gv)/gn, dot(a, b))
	g1 := scale(a, -gn/aa)
	g4 := scale(b, gn/bb)
	fg, hg := dot(f, gv)/(aa*gn), dot(h, gv)/(bb*gn)
	g2 := [3]float64{}
	g3 := [3]float64{}
	for d := 0; d < 3; d++ {
		g2[d] = a[d]*(gn/aa+fg) - b[d]*hg
		g3[d] = b[d]*(hg-gn/bb) - a[d]*fg
	}
	g[0], g[1], g[2], g[3] = g1, g2, g3, g4
	return phi, g
}

// Energy returns the energy [kcal/mol] of the chosen parts at x (Å) and
// adds its gradient into g when g is not nil.
func (m *Molecule) Energy(x []float64, parts int, g []float64) float64 {
	e := 0.0
	if parts&PartBond != 0 {
		for _, b := range m.Bonds {
			d := sub(at(x, b.I), at(x, b.J))
			r := norm(d)
			dr := r - b.R0
			e += 0.5 * b.K * dr * dr
			if g != nil {
				add(g, b.I, d, b.K*dr/r)
				add(g, b.J, d, -b.K*dr/r)
			}
		}
	}
	if parts&PartAngle != 0 {
		for _, a := range m.RawAngles {
			u := sub(at(x, a.I), at(x, a.J))
			v := sub(at(x, a.K), at(x, a.J))
			nu, nv := norm(u), norm(v)
			c := dot(u, v) / (nu * nv)
			c = math.Max(-1, math.Min(1, c))
			th := math.Acos(c)
			dth := th - a.Theta0
			e += 0.5 * a.Kf * dth * dth
			if g != nil {
				s := math.Sqrt(1 - c*c)
				de := a.Kf * dth
				for k := 0; k < 3; k++ {
					gi := (c*u[k]/nu - v[k]/nv) / (nu * s)
					gk := (c*v[k]/nv - u[k]/nu) / (nv * s)
					g[3*a.I+k] += de * gi
					g[3*a.K+k] += de * gk
					g[3*a.J+k] -= de * (gi + gk)
				}
			}
		}
	}
	if parts&PartTorsion != 0 {
		tors := func(idx [4]int, terms []TorsionTerm) {
			phi, dg := Dihedral(at(x, idx[0]), at(x, idx[1]), at(x, idx[2]), at(x, idx[3]))
			de := 0.0
			for _, t := range terms {
				arg := float64(t.N)*phi - t.Phase
				e += t.K * (1 + math.Cos(arg))
				de -= t.K * float64(t.N) * math.Sin(arg)
			}
			if g != nil {
				for q := 0; q < 4; q++ {
					add(g, idx[q], dg[q], de)
				}
			}
		}
		for _, p := range m.Propers {
			tors(p.Atoms, p.Terms)
		}
		for _, p := range m.Impropers {
			tors(p.Atoms, p.Terms)
		}
	}
	if parts&(PartIntra|PartInter|PartStack) != 0 {
		for _, p := range m.Pairs {
			want := PartIntra
			if p.Inter {
				want = PartInter
			} else if p.Stack {
				want = PartStack
			}
			if parts&want == 0 {
				continue
			}
			ai, aj := m.Atoms[p.I], m.Atoms[p.J]
			d := sub(at(x, p.I), at(x, p.J))
			r := norm(d)
			qq := coulombK * ai.Charge * aj.Charge
			sig := 0.5 * (ai.Sigma + aj.Sigma)
			eps := math.Sqrt(ai.Epsilon * aj.Epsilon)
			sc, sl := 1.0, 1.0
			if p.Scale14 {
				sc, sl = m.Coulomb14, m.LJ14
			}
			sr6 := math.Pow(sig/r, 6)
			e += sc*qq/r + sl*4*eps*(sr6*sr6-sr6)
			if g != nil {
				dEdr := -sc*qq/(r*r) + sl*4*eps*(-12*sr6*sr6+6*sr6)/r
				add(g, p.I, d, dEdr/r)
				add(g, p.J, d, -dEdr/r)
			}
		}
	}
	return e
}

// Coords returns the starting coordinates (Å) as one vector.
func (m *Molecule) Coords() []float64 {
	x := make([]float64, 3*len(m.Atoms))
	for i, a := range m.Atoms {
		copy(x[3*i:], a.XYZ[:])
	}
	return x
}

// Hessian returns the second derivatives of the chosen parts at x by
// central differences of the analytic gradient (step h Å), symmetrized.
func (m *Molecule) Hessian(x []float64, parts int, h float64) [][]float64 {
	n := len(x)
	H := make([][]float64, n)
	for i := range H {
		H[i] = make([]float64, n)
	}
	xp := append([]float64(nil), x...)
	gp := make([]float64, n)
	gm := make([]float64, n)
	for j := 0; j < n; j++ {
		for i := range gp {
			gp[i], gm[i] = 0, 0
		}
		xp[j] = x[j] + h
		m.Energy(xp, parts, gp)
		xp[j] = x[j] - h
		m.Energy(xp, parts, gm)
		xp[j] = x[j]
		for i := 0; i < n; i++ {
			H[i][j] = (gp[i] - gm[i]) / (2 * h)
		}
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			s := 0.5 * (H[i][j] + H[j][i])
			H[i][j], H[j][i] = s, s
		}
	}
	return H
}

// Minimize finds a local minimum from x: steepest descent with backtracking
// to get close, then Newton steps on the Hessian with the six rigid-body
// directions removed. It returns the RMS gradient [kcal/mol/Å] at the end.
func (m *Molecule) Minimize(x []float64) (float64, error) {
	n := len(x)
	g := make([]float64, n)
	grad := func(x []float64) float64 {
		for i := range g {
			g[i] = 0
		}
		return m.Energy(x, PartAll, g)
	}
	rms := func() float64 {
		s := 0.0
		for _, v := range g {
			s += v * v
		}
		return math.Sqrt(s / float64(n))
	}
	e := grad(x)
	step := 0.01
	xt := make([]float64, n)
	for it := 0; it < 20000 && rms() > 0.05; it++ {
		gn := 0.0
		for _, v := range g {
			gn += v * v
		}
		gn = math.Sqrt(gn)
		for {
			for i := range x {
				xt[i] = x[i] - step*g[i]/gn
			}
			et := m.Energy(xt, PartAll, nil)
			if et < e {
				copy(x, xt)
				step *= 1.2
				break
			}
			step *= 0.5
			if step < 1e-12 {
				return rms(), fmt.Errorf("%s: descent stalled at rms gradient %.3g", m.Title, rms())
			}
		}
		e = grad(x)
	}
	// Newton refinement.
	for it := 0; it < 50 && rms() > 1e-9; it++ {
		H := m.Hessian(x, PartAll, 1e-4)
		vals, vecs := Eigen(H)
		dx := make([]float64, n)
		for k := range vals {
			if math.Abs(vals[k]) < 1e-6 {
				continue // rigid-body and numerical-zero directions
			}
			c := 0.0
			for i := 0; i < n; i++ {
				c += vecs[i][k] * g[i]
			}
			c /= math.Abs(vals[k]) // |λ|: a saddle direction is descended, not climbed
			for i := 0; i < n; i++ {
				dx[i] -= c * vecs[i][k]
			}
		}
		t := 1.0
		for {
			for i := range x {
				xt[i] = x[i] + t*dx[i]
			}
			et := m.Energy(xt, PartAll, nil)
			if et <= e+1e-12 || t < 1e-6 {
				copy(x, xt)
				break
			}
			t *= 0.5
		}
		e = grad(x)
	}
	return rms(), nil
}

// Eigen returns the eigenvalues (ascending) and eigenvectors (columns) of a
// symmetric matrix by cyclic Jacobi rotations (Golub & Van Loan §8.5).
func Eigen(h [][]float64) ([]float64, [][]float64) {
	n := len(h)
	a := make([][]float64, n)
	v := make([][]float64, n)
	for i := range a {
		a[i] = append([]float64(nil), h[i]...)
		v[i] = make([]float64, n)
		v[i][i] = 1
	}
	scaleA := 0.0
	for i := range a {
		for j := range a {
			scaleA += a[i][j] * a[i][j]
		}
	}
	for sweep := 0; sweep < 100; sweep++ {
		off := 0.0
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off <= 1e-30*scaleA {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if a[p][q] == 0 {
					continue
				}
				th := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := math.Copysign(1, th) / (math.Abs(th) + math.Sqrt(th*th+1))
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p], a[k][q] = c*akp-s*akq, s*akp+c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k], a[q][k] = c*apk-s*aqk, s*apk+c*aqk
				}
				for k := 0; k < n; k++ {
					vkp, vkq := v[k][p], v[k][q]
					v[k][p], v[k][q] = c*vkp-s*vkq, s*vkp+c*vkq
				}
			}
		}
	}
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	for i := 1; i < n; i++ { // insertion sort by eigenvalue
		for j := i; j > 0 && a[idx[j]][idx[j]] < a[idx[j-1]][idx[j-1]]; j-- {
			idx[j], idx[j-1] = idx[j-1], idx[j]
		}
	}
	vals := make([]float64, n)
	vecs := make([][]float64, n)
	for i := range vecs {
		vecs[i] = make([]float64, n)
	}
	for k, c := range idx {
		vals[k] = a[c][c]
		for i := 0; i < n; i++ {
			vecs[i][k] = v[i][c]
		}
	}
	return vals, vecs
}
