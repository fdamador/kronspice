// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kron_schrodinger solves the Schrödinger equation of a chemical
// bond as an electric network, after G. Kron, "Electric circuit models of
// the Schrödinger equation," Phys. Rev. 67, 39–43 (1945): series
// inductors carry the kinetic energy, shunt inductors to ground the
// potential energy, and shunt capacitors the energy eigenvalue; the wave
// function ψ is the junction voltage to ground, and an eigenstate is a
// resonance of the network.
//
// The molecule is H₂⁺, the simplest chemical bond (two protons, one
// electron), and the H atom as its check; both have exact solutions. The
// electron's ground state is axially symmetric (m = 0), so the network is
// laid out on the (ρ, z) half-plane, ρ ≥ 0, with the nuclei on the axis.
// In atomic units the equation, multiplied by ρ so the network is
// symmetric (reciprocal), is
//
//	−½ [∂ρ(ρ ∂ρψ) + ρ ∂z²ψ] + ρ V ψ = E ρ ψ,    V = −Σ Z_a / |r − R_a|
//
// On a grid of spacing h with points at ρ_i = (i + ½)h (so no point sits on
// the axis, and the axis needs no boundary condition: the branch to ρ = 0
// has zero weight) this is A ψ = E B ψ with B = diag(ρ_i). In Kron's
// network: series branches of admittance ρ_face/(2h²) between neighbours,
// a shunt branch ρ_i·V_i at each junction, and a shunt capacitor ρ_i. The
// outer boundary is a short circuit (ψ = 0) far from the nuclei, Kron's
// termination "by a short circuit at a point of zero voltage".
//
// Eigenvalues by Kron's driven network: the network is driven at a trial
// energy σ by sources B·x and its response is the next x (inverse
// iteration; the network's factors are computed once and reused); the
// response converges to the resonance nearest σ, and the Rayleigh quotient
// gives E. As a circuit, the same network is written as a SPICE netlist
// of inductors and capacitors (its dual form, potential shifted by V0 so
// every element is positive, ω² = 1/(E + V0)), and kronspice's AC analysis
// finds the ground state as its highest resonance.
//
// Reference values (exact, non-relativistic, fixed nuclei): H atom
// E = −0.5 Ha; H₂⁺ at R = 2 a₀, electronic energy −1.1026342145 Ha
// (D. R. Bates, K. Ledsham, A. L. Stewart, Phil. Trans. R. Soc. A 246,
// 215, 1953; J. M. Peek, J. Chem. Phys. 43, 3004, 1965); equilibrium
// R_e = 1.997 a₀, dissociation energy D_e = 0.1026 Ha = 2.79 eV.
//
// Water: ten electrons. A one-electron network is exact for H₂⁺ only;
// for H₂O the electrons repel each other, which needs a self-consistent
// field (each electron's network sees the potential of the others, from
// a Poisson network, iterated: Hartree–Fock) on a 3-D grid, since water
// is not axially symmetric. That is the next step and is not done here.
//
//	go run ./examples/kron_schrodinger
package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/solver"
)

// Unit conversions (CODATA 2018).
const (
	Hartree = 27.211386245988 // eV
	Bohr    = 0.529177210903  // Å
)

// Nucleus is a point charge Z on the axis at z.
type Nucleus struct{ Z, At float64 }

// Grid is the (ρ, z) network: nr × nz junctions, spacing h, first junction
// at ρ = h/2 and z = Z0 + h/2.
type Grid struct {
	H      float64
	Nr, Nz int
	Z0     float64
	Nuclei []Nucleus
}

// NewGrid lays a grid of spacing h over ρ ≤ rmax and the z range that
// extends margin beyond the outer nuclei, with every nucleus midway between
// two junctions (its z must be a multiple of h from the first nucleus).
func NewGrid(h, rmax, margin float64, nuc []Nucleus) Grid {
	lo, hi := nuc[0].At, nuc[0].At
	for _, n := range nuc {
		lo, hi = math.Min(lo, n.At), math.Max(hi, n.At)
	}
	k := math.Ceil(margin / h)
	z0 := lo - k*h
	nz := int(math.Round((hi-lo)/h + 2*k))
	return Grid{H: h, Nr: int(math.Ceil(rmax / h)), Nz: nz, Z0: z0, Nuclei: nuc}
}

func (g Grid) n() int            { return g.Nr * g.Nz }
func (g Grid) idx(i, j int) int  { return j*g.Nr + i }
func (g Grid) rho(i int) float64 { return (float64(i) + 0.5) * g.H }
func (g Grid) zAt(j int) float64 { return g.Z0 + (float64(j)+0.5)*g.H }
func (g Grid) potential(i, j int) float64 {
	v := 0.0
	r, z := g.rho(i), g.zAt(j)
	for _, n := range g.Nuclei {
		v -= n.Z / math.Hypot(r, z-n.At)
	}
	return v
}

// Network returns the branches of Kron's network: for each junction its
// shunt potential term ρV and capacitor ρ, and the series admittances to
// its +ρ and +z neighbours (−1: none, the boundary short).
type branch struct {
	a, b int
	y    float64
}

func (g Grid) network() (shunt, cap []float64, series []branch) {
	n := g.n()
	shunt = make([]float64, n)
	cap = make([]float64, n)
	h2 := g.H * g.H
	for j := 0; j < g.Nz; j++ {
		for i := 0; i < g.Nr; i++ {
			k := g.idx(i, j)
			r := g.rho(i)
			cap[k] = r
			shunt[k] = r * g.potential(i, j)
			// Radial branches: weight ρ at the face; the face to ρ = 0
			// has weight 0, and the outer face ends in the short (ψ = 0).
			yOut := (r + g.H/2) / (2 * h2)
			if i+1 < g.Nr {
				series = append(series, branch{k, g.idx(i+1, j), yOut})
			} else {
				series = append(series, branch{k, -1, yOut})
			}
			yz := r / (2 * h2)
			if j+1 < g.Nz {
				series = append(series, branch{k, g.idx(i, j+1), yz})
			} else {
				series = append(series, branch{k, -1, yz})
			}
			if j == 0 {
				series = append(series, branch{k, -1, yz})
			}
		}
	}
	return
}

// Ground finds the lowest state by the driven network at trial energy
// sigma (below the ground state). It returns E [Ha] and ψ.
func (g Grid) Ground(sigma float64) (float64, []float64, error) {
	shunt, cap, series := g.network()
	n := g.n()
	m, err := solver.NewMatrix(n, solver.SparseLU)
	if err != nil {
		return 0, nil, err
	}
	// The network driven at energy σ: A − σB.
	stamp := func(add func(i, j int, v float64)) {
		for k := 0; k < n; k++ {
			add(k, k, shunt[k])
		}
		for _, br := range series {
			add(br.a, br.a, br.y)
			if br.b >= 0 {
				add(br.b, br.b, br.y)
				add(br.a, br.b, -br.y)
				add(br.b, br.a, -br.y)
			}
		}
	}
	stamp(func(i, j int, v float64) { m.Add(i, j, v) })
	for k := 0; k < n; k++ {
		m.Add(k, k, -sigma*cap[k])
	}
	ax := func(x, y []float64) {
		for k := range y {
			y[k] = shunt[k] * x[k]
		}
		for _, br := range series {
			if br.b >= 0 {
				d := br.y * (x[br.a] - x[br.b])
				y[br.a] += d
				y[br.b] -= d
			} else {
				y[br.a] += br.y * x[br.a]
			}
		}
	}
	x := make([]float64, n)
	for j := 0; j < g.Nz; j++ {
		for i := 0; i < g.Nr; i++ {
			// Start from the sum of 1s orbitals on the nuclei.
			s := 0.0
			for _, nu := range g.Nuclei {
				s += math.Exp(-math.Hypot(g.rho(i), g.zAt(j)-nu.At))
			}
			x[g.idx(i, j)] = s
		}
	}
	rhs := make([]float64, n)
	ay := make([]float64, n)
	e, prev := 0.0, math.Inf(1)
	for it := 0; it < 200; it++ {
		for k := range rhs {
			rhs[k] = cap[k] * x[k]
		}
		if err := m.Solve(rhs, x); err != nil {
			return 0, nil, err
		}
		var num, den float64
		ax(x, ay)
		for k := range x {
			num += x[k] * ay[k]
			den += x[k] * cap[k] * x[k]
		}
		e = num / den
		norm := math.Sqrt(den)
		for k := range x {
			x[k] /= norm
		}
		if math.Abs(e-prev) < 1e-13*math.Abs(e) {
			return e, x, nil
		}
		prev = e
	}
	return e, x, fmt.Errorf("driven network did not settle at σ = %g (last E %g)", sigma, e)
}

// Extrapolated returns the h → 0 limit from runs at h, h/2, h/4 assuming
// E(h) = E0 + a·h^p, with p estimated from the three runs (Richardson).
func Extrapolated(e1, e2, e3 float64) (float64, float64) {
	p := math.Log2((e1 - e2) / (e2 - e3))
	return e3 + (e3-e2)/(math.Pow(2, p)-1), p
}

// Netlist writes the network as an LC circuit for kronspice. It is the dual
// of Kron's network (capacitors and inductors exchanged), because the AC
// analysis starts from a DC operating point and Kron's all-inductor
// network has loops of shorts at DC: series capacitors C = y between
// junctions (and to ground at the outer boundary), shunt capacitors
// C = ρ(V + V0) and shunt inductors L = 1/ρ at each junction. Then the
// network's admittance times jω is ω²·(A + V0·B) − B, which vanishes at
// ω² = 1/(E + V0): every eigenstate is a resonance, the ground state the
// highest of them. V0 shifts the potential so every capacitor is
// positive. One junction near the first nucleus is driven by an AC
// current source.
func (g Grid) Netlist(v0 float64) string {
	shunt, cap, series := g.network()
	var s strings.Builder
	fmt.Fprintf(&s, "Kron 1945 Schrodinger network (dual form), %d junctions\n", g.n())
	fmt.Fprintf(&s, "* psi = junction voltage; omega^2 = 1/(E + V0), V0 = %.12g Ha (examples/kron_schrodinger)\n", v0)
	for k := range cap {
		fmt.Fprintf(&s, "CV%d n%d 0 %.15g\nLB%d n%d 0 %.15g\n", k, k, shunt[k]+v0*cap[k], k, k, 1/cap[k])
	}
	for q, br := range series {
		b := "0"
		if br.b >= 0 {
			b = fmt.Sprintf("n%d", br.b)
		}
		fmt.Fprintf(&s, "CS%d n%d %s %.15g\n", q, br.a, b, br.y)
	}
	fmt.Fprintf(&s, "I1 0 n%d AC 1\n", g.idx(0, g.Nz/2))
	// Sweep over bound-state energies E from −1.5 to −0.3 Ha (ω² = 1/(E+V0)).
	f := func(e float64) float64 { return 1 / math.Sqrt(e+v0) / (2 * math.Pi) }
	fmt.Fprintf(&s, ".ac lin 241 %.9g %.9g\n.end\n", f(-0.3), f(-1.5))
	return s.String()
}

// CircuitGround finds the lowest resonance of the netlist with kronspice's
// AC analysis (coarse scan, then a narrowing search for the peak of the
// driving-point response) and returns E = 1/ω² − V0.
func (g Grid) CircuitGround(v0, eLo, eHi float64) (float64, error) {
	d, err := netlist.ParseString(g.Netlist(v0))
	if err != nil {
		return 0, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return 0, err
	}
	node := fmt.Sprintf("n%d", g.idx(0, g.Nz/2))
	resp := func(e float64) (float64, error) {
		f := 1 / math.Sqrt(e+v0) / (2 * math.Pi)
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		v, _ := r.NodeVoltage(node)
		return cmplx.Abs(v[0]), nil
	}
	// The driving-point impedance of a lossless network has a pole at each
	// resonance: scan, then golden-section search on the largest response.
	best, bestE := -1.0, eLo
	const scan = 60
	for k := 0; k <= scan; k++ {
		e := eLo + (eHi-eLo)*float64(k)/scan
		r, err := resp(e)
		if err != nil {
			return 0, err
		}
		if r > best {
			best, bestE = r, e
		}
	}
	w := (eHi - eLo) / scan
	a, b := bestE-w, bestE+w
	gr := (math.Sqrt(5) - 1) / 2
	for it := 0; it < 60 && b-a > 1e-11; it++ {
		c1, c2 := b-gr*(b-a), a+gr*(b-a)
		r1, err := resp(c1)
		if err != nil {
			return 0, err
		}
		r2, err := resp(c2)
		if err != nil {
			return 0, err
		}
		if r1 > r2 {
			b = c2
		} else {
			a = c1
		}
	}
	return (a + b) / 2, nil
}

// CubicMin fits E(R) = c0 + c1 R + c2 R² + c3 R³ by least squares and
// returns the position and value of its minimum inside the data range.
func CubicMin(x, y []float64) (float64, float64) {
	var a [4][5]float64
	for i := range x {
		p := [4]float64{1, x[i], x[i] * x[i], x[i] * x[i] * x[i]}
		for r := 0; r < 4; r++ {
			for c := 0; c < 4; c++ {
				a[r][c] += p[r] * p[c]
			}
			a[r][4] += p[r] * y[i]
		}
	}
	for k := 0; k < 4; k++ { // Gauss–Jordan with partial pivoting
		piv := k
		for r := k + 1; r < 4; r++ {
			if math.Abs(a[r][k]) > math.Abs(a[piv][k]) {
				piv = r
			}
		}
		a[k], a[piv] = a[piv], a[k]
		for r := 0; r < 4; r++ {
			if r != k {
				f := a[r][k] / a[k][k]
				for c := k; c < 5; c++ {
					a[r][c] -= f * a[k][c]
				}
			}
		}
	}
	c := [4]float64{a[0][4] / a[0][0], a[1][4] / a[1][1], a[2][4] / a[2][2], a[3][4] / a[3][3]}
	f := func(r float64) float64 { return c[0] + r*(c[1]+r*(c[2]+r*c[3])) }
	// E′(R) = c1 + 2c2 R + 3c3 R² = 0, the root with E″ > 0.
	disc := 4*c[2]*c[2] - 12*c[3]*c[1]
	best, bestE := x[0], f(x[0])
	if disc >= 0 {
		for _, sgn := range []float64{-1, 1} {
			r := (-2*c[2] + sgn*math.Sqrt(disc)) / (6 * c[3])
			if r >= x[0] && r <= x[len(x)-1] && 2*c[2]+6*c[3]*r > 0 && f(r) < bestE {
				best, bestE = r, f(r)
			}
		}
	}
	return best, bestE
}

// H2Plus returns the grid for H₂⁺ at bond length r [a₀].
func H2Plus(h, r float64) Grid {
	return NewGrid(h, 12, 12, []Nucleus{{1, -r / 2}, {1, r / 2}})
}

// Hydrogen returns the grid for the H atom.
func Hydrogen(h float64) Grid { return NewGrid(h, 14, 14, []Nucleus{{1, 0}}) }

func main() {
	log.SetFlags(0)
	fmt.Println("Kron's Schrödinger network (Phys. Rev. 67, 39, 1945): the electron of a chemical bond as a circuit")
	fmt.Println()
	hs := []float64{0.2, 0.1, 0.05}
	fmt.Println("H atom (check):")
	var eh []float64
	for _, h := range hs {
		g := Hydrogen(h)
		e, _, err := g.Ground(-0.7)
		if err != nil {
			log.Fatal(err)
		}
		eh = append(eh, e)
		fmt.Printf("  h = %.3f a₀  (%6d junctions): E = %.8f Ha\n", h, g.n(), e)
	}
	e0, p := Extrapolated(eh[0], eh[1], eh[2])
	fmt.Printf("  h → 0 (order %.2f): E = %.6f Ha = %.4f eV   (exact −0.5 Ha = −13.6057 eV)\n\n", p, e0, e0*Hartree)

	fmt.Println("H₂⁺ bond (two protons, one electron), electronic energy at R = 2 a₀:")
	var e2 []float64
	for _, h := range hs {
		g := H2Plus(h, 2)
		e, _, err := g.Ground(-1.4)
		if err != nil {
			log.Fatal(err)
		}
		e2 = append(e2, e)
		fmt.Printf("  h = %.3f a₀  (%6d junctions): E = %.8f Ha\n", h, g.n(), e)
	}
	ee, p2 := Extrapolated(e2[0], e2[1], e2[2])
	fmt.Printf("  h → 0 (order %.2f): E = %.6f Ha   (exact −1.1026342 Ha)\n\n", p2, ee)

	fmt.Println("Bond curve E(R) + 1/R (the protons' repulsion), each R extrapolated h → 0 from h = 0.2, 0.1:")
	type pt struct{ r, e float64 }
	var curve []pt
	// Each R from h = 0.2 and 0.1, extrapolated with the order p2 found
	// at R = 2 from three grids (the order is set by the Coulomb cusps,
	// the same at every R).
	// R stays a multiple of the coarse spacing so both protons sit midway
	// between junctions on every grid (otherwise the discretization error
	// jumps from one R to the next).
	for _, r := range []float64{1.6, 1.8, 2.0, 2.2, 2.4, 2.8} {
		var es []float64
		for _, h := range hs[:2] {
			e, _, err := H2Plus(h, r).Ground(-1.5)
			if err != nil {
				log.Fatal(err)
			}
			es = append(es, e)
		}
		ex := es[1] + (es[1]-es[0])/(math.Pow(2, p2)-1)
		tot := ex + 1/r
		curve = append(curve, pt{r, tot})
		fmt.Printf("  R = %.1f a₀ (%.3f Å): total %.6f Ha\n", r, r*Bohr, tot)
	}
	// Least-squares cubic through the points from 1.6 to 2.4 a₀, and its
	// minimum.
	var rs, ts []float64
	for _, p := range curve {
		if p.r <= 2.4 {
			rs, ts = append(rs, p.r), append(ts, p.e)
		}
	}
	re, emin := CubicMin(rs, ts)
	de := -0.5 - emin
	fmt.Printf("  bond length R_e = %.3f a₀ = %.3f Å   (exact 1.997 a₀ = 1.057 Å)\n", re, re*Bohr)
	fmt.Printf("  bond energy D_e = %.4f Ha = %.3f eV   (exact 0.1026 Ha = 2.79 eV)\n\n", de, de*Hartree)

	fmt.Println("The same network as an LC circuit in kronspice (coarse grid, h = 0.4 a₀), ground state by AC analysis:")
	g := H2Plus(0.4, 2)
	eg, _, err := g.Ground(-1.4)
	if err != nil {
		log.Fatal(err)
	}
	v0 := 0.0
	for i := 0; i < g.Nr; i++ {
		for j := 0; j < g.Nz; j++ {
			v0 = math.Max(v0, -g.potential(i, j))
		}
	}
	v0 += 1
	ec, err := g.CircuitGround(v0, eg-0.3, eg+0.3)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  %d junctions: circuit resonance E = %.8f Ha, driven-network eigenvalue %.8f Ha (same network, coarse grid)\n",
		g.n(), ec, eg)
	if err := os.MkdirAll("examples/kron_schrodinger/netlists", 0o755); err == nil {
		_ = os.WriteFile("examples/kron_schrodinger/netlists/h2plus_coarse.cir", []byte(g.Netlist(v0)), 0o644)
		fmt.Println("  netlist: examples/kron_schrodinger/netlists/h2plus_coarse.cir")
	}
	fmt.Println("\nWater has ten interacting electrons: the next step is a self-consistent field (Hartree–Fock) on a 3-D network.")
}
