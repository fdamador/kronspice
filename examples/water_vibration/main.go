// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command water_vibration models the chemical bonds of a water molecule as
// an electric circuit and finds its vibration spectrum with an AC sweep.
//
// Mechanical–electrical analogy (the mobility or "force–current"
// analogy; F. A. Firestone, J. Acoust. Soc. Am. 4, 249, 1933; formalized
// for networks by Branin 1966, as used by Kron for his network analogies):
//
//	velocity of a coordinate   ↔ node voltage
//	force                      ↔ current
//	mass m                     ↔ capacitor C = m to ground
//	spring of stiffness k      ↔ inductor L = 1/k
//	damper                     ↔ resistor to ground
//
// Each of the six in-plane Cartesian coordinates of the three atoms (O,
// H1, H2) is a node with its atom's mass as a capacitor. Each bond
// coordinate (the two O–H stretches and the H–O–H bend) is a spring: the
// rate at which it changes is a linear combination of the atomic
// velocities, s′ = Σ B_kj v_j (E. B. Wilson's B matrix; Wilson, Decius &
// Cross, Molecular Vibrations, McGraw-Hill 1955, Ch. 4), built by a chain
// of voltage-controlled voltage sources; the spring force flows in an
// inductor; and that force acts back on every atom, B_kj·f_k, through
// current-controlled current sources. The E/F pair is an ideal multiport
// transformer — Kron's connection tensor B between the bond frame and the
// atom frame. The force constants form the inductance matrix L = F⁻¹, so
// the stretch–stretch interaction f_rr is a mutual inductance (K element).
//
// Units in the circuit: capacitance in daltons (u), inductance in Å/mdyn,
// so the circuit's angular frequency ω_c = 1/√(LC) maps to the molecule's
// ω = ω_c·√(100 N/m / 1 u) (1 mdyn/Å = 100 N/m), reported in cm⁻¹.
//
// Test of the model: the three force constants are fitted to H₂O alone.
// The same force constants then predict, with no refit, the spectra of
// D₂O, HDO and H₂¹⁸O (isotopes change masses, not bonds — Born–Oppenheimer),
// against measured fundamentals. A harmonic force field fitted to
// fundamentals leaves an error from anharmonicity, larger for H than for D,
// which the run reports.
//
//	go run ./examples/water_vibration [-netlists dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// Atomic masses [u] (CODATA/AME 2016 via NIST "Atomic Weights and Isotopic
// Compositions").
const (
	MassH   = 1.00782503223
	MassD   = 2.01410177812
	MassO16 = 15.99491461957
	MassO18 = 17.99915961286
)

// Equilibrium geometry of water: r_e = 0.958 Å, ∠HOH = 104.4776°
// (Hoy & Bunker 1979, via NIST CCCBDB).
const (
	Re    = 0.958
	Theta = 104.4776 * math.Pi / 180
)

// wavenumber converts the circuit's ω_c (u, Å/mdyn units) to cm⁻¹.
var wavenumber = math.Sqrt(100/1.66053906660e-27) / (2 * math.Pi * 2.99792458e10)

// Molecule is one isotopologue: masses of O, H1, H2.
type Molecule struct {
	Name         string
	MO, MH1, MH2 float64
	// Measured gas-phase fundamentals [cm⁻¹], in the order the modes are
	// reported (ascending here).
	Measured []float64
	Source   string
}

// Isotopologues with measured fundamentals (NIST CCCBDB, experimental
// data: Shimanouchi, Tables of Molecular Vibrational Frequencies; Huber &
// Herzberg 1979). H₂¹⁸O is not in CCCBDB and its band origins have not
// yet been checked against a primary source here, so it is reported as a
// prediction only (Measured nil).
var (
	H2O   = Molecule{"H2O", MassO16, MassH, MassH, []float64{1595, 3657, 3756}, "NIST CCCBDB (Shimanouchi; Huber & Herzberg)"}
	D2O   = Molecule{"D2O", MassO16, MassD, MassD, []float64{1178, 2671, 2788}, "NIST CCCBDB (Shimanouchi)"}
	HDO   = Molecule{"HDO", MassO16, MassH, MassD, []float64{1402, 2727, 3707}, "NIST CCCBDB (Shimanouchi; Huber & Herzberg)"}
	H2O18 = Molecule{"H2-18O", MassO18, MassH, MassH, nil, "no verified measurement yet: prediction only"}
)

// HarmonicH2O are the harmonic frequencies of H₂O (NIST CCCBDB).
var HarmonicH2O = []float64{1649, 3832, 3943}

// ForceField is the valence force field [mdyn/Å]: O–H stretch f_r, the
// stretch–stretch interaction f_rr, and the bend f_a (for the coordinate
// r_e·Δθ, so in mdyn/Å). The stretch–bend interaction is left out, so
// three constants match three H₂O frequencies.
type ForceField struct{ Fr, Frr, Fa float64 }

// geometry returns Cartesian positions (O, H1, H2) in the molecular plane.
func geometry() [3][2]float64 {
	s, c := math.Sin(Theta/2), math.Cos(Theta/2)
	return [3][2]float64{{0, 0}, {Re * s, Re * c}, {-Re * s, Re * c}}
}

// BMatrix is Wilson's B (3 internal × 6 Cartesian: Ox, Oy, H1x, H1y, H2x,
// H2y) for r1, r2 and r_e·θ (Wilson, Decius & Cross §4-1, eqs. for stretch
// and bend).
func BMatrix() [3][6]float64 {
	g := geometry()
	var e [2][2]float64 // unit vectors O→H1, O→H2
	for k := 0; k < 2; k++ {
		dx, dy := g[k+1][0]-g[0][0], g[k+1][1]-g[0][1]
		r := math.Hypot(dx, dy)
		e[k] = [2]float64{dx / r, dy / r}
	}
	var b [3][6]float64
	for k := 0; k < 2; k++ {
		for d := 0; d < 2; d++ {
			b[k][2+2*k+d] = e[k][d]
			b[k][d] = -e[k][d]
		}
	}
	cs, sn := math.Cos(Theta), math.Sin(Theta)
	for d := 0; d < 2; d++ {
		t1 := (cs*e[0][d] - e[1][d]) / (Re * sn)
		t2 := (cs*e[1][d] - e[0][d]) / (Re * sn)
		b[2][2+d] = Re * t1
		b[2][4+d] = Re * t2
		b[2][d] = -Re * (t1 + t2)
	}
	return b
}

// fmat returns the 3×3 force-constant matrix.
func (f ForceField) fmat() [3][3]float64 {
	return [3][3]float64{{f.Fr, f.Frr, 0}, {f.Frr, f.Fr, 0}, {0, 0, f.Fa}}
}

// Modes returns the vibration wavenumbers [cm⁻¹], ascending, from the
// GF eigenproblem: the mass-weighted Cartesian Hessian M^{-1/2} Bᵀ F B
// M^{-1/2} has three zero eigenvalues (two translations, one rotation in
// the plane) and three vibrations.
func Modes(m Molecule, f ForceField) []float64 {
	b := BMatrix()
	fm := f.fmat()
	mass := [6]float64{m.MO, m.MO, m.MH1, m.MH1, m.MH2, m.MH2}
	var h [6][6]float64
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			var s float64
			for k := 0; k < 3; k++ {
				for l := 0; l < 3; l++ {
					s += b[k][i] * fm[k][l] * b[l][j]
				}
			}
			h[i][j] = s / math.Sqrt(mass[i]*mass[j])
		}
	}
	ev := jacobi(h)
	sort.Float64s(ev)
	out := make([]float64, 0, 3)
	for _, l := range ev[3:] {
		out = append(out, math.Sqrt(math.Max(l, 0))*wavenumber)
	}
	return out
}

// jacobi returns the eigenvalues of a symmetric 6×6 matrix (cyclic Jacobi
// rotations; Golub & Van Loan §8.5).
func jacobi(a [6][6]float64) []float64 {
	for sweep := 0; sweep < 100; sweep++ {
		off := 0.0
		for i := 0; i < 6; i++ {
			for j := i + 1; j < 6; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < 6; p++ {
			for q := p + 1; q < 6; q++ {
				if math.Abs(a[p][q]) < 1e-300 {
					continue
				}
				th := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := math.Copysign(1, th) / (math.Abs(th) + math.Sqrt(th*th+1))
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < 6; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p], a[k][q] = c*akp-s*akq, s*akp+c*akq
				}
				for k := 0; k < 6; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k], a[q][k] = c*apk-s*aqk, s*apk+c*aqk
				}
			}
		}
	}
	ev := make([]float64, 6)
	for i := range ev {
		ev[i] = a[i][i]
	}
	return ev
}

// Fit finds the force field that reproduces target wavenumbers (ascending:
// bend, symmetric stretch, antisymmetric stretch) for molecule m, by
// Newton's method on the squared frequencies (which are linear in F).
func Fit(m Molecule, target []float64) (ForceField, error) {
	x := []float64{8, -0.1, 0.7}
	for it := 0; it < 50; it++ {
		ff := ForceField{x[0], x[1], x[2]}
		r := Modes(m, ff)
		res := make([]float64, 3)
		maxr := 0.0
		for i := range res {
			res[i] = r[i]*r[i] - target[i]*target[i]
			maxr = math.Max(maxr, math.Abs(res[i])/(target[i]*target[i]))
		}
		if maxr < 1e-13 {
			return ff, nil
		}
		var jac [3][3]float64
		for j := 0; j < 3; j++ {
			d := 1e-6 * math.Max(1, math.Abs(x[j]))
			xp := append([]float64(nil), x...)
			xp[j] += d
			rp := Modes(m, ForceField{xp[0], xp[1], xp[2]})
			for i := 0; i < 3; i++ {
				jac[i][j] = (rp[i]*rp[i] - r[i]*r[i]) / d
			}
		}
		dx, ok := solve3(jac, res)
		if !ok {
			return ForceField{}, fmt.Errorf("force-field fit: singular Jacobian")
		}
		for j := range x {
			x[j] -= dx[j]
		}
	}
	return ForceField{}, fmt.Errorf("force-field fit did not converge")
}

func solve3(a [3][3]float64, b []float64) ([]float64, bool) {
	det := func(m [3][3]float64) float64 {
		return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
	}
	d := det(a)
	if math.Abs(d) < 1e-300 {
		return nil, false
	}
	x := make([]float64, 3)
	for k := 0; k < 3; k++ {
		m := a
		for i := 0; i < 3; i++ {
			m[i][k] = b[i]
		}
		x[k] = det(m) / d
	}
	return x, true
}

// Netlist writes the molecule as a SPICE netlist for kronspice (and any SPICE
// simulator with E, F, K elements): six velocity nodes with mass
// capacitors, a damper (R) per node so the spectrum lines have a finite
// width and the operating point has a DC path, three spring inductors fed
// by VCVS chains (B·v), current sensors, CCCS back-coupling (Bᵀ·f), the
// stretch–stretch mutual inductance, and an AC force on H1 (both
// directions, so every mode is excited). Frequencies are in circuit units:
// multiply by Wavenumber() for cm⁻¹.
func Netlist(m Molecule, f ForceField, rDamp float64) string {
	b := BMatrix()
	node := []string{"ox", "oy", "h1x", "h1y", "h2x", "h2y"}
	mass := []float64{m.MO, m.MO, m.MH1, m.MH1, m.MH2, m.MH2}
	var s strings.Builder
	fmt.Fprintf(&s, "%s bonds as a circuit\n", m.Name)
	fmt.Fprintf(&s, "* %s vibration as a circuit (examples/water_vibration). Units: C in u, L in Angstrom/mdyn;\n", m.Name)
	fmt.Fprintf(&s, "* circuit frequency f_c [Hz] x %.6f = wavenumber [1/cm]. Force field (mdyn/A): fr %.6f frr %.6f fa %.6f\n",
		2*math.Pi*wavenumber, f.Fr, f.Frr, f.Fa)
	for i, n := range node {
		fmt.Fprintf(&s, "CM_%s %s 0 %.12g\nRD_%s %s 0 %.12g\n", n, n, mass[i], n, n, rDamp)
	}
	// Compliance matrix L = F⁻¹ for the stretch block, 1/fa for the bend.
	det := f.Fr*f.Fr - f.Frr*f.Frr
	l1 := f.Fr / det
	mut := -f.Frr / det
	ls := []float64{l1, l1, 1 / f.Fa}
	coord := []string{"r1", "r2", "bend"}
	for k, c := range coord {
		prev := "0"
		t := 0
		for j, n := range node {
			if math.Abs(b[k][j]) < 1e-15 {
				continue
			}
			next := fmt.Sprintf("%s_s%d", c, t)
			fmt.Fprintf(&s, "E_%s_%s %s %s %s 0 %.15g\n", c, n, next, prev, n, b[k][j])
			prev = next
			t++
		}
		fmt.Fprintf(&s, "VS_%s %s %s_l 0\nL_%s %s_l 0 %.15g\n", c, prev, c, c, c, ls[k])
		for j, n := range node {
			if math.Abs(b[k][j]) < 1e-15 {
				continue
			}
			fmt.Fprintf(&s, "F_%s_%s %s 0 VS_%s %.15g\n", c, n, n, c, b[k][j])
		}
	}
	fmt.Fprintf(&s, "K_rr L_r1 L_r2 %.15g\n", mut/l1)
	s.WriteString("I_fx 0 h1x AC 1\nI_fy 0 h1y AC 1\n")
	// The IR window 500–4500 cm⁻¹ in circuit hertz (Peaks sweeps the same).
	k := 2 * math.Pi * wavenumber
	fmt.Fprintf(&s, ".ac lin 4001 %.9g %.9g\n.end\n", 500/k, 4500/k)
	return s.String()
}

// Spectrum runs kronspice's AC analysis of the netlist from fLo to fHi
// [cm⁻¹] with n points and returns the frequencies and the summed kinetic
// response Σ m|v|² at each.
func Spectrum(src string, m Molecule, fLo, fHi float64, n int) ([]float64, []float64, error) {
	d, err := netlist.ParseString(src)
	if err != nil {
		return nil, nil, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil, nil, err
	}
	k := 2 * math.Pi * wavenumber // circuit Hz → cm⁻¹
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: n, FStart: fLo / k, FStop: fHi / k}, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	mass := map[string]float64{"ox": m.MO, "oy": m.MO, "h1x": m.MH1, "h1y": m.MH1, "h2x": m.MH2, "h2y": m.MH2}
	resp := make([]float64, len(r.Freq))
	for name, mm := range mass {
		v, ok := r.NodeVoltage(name)
		if !ok {
			return nil, nil, fmt.Errorf("node %s missing", name)
		}
		for i := range v {
			a := cmplx.Abs(v[i])
			resp[i] += mm * a * a
		}
	}
	freq := make([]float64, len(r.Freq))
	for i, fc := range r.Freq {
		freq[i] = fc * k
	}
	return freq, resp, nil
}

// Peaks finds each mode's line in the circuit's AC response: a coarse
// sweep to locate the local maxima, then a fine sweep around each and a
// parabola through the top three points.
func Peaks(src string, m Molecule, lo, hi float64) ([]float64, error) {
	f, r, err := Spectrum(src, m, lo, hi, 4001)
	if err != nil {
		return nil, err
	}
	var out []float64
	for i := 1; i+1 < len(r); i++ {
		if r[i] > r[i-1] && r[i] >= r[i+1] {
			w := f[i+1] - f[i]
			ff, rr, err := Spectrum(src, m, f[i]-w, f[i]+w, 401)
			if err != nil {
				return nil, err
			}
			j := 1
			for k := 1; k+1 < len(rr); k++ {
				if rr[k] > rr[j] {
					j = k
				}
			}
			y0, y1, y2 := math.Log(rr[j-1]), math.Log(rr[j]), math.Log(rr[j+1])
			den := y0 - 2*y1 + y2
			off := 0.0
			if den != 0 {
				off = 0.5 * (y0 - y2) / den
			}
			out = append(out, ff[j]+off*(ff[1]-ff[0]))
		}
	}
	return out, nil
}

// Wavenumber is the conversion from the netlist's circuit frequency [Hz]
// to cm⁻¹.
func Wavenumber() float64 { return 2 * math.Pi * wavenumber }

func main() {
	log.SetFlags(0)
	dir := flag.String("netlists", "examples/water_vibration/netlists", "directory for the generated .cir files")
	flag.Parse()
	fmt.Println("Water bonds as a circuit: atoms are mass capacitors, bonds are spring inductors, Wilson's B as E/F transformers")
	fmt.Printf("geometry r_e %.3f Å, angle %.4f° (Hoy & Bunker 1979)\n\n", Re, Theta*180/math.Pi)
	fund, err := Fit(H2O, H2O.Measured)
	if err != nil {
		log.Fatal(err)
	}
	harm, err := Fit(H2O, HarmonicH2O)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("force field fitted to H2O fundamentals: f_r %.4f, f_rr %.4f, f_a %.4f mdyn/Å\n", fund.Fr, fund.Frr, fund.Fa)
	fmt.Printf("force field fitted to H2O harmonic:     f_r %.4f, f_rr %.4f, f_a %.4f mdyn/Å\n\n", harm.Fr, harm.Frr, harm.Fa)
	if err := os.MkdirAll(*dir, 0o755); err != nil {
		log.Fatal(err)
	}
	const rDamp = 1e6 // damper conductance 1e-6: line width ~1e-6 of the frequency
	fmt.Println("molecule   mode        circuit AC peak   GF eigen     measured   error")
	for _, m := range []Molecule{H2O, D2O, HDO, H2O18} {
		src := Netlist(m, fund, rDamp)
		path := filepath.Join(*dir, strings.ToLower(strings.ReplaceAll(m.Name, "-", ""))+".cir")
		if err := os.WriteFile(path, []byte(src), 0o644); err != nil {
			log.Fatal(err)
		}
		pk, err := Peaks(src, m, 500, 4500)
		if err != nil {
			log.Fatal(err)
		}
		ev := Modes(m, fund)
		names := []string{"bend", "stretch", "stretch"}
		for i := range ev {
			p := math.NaN()
			if i < len(pk) {
				p = pk[i]
			}
			if m.Measured == nil {
				fmt.Printf("%-9s  %-10s  %10.2f cm⁻¹  %10.2f        —        —\n", m.Name, names[i], p, ev[i])
				continue
			}
			fmt.Printf("%-9s  %-10s  %10.2f cm⁻¹  %10.2f   %8.1f   %+5.2f %%\n", m.Name, names[i], p, ev[i], m.Measured[i], 100*(ev[i]/m.Measured[i]-1))
		}
		fmt.Printf("           netlist: %s   (measured: %s)\n", path, m.Source)
	}
	fmt.Println("\nH2O is fitted (0 % by construction); D2O, HDO and H2-18O are predictions with the same bonds.")
	fmt.Println("The remaining error is anharmonicity. A measured fundamental sits below the harmonic frequency by")
	fmt.Println("a shift that shrinks with heavier atoms (ω·x ∝ 1/μ against ω ∝ 1/√μ). A harmonic force field fitted")
	fmt.Println("to H2O's fundamentals absorbs H's larger shift, so it comes out soft and under-predicts deuterium")
	fmt.Println("modes by about 1 %; the O–H stretch of HDO, closest to H2O, is within 0.1 %.")
}
