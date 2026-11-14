// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command dirac_qed builds the standard model of the electron — Dirac's
// equation plus quantum electrodynamics — as the benchmark that the
// repository's electron hypotheses (H21, H23, H44, H49, H50) are scored
// against. It is not a hypothesis: it is the standard prediction, computed
// here rather than quoted, so its numbers can be checked line by line.
//
//  1. Dirac's equation, free: the spectrum ±√(p² + m²), and spin ½.
//  2. Dirac's equation in a magnetic field: the Landau levels, and g = 2.
//  3. Dirac's equation in hydrogen: the n = 1, 2 levels by shooting on the
//     radial equations, against the exact formula and the measured fine
//     structure; the Lamb shift is what Dirac misses.
//  4. QED: Schwinger's one-loop integral and Petermann's two-loop closed
//     form computed here, the higher loops, hadronic and weak terms from
//     the literature; a_e against the measured value with α from Rb and Cs.
//  5. The measurement in kronspice: the trapped electron's axial motion as
//     a series LC that notches a detector tank; the trap's clocks read as
//     tank resonances, and g/2 recovered by the Brown–Gabrielse theorem.
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// Alpha is the fine-structure constant used for hydrogen (CODATA 2022).
const Alpha = 7.2973525643e-3

// HzPerMC2 is m c²/h [Hz].
const HzPerMC2 = MeC2 / Hpl

func main() {
	viz := flag.String("viz", "", "write CSVs for the pictures to DIR")
	flag.Parse()

	fmt.Println("dirac_qed: the standard model of the electron (Dirac + QED), built as the benchmark")
	fmt.Println("Units in Parts 1–3: ħ = c = m = 1 (energies in m c², lengths in ƛ = ħ/mc).")

	fmt.Println()
	fmt.Println("Part 1. Dirac's equation, free: H = α·p + β, diagonalized (4×4 Hermitian, by Jacobi)")
	fmt.Println("  p                     eigenvalues                                  √(p² + 1)     helicity on E > 0")
	for _, p := range [][3]float64{{0, 0, 0.001}, {0.3, 0.4, 1.2}, {2, -1, 0.5}, {10, 0, 0}} {
		e, h := HelicityOnPositive(p)
		pn := math.Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
		fmt.Printf("  (%5.3g,%5.3g,%5.3g)   %+.12f %+.12f %+.12f %+.12f   %.12f   %+.6f %+.6f\n", p[0], p[1], p[2], e[0], e[1], e[2], e[3], math.Sqrt(pn*pn+1), h[0], h[1])
	}
	fmt.Printf("  S² = %.12f on every state = s(s + 1) with s = ½.\n", SpinSquared())
	fmt.Println("  Two energies ±√(p² + m²), each twice: particle and antiparticle, spin up and down. Spin ½ is not put in; it")
	fmt.Println("  comes from the 4-component equation that relativity forces.")

	fmt.Println()
	fmt.Println("Part 2. Dirac's equation in a field B along z: Landau levels in the oscillator basis (40 levels)")
	fmt.Println("  b = B/B_crit   E(n=0,↓) − 1        spin splitting            cyclotron spacing         g = 2·spin/cyclotron")
	var landRows []string
	for _, b := range []float64{1e-6, 1e-4, 1e-2, 0.1, 1} {
		lv := LandauDirac(b, 40)
		g, sp, cy := GFromLandau(lv)
		fmt.Printf("  %-12.0e   %+.3e         %.15e   %.15e   %.12f\n", b, lv[0].E-1, sp, cy, g)
		landRows = append(landRows, fmt.Sprintf("%g,%.15g,%.15g,%.15g", b, sp, cy, g))
	}
	fmt.Printf("  B_crit = m²c²/(eħ) = %.4g T. In every field the spin flip costs exactly one cyclotron quantum, so g = 2:\n", Me*Me*9e16/(Qe*1.054571817e-34))
	fmt.Println("  the levels are √(1 + 2b·k), k = n + ½ + s_z, and spin and orbit enter only through k (Dirac's g = 2).")

	fmt.Println()
	fmt.Println("Part 3. Hydrogen from Dirac's radial equations (RK4 shooting in ln r, bisection), Z = 1, infinite nuclear mass")
	fmt.Println("  level     κ    E/mc² (shooting)       E/mc² (exact formula)   binding [eV]")
	type lv struct {
		name string
		n, k int
	}
	levels := []lv{{"1S½", 1, -1}, {"2S½", 2, -1}, {"2P½", 2, 1}, {"2P3/2", 2, -2}}
	E := map[string]float64{}
	for _, l := range levels {
		e := HydrogenLevel(l.n, l.k, Alpha)
		E[l.name] = e
		fmt.Printf("  %-8s  %+d   %.15f      %.15f       %.9f\n", l.name, l.k, e, DiracExact(l.n, l.k, Alpha), (1-e)*MeC2/Qe)
	}
	fs := (E["2P3/2"] - E["2P½"]) * HzPerMC2
	mred := 1 / (1 + 1/1836.152673426)
	fmt.Printf("  2P3/2 − 2P½ fine structure: %.4f GHz (infinite mass); × reduced mass %.6f → %.4f GHz\n", fs/1e9, mred, fs*mred/1e9)
	fmt.Println("  with QED the interval is 10.96905 GHz (the value the n = 2 measurements are fitted with): Dirac alone is 0.23 %")
	fmt.Println("  low, the gap chiefly the anomalous moment. Dirac puts 2S½ and 2P½ at the same energy exactly; measured, 2S½")
	fmt.Println("  is 1057.845(9) MHz higher (Lundeen & Pipkin 1981; Lamb & Retherford 1947) — QED's first success.")
	fmt.Printf("  2S½ − 2P½ from Dirac: %.3e Hz (zero).\n", (E["2S½"]-E["2P½"])*HzPerMC2)

	fmt.Println()
	fmt.Println("Part 4. QED's correction to g: a_e = (g − 2)/2 = Σ A⁽²ⁿ⁾(α/π)ⁿ + hadronic + weak")
	fmt.Printf("  one loop, Feynman-parameter integral (computed): F2(0) = %.13f (α/π)  → Schwinger's α/2π\n", SchwingerF2(0, 400))
	fmt.Printf("  … and its fall with momentum transfer: F2(Q² = m²) = %.6f, F2(Q² = 4m²) = %.6f (α/π)\n", SchwingerF2(1, 400), SchwingerF2(4, 400))
	fmt.Printf("  two loops, Petermann's closed form (computed): A⁽⁴⁾ = %.15f   (literature −0.328478965579194)\n", A4Petermann())
	var qedRows []string
	for _, c := range []struct {
		name string
		ai   float64
	}{{"Rb (Morel 2020)", AlphaInvRb}, {"Cs (Parker 2018)", AlphaInvCs}} {
		terms, tot := AeSeries(c.ai)
		fmt.Printf("  with α⁻¹ = %.9f from %s:\n", c.ai, c.name)
		for _, t := range terms {
			fmt.Printf("    %-46s %+.12e\n", t.Name, t.Value)
			if c.ai == AlphaInvRb {
				qedRows = append(qedRows, fmt.Sprintf("%q,%.15g", t.Name, t.Value))
			}
		}
		d := AeMeasured - tot
		fmt.Printf("    a_e (theory) = %.15e   measured %.14e (±%.2g)   measured − theory = %+.2f × 10⁻¹³\n", tot, AeMeasured, AeMeasuredErr, d*1e13)
		fmt.Printf("    g (theory) = %.14f   g/2 − measured g/2 = %+.2e (relative)\n", GFromAe(tot), d/(1+AeMeasured))
	}
	fmt.Println("  Theory and measurement agree to ~3×10⁻¹³ (Rb) or ~1×10⁻¹² (Cs) in a_e — parts in 10¹⁰ of the anomaly,")
	fmt.Println("  parts in 10¹³ of g. The two α measurements disagree with each other by more than that (5.5σ).")

	fmt.Println()
	fmt.Println("Part 5. The measurement in kronspice (Penning trap, B = 5.3 T, ν̄_z = 114 MHz; Fan et al. 2023)")
	tr := DefaultTrap()
	le, ce := tr.ElectronLC()
	ct, rt := tr.TankRC()
	fmt.Printf("  the electron's axial motion as a series LC (Wineland & Dehmelt): L_e = m d²/e² = %.4g H, C_e = %.4g F (d = %.1f mm)\n", le, ce, tr.D*1e3)
	fmt.Printf("  detector tank: L = %.3g H, C = %.4g F, R = %.4g Ω (Q = %.0f); electron coupled through an ideal 1:%g transformer\n", tr.Lt, ct, rt, tr.Qt, TransformerN)
	fc, w, depth, freqs, with, without, err := tr.AxialDip(5*tr.DipWidth(), 8001)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  kronspice AC: dip at ν_z %+.2e Hz off 114 MHz, depth %.2e, width at half depth %.4f Hz × √3 = %.4f Hz\n", fc-tr.Fz, depth, w, w*math.Sqrt(3))
	fmt.Printf("  expected γ_z/2π = R e²/(2π m d²) = %.4f Hz: the trapped electron is the series LC the theorem says.\n", tr.DipWidth())

	_, aeRb := AeSeries(AlphaInvRb)
	f := tr.Frequencies(GFromAe(aeRb))
	fmt.Printf("  the trap's clocks for g = %.14f: ν_c = %.9g Hz, ν̄_c = %.9g Hz, ν̄_m = %.6g Hz, ν̄_a = %.9g Hz\n", GFromAe(aeRb), f.Nuc, f.NucBar, f.NuM, f.NuA)
	rc, err := ReadTank(f.NucBar, 1e6)
	if err != nil {
		log.Fatal(err)
	}
	ra, _ := ReadTank(f.NuA, 1e5)
	rz, _ := ReadTank(f.NuZ, 1e5)
	rm, _ := ReadTank(f.NuM, 1e4)
	g2 := GFromClocks(rc, rz, rm, ra)
	fmt.Printf("  each read as a tank resonance in kronspice AC (relative errors %.1e, %.1e, %.1e, %.1e)\n", rc/f.NucBar-1, ra/f.NuA-1, rz/f.NuZ-1, rm/f.NuM-1)
	fmt.Printf("  g/2 = (ν̄_c + ν̄_a)/√(ν̄_c² + ν̄_z² + ν̄_m²) = %.15f  (put in %.15f): the chain reads g/2 to %.0e.\n", g2, 1+aeRb, math.Abs(g2-1-aeRb))
	fmt.Println("  The magnetron frequency comes out 43.8 kHz against Fan et al.'s measured 43 kHz.")

	fmt.Println()
	fmt.Println("Benchmark scorecard (for the electron ledger): spin ½ derived (Part 1); g = 2 derived (Part 2); hydrogen fine")
	fmt.Println("structure 0.23 % from Dirac alone, Lamb shift missed (Part 3); a_e within 3×10⁻¹³ (Part 4); charge, mass and α are inputs.")

	if *viz != "" {
		write := func(name, head string, rows []string) {
			fp, err := os.Create(filepath.Join(*viz, name))
			if err != nil {
				log.Fatal(err)
			}
			fmt.Fprintln(fp, head)
			for _, r := range rows {
				fmt.Fprintln(fp, r)
			}
			fp.Close()
		}
		write("landau_g.csv", "b,spin,cyclotron,g", landRows)
		write("qed_terms.csv", "term,value", qedRows)
		var dip []string
		for i := range freqs {
			dip = append(dip, fmt.Sprintf("%.6f,%.6g,%.6g", freqs[i]-tr.Fz, with[i], without[i]))
		}
		write("axial_dip.csv", "df_hz,v_with,v_without", dip)
		lv := LandauDirac(0.05, 40)
		var lrows []string
		for _, l := range lv[:11] {
			lrows = append(lrows, fmt.Sprintf("%.12g,%d,%v", l.E, l.N, l.Up))
		}
		write("landau_levels.csv", "E,n,up", lrows)
		var hrows []string
		for _, l := range levels {
			hrows = append(hrows, fmt.Sprintf("%s,%.15g", l.name, (1-E[l.name])*MeC2/Qe))
		}
		write("hydrogen.csv", "level,binding_eV", hrows)
	}
}
