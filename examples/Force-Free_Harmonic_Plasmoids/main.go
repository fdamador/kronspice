// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command Force-Free_Harmonic_Plasmoids studies J. Nachamkin, "Force-Free
// Time-Harmonic Plasmoids" (PL-TR-92-3044, Phillips Laboratory, 1992) as
// written: his field and Appendix B recomputed, the harmonics his eq. 35
// allows for a plasmoid of each size, a size and energy map over the
// reported range of ball lightning, and kronspice analysis of the
// plasmoid's spectrum, ring-down and charge-density condition.
//
//	go run ./examples/Force-Free_Harmonic_Plasmoids
//	go run ./examples/Force-Free_Harmonic_Plasmoids -viz examples/Force-Free_Harmonic_Plasmoids/viz
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

func main() {
	viz := flag.String("viz", "", "also write the CSVs behind the pictures to this directory")
	flag.Parse()
	if *viz != "" {
		if err := writeViz(*viz); err != nil {
			fmt.Fprintln(os.Stderr, "Force-Free_Harmonic_Plasmoids:", err)
			os.Exit(1)
		}
	}
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "Force-Free_Harmonic_Plasmoids:", err)
		os.Exit(1)
	}
}

func lt(m Mode) string {
	if !m.Found {
		return "not found"
	}
	return fmt.Sprintf("%.3g s", m.Lifetime)
}

func run() error {
	fmt.Println("== 1. The paper's field (eq. 8, l = 1) and Appendix B ==")
	curl, div := CurlCheck(44.93)
	fmt.Printf("∇×F = αF to %.1e, ∇·F to %.1e; ⟨|F|²⟩/max|F|² over the sphere = %.6f\n", curl, div, meanSquare())
	al := J1Zero / 0.1
	fmt.Printf("∫|M|²dV numerically = %.4f/α³ (paper B3: 35.87/α³)\n", IM(al)*al*al*al)
	cases := []AppB{
		AppendixB("as printed (B4 uses 35.87/α²), k = 2686", 35.87/(al*al), 2686),
		AppendixB("as printed, nearest zero", 35.87/(al*al), 0),
		AppendixB("B3 as printed (35.87/α³), nearest zero", IM(al), 0),
	}
	fmt.Println(" case                                      | E0 (V/m)  | ω (B5)     | k0 (1/m) | kR used (zero n)  | f (GHz) | ρ_e (C/m³) | electrons  | v/c  | peak E/air bd")
	for _, c := range cases {
		fmt.Printf(" %-41s | %.4g | %.4g | %8.1f | %8.3f (%4d)   | %7.2f | %10.2f | %.4g | %4.2f | %.0f\n",
			c.Label, c.E0, c.Omega, c.K0, c.K*c.R, nearestIndex(c.K*c.R), c.K*C0/(2*math.Pi)/1e9, c.RhoE, c.Electrons, c.Beta, c.PeakOverBreakdown)
	}
	e := cases[0]
	fmt.Println("paper: E0 3.566e9, ω 7.920e11, k 2642 → 2686, ρ_e −32.63, 8.532e17 electrons")
	fmt.Printf("eq. 49 = cold-plasma dispersion ω² = ωp² + α²c² to %.0e; ionized fraction %.2e (paper: one to two molecules in 10⁵)\n", e.DispersionResidue, e.IonFrac)

	fmt.Println("\n   his scaling at fixed R/λ and ξ (eq. 86: energy ∝ R), from the printed example:")
	for _, r := range []float64{1e-3, 50e-6} {
		x := e.K * e.R
		ne := Eps0 * Me * (x*x - J1Zero*J1Zero) * C0 * C0 / (r * r) / (Qe * Qe)
		vol := 4.0 / 3 * math.Pi * r * r * r
		fmt.Printf("   R = %-7s W = %8.4g J, λ = %.3g µm, electrons %.3g, density %.2f × Loschmidt\n",
			fmtLen(r), 1e6*r/0.1, 2*math.Pi*r/x*1e6, ne*vol, ne/Loschmidt)
	}
	fmt.Println("   paper: 1 mm → about 15,000 J; 50 µm → 500 J, λ = 1.17 µm, 4.3×10¹³ electrons (as printed), ≈ 40 × sea-level density")

	fmt.Println("\n== 2. Spatial harmonics: minimum radius αR for order l (first zero of j_l) ==")
	for l := 1; l <= 4; l++ {
		fmt.Printf("l = %d: αR = %.4f\n", l, FirstZero(l))
	}

	R := 0.12 // average reported ball-lightning diameter 24 cm
	fmt.Printf("\n== 3. Temporal harmonics allowed by eq. 35 for R = %.0f cm (average reported diameter 24 cm) ==\n", R*100)
	fmt.Printf("critical frequency αc/2π = %.4f GHz; allowed kR = zeros of j1 above αR\n", J1Zero*C0/R/(2*math.Pi)/1e9)
	fmt.Println(" n | kR      | f (GHz)  | k/α   | electrons /m³ | ionized fraction | TE lifetime | TM lifetime")
	z := J1Zeros(40)
	for i := 1; i < 8 && i < len(z); i++ {
		d := DesignAt(R, 1e5, z[i])
		fmt.Printf("%2d | %7.3f | %8.4f | %5.2f | %13.4g | %16.3g | %11s | %s\n", i+1, z[i], d.F/1e9, z[i]/J1Zero, d.Ne, d.IonFrac, sec(d.TauTE), sec(d.TauTM))
	}

	fmt.Println("\n== 4. Size and energy map over the reported range (diameter 2–100 cm, energy 10 kJ–1 MJ) ==")
	fmt.Println("(a) the paper's recipe, C = 0.01 (Appendix B):")
	fmt.Println(" R (cm) | W (J)  | kR (n)        | f (GHz)  | v/c    | ionized   | E0/air bd | TE life   | TM life")
	sizes := []float64{0.01, 0.05, 0.10, 0.12, 0.20, 0.25, 0.50}
	energies := []float64{1e4, 1e6}
	for _, r := range sizes {
		for _, w := range energies {
			d := DesignFromC(r, w, 0.01)
			fmt.Printf(" %6.0f | %6.0e | %7.2f (%4d) | %8.3f | %6.3f | %9.3g | %9.0f | %9s | %s\n", r*100, w, d.X, d.N, d.F/1e9, d.Beta, d.IonFrac, d.EOverBreakdown, sec(d.TauTE), sec(d.TauTM))
		}
	}
	fmt.Println("(b) the lowest harmonic that keeps the electrons at v ≤ 0.1 c (the paper's non-relativistic assumption):")
	fmt.Println(" R (cm) | W (J)  | kR (n)            | f (THz)  | C        | ionized   | E0/air bd | TE life   | TM life")
	var allMin, allMax float64 = math.Inf(1), 0
	for _, r := range sizes {
		for _, w := range energies {
			d := DesignFromBeta(r, w, 0.1)
			fmt.Printf(" %6.0f | %6.0e | %9.1f (%6d) | %8.3f | %8.2e | %9.3g | %9.0f | %9s | %s\n", r*100, w, d.X, d.N, d.F/1e12, d.C, d.IonFrac, d.EOverBreakdown, sec(d.TauTE), sec(d.TauTM))
			for _, t := range []float64{d.TauTE, d.TauTM} {
				allMin = math.Min(allMin, t)
				allMax = math.Max(allMax, t)
			}
		}
	}
	fmt.Printf("radiative lifetimes over the map (b): %.3g – %.3g s; reported ball lightning: peak 2–5 s, mean ≈ 10 s\n", allMin, allMax)
	fmt.Println("(c) what a free plasmoid needs to last as long as ball lightning (TE lifetime; independent of the energy):")
	fmt.Printf(" fully stripped sea-level air holds %.1f × Loschmidt free electrons\n", AirElectronsPerMolecule)
	fmt.Println(" R (cm) | τ (s) | kR        | f (THz)  | λ (µm) | free electrons ÷ Loschmidt")
	for _, r := range []float64{0.05, 0.12, 0.25, 0.50} {
		for _, tau := range []float64{1, 3, 10} {
			y, f, ne := ForLifetime(r, tau)
			fmt.Printf(" %6.0f | %5.0f | %9.4g | %8.1f | %6.3f | %.3g\n", r*100, tau, y, f/1e12, C0/f*1e6, ne/Loschmidt)
		}
	}

	fmt.Println("lifetime check, root finder vs large-kR closed form (R = 10 cm, plasma at eq. 49):")
	for _, y := range []float64{268.6, 1000, 3000} {
		ne := Eps0 * Me * (y*y - J1Zero*J1Zero) * C0 * C0 / (0.01 * Qe * Qe)
		s := Sphere{R: 0.1, Wp: math.Sqrt(ne * Qe * Qe / (Me * Eps0))}
		te, tm := s.FindMode(true), s.FindMode(false)
		fmt.Printf("   kR = %6.1f: TE %.4g s vs %.4g s, TM %s vs %.4g s\n", y, te.Lifetime, TauTE(0.1, y), lt(tm), TauTM(0.1, y))
	}

	fmt.Println("\n== 5. kronspice: the Appendix B plasmoid's spectrum (4 TE + 4 TM l = 1 overtones, Foster chain, AC) ==")
	modes := Spectrum(e.R, e.Wp, 4)
	peaks, err := SpectrumAC(modes)
	if err != nil {
		return err
	}
	fmt.Println(" mode | interior x | f mode (GHz)  | f kronspice   | Q mode     | Q kronspice | f err   | Q err")
	for i, p := range peaks {
		kind := "TE"
		if !p.Mode.TE {
			kind = "TM"
		}
		x := []float64{}
		if p.Mode.TE {
			x = J1Zeros(20)
		} else {
			x = TMZeros(4)
		}
		fmt.Printf(" %s%d  | %10.4f | %13.8f | %13.8f | %10.4g | %11.4g | %.1e | %.1e\n", kind, i%4+1, x[i%4], p.FMode/1e9, p.FPeak/1e9, p.QMode, p.QKron, p.FErr, p.QErr)
	}
	fmt.Println("(TE1 is Nachamkin's mode: αR = 4.4934, his 128 GHz)")
	rd, err := RunRingDown(modes[0], 400)
	if err != nil {
		return err
	}
	fmt.Printf("ring-down of TE1 (400 periods, normalized ω = 1): Q from energy decay %.5g vs %.5g; ledger heat %.4e = energy lost %.4e (mismatch %.1e of stored)\n",
		rd.QFit, rd.Q, rd.Heat, rd.Lost, rd.Mismatch)
	fmt.Printf("TE1 energy lifetime Q/ω = %.3g s at 128 GHz\n", modes[0].Lifetime)

	fmt.Println("\n== 6. kronspice: eq. 49 as a plasma-loaded line (Kron field network, shorted half wave π/α) ==")
	for _, c := range []struct {
		r float64
		n int
	}{{e.K / e.Alpha, 200}, {1.5, 50}, {1.5, 200}} {
		ld, err := RunLadder(c.r, c.n)
		if err != nil {
			return err
		}
		fmt.Printf("k/α = %7.3f, N = %3d: resonance ÷ αc = %.8f (kronspice), %.8f (lattice exact), %.8f (eq. 49); kronspice − lattice %.1e, lattice − eq. 49 %.1e\n",
			c.r, c.n, ld.WKron, ld.WLattice, ld.WPaper, math.Abs(ld.WKron/ld.WLattice-1), math.Abs(ld.WLattice/ld.WPaper-1))
	}
	return nil
}

func sec(t float64) string { return fmt.Sprintf("%.3g s", t) }

func nearestIndex(x float64) int {
	z := J1Zeros(x + 4)
	best := 0
	for i := range z {
		if math.Abs(z[i]-x) < math.Abs(z[best]-x) {
			best = i
		}
	}
	return best + 1
}

func fmtLen(r float64) string {
	if r >= 1e-3 {
		return fmt.Sprintf("%.0f mm", r*1e3)
	}
	return fmt.Sprintf("%.0f µm", r*1e6)
}
