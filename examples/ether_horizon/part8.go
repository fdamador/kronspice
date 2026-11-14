// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strconv"

	"github.com/fdamador/kronspice/output"
)

// FountainDT is the streamline tracing step (units of a/|u|).
const FountainDT = 2e-4

func part8(vizDir string) {
	fmt.Println()
	fmt.Println("Part 8a. The fountain: in round the equator, swirling, out at the poles (steady Euler, Bragg–Hawthorne)")
	fmt.Println("  swirl C(ψ) = √(C0² + k²ψ²): Δ*ψ + k²ψ = 0 exactly; C0 is a line vortex on the axis (the void's swirl)")
	fmt.Println("  l  flow                               ka       cell centre (r, θ)   residual   net flux r=½   cell flux 2πψ_O")
	for _, l := range []int{1, 2} {
		f := NewFountain(l, 0)
		ro, to := f.OPoint()
		name := "one cell: in at the south pole, out at the north"
		if l == 2 {
			name = "two cells: in at the equator, out at both poles"
		}
		fmt.Printf("  %d  %-33s  %.4f   (%.4f a, %5.2f°)   %.1e    %.1e        %.4f\n", l, name, f.K, ro, to*180/math.Pi, f.Residual(0.6, 1.0), f.SphereFlux(0.5, 4000), f.CellFlux())
	}
	fmt.Println("  Winding q (toroidal turns per poloidal turn) of the streamline tori, by tracing, against the axis circulation β = C0/(kψ_O):")
	xs := []float64{0.05, 0.3, 0.6, 0.9}
	for _, l := range []int{1, 2} {
		for _, beta := range []float64{0, 1, 2} {
			f := NewFountain(l, beta*NewFountain(l, 0).K*NewFountain(l, 0).PsiO())
			fmt.Printf("    l = %d, β = %.0f:", l, beta)
			for _, s := range f.QProfile(xs, FountainDT) {
				fmt.Printf("  q = %.3f (tube/ring %.2f)", s.Q, s.Minor/s.Major)
			}
			fmt.Println()
		}
	}
	fmt.Println("  Without the axis vortex (β = 0) every surface winds less than once: no 2:1 surface. The void's")
	fmt.Println("  swirl is what makes one. Fitting β so the q = 2 surface has the photon torus's tube/ring = ½:")
	for _, l := range []int{1, 2} {
		beta, s, ok := FitTorus(l, 2, 0.5, FountainDT)
		if !ok {
			log.Fatal("no q = 2 surface")
		}
		a := RingRadius() / s.Major
		where := "one ring, in the equatorial plane"
		if l == 2 {
			ro, to := NewFountain(l, 0).OPoint()
			where = fmt.Sprintf("two rings, at z = ±%.3f a (±%.0f fm)", ro*math.Cos(to), ro*math.Cos(to)*a*1e15)
		}
		fmt.Printf("    l = %d: β = %.4f, q = %.4f, tube/ring = %.4f, loop closes to %.0e a; scaled to R = %.0f fm the sphere is a = %.0f fm; %s\n",
			l, beta, s.Q, s.Minor/s.Major, s.Closure, RingRadius()*1e15, a*1e15, where)
		if vizDir != "" && l == 2 {
			writeFountainViz(vizDir, NewFountain(l, beta*NewFountain(l, 0).K*NewFountain(l, 0).PsiO()), s)
		}
		if vizDir != "" && l == 1 {
			writeFountainVizL1(vizDir, NewFountain(l, beta*NewFountain(l, 0).K*NewFountain(l, 0).PsiO()), s)
		}
	}
	f2 := NewFountain(2, 1.6*NewFountain(2, 0).K*NewFountain(2, 0).PsiO())
	fmt.Printf("  Books: nothing is destroyed (net flux through every sphere 0); each cell carries %.4f; angular momentum %.4f,\n", f2.CellFlux(), f2.AngularMomentum(400))
	fmt.Printf("  energy outside a void of radius 0.01 a: %.4f, 0.001 a: %.4f (the axis vortex's ln(1/ξ)) — units ρ = a = A = 1.\n", f2.Energy(0.01, 300), f2.Energy(0.001, 300))
	fmt.Println("  The 2:1 surface is a fit (one knob, one condition), not a prediction. The question's flow (l = 2) gives")
	fmt.Println("  two tori, above and below the equator; the photon electron's single equatorial ring sits on l = 1.")

	fmt.Println()
	fmt.Println("Part 8b. The spinning core as a battery (membrane circuit): EMF ∝ Ω_H, horizon R_H inside, load R_L outside")
	fmt.Println("  R_L/R_H   Ω_F/Ω_H   P_load   P_horizon   P_spin   ledger   Γ_load    Γ_core (Z = −R_L)   Γ_core·Γ_load")
	ratios := []float64{0.1, 0.5, 1, 2, 10}
	pts, err := SpinSweep(ratios)
	if err != nil {
		log.Fatal(err)
	}
	for _, p := range pts {
		gc := fmt.Sprintf("%-8.4f", real(p.GammaCore))
		prod := fmt.Sprintf("%.6f", real(p.GammaCore*p.GammaLoad))
		if math.IsInf(real(p.GammaCore), 0) {
			gc, prod = "∞       ", "— (0·∞, matched)"
		}
		fmt.Printf("  %-8.2f  %-8.4f  %-7.4f  %-10.4f  %-7.4f  %-7.0g  %-8.4f  %s            %s\n",
			p.RLoverRH, p.OmegaF, p.PLoad, p.PHorizon, p.PSpin, p.Residual, real(p.GammaLoad), gc, prod)
	}
	cs, err := CoreSmith()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  The core seen by kronspice's network-analyser sweep (sources off): Z = %.6f R_H, Γ = %.1e, %s.\n",
		real(cs.Points[0].Z), real(cs.Points[0].Gamma), cs.Points[0].Region)
	fmt.Println("  So the horizon itself is +R_H (it absorbs). At the operating point the same core looks like −R_L to")
	fmt.Println("  its load: outside the unit circle, paid by the spin column. Maximum power at the match, Ω_F = Ω_H/2.")
	fmt.Printf("  Superradiance for the ring (m = ½, from the 720° double loop): the core must spin faster than Ω = 2ω_C = %.3g rad/s;\n", SuperradianceOmega())
	fmt.Printf("  its rim speed is c·r_c/R, below c for any core smaller than the ring (r_c = R/4: %.2f c).\n", RimSpeed(RingRadius()/4)/C0)
	fmt.Printf("  But the ring leaks (photon_electron §3: Q = %.2f), and a spin of ħ/2 at that rate holds E_rot = ħΩ/4 = m_e c²/2:\n", RingQ())
	fmt.Printf("  through a matched load it keeps the ring alive for %.2g s — about %.2f of a period — then the spin is gone.\n",
		SpinLifetime(SuperradianceOmega()), SpinLifetime(SuperradianceOmega())*OmegaC()/(2*math.Pi))
	fmt.Println("  The electron is stable for > 6.6×10²⁸ years: the spin would need resupplying, from a column not yet in the model.")
	if vizDir != "" {
		writeSpinSmith(vizDir)
	}
}

func writeCSV(path string, header []string, rows [][]float64) {
	f, err := os.Create(path)
	if err != nil {
		log.Fatal(err)
	}
	w := csv.NewWriter(f)
	w.Write(header)
	for _, r := range rows {
		rec := make([]string, len(r))
		for i, v := range r {
			rec[i] = strconv.FormatFloat(v, 'g', 7, 64)
		}
		w.Write(rec)
	}
	w.Flush()
	if err := f.Close(); err != nil {
		log.Fatal(err)
	}
}

func fountainGrid(f Fountain) [][]float64 {
	var rows [][]float64
	for i := 0; i <= 120; i++ {
		s := float64(i) / 120
		for j := 0; j <= 240; j++ {
			z := -1 + 2*float64(j)/240
			r := math.Hypot(s, z)
			psi := math.NaN()
			if r <= 1 && r > 0 && s > 0 {
				psi = f.Psi(r, math.Atan2(s, z))
			}
			rows = append(rows, []float64{s, z, psi})
		}
	}
	return rows
}

func surfaceRows(s Surface) [][]float64 {
	var rows [][]float64
	for i := 0; i < len(s.S); i += 20 {
		rows = append(rows, []float64{s.S[i], s.Z[i]})
	}
	return rows
}

func writeFountainViz(dir string, f Fountain, s Surface) {
	writeCSV(filepath.Join(dir, "fountain_l2_psi.csv"), []string{"s", "z", "psi"}, fountainGrid(f))
	writeCSV(filepath.Join(dir, "fountain_l2_q2.csv"), []string{"s", "z"}, surfaceRows(s))
}

func writeFountainVizL1(dir string, f Fountain, s Surface) {
	writeCSV(filepath.Join(dir, "fountain_l1_psi.csv"), []string{"s", "z", "psi"}, fountainGrid(f))
	writeCSV(filepath.Join(dir, "fountain_l1_q2.csv"), []string{"s", "z"}, surfaceRows(s))
}

func writeSpinSmith(dir string) {
	var load, core []complex128
	for i := 0; i <= 200; i++ {
		x := math.Pow(10, -1.5+3*float64(i)/200)
		load = append(load, output.SmithGamma(complex(x, 0), 1))
		core = append(core, output.SmithGamma(complex(-x, 0), 1))
	}
	cs, err := CoreSmith()
	if err != nil {
		log.Fatal(err)
	}
	f, err := os.Create(filepath.Join(dir, "spin_core_smith.svg"))
	if err != nil {
		log.Fatal(err)
	}
	err = output.WriteSmithSVG(f, []output.SmithTrace{
		{Name: "load: +R_L", Gamma: load, Color: "#1f6fd1"},
		{Name: "core running: −R_L", Gamma: core, Color: "#d1361f"},
		{Name: "core, sources off: +R_H", Gamma: []complex128{cs.Points[0].Gamma, cs.Points[0].Gamma}, Color: "#2a9d3c"},
	}, output.SmithChartOptions{Title: "H49 Part 8b: spinning core and its load, Z0 = R_H (R_L/R_H from 0.03 to 30)", Extent: 4})
	if err != nil {
		log.Fatal(err)
	}
	if err := f.Close(); err != nil {
		log.Fatal(err)
	}
}
