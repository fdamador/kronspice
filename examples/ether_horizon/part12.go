// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// DonutShapes are the tube/ring ratios a/R printed.
var DonutShapes = []float64{0.02, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0}

// DonutPressures are the ether pressures Π = p∞ħc/K² compared.
var DonutPressures = []float64{0, 1e-4, 1e-2, 1}

// invDCache memoizes TorusInvD at the resolution used here.
var invDCache = map[float64]float64{}

// WallInvD is ⟨1/d⟩ for charge on the core wall at aspect t.
func WallInvD(t float64) float64 {
	if v, ok := invDCache[t]; ok {
		return v
	}
	v := TorusInvD(t, 300)
	invDCache[t] = v
	return v
}

// ShapeGrid is the fine t grid for the shape scan and the graph: 0.01 … 1.
func ShapeGrid() []float64 {
	var ts []float64
	for i := 0; i <= 60; i++ {
		ts = append(ts, 0.01*math.Pow(100, float64(i)/60))
	}
	return ts
}

func part12(vizDir string) {
	fmt.Println()
	fmt.Println("Part 12a. The donut void: the double-loop photon in the hollow core of an ether vortex ring")
	fmt.Println("  the wave turns over once per circuit (k_z = 1/2R), so it carries ħ/2 for every R and a; the core's")
	fmt.Println("  TE₁₁ cutoff adds a transverse x/a. g = Mc²/E_L still holds (the group velocity cancels).")
	fmt.Println("  a/R     light energy / free double-loop photon   share that goes round the ring")
	for _, t := range DonutShapes {
		d := DonutAt(t, 0, 0)
		fmt.Printf("  %-6.2f  %-40.2f  %.5f\n", t, math.Sqrt(XGuide*XGuide/(t*t)+0.25)/0.5, d.LongFrac)
	}
	fmt.Println("  The double-loop photon's wavelength (4πR) is far longer than any core can guide: most of its")
	fmt.Println("  energy is held across the channel, not carried round it.")

	fmt.Println()
	fmt.Println("Part 12b. The size balance at fixed shape: e(λ) = A/λ + B·λ + C·λ³ (light + charge, swirl, pressure)")
	fmt.Println("  the swirl's energy K·R·(ln 8R/a − 2) is linear in R at fixed shape: a constant-force string.")
	fmt.Println("  Virial: E_V = E_L + U − 3E_P, so g = 2(1 + U/E_L) − 2E_P/E_L.")
	fmt.Println("  Π       a/R    R/ƛ_C     swirl share  pressure share  g          stable in size")
	for _, pi := range DonutPressures {
		for _, t := range []float64{0.1, 0.4, 1.0} {
			d := DonutAt(t, pi, WallInvD(t))
			st := "yes"
			if d.D2 <= 0 {
				st = "no"
			}
			fmt.Printf("  %-6.0e  %-5.2f  %-8.3f  %-11.4f  %-14.4f  %-9.6f  %s\n", pi, t, d.RadiusOverLbar(), d.EV/d.M, d.EP/d.M, d.G, st)
		}
	}
	d0 := DonutAt(0.1, 0, 0)
	fmt.Printf("  With no pressure and no charge g = %.12f: the vortex ring is the string Part 11 needed.\n", d0.G)
	fmt.Println("  Pressure always pulls g below 2; the charge's static field pushes it above.")

	fmt.Println()
	fmt.Println("Part 12c. Does the ether hold the donut's shape? Total energy along a/R (size relaxed at each shape)")
	ts := ShapeGrid()
	fmt.Println("  Π       M at a/R = 0.01, 0.1, 1   [E0]          falls all the way to a/R = 1?")
	for _, pi := range DonutPressures {
		ds, fat := ShapeScan(ts, pi, WallInvD)
		fmt.Printf("  %-6.0e  %-9.4f %-9.4f %-9.4f       %v\n", pi, ds[0].M, ds[30].M, ds[len(ds)-1].M, fat)
	}
	fmt.Println("  No: for every pressure the mass falls as the core swells. Both the light (x/a) and the swirl")
	fmt.Println("  (ln 8R/a) push the wall out; only the pressure pushes in, and it loses. The donut fattens to a")
	fmt.Println("  horn torus (a = R), where the thin-ring formulas end; its end state needs the full ether solve.")

	fmt.Println()
	fmt.Println("Part 12d. The anomaly with no pressure: (g − 2)/2 = U/E_L, charge on the core wall")
	fmt.Println("  a/R     ⟨1/d⟩·R   (g − 2)/2      / measured")
	for _, t := range DonutShapes {
		d := DonutAt(t, 0, WallInvD(t))
		a := d.G/2 - 1
		fmt.Printf("  %-6.2f  %-8.4f  %-13.6e  %.4f\n", t, WallInvD(t), a, a/AeMeasured)
	}
	dh := DonutAt(1, 0, WallInvD(1))
	fmt.Printf("  It grows as the donut fattens and peaks at the horn torus, %.6f = %.1f%% of the measured\n", dh.G/2-1, 100*(dh.G/2-1)/AeMeasured)
	fmt.Println("  0.00115965 — the same direction as the swelling of 12c, and with no number fitted, but a torus")
	fmt.Println("  cannot get closer, and any ether pressure takes some back.")
	fmt.Printf("  Size: the horn donut scaled to m_e c² has R = %.2f ƛ_C = %.3g m — still Compton-sized.\n", dh.RadiusOverLbar(), dh.RadiusOverLbar()*Hbar/(Me*C0))

	fmt.Println()
	fmt.Println("Part 12e. kronspice: the guided crossed ring (shunt inductors give the cutoff)")
	fz := GuideRing{N: 200, L: 1e-6, Cd: 1e-9, Ld: 1e30}.Analytic()
	gr := GuideRing{N: 200, L: 1e-6, Cd: 1e-9}
	gr.Ld = 1 / (math.Pow(2*math.Pi*3*fz, 2) * gr.Cd) // cutoff at 3 × the ring's own note
	pk, err := gr.Peak()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  ring note %.6g Hz, cutoff %.6g Hz: kronspice %.6g Hz, √(f_c² + f_z²) = %.6g Hz (ladder %.6g)\n",
		fz, gr.Cutoff(), pk, math.Hypot(gr.Cutoff(), fz), gr.Analytic())

	if vizDir != "" {
		if err := writeDonutCSV(filepath.Join(vizDir, "donut_scan.csv"), ts); err != nil {
			log.Fatal(err)
		}
	}
}

func writeDonutCSV(path string, ts []float64) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "pi,t,lambda,M,EL,EV,EP,EU,g,R_over_lbar,long_frac")
	for _, pi := range DonutPressures {
		for _, t := range ts {
			d := DonutAt(t, pi, WallInvD(t))
			fmt.Fprintf(f, "%g,%.6g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.10g,%.6g,%.6g\n",
				pi, t, d.Lambda, d.M, d.EL, d.EV, d.EP, d.EU, d.G, d.RadiusOverLbar(), d.LongFrac)
		}
	}
	return nil
}
