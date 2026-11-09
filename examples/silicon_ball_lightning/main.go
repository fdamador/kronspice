// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command silicon_ball_lightning tests the chemical model of ball lightning
// (Abrahamson & Dinniss, Nature 403, 519, 2000) against the first
// spectrum of a natural ball (Cen, Yuan & Xue, PRL 112, 035001, 2014,
// reference/PhysRevLett.112.035001.pdf) and the laboratory silicon balls,
// and asks whether a welding arc can make them: the energy each element
// stores per gram, the time silicon nanoparticles take to burn
// (Deal–Grove), a welder's arc at Paiva et al.'s operating point in
// kronspice, and the 35 kV power line beside the Qinghai ball.
//
//	go run ./examples/silicon_ball_lightning
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

func main() {
	viz := flag.String("viz", "", "also write the CSVs behind the figure to this directory")
	flag.Parse()
	if *viz != "" {
		if err := writeViz(*viz); err != nil {
			fmt.Fprintln(os.Stderr, "silicon_ball_lightning:", err)
			os.Exit(1)
		}
	}
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "silicon_ball_lightning:", err)
		os.Exit(1)
	}
}

func run() error {
	fmt.Println("== 1. Observations ==")
	for _, o := range Observations {
		life := "—"
		if o.LifetimeS[1] > 0 {
			life = fmt.Sprintf("%g–%g s", o.LifetimeS[0], o.LifetimeS[1])
		}
		fmt.Printf("%-55s | %-34s | diameter %g–%g m | lifetime %s\n   %s\n", o.Source, o.Kind, o.DiameterM[0], o.DiameterM[1], life, o.Note)
	}

	fmt.Println("\n== 2. The chemical store: heat on burning in air vs the heat to make the vapor ==")
	fmt.Println(" element → oxide   | burns (kJ/g) | vaporize (kJ/g) | return | mass for 1 kJ (mg)")
	for _, e := range Elements {
		fmt.Printf(" %-2s → %-12s | %12.2f | %15.2f | %6.2f | %8.1f\n", e.Name, e.Oxide, e.OxidationHeat()/1e3, e.VaporizationCost()/1e3, e.Return(), e.MassFor(1e3)*1e3)
	}
	fmt.Println(" reduced silicon species (Abrahamson & Dinniss):")
	for _, s := range Species {
		fmt.Printf("   %-18s %6.2f kJ/g\n", s.Name, s.HPerMol/s.M)
	}
	fmt.Println(" the 'average' ball (100 W for 10 s ≈ 1 kJ, an upper bound if all heat leaves as light and heat):")
	si := Elements[0]
	fmt.Printf("   %.0f mg of silicon; soil holds it as SiO2, so the strike must first reduce it, paying the same %.1f kJ/g (and the vapor's %.1f kJ/g)\n",
		si.MassFor(1e3)*1e3, si.OxidationHeat()/1e3, si.VaporizationCost()/1e3)
	fmt.Println("   oxides (SiO2 sand, TiO2 rutile flux) store nothing: only reduced metal burns")

	fmt.Println("\n== 3. How long silicon nanoparticles take to burn (Deal–Grove, dry O2, scaled by pO2) ==")
	radii := []float64{5e-9, 10e-9, 25e-9, 50e-9, 100e-9}
	temps := []float64{1100, 1200, 1300, 1400, 1500}
	for _, p := range []struct {
		name string
		po2  float64
	}{{"sea-level air (pO2 0.21 atm)", 0.21}, {"Qinghai, 2530 m (pO2 0.16 atm)", 0.16}} {
		fmt.Printf(" %s — burn-through time (s):\n   r (nm) |", p.name)
		for _, T := range temps {
			fmt.Printf(" %7.0f K |", T)
		}
		fmt.Println()
		for _, r := range radii {
			fmt.Printf("   %6.0f |", r*1e9)
			for _, T := range temps {
				fmt.Printf(" %9.3g |", BurnTime(r, T, p.po2))
			}
			fmt.Println()
		}
	}
	fmt.Println(" (Deal–Grove is fitted at 973–1473 K on planar wafers; values above 1473 K and at nanometre radii are extrapolations)")
	fmt.Println(" at 1100–1500 K Deal–Grove burns 5–100 nm particles in minutes to days, far slower than the observed seconds.")
	fmt.Println(" temperature at which a particle burns through in the observed lifetime (sea-level air; extrapolated beyond 1473 K):")
	for _, r := range []float64{10e-9, 25e-9, 50e-9} {
		fmt.Printf("   r = %3.0f nm:", r*1e9)
		for _, life := range []float64{1.64, 3.5, 10} {
			if T := tempFor(r, life, 0.21); T > 4990 {
				fmt.Printf("  %.3g s → above 5000 K;", life)
			} else {
				fmt.Printf("  %.3g s → %.0f K;", life, T)
			}
		}
		fmt.Println()
	}
	fmt.Println(" (silicon melts at 1687 K; Stephan & Massey estimate about 3140 K for their silicon balls)")

	fmt.Println("\n== 4. kronspice: an AC stick welder at Paiva et al.'s operating point ==")
	w, err := TuneWeld(70, 22.5, 60, 120)
	if err != nil {
		return err
	}
	fmt.Printf("open-circuit 70 V rms, arc 22.5 V, 60 Hz: leakage L = %.3f mH gives %.1f A rms; arc voltage %.1f V rms; arc power %.0f W; ledger %.1e\n",
		w.L*1e3, w.Irms, w.Vrms, w.Power, w.LedgerMismatch)
	fmt.Println("(Paiva et al.: 20–25 V, 100–140 A, 60 Hz)")
	fmt.Println(" one second of arc, if a fraction f of its energy vaporizes the electrode metal:")
	fmt.Println("   metal | f     | vapor (mg) | chemical store (J) | 100 W ball lasts (s)")
	for _, name := range []string{"Si", "Fe", "Mn"} {
		e := byName(name)
		for _, f := range []float64{0.01, 0.1} {
			m := f * w.Power / e.VaporizationCost() // g
			store := m * e.OxidationHeat()
			fmt.Printf("   %-5s | %4.0f%% | %10.1f | %18.1f | %.2f\n", name, f*100, m*1e3, store, store/100)
		}
	}

	fmt.Println("\n== 5. kronspice: the 35 kV line 20 m from the Qinghai ball ==")
	for _, a := range []float64{0.15, 2.5} {
		l, err := RunLine(35e3/math.Sqrt(3), 50, a, 20)
		if err != nil {
			return err
		}
		fmt.Printf("ball radius %.2f m: C_m %.3g F, C_ball %.3g F; the most power the line can drive into it %.3g W (at R = %.3g Ω); power pulses at %.2f Hz\n",
			a, l.Cm, l.Cb, l.BestPower, l.BestR, l.RippleHz)
	}
	fmt.Println("observed: light oscillating at 99.4 Hz (Cen et al.); the ball shines ~ 100 W")
	return nil
}

func byName(n string) Element {
	for _, e := range Elements {
		if e.Name == n {
			return e
		}
	}
	return Element{}
}

// tempFor is the temperature at which a particle of radius r burns through
// in time t (bisection; the burn time falls as the temperature rises).
func tempFor(r, t, po2 float64) float64 {
	lo, hi := 800.0, 5000.0
	for i := 0; i < 200; i++ {
		m := 0.5 * (lo + hi)
		if BurnTime(r, m, po2) > t {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi)
}
