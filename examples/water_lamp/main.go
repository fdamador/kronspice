// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command water_lamp puts water in place of the neon lamp of the classic
// flashing-lamp (relaxation) oscillator, three ways, beside the neon lamp
// itself, and then lets each one drive the Tesla bifilar coil of
// examples/pulsed_bifilar.
//
// The circuit: a supply Vs charges a capacitor C through R; the lamp sits
// across C. When the lamp strikes it discharges C down to its maintaining
// voltage and goes out; then C charges again. The flash rate is set by
// R·C and the lamp's strike and maintaining voltages.
//
// A gas lamp strikes only when a free electron is there to start the
// avalanche. Neon lamps were made reliable with a trace of a radioactive
// gas (krypton-85, a beta emitter, TekWiki "Neon Lamps"); without it, in
// darkness, the lamp waits for a cosmic-ray ionization (the "dark
// effect"). Here the wait is modelled as a random process: free electrons
// appear at a declared rate S, and each one starts a discharge with
// probability 1 − Vb/V once the voltage V is past breakdown Vb (declared
// form). All randomness comes from one seeded generator (CLAUDE.md §4).
//
// The four lamps (numbers cited or declared in lamps.go):
//
//  1. Neon NE-2 (reference): strikes at 76 V, maintains 58 V.
//  2. Water-vapour glow tube: a sealed tube of low-pressure water vapour,
//     at the measured Paschen minimum of water vapour (480 V at 0.6 Torr·cm,
//     Škoro et al., Phys. Rev. E 84, 055401, 2011). No radioactive starter;
//     first electrons from cosmic rays in the dark, or from an ultraviolet
//     LED (declared rate).
//  3. Spark gap in liquid water: two electrodes 1 mm apart in pure water.
//     Breakdown field from the Sandia water relation E·τ^(1/3) = 0.135 ± 12 %
//     (MV/cm, µs; Stygar et al., Phys. Rev. ST Accel. Beams 9, 070401,
//     2006), measured for τ = 0.011–0.6 µs: the slow charging here is far
//     outside that range, so its breakdown voltage is an extrapolation and is
//     labelled so. Pure water conducts (0.055 µS/cm), so a leakage current
//     flows all the time and splits water (electrolysis, chemical channel).
//  4. Air gap filled with water mist: Paschen's law for air at 1 atm, 1 mm,
//     raised by the factor clean fog gave Liu et al. 2018 for rod–plane gaps
//     (29.4/21.7 kV at 5 cm; extrapolated to 1 mm and labelled so), with the
//     hold-off recovery of a spark gap (dry value declared from Bhoomireddy
//     2016; the faster recovery with mist is hypothesis H5, a parameter).
//
// Part 2 drives the bifilar coil: the lamp's discharge current flows
// through the 5-turn primary of examples/pulsed_bifilar, coupled to the
// bifilar solenoid's first mode (4.853 mH, f1 = 41.6 kHz, Q ≈ 83,
// primary mutual 28.1 µH). The reduced coil model is checked against a
// kronspice transient of the same network.
//
//	go run ./examples/water_lamp [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

func main() {
	viz := flag.String("viz", "", "directory for traces (viz/render.py)")
	flag.Parse()

	fmt.Println("Part 1 — four lamps in the same flashing-lamp circuit (lamps.go)")
	var all []Run
	for _, c := range Cases() {
		r, err := c.Run(c.Seconds, 7)
		if err != nil {
			log.Fatal(err)
		}
		all = append(all, r)
		printRun(r)
		if *viz != "" {
			if err := writeTrace(*viz, c.Tag, r); err != nil {
				log.Fatal(err)
			}
		}
	}
	fmt.Println()
	fmt.Printf("   %-44s %10s %18s %20s %14s\n", "lamp", "supply", "strike [V]", "flash interval", "energy/flash")
	for _, r := range all {
		iv := "none in run"
		if r.Flashes > 1 {
			iv = fmt.Sprintf("%.3g s ± %.0f %%", r.MeanInterval, 100*r.JitterInterval/r.MeanInterval)
		}
		fmt.Printf("   %-44s %8.0f V %11.0f ± %3.0f %20s %11.3g J\n", r.Case.Name, r.Case.Vs, r.MeanStrike, r.SDStrike, iv, r.EnergyPerFlash)
	}

	fmt.Println("\nPart 2 — each lamp drives the bifilar coil through the 5-turn primary (coil.go)")
	coil := DefaultCoil()
	tw, err := coil.Twin()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   coil: %.4g mH, first mode %.4g kHz, Q %.0f, primary %.3g µH, mutual %.3g µH (k = %.3f)\n",
		coil.Ls*1e3, coil.F1/1e3, coil.Q, coil.Lp*1e6, coil.M*1e6, coil.M/math.Sqrt(coil.Lp*coil.Ls))
	fmt.Printf("   reduced coil model vs kronspice transient for the same pulse: largest difference %.2g of the peak\n", tw)
	fmt.Printf("   %-44s %14s %18s %16s %16s\n", "lamp", "flashes", "coil peak [V]", "energy to coil", "of flash energy")
	for _, c := range Cases() {
		if !c.DriveCoil {
			continue
		}
		r, err := c.RunCoil(coil, c.CoilSeconds, 7)
		if err != nil {
			log.Fatal(err)
		}
		share := 0.0
		if r.LampEnergy+r.CoilEnergy > 0 {
			share = r.CoilEnergy / r.FlashEnergy
		}
		fmt.Printf("   %-44s %14d %18.4g %13.3g J %15.1f %%\n", c.Name, r.Flashes, r.CoilPeak, r.CoilEnergy, 100*share)
		fmt.Printf("      ledger: supply %.4g J = R %.4g + lamp %.4g + leakage %.4g + coil %.4g + stored %.4g J (residual %.1e)\n",
			r.Supply, r.RHeat, r.LampEnergy, r.Leak, r.CoilEnergy, r.Stored, r.Residual)
		if *viz != "" {
			if err := writeCoil(*viz, c.Tag, r); err != nil {
				log.Fatal(err)
			}
		}
	}
	if err := part3and4(*viz); err != nil {
		log.Fatal(err)
	}
	fmt.Println("\nEvery joule comes from the supply; the water adds no energy. Without a radioactive starter a gas lamp in the dark")
	fmt.Println("waits for cosmic rays; an ultraviolet LED (or room light on a low-work-function cathode) does the starter's job.")
}

func printRun(r Run) {
	c := r.Case
	fmt.Printf("\n   == %s ==\n", c.Name)
	fmt.Printf("      %s\n", c.Note)
	fmt.Printf("      supply %.0f V, R %.3g Ω, C %.3g F; breakdown %s, maintains %.0f V; seed electrons %s\n",
		c.Vs, c.R, c.C, c.BreakdownText(), c.Lamp.Vm, c.SeedText())
	fmt.Printf("      %.0f s simulated: %d flashes", r.T, r.Flashes)
	if r.Flashes > 0 {
		fmt.Printf("; strike at %.0f ± %.0f V", r.MeanStrike, r.SDStrike)
	}
	if r.Flashes > 1 {
		fmt.Printf("; every %.3g s ± %.0f %%", r.MeanInterval, 100*r.JitterInterval/r.MeanInterval)
	}
	if c.Lamp.Seeds > 0 {
		fmt.Printf("\n      expected wait for a starting electron once past breakdown at the supply voltage: %.3g s", 1/(c.Lamp.Seeds*(1-c.Lamp.Vb/c.Vs)))
	}
	fmt.Println()
	fmt.Printf("      ledger: supply %.4g J = heat in R %.4g + lamp (light + heat) %.4g + leakage %.4g + stored %.4g J (residual %.1e)\n",
		r.Supply, r.RHeat, r.LampEnergy, r.Leak, r.Stored, r.Residual)
	if r.Chemical > 0 {
		fmt.Printf("      leakage through the water splits it: %.3g J into hydrogen and oxygen (1.23 V × charge), %.3g mol H₂, rest heat\n", r.Chemical, r.Hydrogen)
	}
}

func writeTrace(dir, tag string, r Run) error {
	f, err := os.Create(filepath.Join(dir, tag+"_trace.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "t_s,v_cap")
	for i := range r.TraceT {
		fmt.Fprintf(f, "%g,%g\n", r.TraceT[i], r.TraceV[i])
	}
	g, err := os.Create(filepath.Join(dir, tag+"_flashes.csv"))
	if err != nil {
		return err
	}
	defer g.Close()
	fmt.Fprintln(g, "t_s")
	for _, t := range r.StrikeTimes {
		fmt.Fprintf(g, "%g\n", t)
	}
	return nil
}

func writeCoil(dir, tag string, r CoilRun) error {
	f, err := os.Create(filepath.Join(dir, tag+"_coil.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "t_ms,v_coil")
	for i := range r.TraceT {
		fmt.Fprintf(f, "%g,%g\n", r.TraceT[i]*1e3, r.TraceV[i])
	}
	return nil
}
