// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command magnetic_gas asks whether a gas can be made magnetic by tearing
// an electron off it, compressed to 500 psi, sparked, and used to make
// pulses: argon first, with hydrogen, oxygen and oxyhydrogen (HHO) beside
// it.
//
// Parts:
//
//  1. The gases (gas.go): measured Paschen minima, ionization energies,
//     magnetism of the neutral molecule and of its ion.
//
//  2. Breakdown (gas.go): the sparking voltage of each gas at 1 atm and at
//     500 psi for gaps of 0.05–1 mm, from the measured Paschen minimum,
//     with a check against air and Raizer's coefficients.
//
//  3. Magnetism (magnet.go): neutral susceptibility at 1 atm and 500 psi,
//     the pull of a strong magnet against the gas's weight, and the
//     susceptibility after one electron is removed (ions and free
//     electrons, Curie and Landau).
//
//  4. Compression (magnet.go): magnetic pressure B²/2μ₀, the field that
//     holds the 500 psi fill cold and as a plasma, the heating of a fast
//     compression, the energy to ionize the whole fill, and for HHO the
//     chemical energy it holds.
//
//  5. Marx generator (circuits.go): four stages charged in parallel and
//     fired in series through gas gaps at 500 psi, in kronspice, charging
//     and firing phases each with its energy ledger.
//
//  6. ECE2 beside the standard prediction (H48): the argon Marx run plain,
//     with the ECE2 torsion term and with the ECE2 scalar term.
//
//  7. Plasma varactor charge pump (circuits.go): a gas gap whose
//     capacitance a plasma column changes, steered by two diodes; the work
//     that changes the capacitance is on the ledger.
//
//     go run ./examples/magnetic_gas [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// eng prints a value with an SI prefix.
func eng(v float64, unit string) string {
	if v == 0 || math.IsNaN(v) || math.IsInf(v, 0) {
		return fmt.Sprintf("%g %s", v, unit)
	}
	pre := []struct {
		f float64
		p string
	}{{1e12, "T"}, {1e9, "G"}, {1e6, "M"}, {1e3, "k"}, {1, ""}, {1e-3, "m"}, {1e-6, "µ"}, {1e-9, "n"}, {1e-12, "p"}, {1e-15, "f"}}
	a := math.Abs(v)
	for _, x := range pre {
		if a >= x.f*0.9995 {
			return fmt.Sprintf("%.3g %s%s", v/x.f, x.p, unit)
		}
	}
	return fmt.Sprintf("%.3g %s", v, unit)
}

func main() {
	viz := flag.String("viz", "", "directory for CSV picture data")
	flag.Parse()
	if *viz != "" {
		if err := os.MkdirAll(*viz, 0o755); err != nil {
			log.Fatal(err)
		}
	}
	p500 := 500 * psi
	fmt.Printf("500 psi = %.4g MPa = %.2f atm = %.0f Torr (taken as absolute pressure; 500 psig would be %.1f atm)\n\n",
		p500/1e6, p500/atm, p500/torr, (p500+atm)/atm)

	fmt.Println("Part 1 — the gases")
	for _, g := range Gases {
		fmt.Println("  " + g.Describe())
	}
	fmt.Println("  HHO is 2 H₂ : 1 O₂. Its Paschen numbers are mole-weighted from H₂ and O₂ (an ESTIMATE: no measured curve found).")
	fmt.Println("  HHO is fuel and oxidizer together: any spark in it ignites it (Part 4 gives the energy it holds).")

	fmt.Println("\nPart 2 — breakdown voltage (uniform field; Townsend through the measured Paschen minimum)")
	gaps := []float64{0.05e-3, 0.1e-3, 0.2e-3, 0.5e-3, 1e-3}
	fmt.Printf("  %-10s", "gap")
	for _, g := range Gases {
		fmt.Printf(" %13s", g.Label()+" 1 atm")
	}
	for _, g := range Gases {
		fmt.Printf(" %15s", g.Label()+" 500 psi")
	}
	fmt.Println()
	for _, d := range gaps {
		fmt.Printf("  %-10s", fmt.Sprintf("%.2f mm", d*1e3))
		for _, g := range Gases {
			fmt.Printf(" %13s", eng(g.Breakdown(atm, d), "V"))
		}
		for _, g := range Gases {
			fmt.Printf(" %15s", eng(g.Breakdown(p500, d), "V"))
		}
		fmt.Println()
	}
	air := Air.Breakdown(atm, 0.01)
	fmt.Printf("  check, air 1 cm at 1 atm: %s by this rule against about 30 kV measured: it reads %.2f× high at large pd.\n", eng(air, "V"), air/30e3)
	fmt.Println("  So read the 500 psi column as an upper bound; above about 10 bar measured gaps also fall below Paschen's law.")
	for _, g := range []Gas{Argon, Hydrogen} {
		fmt.Printf("  cross-check %s, 0.2 mm at 500 psi: %s (this rule) vs %s (Raizer A = %g, B = %g, γ = 0.01)\n",
			g.Formula, eng(g.Breakdown(p500, 0.2e-3), "V"), eng(g.BreakdownRaizer(p500, 0.2e-3), "V"), g.A, g.B)
	}
	if *viz != "" {
		writePaschen(filepath.Join(*viz, "paschen.csv"))
	}

	fmt.Println("\nPart 3 — magnetism (SI volume susceptibility χ; + is pulled into a field, − pushed out)")
	fmt.Printf("  %-12s %12s %12s %16s %16s\n", "gas", "χ 1 atm", "χ 500 psi", "pull/weight 1 atm", "pull/weight 500 psi")
	for _, g := range Gases {
		fmt.Printf("  %-12s %12.3g %12.3g %16.3g %16.3g\n", g.Label(), g.ChiNeutral(atm, Troom), g.ChiNeutral(p500, Troom),
			g.ForceOverWeight(atm, Troom, 100), g.ForceOverWeight(p500, Troom, 100))
	}
	fmt.Println("  (pull/weight: magnetic force over the gas's own weight in a 100 T²/m gradient, a strong lab magnet; ratio is independent of pressure)")
	fmt.Println("  Removing one electron, sealed 500 psi fill, at 1 eV (11,605 K) — ions + free electrons against the neutral gas:")
	fmt.Printf("  %-12s %-5s %10s %12s %12s %12s %12s\n", "gas", "ion", "ionized", "χ neutral", "χ ions", "χ electrons", "χ total")
	T1 := qe / kB
	for _, g := range Gases {
		for _, f := range []float64{0.01, 1} {
			z := g.Ionized(p500, f, T1)
			fmt.Printf("  %-12s %-5s %9.0f%% %12.3g %12.3g %12.3g %12.3g\n", g.Label(), g.Ion, 100*f, z.Neutral, z.Ions, z.Electrons, z.Total)
		}
	}
	z0, z1 := Argon.ChiNeutral(p500, Troom), Argon.Ionized(p500, 1, Troom).Total
	fmt.Printf("  argon at room temperature, fully ionized (an upper bound): χ %.3g → %.3g, ×%.0f — argon does turn paramagnetic.\n", z0, z1, math.Abs(z1/z0))
	o0, o1 := Oxygen.ChiNeutral(p500, Troom), Oxygen.Ionized(p500, 1, Troom).Total
	fmt.Printf("  oxygen goes the other way: O₂ has 2 unpaired electrons, O₂⁺ has 1, so χ %.3g → %.3g at room temperature.\n", o0, o1)

	fmt.Println("\nPart 4 — compression")
	for _, B := range []float64{1, 3, 10, 30, 100} {
		fmt.Printf("  B = %5.0f T: magnetic pressure %s = %.4g atm = %.4g psi\n", B, eng(MagneticPressure(B), "Pa"), MagneticPressure(B)/atm, MagneticPressure(B)/psi)
	}
	fmt.Printf("  field to hold the 500 psi fill at room temperature (β = 1): %.2f T\n", FieldToHold(p500))
	for _, f := range []float64{0.01, 1} {
		pp := PlasmaPressure(p500, f, T1)
		fmt.Printf("  same fill heated to 1 eV, %3.0f%% ionized: %s (%.0f atm), field to hold it %.1f T\n", 100*f, eng(pp, "Pa"), pp/atm, FieldToHold(pp))
	}
	fmt.Println("  Fast one-stroke compression from 1 atm, 20 °C, to 500 psi heats the gas (adiabatic):")
	for _, g := range Gases {
		T2 := g.AdiabaticT(Troom, atm, p500)
		note := ""
		if g.Fuel {
			note = "  ← near hydrogen's auto-ignition (about 560–585 °C, declared): compress HHO slowly and cooled, or not at all"
		}
		fmt.Printf("  %-12s γ %.3f → %4.0f K (%4.0f °C)%s\n", g.Label(), g.Gamma, T2, T2-273.15, note)
	}
	V1 := 1e-6
	for _, g := range Gases {
		fmt.Printf("  energy to singly ionize all of 1 cm³ of %s at 500 psi: %s\n", g.Label(), eng(g.IonizeEnergy(p500, V1), "J"))
	}
	fmt.Printf("  chemical energy in 1 cm³ of HHO at 500 psi (burned to vapour): %s — a spark ignites it; detonation multiplies the fill pressure many times over (not computed here)\n", eng(HHOChemical(p500, V1), "J"))

	fmt.Println("\nPart 5 — Marx generator: 4 stages × 10 nF, 0.2 mm gaps at 500 psi, charged to 90 % of breakdown, 1 µH, 50 Ω load (kronspice)")
	fmt.Printf("  %-12s %10s %10s %10s %10s %7s %9s %10s %10s %10s %9s %9s\n", "gas", "breakdown", "charge to", "stage min", "output pk", "× Vc", "rise", "stored", "to load", "in arcs", "books ch.", "books fire")
	var argonMarx MarxResult
	for _, g := range Gases {
		r, err := DefaultMarx(g).Run()
		if err != nil {
			log.Fatal(err)
		}
		if g.Formula == Argon.Formula {
			argonMarx = r
		}
		smin := math.Inf(1)
		for _, v := range r.StageV {
			smin = math.Min(smin, v)
		}
		fmt.Printf("  %-12s %10s %10s %10s %10s %7.2f %9s %10s %10s %10s %9.1e %9.1e\n", g.Label(), eng(r.Vb, "V"), eng(r.Vc, "V"), eng(smin, "V"),
			eng(r.Peak, "V"), r.Gain, eng(r.Rise, "s"), eng(r.Stored, "J"), eng(r.ELoad, "J"), eng(r.EGaps, "J"), r.ChargeMismatch, r.Mismatch)
		if *viz != "" {
			writeMarx(filepath.Join(*viz, "marx_"+strings.ToLower(strings.Fields(g.Name)[0])+".csv"), DefaultMarx(g))
		}
	}
	fmt.Println("  (books: the ledger's integration mismatch, relative, charging and firing; point residuals are 0 to round-off.)")
	fmt.Printf("  The peak is %.2f·Vc, not 4·Vc: the erected bank (C/4 = 2.5 nF) drains into the load in RC = 130 ns while the 1 µH loop holds the\n", argonMarx.Gain)
	fmt.Println("  rise back (L/R = 20 ns); with no loop inductance it would reach 4·Vc·50/52 = 3.85·Vc. Charging through Rs heats it by half of what")
	fmt.Println("  the supply gives, as for any RC charge. Oxygen holds the most, but sparking pure O₂ at 500 psi burns electrode metal: an oxygen")
	fmt.Println("  fire hazard. HHO's row is a shot into fuel and oxidizer: the gap detonates the fill (Part 4).")
	fmt.Printf("  argon: supply gave %s while charging, %s went to heat, %s stored; firing: %s to the load + %s in the arcs, integration mismatch %.1e of the stored energy\n",
		eng(argonMarx.ChargeIn, "J"), eng(argonMarx.ChargeHeat, "J"), eng(argonMarx.Stored, "J"), eng(argonMarx.ELoad, "J"), eng(argonMarx.EGaps, "J"), argonMarx.Mismatch)

	fmt.Println("\nPart 6 — ECE2 (H48, hypothetical) beside the standard prediction, argon Marx")
	names := []string{"standard (Maxwell)", "ECE2 torsion, ω₀ = 8.339e-8 /m", "ECE2 scalar, κ = 0.05"}
	var base MarxResult
	for i, m := range MarxECE(Argon) {
		r, err := m.Run()
		if err != nil {
			log.Fatal(err)
		}
		if i == 0 {
			base = r
		}
		fmt.Printf("  %-32s output peak %s (%+.3g %%), stored %s, heat %s, model-supplied %s, books %.1e\n", names[i], eng(r.Peak, "V"),
			100*(r.Peak/base.Peak-1), eng(r.StoredStart, "J"), eng(r.EHeat, "J"), eng(r.EModel, "J"), r.Mismatch)
	}
	fmt.Println("  Torsion: a series resistance −2cω₀L = −50 µΩ in the 1 µH loop, a 7e-7 change. Scalar: every terminal reads 1/(1 − 2κ) = 1.111× the")
	fmt.Println("  body voltage, so the probe sees 11 % more for the same stored charge, energy and heat; it delivers no energy (model-supplied ≈ 0).")
	fmt.Println("  ECE2 power is booked as model-supplied, never inflow. A probe and a calorimeter on one bench shot, against the standard row, decide.")

	fmt.Println("\nPart 7 — plasma varactor charge pump: 100 V in, 10 nF gap, plasma fills 2/3 of it (C ×3), 10 kHz, 100 nF into 1 MΩ")
	pr, err := DefaultPump(Argon).Run()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  argon: ideal r·Vin %s, simulated %s after %d cycles\n", eng(pr.Ideal, "V"), eng(pr.Vout, "V"), pr.Cycles)
	fmt.Printf("  power over the last 10 %%: from the 100 V supply %s, through the varactor's port (paid by the plasma supply) %s, into the load %s\n",
		eng(pr.PIn, "W"), eng(pr.PPlasma, "W"), eng(pr.PLoad, "W"))
	fmt.Printf("  ledger worst residual %.1e, integration mismatch %.1e of the energy in\n", pr.Residual, pr.Mismatch)
	fmt.Println("  What each gas's 1 mm varactor gap holds at 500 psi with the plasma off (the pump's ceiling), and what pressure alone would do:")
	for _, g := range Gases {
		q := DefaultPump(g)
		fmt.Printf("  %-12s holds %s; compressing the neutral gas 1 atm → 500 psi raises C only ×%.4f\n", g.Label(), eng(g.Breakdown(q.P, q.Gap), "V"),
			(1+g.EpsMinus1*q.P/atm)/(1+g.EpsMinus1))
	}
	if *viz != "" {
		writePump(filepath.Join(*viz, "pump_argon.csv"), DefaultPump(Argon))
	}
}
