// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command forward_casimir works through R. L. Forward, "Extracting
// electrical energy from the vacuum by cohesion of charged foliated
// conductors," Phys. Rev. B 30, 1700 (1984): the "vacuum-fluctuation
// battery".
//
// Charged aluminium leaves are drawn together by the Casimir force while a
// servo-controlled supply holds their electrostatic repulsion just below
// it, so the closing leaves push charge back into the supply.
//
//  1. Casimir energy and pressure: Forward's ideal-mirror formulas and the
//     Lifshitz result for real aluminium.
//
//  2. One leaf pair as a Kron electromechanical system (charge q and gap a
//     as coordinates), closed under Forward's servo, full energy ledger.
//
//  3. The electrical side in kronspice (a capacitor varying with the gap).
//
//  4. Reopening ("recharging", Forward p. 1702): the whole cycle.
//
//  5. Scale, and what it means for electron pumps.
//
//     go run ./examples/forward_casimir [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"path/filepath"
)

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()

	fmt.Println("1. Casimir energy and pressure (per m²) — Forward eqs. 1–2, and real aluminium (Lifshitz, plasma model, λp = 107 nm)")
	fmt.Printf("   %8s %14s %14s %12s %14s %14s\n", "gap", "energy ideal", "energy Al", "Al ÷ ideal", "pressure ideal", "pressure Al")
	al := NewTable(Aluminium, 5e-9, 10e-6, 400)
	for _, a := range []float64{10e-9, 20e-9, 100e-9, 1e-6} {
		fmt.Printf("   %6.0f nm %11.3e J %11.3e J %12.3f %11.3e Pa %11.3e Pa\n", a*1e9, PerfectEnergy(a), al.Energy(a), al.Energy(a)/PerfectEnergy(a), PerfectPressure(a), al.Pressure(a))
	}
	fmt.Println("   Real metal is not a perfect mirror above its plasma frequency: at 20 nm aluminium gives a fifth of Forward's")
	fmt.Println("   ideal energy. (Forward eq. 4: below λ/2π the force turns to the van der Waals form; the plasma model includes it.)")

	fmt.Println("\n2. One leaf pair closing from 1 µm to 20 nm under Forward's servo (100 nm Al leaves, 1 mm/s; Kron's (q, a))")
	var strokes []Stroke
	var leaves []Leaf
	for _, m := range []Material{Perfect, Aluminium} {
		l := DefaultLeaf(m)
		c := l.Run(l.A0, l.AMin)
		strokes, leaves = append(strokes, c), append(leaves, l)
		fmt.Printf("   %s leaves:\n", m.Name)
		fmt.Printf("      Casimir work %.4e J/m²; servo voltage rises to %.1f V\n", c.Casimir, c.PeakVolt)
		fmt.Printf("      during the stroke the supply gives %.4e J/m² and the field stores %.4e J/m²\n", -c.ToSupply, c.Stored)
		fmt.Printf("      after discharging at the stop the supply has gained %.4e J/m² (heat %.1e, motion %.1e, residual %.1e)\n", c.Net, c.Heat, c.Kinetic, c.Residual)
		fmt.Printf("      the supply sends %.3e C/m² into the leaves during the stroke; at discharge all %.2e electrons/m² come back,\n", c.ChargeIn, c.Electrons)
		fmt.Printf("      each carrying %.3f eV of Casimir energy out (net charge moved over the whole operation: zero)\n", c.PerElectron)
	}
	fmt.Println("   The Casimir pull grows as 1/a⁴, so the servo must raise the voltage steeply to hold it: the field stores")
	fmt.Println("   hundreds of times the Casimir work. Forward's gain appears only after the leaves are stopped and discharged:")
	fmt.Println("   the supply then ends up with the Casimir work, less the servo's losses.")

	fmt.Println("\n3. The electrical side in kronspice (capacitor C(t) = 4ε0/(L − a(t)) from the stroke, servo voltage, 1 µΩ·m²)")
	for i, l := range leaves {
		tw, err := KronTwin(l, strokes[i])
		if err != nil {
			log.Fatal(err)
		}
		c := strokes[i]
		fmt.Printf("   %-15s supply gives %.5e J (own %.5e), sends %.5e C (own %.5e); leaves' work on the capacitor %.5e J (Casimir %.5e); ledger mismatch %.0e\n",
			l.Mat.Name+":", -tw.ToSupply, -c.ToSupply, -tw.Charge, c.ChargeIn, tw.PumpWork, c.Casimir, tw.Mismatch)
	}

	fmt.Println("\n4. Reopening the leaves to 1 µm (Forward's \"recharge\") and the whole cycle")
	var opens []Stroke
	for i, l := range leaves {
		o := l.Run(l.AMin, l.A0)
		opens = append(opens, o)
		c := strokes[i]
		fmt.Printf("   %-15s close +%.4e J/m², reopen %.4e J/m², cycle %.3e J/m² — the supply pays for the reopening\n",
			l.Mat.Name+":", c.Net, o.Net, c.Net+o.Net)
	}
	fmt.Println("   The Casimir force is conservative: closing gives back what opening costs, and the servo's losses make every")
	fmt.Println("   cycle a net loss. It is a one-shot release of stored surface energy, like a spring, not a continuous source.")

	fmt.Println("\n5. Scale")
	c := strokes[1]
	m := 2 * Aluminium.Density * leaves[1].T
	fmt.Printf("   One stroke of real Al: %.2e J per m² of gap, %.3f J per kg of leaf metal (%.0f nm leaves, 2 per gap).\n", c.Net, c.Net/m, leaves[1].T*1e9)
	fmt.Println("   For comparison: a lithium-ion cell stores about 900,000 J/kg. Even pressing the leaves into full contact cannot")
	fmt.Println("   release more than their surface energy, about 1.1 J/m² per surface for aluminium (Tyson & Miller 1977) —")
	fmt.Println("   the same energy that had to be spent to make the separate leaves.")
	fmt.Printf("   As an electron pump: one stroke is a variable-capacitor charge pump driven by the Casimir force. It lifts\n")
	fmt.Printf("   %.1e electrons per m² by %.3f eV each, once; reopening costs that back. A continuous pump needs a continuous\n", c.Electrons, c.PerElectron)
	fmt.Println("   source of work, as the pumps of photon_electron do.")
	fmt.Println("\nNo energy is created: the ledger closes; what the supply gains is the Casimir (surface) energy of the leaves.")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "casimir.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(f, "gap_nm,ideal,al")
	for a := 8e-9; a <= 2e-6; a *= 1.05 { // inside the table (5 nm – 10 µm)
		fmt.Fprintf(f, "%g,%g,%g\n", a*1e9, -PerfectPressure(a), -al.Pressure(a))
	}
	f.Close()
	g, err := os.Create(filepath.Join(*viz, "stroke.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(g, "case,gap_nm,volt,charge")
	for i, s := range append(strokes, opens...) {
		name := []string{"close, perfect", "close, Al", "reopen, perfect", "reopen, Al"}[i]
		for k := 0; k < len(s.Gap); k += 5 {
			fmt.Fprintf(g, "%q,%g,%g,%g\n", name, s.Gap[k]*1e9, s.Volt[k], s.Charge[k])
		}
	}
	g.Close()
}
