// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command sparkgap_diodebridge combines water diodes (aluminium electrolytic
// cells, see examples/sparkgap_diode) and spark gaps into bridges and
// related circuits for different jobs:
//
//	A. Mains DC supply — four lead–aluminium cells (Pollak's bridge).
//	B. Surge-protected supply — the same bridge behind an air-gap arrester
//	   and an aluminium–aluminium cell clamp, hit by a 4 kV surge.
//	C. Spark pulser — a water spark-gap bridge charges a capacitor that an
//	   air gap dumps into a load; three output designs.
//	D. Voltage doubler — Greinacher's circuit with water cells.
//	E. A varactor with memory — the aluminium oxide film as a tuning
//	   capacitor set by its formation voltage.
//	F. Frequency doubler — a bridge's output repeats twice per cycle.
//
// Circuits with only cells run in kronspice; circuits with spark gaps (a
// switch with memory) run in this example's small nodal solver (sim.go),
// checked against kronspice on the supply. Every run has an energy ledger.
//
//	go run ./examples/sparkgap_diodebridge [-viz dir]
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
	must := func(err error) {
		if err != nil {
			log.Fatal(err)
		}
	}

	fmt.Println("A. Mains DC supply: four lead–aluminium water cells (100 V peak, 60 Hz, 2,200 µF, 50 Ω)")
	s, l, err := RunSupply()
	must(err)
	fmt.Printf("   output %.1f V DC, ripple %.1f V, efficiency %.0f %%\n", s.Vdc, s.Ripple, 100*s.Efficiency)
	fmt.Printf("   twin: kronspice gives %.2f V (difference %.2g %%); ledger residual %.2g J of %.1f J supplied\n",
		s.KronVdc, 100*(s.Vdc-s.KronVdc)/s.KronVdc, l.Residual, l.Supplied)

	fmt.Println("\nB. Surge-protected supply: 4 kV, 1.2/50 µs surge from a 2 Ω source on the mains line, at the mains peak")
	var surges []SurgeCase
	for _, p := range []bool{false, true} {
		sc, err := RunSurge(p, 4000)
		must(err)
		surges = append(surges, sc)
		fmt.Printf("   %s\n", sc.Name)
		fmt.Printf("      largest voltage across a bridge cell %.0f V (its oxide sparks over at 140 V); extra energy into the bridge cells %.2f J\n", sc.PeakArm, sc.BridgeSpark)
		fmt.Printf("      arrester %.1f J (%d strike), clamp %.3f J, DC output peak %.0f V; ledger residual %.1g J\n", sc.Arrester, sc.Strikes, sc.Clamp, sc.PeakOut, sc.Ledger.Residual)
	}
	fmt.Printf("   The 0.2 mm air gap (breaks down at %.0f V, Paschen) takes the surge; the aluminium–aluminium cell, which blocks both\n", AirGap(0.2e-3).Vpos())
	fmt.Println("   ways up to 140 V, clamps what gets past the line inductance. The bridge cells stay below their sparking voltage.")

	fmt.Println("\nC. Spark pulser: 6 kV, 60 Hz transformer (1 kΩ) → water spark-gap bridge (ratio 3) → 0.1 µF → air gap → 1 µH + load")
	var pulsers []Pulser
	for _, d := range PulserDesigns {
		r, err := RunPulser(d, 6)
		must(err)
		pulsers = append(pulsers, r)
		fmt.Printf("   %s\n", r.Name)
		fmt.Printf("      %d pulses/s, peak %.0f A, %.2f J per pulse to the load; load %.0f W of %.0f W from the mains (%.1f %%)\n",
			int(r.Rate+0.5), r.PeakAmps, r.PerPulse, r.LoadPower, r.InPower, 100*r.LoadPower/r.InPower)
		fmt.Printf("      output gap conducting %.2g s of 0.1 s; bridge gaps' backward sparks %d (%.2f J); ledger residual %.1g J\n", r.OnTime, r.BackStrikes, r.BackEnergy, r.Ledger.Residual)
	}
	fmt.Println("   A self-breaking output gap fails both ways: with a ringing load the capacitor swings negative and the bridge gaps")
	fmt.Println("   re-strike, shorting the transformer; with a damped load the mains keeps the arc lit (power follow). Firing the")
	fmt.Println("   output gap at each zero crossing of the mains — when the transformer cannot follow — works: about half the mains")
	fmt.Println("   power reaches the load. Its backward sparks are small: the transformer and the full capacitor nearly balance.")

	fmt.Println("\nD. Voltage doubler (Greinacher) of lead–aluminium cells: 100 V peak, 60 Hz, 2,200 µF, 200 Ω (kronspice)")
	var doublers []Doubler
	for _, m := range []int{1, 2} {
		d, err := RunDoubler(m)
		must(err)
		doublers = append(doublers, d)
		fmt.Printf("   %d cell(s) per position: %.0f V DC, efficiency %.0f %%, reverse peak %.0f V; ledger mismatch %.1g\n", m, d.Vdc, 100*d.Efficiency, d.ReversePeak, d.Mismatch)
	}
	fmt.Println("   Each position sees twice the peak in reverse, so one cell sparks over; two in series hold.")

	fmt.Println("\nE. A varactor with memory: the aluminium oxide film tunes a 1 mH tank (kronspice AC)")
	var vars []Varactor
	for _, vf := range []float64{10, 20, 40, 70, 140} {
		v, err := RunVaractor(vf, vf/2)
		must(err)
		vars = append(vars, v)
		fmt.Printf("   formed at %3.0f V: film %6.2f µF → resonance %6.0f Hz (formula %6.0f Hz), Q %.0f\n", v.Vform, v.FilmC*1e6, v.F0, v.F0Formula, v.Q)
	}
	fwd, err := RunVaractor(140, -5)
	must(err)
	fmt.Printf("   biased forward (−5 V): resonance %s — the conducting cell shorts the tank\n", map[bool]string{true: "none", false: fmt.Sprintf("%.0f Hz", fwd.F0)}[fwd.F0 == 0])
	fmt.Println("   The film grows about 1.3 nm per forming volt, so its capacitance falls as 1/V(form) and the resonance rises as")
	fmt.Println("   √V(form). Unlike a semiconductor varactor the capacitance does not follow the bias: below the forming voltage it")
	fmt.Println("   stays fixed. It is set by the highest voltage the film has seen — a capacitor with memory. Re-forming higher")
	fmt.Println("   thickens the film; thinning it again takes time (the film slowly dissolves) and is not modelled.")

	fmt.Println("\nF. Frequency doubler: a bridge's output repeats at twice the mains frequency (amplitudes ÷ source peak)")
	hs, err := RunDoubling()
	must(err)
	hs = append([]Harmonics{IdealDoubling()}, hs...)
	fmt.Printf("   %-36s %8s %8s %8s %12s\n", "bridge", "DC", "2f", "4f", "2f share of AC")
	for _, h := range hs {
		share := fmt.Sprintf("%.0f %%", 100*h.H2Fraction)
		fmt.Printf("   %-36s %8.3f %8.3f %8.3f %12s\n", h.Name, h.DC, h.H2, h.H4, share)
	}
	fmt.Println("   60 Hz in, 120 Hz out: a filter tuned to 2f picks it out. The spark-gap bridge's jumps put more into 2f and")
	fmt.Println("   higher harmonics than a smooth bridge.")
	fmt.Println("\nThe cells and gaps add no energy: every joule at a load comes from the mains or the surge source; the rest is heat.")

	if *viz == "" {
		return
	}
	write := func(name, head string, rows func(f *os.File)) {
		f, err := os.Create(filepath.Join(*viz, name))
		must(err)
		defer f.Close()
		fmt.Fprintln(f, head)
		rows(f)
	}
	write("surge.csv", "case,t_us,v_arm", func(f *os.File) {
		for _, sc := range surges {
			for i := range sc.Time {
				fmt.Fprintf(f, "%q,%g,%g\n", sc.Name, sc.Time[i]*1e6, sc.VArm[i])
			}
		}
	})
	write("pulser.csv", "design,t_ms,v_store,i_out", func(f *os.File) {
		for _, r := range pulsers {
			for i := range r.Time {
				fmt.Fprintf(f, "%q,%g,%g,%g\n", r.Name, r.Time[i]*1e3, r.VStore[i], r.IOut[i])
			}
		}
	})
	write("varactor.csv", "vform,film_uF,f0,f0_formula,q", func(f *os.File) {
		for _, v := range vars {
			fmt.Fprintf(f, "%g,%g,%g,%g,%g\n", v.Vform, v.FilmC*1e6, v.F0, v.F0Formula, v.Q)
		}
	})
	write("doubling.csv", "bridge,dc,h2,h4", func(f *os.File) {
		for _, h := range hs {
			fmt.Fprintf(f, "%q,%g,%g,%g\n", h.Name, h.DC, h.H2, h.H4)
		}
	})
}
