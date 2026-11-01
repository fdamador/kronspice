// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "fmt"

func part3and4(viz string) error {
	coil := DefaultCoil()
	fmt.Println("\nPart 3 — tune the flash circuit to the coil (tuned.go)")
	fmt.Printf("   tuning Lp·Cp = Ls·Cs: with the 5-turn primary Cp = %.3g µF; coupling k = %.3f, full transfer after 1/(2k·f1) = %.3g µs\n",
		coil.TunedC()*1e6, coil.K(), coil.Transfer()*1e6)
	rows, err := TunedComparison()
	if err != nil {
		return err
	}
	fmt.Printf("   %-56s %8s %12s %14s %10s %12s %12s\n", "", "flashes", "spark [A]", "to coil/flash", "of flash", "of supply", "coil peak")
	for _, r := range rows {
		fmt.Printf("   %-56s %8d %12.4g %12.3g J %9.1f %% %11.2f %% %10.4g V\n", r.Label, r.Run.Flashes, r.Run.PeakIp, r.PerFlash, 100*r.Share, 100*r.OfSupply, r.Run.CoilPeak)
	}
	fmt.Println("   (charging a capacitor from a fixed supply through a resistor loses half the energy in the resistor, so 'of supply' is at most 50 %)")
	ev, err := EvolveTuning(24, 12, 1)
	if err != nil {
		return err
	}
	b := ev.Best
	fmt.Printf("   genetic algorithm (24 × 12, seed 1) on the mist gap: capacitor %.3g µF, %d primary turns, quench after %.2f of the transfer time\n",
		b.C*1e6, b.Turns, ev.Quench)
	fmt.Printf("      → %.1f %% of each flash's energy reaches the coil (%.3g J per flash, spark %.0f A, coil peak %.4g V); %.1f %% of the supply energy; ledger %.1e\n",
		100*b.Share, b.PerFlash, b.Run.PeakIp, b.Run.CoilPeak, 100*b.OfSupply, b.Run.Residual)
	fmt.Println("      (20 turns is the top of the allowed range: a longer primary would raise its impedance further)")
	fmt.Printf("      Caution: in a bifilar coil neighbouring turns carry about half the coil voltage, here about %.0f kV, far above what\n", b.Run.CoilPeak/2/1e3)
	fmt.Println("      enamelled wire holds (a few kV, declared); this coil would need much thicker insulation or a smaller drive (§3e rule).")

	fmt.Println("\nPart 4 — plates in parallel on the spark side (parallel.go)")
	g, err := ParallelGaps(b.C, b.Turns, ev.Quench)
	if err != nil {
		return err
	}
	fmt.Println("   (a) N mist gaps in parallel across the evolved capacitor: only the first to break down carries the flash")
	fmt.Printf("   %6s %18s %20s %12s %16s %12s\n", "gaps", "strike [V]", "flash interval", "spark [A]", "to coil/flash", "of flash")
	for _, r := range g {
		fmt.Printf("   %6d %11.0f ± %4.0f %13.3g s ± %3.0f%% %12.0f %14.3g J %10.1f %%\n", r.N, r.MeanStrike, r.SDStrike, r.MeanInterval, 100*r.Jitter/r.MeanInterval, r.PeakIp, r.ToCoil, 100*r.Share)
	}
	fmt.Println("   More gaps give more starting electrons, so a gap also re-strikes before the others have recovered: smaller,")
	fmt.Println("   more frequent flashes; the spark current does not add.")
	p, err := ParallelPlates(b.C, ev.Quench)
	if err != nil {
		return err
	}
	fmt.Println("   (b) n capacitor sections in parallel, primary turns re-chosen to stay tuned")
	fmt.Printf("   %6s %10s %8s %10s %12s %14s %16s %10s\n", "plates", "C [µF]", "turns", "detuning", "spark [A]", "flash energy", "to coil/flash", "of flash")
	for _, r := range p {
		fmt.Printf("   %6d %10.3g %8d %9.1f%% %12.0f %12.3g J %14.3g J %8.1f %%\n", r.N, r.C*1e6, r.Turns, 100*r.Detuning, r.PeakIp, r.PerFlash, r.ToCoil, 100*r.Share)
	}
	fmt.Println("   More plates raise the energy per flash (×n) and the spark current nearly as much (the retuned primary has less")
	fmt.Println("   inductance); fewer primary turns couple less,")
	fmt.Println("   so the share reaching the coil falls — the coil and primary would be redesigned together for a bigger bank.")
	if viz == "" {
		return nil
	}
	c, cl := tuneDesign([]float64{b.C, float64(b.Turns), ev.Quench})
	r, err := c.RunCoil(cl, 6*c.R*c.C, 7)
	if err != nil {
		return err
	}
	return writeCoil(viz, "evolved", r)
}
