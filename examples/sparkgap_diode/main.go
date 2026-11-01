// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command sparkgap_diode asks whether a water gap can act as a diode
// because of its two electrode materials, and builds diodes and diode
// bridges from it.
//
// It can, and it was done: a valve metal (aluminium) in water grows an
// insulating anodic oxide when it is the positive electrode and passes
// current freely when it is the negative one. Paired with a non-valve metal
// (lead, iron) it is the aluminium electrolytic rectifier of Pollak (1896)
// and Nodon (1899), rated about 140 V per cell at about 75–80 % efficiency
// (Encyclopaedia Britannica 1911, "Transformers"). When the reverse voltage
// exceeds the oxide's sparking voltage the oxide breaks down in sparks:
// the cell is then a spark gap in reverse, like a Zener diode.
//
// At thousands of volts the water itself breaks down in a spark, and what
// makes a spark one-way there is the electrodes' shape: a sharp positive
// point breaks down at a lower voltage than a sharp negative one. That is
// the second device, the polarity spark gap.
//
// Parts:
//
//  1. I–V of electrode pairs in water (cell.go), each curve checked against
//     a kronspice DC sweep of the same element.
//
//  2. Half-wave and bridge rectifiers of aluminium–lead cells, solved in
//     kronspice (.TRAN) with the energy ledger; 60 Hz and 1 kHz.
//
//  3. Polarity spark gaps as rectifiers (gap.go, gapsim.go): when they work
//     and when they break down backwards.
//
//  4. Pulse test (pulse.go): the same 10 J through each device; peak amps,
//     rise time and the energy that reaches the load.
//
//  5. Frequency sweep (sweep.go): every device in the same bridge from
//     50 Hz to 100 kHz, as amps per volt against frequency.
//
//     go run ./examples/sparkgap_diode [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// duration prints seconds as a readable span.
func duration(s float64) string {
	switch {
	case s < 120:
		return fmt.Sprintf("%.0f s", s)
	case s < 7200:
		return fmt.Sprintf("%.0f min", s/60)
	case s < 2*86400:
		return fmt.Sprintf("%.1f h", s/3600)
	}
	return fmt.Sprintf("%.1f days", s/86400)
}

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()

	fmt.Println("Part 1 — electrode pairs in water (cell.go); current from the first electrode to the second")
	pairs := []Cell{DefaultCell(Pb, Al), DefaultCell(Fe, Al), DefaultCell(CuAgTips, Al), DefaultCell(FeAgCoat, Al),
		DefaultCell(Fe, Cu), DefaultCell(Pb, Fe), DefaultCell(Al, Al)}
	fmt.Printf("   %-38s %10s %10s %10s %10s %10s %14s %12s\n", "pair", "+5 V", "−5 V", "+100 V", "−100 V", "−200 V", "fwd/rev 100 V", "kronspice")
	for _, c := range pairs {
		tw, err := IVTwin(c, 200)
		if err != nil {
			log.Fatal(err)
		}
		ratio := "—"
		if r := math.Abs(c.Current(100) / c.Current(-100)); r > 1.5 {
			ratio = fmt.Sprintf("%.3g", r)
		}
		fmt.Printf("   %-38s %9.3g A %9.3g A %9.3g A %9.3g A %9.3g A %14s %12.1e\n", c.Name(), c.Current(5), c.Current(-5), c.Current(100), c.Current(-100), c.Current(-200), ratio, tw)
	}
	c := DefaultCell(Pb, Al)
	fmt.Println("   " + c.Describe())
	fmt.Println("   A valve metal against an ordinary metal rectifies; two ordinary metals do not (they make a weak battery:")
	for _, p := range pairs[4:6] {
		fmt.Printf("   %s %.3f V open circuit);", p.Name(), p.Galvanic())
	}
	fmt.Println("\n   two valve metals block both ways. Past the oxide's sparking voltage the blocked direction sparks over.")

	fmt.Println("\nPart 2 — aluminium–lead water cells as rectifiers (kronspice .TRAN, ledger on)")
	fmt.Printf("   %-50s %9s %9s %10s %12s %12s %10s\n", "circuit", "Vdc", "ripple", "load", "efficiency", "reverse pk", "ledger")
	type rr struct {
		r   Rectifier
		res RectResult
	}
	var runs []rr
	for _, r := range []Rectifier{
		{Name: "half-wave, 1 cell, 100 V peak, 60 Hz", Arm: Arm{c, 1}, Vp: 100, F: 60},
		{Name: "half-wave, 2 cells in series, 100 V peak, 60 Hz", Arm: Arm{c, 2}, Vp: 100, F: 60},
		{Name: "bridge, 1 cell per arm, 100 V peak, 60 Hz", Bridge: true, Arm: Arm{c, 1}, Vp: 100, F: 60},
		{Name: "bridge, 1 cell per arm, 170 V peak (120 V rms), 60 Hz", Bridge: true, Arm: Arm{c, 1}, Vp: 170, F: 60},
		{Name: "bridge, 2 cells per arm, 170 V peak, 60 Hz", Bridge: true, Arm: Arm{c, 2}, Vp: 170, F: 60},
		{Name: "bridge, 1 cell per arm, 100 V peak, 1 kHz", Bridge: true, Arm: Arm{c, 1}, Vp: 100, F: 1000},
		{Name: "bridge, 1 cell per arm, 100 V peak, 20 kHz", Bridge: true, Arm: Arm{c, 1}, Vp: 100, F: 20000},
	} {
		r.Rs, r.Cload, r.Rload, r.Cycles, r.Skip = 0.5, 2200e-6, 50, 30, 20
		if r.F > 100 {
			r.Cycles, r.Skip = int(r.F*0.3), int(r.F*0.3)-50 // 0.3 s: the filter (RC 0.11 s) settles
		}
		res, err := r.Run()
		if err != nil {
			log.Fatal(err)
		}
		runs = append(runs, rr{r, res})
		fmt.Printf("   %-50s %7.1f V %7.2f V %8.1f W %11.1f %% %10.0f V %10.1e\n", r.Name, res.Vdc, res.Ripple, res.Pload, 100*res.Efficiency, res.ReversePeak, res.Mismatch/res.EnergyIn)
	}
	fmt.Println("   Silver counter-electrodes against aluminium (bridge, 1 cell per arm, 100 V peak, 60 Hz):")
	for _, e := range []Electrode{CuAgTips, FeAgCoat} {
		cell := DefaultCell(e, Al)
		r := Rectifier{Name: "bridge, " + cell.Name(), Bridge: true, Arm: Arm{cell, 1}, Vp: 100, F: 60, Rs: 0.5, Cload: 2200e-6, Rload: 50, Cycles: 30, Skip: 20}
		res, err := r.Run()
		if err != nil {
			log.Fatal(err)
		}
		iCell := res.Idc / 2 // each arm carries the load current half the time
		m, q, sec := cell.SilverLife(e, iCell)
		fmt.Printf("   %-50s %7.1f V %7.2f V %8.1f W %11.1f %% %10.0f V %10.1e\n", r.Name, res.Vdc, res.Ripple, res.Pload, 100*res.Efficiency, res.ReversePeak, res.Mismatch/res.EnergyIn)
		fmt.Printf("      the silver is the positive electrode when the cell conducts (%.2f A on average): %.3g g of silver passes %.3g C,\n", iCell, m*1e3, q)
		fmt.Printf("      so if it all dissolved it would last %s (worst case; then the %s underneath carries the current and corrodes)\n", duration(sec), e.Base)
	}
	fmt.Println("   Efficiency = load power ÷ source power over the last cycles (Pollak's cells: about 75 %, Nodon's 75–80 %, EB 1911).")
	fmt.Println("   A half-wave rectifier with a filter capacitor sees twice the peak in reverse, so one cell sparks over; two in series")
	fmt.Printf("   hold. The oxide film is a %.2g µF capacitor: at high frequency it passes current both ways (the 1911 article's\n", c.FilmC()*1e6)
	fmt.Println("   warning); compare the 1 kHz and 20 kHz rows.")

	fmt.Println("\nPart 3 — polarity spark gaps in water as rectifiers (gap.go, gapsim.go)")
	g := DefaultGap()
	fmt.Println("   " + g.Describe())
	fmt.Println("   A one-way spark needs: forward voltage above +Vb, reverse below −Vb·ratio. With a filter capacitor a half-wave")
	fmt.Println("   rectifier sees 2·Vpeak in reverse, and a bridge needs 2·Vb forward (two gaps in series) while each blocked gap sees")
	fmt.Println("   the full output: both need a polarity ratio above 2.")
	fmt.Printf("   %-52s %6s %9s %9s %10s %12s %14s\n", "circuit", "ratio", "Vdc", "ripple", "efficiency", "reverse pk", "strikes f/b")
	var gapRuns []RectResult
	for _, cse := range []struct {
		name  string
		ratio float64
		r     Rectifier
	}{
		{"half-wave, resistive load, 4 kV peak", 2, Rectifier{Vp: 4000}},
		{"half-wave, filter capacitor, 4 kV peak", 2, Rectifier{Vp: 4000, Cload: 0.1e-6}},
		{"bridge, filter capacitor, 4 kV peak", 2, Rectifier{Bridge: true, Vp: 4000, Cload: 0.1e-6}},
		{"bridge, filter capacitor, 6 kV peak", 2, Rectifier{Bridge: true, Vp: 6000, Cload: 0.1e-6}},
		{"bridge, filter capacitor, 6 kV peak", 3, Rectifier{Bridge: true, Vp: 6000, Cload: 0.1e-6}},
	} {
		gg := g
		gg.NegRatio = cse.ratio
		r := cse.r
		r.F, r.Rs, r.Rload, r.Cycles, r.Skip = 60, 1000, 1e6, 20, 10
		res := r.RunGaps(gg)
		gapRuns = append(gapRuns, res)
		verdict := "rectifies"
		switch {
		case res.Forward == 0:
			verdict = "never strikes"
		case res.Backward > 0:
			verdict = "breaks down backwards"
		}
		fmt.Printf("   %-52s %6.1f %7.0f V %7.0f V %9.1f %% %10.0f V %7d/%-6d %s\n", cse.name, cse.ratio, res.Vdc, res.Ripple, 100*res.Efficiency, res.ReversePeak, res.Forward, res.Backward, verdict)
	}
	fmt.Printf("   A spark gap is not a low-drop diode: it conducts only once the voltage passes its breakdown (%.1f kV here), so the\n", g.Vpos()/1e3)
	fmt.Println("   output sags until the source exceeds the load voltage by the gaps' breakdown, then jumps back up.")
	fmt.Println("   The ratio of 2 is declared (J. C. Martin's constants, to be checked) and the breakdown voltages are extrapolated")
	fmt.Println("   to 60 Hz from microsecond data: a bench measurement of both decides whether a water spark-gap bridge can work.")
	fmt.Println("\nPart 4 — pulse test (pulse.go): 10 J from a charged capacitor through the device, 1 µH and a 0.1 Ω load")
	pulses, err := Pulses()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   %-46s %8s %10s %10s %10s %10s %10s %10s\n", "device", "charge", "C", "peak", "rise", "to load", "in device", "left in C")
	for _, p := range pulses {
		fmt.Printf("   %-46s %6.0f V %8.3g F %8.0f A %7.3g µs %8.2f J %8.2f J %8.2f J\n", p.Name, p.V0, p.C, p.Peak, p.Rise*1e6, p.ToLoad, p.InDevice, p.Left)
	}
	fmt.Println("   Spark gaps give the highest and fastest pulses (kiloamps in a fraction of a microsecond) because they run at")
	fmt.Println("   kilovolts; the water cells at 100 V give hundreds of amps in microseconds, limited by the water's resistance.")
	fmt.Println("   A gap's arc goes out at the first current zero, so energy left ringing in the capacitor stays there (it rings")
	fmt.Println("   again once the gap re-strikes). The water gap's 40 kV is Stygar's relation at 1 µs, inside its measured range.")

	fmt.Println("\nPart 5 — frequency sweep (sweep.go): bridge, source 0.1 Ω, load 1 Ω, no filter; amps per volt = mean load current ÷ peak")
	sweep, err := Sweep()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   perfect diode bridge: %.3f A/V at every frequency\n", IdealAperV())
	fmt.Printf("   %-62s", "device (peak volts)")
	for _, f := range SweepFreqs {
		fmt.Printf(" %8s", hz(f))
	}
	fmt.Println()
	for _, d := range sweep {
		fmt.Printf("   %-62s", d.Name)
		for _, r := range d.Rows {
			fmt.Printf(" %8.3f", r.AperV)
		}
		fmt.Println("  A/V")
		fmt.Printf("   %-62s", "      efficiency")
		for _, r := range d.Rows {
			fmt.Printf(" %7.0f%%", 100*r.Efficiency)
		}
		fmt.Println()
	}
	fmt.Println("   Cells: steady to about 3 kHz; above it the oxide film (a capacitor) passes current both ways, so the mean")
	fmt.Println("   current holds up but the efficiency collapses (34 % at 100 kHz). Silver tips on 10 % of the plate carry less.")
	fmt.Println("   Water gaps: at a fixed 6 kV the bridge works to 100 Hz and then never strikes, because breakdown rises with")
	fmt.Println("   frequency (E·τ^⅓, extrapolated). Raising the voltage to the middle of the working window keeps it striking, but")
	fmt.Println("   above about 100 Hz the gaps that just carried an arc have not recovered (1 ms recovery, declared; never below")
	fmt.Println("   water vapour's 480 V Paschen minimum) and break down backwards: the bridge shorts. Spark gaps are pulse")
	fmt.Println("   switches, not high-frequency rectifiers.")

	fmt.Println("\nThe water and the metals add no energy: every watt at the load comes from the AC source, and the rest is heat")
	fmt.Println("(and, for the cells, some electrolysis).")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "iv.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprint(f, "v")
	for _, p := range pairs {
		fmt.Fprintf(f, ",%s", p.Name())
	}
	fmt.Fprintln(f)
	for v := -200.0; v <= 200; v += 1 {
		fmt.Fprintf(f, "%g", v)
		for _, p := range pairs {
			fmt.Fprintf(f, ",%g", p.Current(v))
		}
		fmt.Fprintln(f)
	}
	f.Close()
	write := func(name string, r RectResult, per float64) {
		f, err := os.Create(filepath.Join(*viz, name))
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		fmt.Fprintln(f, "t_ms,vout")
		step := len(r.Time)/3000 + 1
		for i := 0; i < len(r.Time); i += step {
			fmt.Fprintf(f, "%g,%g\n", r.Time[i]*1e3, r.Vout[i])
		}
	}
	write("bridge_cells.csv", runs[2].res, 1/60.0)
	write("halfwave_cells.csv", runs[0].res, 1/60.0)
	write("gap_resistive.csv", gapRuns[0], 1/60.0)
	write("gap_bridge_ratio3.csv", gapRuns[4], 1/60.0)
	sf, err := os.Create(filepath.Join(*viz, "frequency_trend.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(sf, "device,f_hz,vpeak,a_per_v,efficiency")
	for _, d := range sweep {
		for _, r := range d.Rows {
			fmt.Fprintf(sf, "%q,%g,%g,%g,%g\n", d.Name, r.F, r.Vp, r.AperV, r.Efficiency)
		}
	}
	sf.Close()
	pf, err := os.Create(filepath.Join(*viz, "pulses.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(pf, "device,t_us,amps")
	for _, p := range pulses {
		n := 0
		for n < len(p.Time) && p.Time[n] <= 50e-6 {
			n++
		}
		step := n/3000 + 1
		for i := 0; i < n; i += step {
			fmt.Fprintf(pf, "%q,%g,%g\n", p.Name, p.Time[i]*1e6, p.Amps[i])
		}
	}
	pf.Close()
}

// hz prints a frequency briefly.
func hz(f float64) string {
	if f >= 1e3 {
		return fmt.Sprintf("%gk", f/1e3)
	}
	return fmt.Sprintf("%g", f)
}
