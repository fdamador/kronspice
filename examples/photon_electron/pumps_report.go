// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
)

// part4 runs and prints the three electron pumps; with viz set it writes
// their curves for viz/render.py.
func part4(viz string) error {
	fmt.Println("\nPart 4 — three electron pumps that work (pumps.go)")
	fmt.Printf("   The model's own frequency, m_e c²/h = %.3g Hz, is about 10⁸ times beyond any circuit,\n", Me*C0*C0/(2*math.Pi*Hbar))
	fmt.Println("   so these pumps run at frequencies circuits reach. Electrons come from metal or ground; energy from a declared source.")

	// 1
	p := DefaultSEPump()
	fmt.Println("\n   1. Single-electron pump (Pothier et al. 1992; orthodox tunnelling theory)")
	fmt.Printf("      3 junctions of %.0f aF, %.0f kΩ; gates %.0f aF at 90°; %.0f mK; charging energy e²/2C = %.3g meV (%.2g K)\n",
		p.C[0]*1e18, p.R[0]/1e3, p.Cg[0]*1e18, p.T*1e3, E/(2*(p.C[0]+p.C[1]+p.Cg[0]))*1e3, E*E/(2*(p.C[0]+p.C[1]+p.Cg[0]))/kB)
	fmt.Printf("      %10s %22s %16s %26s %14s\n", "drive", "electrons per cycle", "current", "heat per electron (gates)", "ledger")
	var seRows []SEResult
	for _, f := range []float64{1e6, 1e7, 3e7, 1e8, 3e8} {
		r := p.Pump(f, 2, 4)
		seRows = append(seRows, r)
		fmt.Printf("      %7.0f MHz %22.9f %13.4g pA %21.3g µeV %14.1e\n", f/1e6, r.PerCycle, math.Abs(r.Current)*1e12, r.Heat/E*1e6, r.Residual)
	}
	pb := p
	pb.V = 20e-6 // left lead higher: moving electrons left → right is uphill
	rb := pb.Pump(1e7, 2, 4)
	fmt.Printf("      uphill against 20 µV at 10 MHz: %.6f electrons per cycle; the gates supply %.1f µeV per cycle =\n", rb.PerCycle, rb.GateWork/E*1e6)
	fmt.Printf("      %.1f µeV lifted into the bias source + %.1f µeV heat (ledger %.1e)\n", -rb.BiasWork/E*1e6, rb.Heat/E*1e6, rb.Residual)
	fmt.Println("      Exactly one electron per cycle at low drive (I = e·f); errors grow when the drive outruns the tunnelling (RC ≈ 40 ps).")

	// 2
	m := DefaultMultipactor()
	L, C, _, Is := m.Tank()
	fmt.Println("\n   2. Multipactor pump in an RF resonator (Farnsworth US 2,091,439; Kishek & Lau 1998)")
	fmt.Printf("      %.0f mm copper gap, %.0f cm², resonator %.0f GHz (C = %.3g pF, L = %.3g nH, Q = %.0f); drive %.3g mA for %.0f V unloaded\n",
		m.Gap*1e3, m.Area*1e4, m.F/1e9, C*1e12, L*1e9, m.Q, Is*1e3, m.V0)
	tw, err := m.TankTwin(5)
	if err != nil {
		return err
	}
	fmt.Printf("      resonator alone, own integrator vs kronspice transient: largest difference %.2g of V0\n", tw)
	r := m.Run(2000, true)
	off := m
	off.V0 = 300
	ro := off.Run(40, true)
	peak := 0.0
	for _, s := range r.Samples {
		peak = math.Max(peak, s.Electrons)
	}
	fmt.Printf("      in step (%.0f V, seed %.0f electrons): grows to %.3g electrons in the gap, then pulls the gap voltage down to %.0f V,\n",
		m.V0, m.Seed, peak, minGap(r.Samples))
	fmt.Printf("      stops multiplying; the resonator recovers and a stray electron (one per %d ns while the gap is empty, declared) restarts it:\n", m.SeedEvery)
	fmt.Printf("      bursts, as Kishek & Lau found (%.0f strays used)\n", r.Strays)
	fmt.Printf("      out of step (300 V): %.3g electrons left after 40 periods — no pumping.\n", ro.Samples[len(ro.Samples)-1].Electrons)
	fmt.Printf("      over %d periods (%.0f ns): %.3g electrons collected, average collector current %.3g µA (second 1000 ns)\n", len(r.Samples), float64(len(r.Samples))/m.F*1e9, r.Collected, r.Current*1e6)
	fmt.Printf("      ledger [nJ]: drive %.4g = resonator heat %.4g + impacts on the plates %.4g + collector %.4g + stored change %.4g; residual %.1e\n",
		r.SourceWork*1e9, r.TankHeat*1e9, r.WallHeat*1e9, r.Collector*1e9, r.Stored*1e9, r.Residual)
	fmt.Println("      (the stored change is negative: the bursts also drain the energy the resonator held at the start)")
	fmt.Printf("      energy per collected electron: %.0f eV from the RF drive\n", r.SourceWork/r.Collected/E)

	// 3
	g := DefaultGroundPump()
	gr, err := g.Run()
	if err != nil {
		return err
	}
	fmt.Println("\n   3. Ground charge pump with a resonant drive (§3j; Cockcroft & Walton 1932)")
	fmt.Printf("      %.0f V at %.0f kHz into a series resonator (Q = %.0f) → %d-stage negative ladder → %.0f nF bank, %.0f MΩ load; %d cycles\n",
		g.Vs, g.F/1e3, g.Q, g.Stages, g.Cbank*1e9, g.Rload/1e6, g.Cycles)
	fmt.Printf("      resonator step-up: %.0f V → %.1f V peak (×%.1f); bank %.1f V against the ideal −2·n·peak = %.1f V\n", g.Vs, gr.TankPeak, gr.TankPeak/g.Vs, gr.Bank, gr.Ideal)
	fmt.Printf("      electrons on the bank's negative plate %.3g (plus %.3g through the load), lent by the earth; the bank's earthed plate\n",
		-gr.BankCharge/E, -gr.LoadCharge/E)
	fmt.Println("      returns the same number, so the earth's net charge is unchanged.")
	fmt.Printf("      ledger: source %.4g J = heat %.4g J + stored %.4g J (bank %.3g J) + time-stepping error %.1e J\n", gr.EnergyIn, gr.Heat, gr.Stored, gr.BankEnergy, gr.Mismatch)

	if viz == "" {
		return nil
	}
	f, err := os.Create(filepath.Join(viz, "pump_single.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "f_hz,per_cycle")
	for _, s := range seRows {
		fmt.Fprintf(f, "%g,%.12g\n", s.F, s.PerCycle)
	}
	f.Close()
	f, err = os.Create(filepath.Join(viz, "pump_multipactor.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "period,electrons,gap_v")
	for _, s := range r.Samples {
		fmt.Fprintf(f, "%d,%g,%g\n", s.Period, s.Electrons, s.GapV)
	}
	f.Close()
	f, err = os.Create(filepath.Join(viz, "pump_ground.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "time_ms,bank_v")
	step := len(gr.Time)/2000 + 1
	for i := 0; i < len(gr.Time); i += step {
		fmt.Fprintf(f, "%g,%g\n", gr.Time[i]*1e3, gr.BankV[i])
	}
	return f.Close()
}

func minGap(s []MPSample) float64 {
	m := math.Inf(1)
	for _, x := range s {
		m = math.Min(m, x.GapV)
	}
	return m
}

// part5 evolves the magnetron-driven multipactor pump and prints it.
func part5(viz string) error {
	fmt.Println("\nPart 5 — the multipactor pump driven by a standard microwave-oven magnetron, evolved (magnetron.go)")
	mg := LG2M214
	fmt.Printf("   source: %s, %.0f MHz, %.0f W RF; declared efficiency %.0f %% + %.1f W heater → %.0f W from the mains\n",
		mg.Name, mg.F/1e6, mg.Pout, 100*mg.Efficiency, mg.Filament, mg.Mains())
	fmt.Println("   chain: magnetron → circulator → coupling port → copper re-entrant cavity (Q0 = 2000) whose gap multipacts → grid → collector;")
	fmt.Println("   power the cavity does not take goes to the circulator's water load. Genes: gap, gap area (side ≤ λ/10), coupling β,")
	fmt.Println("   grid transparency, surface. Fitness: mains watts per collected µA. Limits: Kilpatrick field, gap ≤ 10 kV (non-relativistic).")
	o := DefaultEvolve()
	ev, err := Evolve(o)
	if err != nil {
		return err
	}
	fmt.Printf("   %d designs × %d generations (seed %d); search at ≥ %d steps per period over %d periods, results re-run at 800 steps, 2000 periods\n",
		o.Population, o.Generations, o.Seed, o.Steps, o.Periods)
	fmt.Printf("   %-30s %16s %16s\n", "", "best of gen. 0", "evolved")
	row := func(name, f string, a, b float64) {
		fmt.Printf("   %-30s %16s %16s\n", name, fmt.Sprintf(f, a), fmt.Sprintf(f, b))
	}
	a, b := ev.First, ev.Best
	row("gap [mm]", "%.3f", a.Design.Gap*1e3, b.Design.Gap*1e3)
	row("gap area [cm²]", "%.3f", a.Design.Area*1e4, b.Design.Area*1e4)
	row("coupling β", "%.3g", a.Design.Beta, b.Design.Beta)
	row("grid transparency", "%.2f", a.Design.Grid, b.Design.Grid)
	fmt.Printf("   %-30s %16s %16s\n", "surface", a.Design.Surface.Name, b.Design.Surface.Name)
	row("gap amplitude, unloaded [V]", "%.0f", a.V0, b.V0)
	row("peak field / Kilpatrick", "%.3f", a.PeakField/a.Kilpatrick, b.PeakField/b.Kilpatrick)
	row("collector current [µA]", "%.4g", a.Current*1e6, b.Current*1e6)
	row("mains W per µA", "%.4g", a.Mains/(a.Current*1e6), b.Mains/(b.Current*1e6))
	row("electrons per mains joule", "%.3g", a.ElectronsPerJoule, b.ElectronsPerJoule)
	row("into the cavity [W]", "%.3g", a.Delivered, b.Delivered)
	row("  cavity wall heat [W]", "%.3g", a.CavityHeat, b.CavityHeat)
	row("  electrons: plates + collector [W]", "%.3g", a.Electrons, b.Electrons)
	row("to the water load [W]", "%.1f", a.Dummy, b.Dummy)
	row("ledger residual", "%.1e", a.Residual, b.Residual)
	fmt.Printf("   The collector current varies by about ±10 %% with step size and run length (burst timing); the order of magnitude is stable.\n")
	fmt.Printf("   The gap needs only %.2g W, so with weak coupling %.1f W of the 900 W goes to the water load: one magnetron could feed\n", b.Delivered, b.Dummy)
	fmt.Printf("   about %.0f such gaps through a power divider (not modelled). Per watt into the cavity: %.3g µA.\n", mg.Pout/b.Delivered, b.Current*1e6/b.Delivered)
	if viz == "" {
		return nil
	}
	f, err := os.Create(filepath.Join(viz, "magnetron_history.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "generation,best_w_per_ua")
	for _, h := range ev.History {
		fmt.Fprintf(f, "%d,%g\n", h.Generation, h.Best)
	}
	f.Close()
	f, err = os.Create(filepath.Join(viz, "magnetron_best.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "ns,electrons,gap_v")
	for _, s := range b.Run.Samples {
		fmt.Fprintf(f, "%g,%g,%g\n", float64(s.Period)/b.Design.Mag.F*1e9, s.Electrons, s.GapV)
	}
	return f.Close()
}

// part6 scales the magnetron-driven multipactor pump to many gaps.
func part6(viz string) error {
	fmt.Println("\nPart 6 — scaling up: many multipactor gaps sharing one magnetron (magnetron.go, gene \"gaps\")")
	fmt.Println("   N identical gaps take Pout/N each through an ideal isolated power divider (declared lossless); the genetic")
	fmt.Println("   algorithm now also picks N (1 to 1,000,000). Electrons multiply in proportion to the copper area in the field, so")
	fmt.Println("   current scales with total gap area and with microwave power; the fitness is still mains watts per collected µA.")
	type res struct {
		m  Magnetron
		ev Evolution
	}
	var out []res
	for _, m := range []Magnetron{LG2M214, StellantCWM100} {
		o := DefaultEvolve()
		o.Mag, o.MaxGaps, o.Population, o.Generations = m, 1000000, 40, 50
		ev, err := Evolve(o)
		if err != nil {
			return err
		}
		out = append(out, res{m, ev})
	}
	fmt.Printf("   %-30s %20s %26s\n", "", out[0].m.Name, out[1].m.Name)
	row := func(name string, f func(r MagResult) string) {
		fmt.Printf("   %-30s %20s %26s\n", name, f(out[0].ev.Best), f(out[1].ev.Best))
	}
	row("frequency, RF power", func(r MagResult) string {
		return fmt.Sprintf("%.0f MHz, %.0f kW", r.Design.Mag.F/1e6, r.Design.Mag.Pout/1e3)
	})
	row("gaps in parallel", func(r MagResult) string { return fmt.Sprintf("%d", r.N) })
	row("gap; area per gap", func(r MagResult) string {
		return fmt.Sprintf("%.2f mm; %.2f cm²", r.Design.Gap*1e3, r.Design.Area*1e4)
	})
	row("total gap area", func(r MagResult) string { return fmt.Sprintf("%.3g m²", float64(r.N)*r.Design.Area) })
	row("gap voltage (unloaded)", func(r MagResult) string { return fmt.Sprintf("%.0f V", r.V0) })
	row("collected current", func(r MagResult) string { return fmt.Sprintf("%.3g A", r.Current) })
	row("per gap", func(r MagResult) string { return fmt.Sprintf("%.3g µA", r.Current/float64(r.N)*1e6) })
	row("into the gaps / water load", func(r MagResult) string {
		return fmt.Sprintf("%.3g / %.3g kW", r.Delivered/1e3, r.Dummy/1e3)
	})
	row("mains", func(r MagResult) string { return fmt.Sprintf("%.3g kW", r.Mains/1e3) })
	row("amps per kW of mains", func(r MagResult) string { return fmt.Sprintf("%.3g", r.Current/(r.Mains/1e3)) })
	row("amps per m² of gap", func(r MagResult) string { return fmt.Sprintf("%.3g", r.Current/(float64(r.N)*r.Design.Area)) })
	row("ledger residual (one gap)", func(r MagResult) string { return fmt.Sprintf("%.1e", r.Residual) })
	best := out[0].ev.Best
	if a, b := out[1].ev.Best, out[0].ev.Best; a.Current/a.Mains > b.Current/b.Mains {
		best = a
	}
	aPerW, aPerM2 := best.Current/best.Mains, best.Current/(float64(best.N)*best.Design.Area)
	fmt.Println("   Scaling the better of the two designs (same gaps, more of them, more magnetrons):")
	for _, target := range []float64{1, 10, 100, 300} {
		fmt.Printf("      %5.0f A needs about %8.3g kW from the mains and %8.3g m² of gap plates\n", target, target/aPerW/1e3, target/aPerM2)
	}
	fmt.Println("   The searches were still improving when they stopped, so these are lower bounds on what the design space allows.")
	fmt.Println("   Feeding thousands of gaps needs a real power-distribution structure (a large multi-gap cavity rather than a divider),")
	fmt.Println("   which is not modelled. The electrons come from the copper; every watt comes from the magnetron.")
	if viz == "" {
		return nil
	}
	f, err := os.Create(filepath.Join(viz, "array_history.csv"))
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "generation,lg_w_per_ua,stellant_w_per_ua")
	for i := range out[0].ev.History {
		s := out[1].ev.History[min(i, len(out[1].ev.History)-1)]
		fmt.Fprintf(f, "%d,%g,%g\n", i, out[0].ev.History[i].Best, s.Best)
	}
	return f.Close()
}
