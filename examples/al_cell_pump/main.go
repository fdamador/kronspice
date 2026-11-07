// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command al_cell_pump: the resonant-driven negative Cockcroft–Walton
// ground pump of examples/photon_electron with its rectifiers made of
// aluminium-water cells (the valve cell of examples/sparkgap_diode, built
// from plates: cells in series for voltage, plates in parallel for
// current), evolved by NSGA-II for "catalytic efficiency" — how close a
// design comes to spending only the e·|V| that lifting each electron
// costs — together with output voltage and current. Register entry H38.
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
)

func main() {
	pop := flag.Int("pop", 40, "evolution population")
	gens := flag.Int("gens", 25, "generations")
	seed := flag.Uint64("seed", 1, "seed")
	data := flag.String("data", "", "write figure data (CSV) to this directory")
	study := flag.String("study", "all", "plates, pipes or all")
	flag.Parse()
	var err error
	if *study != "pipes" {
		err = run(*pop, *gens, *seed, *data)
	}
	if err == nil && *study != "plates" {
		err = runPipes(*pop, *gens, *seed, *data)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "al_cell_pump:", err)
		os.Exit(1)
	}
}

func line(name string, p Pump, r Result) string {
	return fmt.Sprintf("  %-34s %8.0f V %9.3g mA %9.3g W %8.4f | coil %5.1f %%  forward %5.1f %%  leakage %5.1f %%  rest %5.1f %%",
		name, r.Out, r.Iout*1e3, r.Pin, r.Eff,
		100*r.LossCoil/r.Pin, 100*r.LossForward/r.Pin, 100*r.LossLeak/r.Pin, 100*r.Rest/r.Pin)
}

func run(pop, gens int, seed uint64, data string) error {
	fmt.Println("H38 — aluminium-cell ground pump, evolved for catalytic efficiency")
	fmt.Println("  catalytic efficiency η = power into the load ÷ power from the source, in steady state: lifting an")
	fmt.Println("  electron through |V| costs e·|V| whatever the circuit, so η = 1 means nothing more is spent")
	fmt.Println()
	base := DefaultPump()
	fmt.Println("  baseline (photon_electron's ground pump: 10 V at 100 kHz, 1 µF tank, 3 stages, 10 MΩ load; coil Q 50)")
	fmt.Println("                                         bank       load       input    η       where the input goes")
	var baseRuns []Result
	for _, c := range []struct {
		name  string
		diode bool
	}{{"silicon diodes", true}, {"aluminium stacks (2 × 1 × 1000 mm²)", false}} {
		p := base
		p.Diode = c.diode
		r, err := p.Run(analysis.Options{})
		if err != nil && err != ErrNoPump {
			return err
		}
		baseRuns = append(baseRuns, r)
		fmt.Println(line(c.name, p, r))
	}
	fmt.Println("  the 1 µF tank circulates hundreds of amps to hold its swing: almost all the drive heats the coil;")
	fmt.Printf("  with aluminium stacks at 100 kHz each oxide film (%.3g µF) all but shorts its rectifier\n", base.Cell.C0()*1e6)
	fmt.Println("  (the high-frequency failure the 1911 Britannica notes).")
	fmt.Println()

	fmt.Printf("  evolving (NSGA-II, %d × %d, seed %d): 11 genes — drive frequency and amplitude, tank capacitor, stages,\n", pop, gens, seed)
	fmt.Println("  ladder capacitors, load, and the stacks' cells in series, plates in parallel, plate area, gap, forming voltage")
	ev, err := EvolvePump(pop, gens, seed, func(s string) { fmt.Println("   ", s) })
	if err != nil {
		return err
	}
	fmt.Printf("  Pareto designs (%d), by efficiency; each re-measured with silicon diodes and with the 2 %% varactor film:\n", len(ev))
	fmt.Println("     η      bank V   load mA   input W  | drive         tank     stages ladder   load    | stack                                          | η Si    η κ=2%")
	for _, e := range ev {
		p, r := e.Pump, e.R
		si := p
		si.Diode = true
		rs, _ := si.Run(analysis.Options{})
		kv := p
		kv.Cell.Kappa = kappaDef
		rk, _ := kv.Run(analysis.Options{})
		fmt.Printf("  %6.3f %8.0f %9.3g %9.3g | %5.0f Hz %4.1f V %6.3g F  %d  %7.3g F %7.3g Ω | %-46s | %.3f   %.3f\n",
			r.Eff, r.Out, r.Iout*1e3, r.Pin, p.F, p.Vs, p.Ct, p.Stages, p.Cstage, p.Rload, p.Cell.Short()+fmt.Sprintf(", %.0f V", p.Cell.Vform), rs.Eff, rk.Eff)
	}
	if len(ev) > 0 {
		b := ev[0]
		fmt.Println()
		fmt.Println("  the most efficient design, where its input goes:")
		fmt.Println(line("evolved", b.Pump, b.R))
		fmt.Printf("    %s; tank swing %.0f V, largest reverse %.0f %% of the oxide's\n", b.Pump.Cell.Describe(), b.R.TankPeak, 100*b.R.RevRatio)
		fmt.Printf("    electrons delivered per joule: %.3g (the ideal for %.0f V is %.3g)\n", b.R.ElPerJ, b.R.Out, b.R.Ideal)
		fw := float64(b.Pump.Cell.Ns) * vdCell * 2 * float64(b.Pump.Stages)
		fmt.Printf("    floor from the electrolysis onset alone: each electron crosses 2n = %d stacks of %d cells, %.0f V of\n",
			2*b.Pump.Stages, b.Pump.Cell.Ns, fw)
		fmt.Printf("    %.1f V each, so η ≤ |V|/(|V| + %.0f V) = %.3f\n", vdCell, fw, math.Abs(b.R.Out)/(math.Abs(b.R.Out)+fw))
	}
	fmt.Println()
	fmt.Println("  the parametric route (pumping the tank at 2f through the oxide film's C(v)):")
	fmt.Printf("    a capacitance swinging by m = ΔC/C at 2f grows oscillation at f only if m > 2/Q (Mathieu threshold);\n")
	fmt.Printf("    the film's κ = %.0f %% gives m ≤ %.0f %% even swung to its forming voltage, so it needs Q > %.0f; the\n",
		100*kappaDef, 100*kappaDef, 2/kappaDef)
	fmt.Printf("    coil's Q of %.0f needs m > %.0f %%: the oxide is too linear to pump parametrically (a varactor diode,\n", base.QL, 200/base.QL)
	fmt.Println("    whose C changes several-fold, could)")
	fmt.Println()
	if data != "" {
		if err := writeData(data, baseRuns, ev); err != nil {
			return err
		}
	}
	return nil
}

func writeData(dir string, base []Result, ev []Evolved) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, "pareto.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "eff,out_V,iout_A,pin_W,coil,forward,leak,rest,F,Vs,stages,ns,np,area,gap,vform")
	for _, e := range ev {
		r, p := e.R, e.Pump
		fmt.Fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d,%d,%g,%g,%g\n", r.Eff, r.Out, r.Iout, r.Pin,
			r.LossCoil, r.LossForward, r.LossLeak, r.Rest, p.F, p.Vs, p.Stages, p.Cell.Ns, p.Cell.Np, p.Cell.Area, p.Cell.Gap, p.Cell.Vform)
	}
	g, err := os.Create(filepath.Join(dir, "baseline.csv"))
	if err != nil {
		return err
	}
	defer g.Close()
	fmt.Fprintln(g, "name,eff,out_V,pin_W,coil,forward,leak,rest,pout")
	for i, r := range base {
		name := []string{"silicon diodes", "aluminium stacks"}[i]
		fmt.Fprintf(g, "%s,%g,%g,%g,%g,%g,%g,%g,%g\n", name, r.Eff, r.Out, r.Pin, r.LossCoil, r.LossForward, r.LossLeak, r.Rest, r.Pout)
	}
	if len(ev) > 0 {
		h, err := os.Create(filepath.Join(dir, "best_bank.csv"))
		if err != nil {
			return err
		}
		defer h.Close()
		fmt.Fprintln(h, "t,bank")
		r := ev[0].R
		step := len(r.Time)/2000 + 1
		for i := 0; i < len(r.Time); i += step {
			fmt.Fprintf(h, "%g,%g\n", r.Time[i], r.Bank[i])
		}
	}
	return nil
}
