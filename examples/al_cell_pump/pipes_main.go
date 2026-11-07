// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
)

// Seal study constants: a pipe sealed at the bottom (its top stands above
// the water line, so only the bottom can leak), the seal's width along
// the leak path, and the clearance of a pipe merely standing on the floor.
const (
	sealWidth = 5e-3   // m
	looseGap  = 0.5e-3 // m
)

// filmFor is the water-film thickness [m] that leaves seal resistance R
// under a seal of width sealWidth around a pipe of radius r.
func filmFor(R, r float64) float64 { return rhoWater * sealWidth / (2 * math.Pi * r * R) }

// SealPoint is one seal resistance and the pump's efficiency with it.
type SealPoint struct {
	Rseal, Eff, Seal, Out, Rev float64
}

// SealSweep runs a series-pipe design with leak resistances Rs.
func SealSweep(p Pump, rs []float64) []SealPoint {
	var out []SealPoint
	for _, R := range rs {
		q := p
		pp := *p.Pipes
		pp.Rseal = R
		q.Pipes = &pp
		r, err := q.Run(analysis.Options{})
		if err != nil && err != ErrNoPump {
			out = append(out, SealPoint{R, 0, 0, 0, 0})
			continue
		}
		out = append(out, SealPoint{R, r.Eff, r.LossSeal / math.Max(r.Pin, 1e-300), r.Out, r.RevRatio})
	}
	return out
}

func runPipes(pop, gens int, seed uint64, data string) error {
	fmt.Println("H38b — one tank per rectifier: concentric pipes")
	fmt.Println("  rod at the centre, pipes around it, the tank wall outside; each annulus of water is one cell.")
	fmt.Println("  parallel: pipes alternate aluminium and iron, all aluminium on one bus bar and all iron on the other —")
	fmt.Println("            every annulus is the same valve between the buses; ends need no seals (neighbouring annuli")
	fmt.Println("            sit at the same voltage); the tank holds one forming voltage")
	fmt.Println("  series:   clad pipes (aluminium outside, iron inside) carry the current across every annulus in turn;")
	fmt.Println("            the tank holds N forming voltages, but each pipe's ends must be sealed (water around an end")
	fmt.Println("            shorts a cell); evolved with perfect seals, then swept. Plain, the inner annulus (least")
	fmt.Println("            aluminium, least film) takes most of the reverse voltage; graded, the outer pipes' aluminium")
	fmt.Println("            is masked to give every annulus the same area, and they share it equally")
	fmt.Println()
	fmt.Printf("  evolving (NSGA-II, %d × %d, seed %d): drive, tank, stages, ladder, load, and the tank's wiring, annuli,\n", pop, gens, seed)
	fmt.Printf("  rod radius, water gap, wetted length, forming voltage (pipe walls %.1f mm)\n", wallDef*1e3)
	ev, err := EvolvePipes(pop, gens, seed, func(s string) { fmt.Println("   ", s) })
	if err != nil {
		return err
	}
	var ampsRows []string
	fmt.Printf("  Pareto designs (%d), by efficiency:\n", len(ev))
	fmt.Println("     η      bank V   load mA   input W  | drive          stages | tank (one per rectifier)                                              | worst oxide")
	var bestOne, bestPar, bestSer *Evolved
	for i := range ev {
		e := &ev[i]
		p, r := e.Pump, e.R
		fmt.Printf("  %6.3f %8.0f %9.3g %9.3g | %5.0f Hz %4.1f V  %d   | %-69s | %3.0f %%\n",
			r.Eff, r.Out, r.Iout*1e3, r.Pin, p.F, p.Vs, p.Stages, p.Pipes.Short()+fmt.Sprintf(", %.0f V", p.Pipes.Vform), 100*r.RevRatio)
		switch {
		case p.Pipes.N == 1:
			if bestOne == nil {
				bestOne = e
			}
		case p.Pipes.Series:
			if bestSer == nil {
				bestSer = e
			}
		default:
			if bestPar == nil {
				bestPar = e
			}
		}
	}
	count := map[string]int{}
	for _, e := range ev {
		count[e.Pump.Pipes.kind()]++
	}
	fmt.Printf("  on the front: %d single-cell tanks (one annulus: the wiring makes no difference), %d parallel, %d series\n",
		count["single cell"], count["parallel"], count["series"])
	fmt.Println()
	for _, b := range []*Evolved{bestOne, bestPar, bestSer} {
		if b == nil {
			continue
		}
		fmt.Printf("  most efficient %s tank:\n", b.Pump.Pipes.kind())
		fmt.Println(line("evolved", b.Pump, b.R))
		fmt.Printf("    %s\n", b.Pump.Pipes.Describe())
		fmt.Printf("    %d tanks (2 per stage); tank swing %.0f V; %.3g electrons per joule (ideal %.3g)\n",
			2*b.Pump.Stages, b.R.TankPeak, b.R.ElPerJ, b.R.Ideal)
		if b.Pump.Pipes.Series && b.Pump.Pipes.N > 1 {
			fmt.Print("    largest reverse on each annulus, inner to outer, % of its forming voltage:")
			for _, x := range b.R.CellRev {
				fmt.Printf(" %.0f", 100*x)
			}
			fmt.Println()
			if !b.Pump.Pipes.Graded {
				fmt.Println("    (the inner annuli have the least aluminium, so the least film capacitance: they take the most)")
			}
		}
	}
	if bestOne != nil || bestPar != nil {
		b := bestOne
		if b == nil {
			b = bestPar
		}
		fmt.Println()
		fmt.Println("  parallel annuli for current: the most efficient single cell rewired as a parallel tank of N annuli")
		fmt.Println("  (same rod, gap, length and forming voltage), at its own load and at a tenth and a hundredth of it")
		fmt.Println("      load      N   tank Ø     aluminium      η      bank V   load mA   forward  coil")
		for _, div := range []float64{1, 10, 100} {
			for _, n := range []int{1, 2, 4, 6} {
				p := b.Pump
				q := *p.Pipes
				q.Series, q.N = false, n
				p.Pipes = &q
				p.Rload /= div
				r, err := p.Run(analysis.Options{})
				if err != nil && err != ErrNoPump {
					return err
				}
				fmt.Printf("    %8.3g Ω  %d  %5.0f mm  %8.3g mm²  %6.3f %8.0f %9.3g   %5.1f %%  %4.1f %%\n", p.Rload, n, 2e3*q.TankRadius(), q.Area()*1e6,
					r.Eff, r.Out, r.Iout*1e3, 100*r.LossForward/math.Max(r.Pin, 1e-300), 100*r.LossCoil/math.Max(r.Pin, 1e-300))
				ampsRows = append(ampsRows, fmt.Sprintf("%g,%d,%g,%g,%g,%g", p.Rload, n, q.Area(), r.Eff, r.Out, r.Iout))
			}
		}
	}
	// Seals matter only for a series tank whose reverse voltage needs more
	// than one annulus (N × the worst annulus's share > 1): with less, a
	// leak that shorts the outer annuli leaves one annulus able to hold it.
	var needs *Evolved
	for i := range ev {
		q := ev[i].Pump.Pipes
		if q.Series && float64(q.N)*ev[i].R.RevRatio > 1 {
			needs = &ev[i]
			break
		}
	}
	if bestSer != nil {
		fmt.Println()
		fmt.Printf("  the most efficient series tank's annuli take %.0f %% each: %.2f forming voltages in all, so one annulus\n",
			100*bestSer.R.RevRatio, float64(bestSer.Pump.Pipes.N)*bestSer.R.RevRatio)
		if float64(bestSer.Pump.Pipes.N)*bestSer.R.RevRatio <= 1 {
			fmt.Println("  could hold it alone: its seals do not matter (a leak shorting the outer annuli even saves their forward drops)")
		}
	}
	bestSer = needs
	var sweep []SealPoint
	if bestSer != nil {
		fmt.Println()
		p := bestSer.Pump
		rs := []float64{}
		for e := 0.0; e <= 10; e += 0.5 {
			rs = append(rs, math.Pow(10, e))
		}
		sweep = SealSweep(p, rs)
		a0, _ := p.Pipes.radii(0)
		rMid := a0 + float64(p.Pipes.N/2)*(p.Pipes.Gap+p.Pipes.Wall)
		rLoose := rhoWater * p.Pipes.Wall / (2 * math.Pi * rMid * looseGap)
		fmt.Printf("  seals on the most efficient series tank that needs its annuli (%s, %.0f V;\n", bestSer.Pump.Pipes.Short(), bestSer.Pump.Pipes.Vform)
		fmt.Printf("  η %.3f at %.0f V with perfect seals, %.2f forming voltages across it): leak resistance between\n",
			bestSer.R.Eff, bestSer.R.Out, float64(bestSer.Pump.Pipes.N)*bestSer.R.RevRatio)
		fmt.Printf("  neighbouring annuli, and the water film it means under a %.0f mm-wide bottom seal around a %.0f mm pipe\n", sealWidth*1e3, rMid*1e3)
		fmt.Println("        Rseal       film         η      bank V   seal loss   worst oxide")
		for _, s := range sweep {
			fmt.Printf("    %9.3g Ω %9.3g m  %7.3f %8.0f   %5.1f %%     %4.0f %%\n", s.Rseal, filmFor(s.Rseal, rMid), s.Eff, s.Out, 100*s.Seal, 100*s.Rev)
		}
		// The seal is good enough from the lowest Rseal above which every
		// point keeps the efficiency and no oxide passes 90 % (a leak
		// shorts the outer annuli and loads the rest past their forming
		// voltage).
		need := func(frac float64) float64 {
			ok := math.Inf(1)
			for i := len(sweep) - 1; i >= 0; i-- {
				s := sweep[i]
				if s.Eff < frac*bestSer.R.Eff || s.Rev > 0.9 {
					break
				}
				ok = s.Rseal
			}
			return ok
		}
		r99, r90 := need(0.99), need(0.9)
		fmt.Printf("    to keep 99 %% of the sealed η: Rseal ≥ %.3g Ω (film ≤ %.2g m); 90 %%: ≥ %.3g Ω (film ≤ %.2g m)\n",
			r99, filmFor(r99, rMid), r90, filmFor(r90, rMid))
		loose := SealSweep(p, []float64{rLoose})[0]
		fmt.Printf("    pipes merely standing on the floor (%.1f mm clearance under a %.1f mm wall): Rseal ≈ %.3g Ω, η %.3f, worst oxide %.0f %%\n",
			looseGap*1e3, p.Pipes.Wall*1e3, rLoose, loose.Eff, 100*loose.Rev)
		fmt.Println("    the seal must be bonded (epoxy, potting): a gasket's wet film is microns thick")
	}
	fmt.Println()
	if data != "" {
		return writePipes(data, ev, bestSer, sweep, ampsRows)
	}
	return nil
}

func writePipes(dir string, ev []Evolved, ser *Evolved, sweep []SealPoint, amps []string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, "pipes_pareto.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "eff,out_V,iout_A,pin_W,coil,forward,leak,seal,rest,F,Vs,stages,series,graded,annuli,rod,gap,length,vform,rev")
	for _, e := range ev {
		r, p := e.R, e.Pump
		q := p.Pipes
		fmt.Fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%d,%d,%d,%d,%g,%g,%g,%g,%g\n", r.Eff, r.Out, r.Iout, r.Pin,
			r.LossCoil, r.LossForward, r.LossLeak, r.LossSeal, r.Rest, p.F, p.Vs, p.Stages, int(b2f(q.Series)), int(b2f(q.Graded)), q.N, q.R0, q.Gap, q.Len, q.Vform, r.RevRatio)
	}
	a, err := os.Create(filepath.Join(dir, "pipes_amps.csv"))
	if err != nil {
		return err
	}
	defer a.Close()
	fmt.Fprintln(a, "rload,n,area_m2,eff,out_V,iout_A")
	for _, row := range amps {
		fmt.Fprintln(a, row)
	}
	if ser != nil {
		g, err := os.Create(filepath.Join(dir, "seal_sweep.csv"))
		if err != nil {
			return err
		}
		defer g.Close()
		a0, _ := ser.Pump.Pipes.radii(0)
		rMid := a0 + float64(ser.Pump.Pipes.N/2)*(ser.Pump.Pipes.Gap+ser.Pump.Pipes.Wall)
		fmt.Fprintln(g, "rseal,film_m,eff,seal,out_V,rev")
		for _, s := range sweep {
			fmt.Fprintf(g, "%g,%g,%g,%g,%g,%g\n", s.Rseal, filmFor(s.Rseal, rMid), s.Eff, s.Seal, s.Out, s.Rev)
		}
		h, err := os.Create(filepath.Join(dir, "pipes_cells.csv"))
		if err != nil {
			return err
		}
		defer h.Close()
		fmt.Fprintln(h, "annulus,rev,area_m2")
		for j, x := range ser.R.CellRev {
			_, _, a := ser.Pump.Pipes.Cell(j)
			fmt.Fprintf(h, "%d,%g,%g\n", j, x, a)
		}
	}
	return nil
}
