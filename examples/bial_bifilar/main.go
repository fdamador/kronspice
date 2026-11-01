// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command bial_bifilar is the project's first experiment: a bismuth–
// aluminium capacitor element with a bifilar coil in parallel, self-driven
// (the Bi–Al element is the only energy source), evolved for two competing
// goals: a high quality factor Q at the tank's resonance, and a high stored
// energy per volume. Both element variants are run and compared:
//
//   - dry: Seebeck EMF of the Bi–Al leads in series with the plate
//     capacitance (energy from heat, thermal port);
//   - wet: galvanic Al|electrolyte|Bi cell with a Randles electrode model
//     (energy from the aluminium plate, chemical port).
//
// The trade-off is traced with a sweep of weighted-sum GA runs and the
// non-dominated designs are reported (NSGA-II arrives in M9). Every design
// on a front is also simulated with an ordinary coil of the same geometry,
// and its resonance is related to the strongest emission line of copper,
// aluminium or bismuth as a whole-number subharmonic.
//
//	go run ./examples/bial_bifilar            # full run, report in runs/
//	go run ./examples/bial_bifilar -quick     # small populations
package main

import (
	"context"
	"encoding/csv"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/material"
)

// variants are the element versions the experiment compares.
var variants = []Variant{Dry, Wet, Thermopile}

// energyRef sets the zero of the energy objective: error = log10(ref/E_v),
// 0 only at 1 MJ/m³, far above anything these elements reach.
const energyRef = 1e6

func space(v Variant) evo.Space {
	gap := evo.Gene{Name: "gap_m", Min: 1e-5, Max: 2e-3}
	eps := evo.Gene{Name: "dielectric_eps", Min: 2, Max: 100}
	if v == Wet {
		gap = evo.Gene{Name: "gap_m", Min: 5e-5, Max: 5e-3}
		eps = evo.Gene{Name: "dielectric_eps", Min: 1, Max: 1} // unused: the electrolyte sets ε
	}
	sp := evo.Space{
		{Name: "turns", Min: 10, Max: 2000, Integer: true},
		{Name: "coil_diameter_m", Min: 5e-3, Max: 0.15},
		{Name: "wire_diameter_m", Min: 1e-4, Max: 2e-3},
		{Name: "pitch_ratio", Min: 1.05, Max: 3},
		{Name: "insulation_eps", Min: 2, Max: 6},
		{Name: "plate_side_m", Min: 5e-3, Max: 0.15},
		gap,
		eps,
	}
	if v == Thermopile {
		sp = append(sp,
			evo.Gene{Name: "pairs", Min: 1, Max: 5000, Integer: true},
			evo.Gene{Name: "leg_diameter_m", Min: 5e-5, Max: 2e-3},
			evo.Gene{Name: "leg_length_m", Min: 1e-3, Max: 0.05})
	}
	return sp
}

func measured(e *evo.Eval, v Variant, a Assumptions) (Result, error) {
	r, err := e.Memo("measure", func() (any, error) {
		return Measure(designFrom(v, a, true, e.Values), e.AnalysisOptions())
	})
	if err != nil {
		return Result{}, err
	}
	return r.(Result), nil
}

func problem(v Variant, a Assumptions, wQ, wE float64) evo.Problem {
	return evo.Problem{
		Space: space(v),
		Build: func(x []float64) (*circuit.Circuit, error) {
			d := designFrom(v, a, true, x)
			return d.Build(d.Coil.Resistance(0), false)
		},
		Objectives: []evo.Objective{
			{Name: "1/Q", Weight: wQ, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measured(e, v, a)
				if err != nil {
					return 0, err
				}
				if !(r.Q > 0) {
					return 0, ErrNoResonance
				}
				return 1 / r.Q, nil
			}},
			{Name: "log10(1 MJ/m³ ÷ stored energy per volume)", Weight: wE, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measured(e, v, a)
				if err != nil {
					return 0, err
				}
				if !(r.EnergyDensity > 0) {
					return 0, analysis.ErrNoMeasurement
				}
				return math.Max(0, math.Log10(energyRef/r.EnergyDensity)), nil
			}},
			{Name: "distance to a subharmonic of the line", Weight: 0.01, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measured(e, v, a)
				if err != nil {
					return 0, err
				}
				return math.Min(1, 2*r.SubN*math.Abs(r.SubOffset)), nil
			}},
		},
	}
}

// point is one design on (or near) a Pareto front.
type point struct {
	X        []float64
	Bifilar  Result
	Ordinary Result
	OrdErr   string
}

type sweepResult struct {
	Variant  Variant
	Front    []point
	Evals    int
	Failed   int
	Failures map[string]int // by class, over every evaluated individual
}

func sweep(ctx context.Context, v Variant, a Assumptions, pop, gens int, seed uint64, logf func(string, ...any)) (sweepResult, error) {
	weights := [][2]float64{{1, 0}, {1, 0.1}, {1, 1}, {0.1, 1}, {0, 1}}
	sr := sweepResult{Variant: v, Failures: map[string]int{}}
	type cand struct {
		x     []float64
		terms []float64
	}
	var cands []cand
	for i, w := range weights {
		p := problem(v, a, w[0], w[1])
		opt := evo.Options{Population: pop, Generations: gens, Seed: seed + uint64(i), Sigma: 0.15, Elite: 2,
			Timeout: 20 * time.Second}
		counted := 0 // individuals of the population already tallied (elites are re-listed)
		st, err := evo.Run(ctx, p, opt, nil, func(st *evo.State) error {
			start := 0
			if counted > 0 {
				start = opt.Elite
			}
			for _, ind := range st.Population[start:] {
				if ind.Failure != "" {
					sr.Failures[ind.Failure+": "+failureReason(ind.Error)]++
				}
			}
			counted++
			return nil
		})
		if err != nil {
			return sr, err
		}
		sr.Evals += st.Evaluations
		for _, h := range st.History {
			sr.Failed += h.Failed
		}
		for _, ind := range append(st.Population, st.Best) {
			if ind.Failure == "" {
				cands = append(cands, cand{ind.Values, ind.Terms})
			}
		}
		logf("  %s, weights Q %g / energy %g: best fitness %.4g", v, w[0], w[1], st.Best.Fitness)
	}
	pts := make([][]float64, len(cands))
	for i, c := range cands {
		pts[i] = []float64{c.terms[0], c.terms[1]}
	}
	seen := map[string]bool{}
	for _, i := range evo.NonDominated(pts) {
		key := fmt.Sprint(cands[i].x)
		if seen[key] {
			continue
		}
		seen[key] = true
		x := cands[i].x
		b, err := Measure(designFrom(v, a, true, x), analysis.Options{})
		if err != nil {
			continue
		}
		pt := point{X: x, Bifilar: b}
		if o, err := Measure(designFrom(v, a, false, x), analysis.Options{}); err == nil {
			pt.Ordinary = o
		} else {
			pt.OrdErr = err.Error()
		}
		sr.Front = append(sr.Front, pt)
	}
	sort.Slice(sr.Front, func(i, j int) bool { return sr.Front[i].Bifilar.Q > sr.Front[j].Bifilar.Q })
	return sr, nil
}

// ringDown simulates the EMF switching on and returns the transient ledger
// totals and the ringing frequency seen in v(a).
type ringDown struct {
	Stop, EnergyIn, Dissipated, StoredEnd, Mismatch, Residual, Supplied float64
	RingHz                                                              float64
	Steps                                                               int
}

func simulateStep(d Design, r Result) (ringDown, error) {
	var rd ringDown
	c, err := d.Build(r.Rcoil, true)
	if err != nil {
		return rd, err
	}
	// Long enough to see the ringing decay (≈ 2Q/ω0 per e-fold), at most
	// 2,000 periods; 40 points per period.
	periods := math.Min(2000, math.Max(20, 3*r.Q))
	stop := periods / r.F0
	spec := analysis.TranSpec{Step: 1 / (40 * r.F0), Stop: stop, MaxStep: 1 / (40 * r.F0)}
	wf, err := analysis.Tran(c, spec, analysis.TranOptions{MaxSteps: 400000})
	if err != nil {
		return rd, err
	}
	l := wf.Ledger
	rd = ringDown{Stop: stop, EnergyIn: l.EnergyIn, Dissipated: l.EnergyDissipated, StoredEnd: l.StoredEnd,
		Mismatch: l.IntegrationMismatch, Residual: l.MaxRelativeResidual, Supplied: l.EnergyModelSupplied,
		Steps: wf.Stats.Accepted}
	if v, ok := wf.NodeVoltage("a"); ok {
		// Mean frequency from rising crossings of the final value, counted
		// only while the ringing is above 1 % of its first peak: below that,
		// trapezoidal integration leaves a tiny point-to-point oscillation of
		// the circuit's fastest (stiff) RC mode — a known artefact of the
		// method (Gear-2 is planned), far inside the error tolerance, that
		// would be miscounted as ringing.
		final := v[len(v)-1]
		peak := 0.0
		for _, x := range v {
			peak = math.Max(peak, math.Abs(x-final))
		}
		var ts []float64
		for k := 1; k < len(v); k++ {
			if math.Abs(v[k]-final) > 0.01*peak || math.Abs(v[k-1]-final) > 0.01*peak {
				if v[k-1] < final && v[k] >= final {
					ts = append(ts, wf.Time[k-1]+(wf.Time[k]-wf.Time[k-1])*(final-v[k-1])/(v[k]-v[k-1]))
				}
			}
			if wf.Time[k] > 50/r.F0 {
				break
			}
		}
		if len(ts) >= 3 {
			rd.RingHz = float64(len(ts)-1) / (ts[len(ts)-1] - ts[0])
		}
	}
	return rd, nil
}

func eng(x float64, unit string) string {
	if x == 0 || math.IsNaN(x) || math.IsInf(x, 0) {
		return fmt.Sprintf("%g %s", x, unit)
	}
	prefixes := []struct {
		s string
		f float64
	}{{"G", 1e9}, {"M", 1e6}, {"k", 1e3}, {"", 1}, {"m", 1e-3}, {"µ", 1e-6}, {"n", 1e-9}, {"p", 1e-12}, {"f", 1e-15}, {"a", 1e-18}}
	a := math.Abs(x)
	for _, p := range prefixes {
		if a >= p.f*0.9999995 {
			return fmt.Sprintf("%.4g %s%s", x/p.f, p.s, unit)
		}
	}
	return fmt.Sprintf("%.3e %s", x, unit)
}

func main() {
	quick := flag.Bool("quick", false, "small populations (a few seconds)")
	out := flag.String("out", "runs", "base directory for the report")
	line := flag.String("line", "Cu", "spectral line for the subharmonic: Cu, Al or Bi")
	dT := flag.Float64("deltaT", 100, "dry variant: plate temperature above the coil connections [K]")
	emf := flag.Float64("emf-factor", 1, "wet variant: fraction of the standard cell EMF reached")
	seed := flag.Uint64("seed", 1, "random seed")
	netlists := flag.String("netlists", "", "write SPICE netlists of the best designs from the front CSVs in this directory, then exit")
	flag.Parse()
	a := DefaultAssumptions()
	a.DeltaT, a.EMFFactor = *dT, *emf
	switch strings.ToLower(*line) {
	case "cu", "copper":
		a.Line = material.LineCopper
	case "al", "aluminium", "aluminum":
		a.Line = material.LineAluminium
	case "bi", "bismuth":
		a.Line = material.LineBismuth
	default:
		fmt.Fprintln(os.Stderr, "unknown -line; use Cu, Al or Bi")
		os.Exit(2)
	}
	if *netlists != "" {
		files, err := writeNetlists(*netlists, a)
		for _, f := range files {
			fmt.Println("wrote", f)
		}
		if err != nil {
			fmt.Fprintln(os.Stderr, "bial_bifilar:", err)
			os.Exit(1)
		}
		return
	}
	pop, gens := 40, 30
	if *quick {
		pop, gens = 12, 4
	}
	dir := filepath.Join(*out, time.Now().Format("20060102-150405")+"-bial_bifilar")
	if err := run(context.Background(), a, pop, gens, *seed, dir, func(f string, x ...any) { fmt.Printf(f+"\n", x...) }); err != nil {
		fmt.Fprintln(os.Stderr, "bial_bifilar:", err)
		os.Exit(1)
	}
}

func run(ctx context.Context, a Assumptions, pop, gens int, seed uint64, dir string, logf func(string, ...any)) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	var results []sweepResult
	for _, v := range variants {
		logf("evolving the %s variant (%d runs × %d × %d)", v, 5, pop, gens)
		sr, err := sweep(ctx, v, a, pop, gens, seed, logf)
		if err != nil {
			return err
		}
		results = append(results, sr)
		if err := writeCSV(filepath.Join(dir, fmt.Sprintf("front_%s.csv", v.short())), sr); err != nil {
			return err
		}
	}
	rep, err := report(a, results, pop, gens, seed)
	if err != nil {
		return err
	}
	path := filepath.Join(dir, "report.md")
	if err := os.WriteFile(path, []byte(rep), 0o644); err != nil {
		return err
	}
	logf("report: %s", path)
	return nil
}

func writeCSV(path string, sr sweepResult) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	w := csv.NewWriter(f)
	hdr := geneNames(sr.Variant)
	hdr = append(hdr, "f0_hz", "q", "energy_density_j_m3", "stored_j", "volume_m3", "L_h", "cself_f", "rcoil_ohm",
		"emf_v", "source_power_w", "heat_w", "subharmonic_n", "subharmonic_offset",
		"ordinary_f0_hz", "ordinary_q", "ordinary_energy_density_j_m3")
	w.Write(hdr)
	ff := func(x float64) string { return strconv.FormatFloat(x, 'g', 10, 64) }
	for _, p := range sr.Front {
		row := []string{}
		for _, x := range p.X {
			row = append(row, ff(x))
		}
		b, o := p.Bifilar, p.Ordinary
		row = append(row, ff(b.F0), ff(b.Q), ff(b.EnergyDensity), ff(b.StoredJ), ff(b.Volume), ff(b.L), ff(b.Cself),
			ff(b.Rcoil), ff(b.Element.EMF), ff(b.SourcePower), ff(b.Heat), ff(b.SubN), ff(b.SubOffset),
			ff(o.F0), ff(o.Q), ff(o.EnergyDensity))
		w.Write(row)
	}
	w.Flush()
	if err := w.Error(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

func report(a Assumptions, rs []sweepResult, pop, gens int, seed uint64) (string, error) {
	var b strings.Builder
	p := func(f string, x ...any) { fmt.Fprintf(&b, f, x...) }
	worstResidual, maxSupplied, worstMismatchRel := 0.0, 0.0, 0.0
	transients := 0
	p("# Bi–Al capacitor element with a bifilar coil: first experiment\n\n")
	p("Generated %s by `examples/bial_bifilar` (seed %d; 5 weighted-sum runs per variant, population %d, %d generations).\n\n",
		time.Now().UTC().Format("2006-01-02 15:04 UTC"), seed, pop, gens)
	p("## Set-up\n\n")
	p("A coil wound bifilar (Tesla, US 512,340) in parallel with a bismuth–aluminium element, self-driven: the element is the only energy source. ")
	p("Evolution trades the tank's quality factor Q at resonance against the energy stored per volume (coil + element) once the element's EMF is on. ")
	p("Each non-dominated design is also simulated with an ordinary coil of the same geometry.\n\n")
	p("| Assumption | Value | Source |\n|---|---|---|\n")
	p("| Bi–Al Seebeck coefficient | %.1f µV/K (S_Al − S_Bi, relative to Pt: +3.5 and −72 µV/K) | %s |\n",
		(material.Aluminium.Seebeck-material.Bismuth.Seebeck)*1e6, "Electronics Cooling (2006), Seebeck coefficient table")
	p("| Dry: plates above coil connections | ΔT = %g K → EMF %s | assumption |\n", a.DeltaT,
		eng(material.ThermocoupleEMF(material.Bismuth, material.Aluminium, a.DeltaT), "V"))
	p("| Wet: cell EMF | %s (%g × standard 1.970 V: Al³⁺/Al −1.662 V, Bi³⁺/Bi +0.308 V) | CRC electrochemical series; real cells are lower |\n",
		eng(a.EMFFactor*material.CellEMF(material.Aluminium, material.Bismuth), "V"), a.EMFFactor)
	p("| Wet: electrolyte | ρ = %g Ω·m (sea water), ε = %g | UNESCO 1981 |\n", a.ElectrolyteRho, a.WaterEps)
	p("| Wet: double layer, charge transfer | %g µF/cm², %g Ω·cm² per electrode | Bard & Faulkner (typical); assumed |\n",
		a.DoubleLayer*100, a.ChargeTransfer*1e4)
	p("| Plates | %g mm Bi + %g mm Al | assumption |\n", a.PlateThickness*1e3, a.PlateThickness*1e3)
	p("| Coil | copper; Wheeler inductance; turn-to-turn capacitance (two-wire formula), V/N per gap ordinary, V/2 bifilar; skin effect at f0 | Wheeler 1928; Tesla US 512,340; Terman 1943 |\n")
	p("| Spectral line | %s %.3f nm (air), %.4g THz | NIST Handbook of Basic Atomic Spectroscopic Data |\n\n",
		a.Line.Element, a.Line.AirWavelength*1e9, a.Line.Frequency()/1e12)

	for _, sr := range rs {
		p("## %s variant\n\n", strings.ToUpper(sr.Variant.String()[:1])+sr.Variant.String()[1:])
		p("%d evaluations, %d penalized", sr.Evals, sr.Failed)
		if len(sr.Failures) > 0 {
			keys := make([]string, 0, len(sr.Failures))
			for k := range sr.Failures {
				keys = append(keys, k)
			}
			sort.Slice(keys, func(i, j int) bool { return sr.Failures[keys[i]] > sr.Failures[keys[j]] })
			parts := []string{}
			for _, k := range keys {
				parts = append(parts, fmt.Sprintf("%d × %s", sr.Failures[k], k))
			}
			p(" (%s)", strings.Join(parts, "; "))
		}
		p(". ")
		p("%d non-dominated designs (highest Q first):\n\n", len(sr.Front))
		pileCols, pileSep := "", ""
		if sr.Variant == Thermopile {
			pileCols, pileSep = " pairs | leg Ø × length |", "---|---|"
		}
		p("| # | turns | coil Ø | wire Ø | pitch | plate | gap |%s f0 | Q | stored | per volume | ordinary coil: f0 / Q / per volume |\n", pileCols)
		p("|---|---|---|---|---|---|---|%s---|---|---|---|---|\n", pileSep)
		for i, pt := range sr.Front {
			x, bi, o := pt.X, pt.Bifilar, pt.Ordinary
			ord := pt.OrdErr
			if ord == "" {
				ord = fmt.Sprintf("%s / %.3g / %s", eng(o.F0, "Hz"), o.Q, eng(o.EnergyDensity, "J/m³"))
			}
			pile := ""
			if sr.Variant == Thermopile {
				pile = fmt.Sprintf(" %g | %s × %s |", x[8], eng(x[9], "m"), eng(x[10], "m"))
			}
			p("| %d | %g | %s | %s | %.2f×d | %s | %s |%s %s | %.3g | %s | %s | %s |\n", i+1, x[0], eng(x[1], "m"), eng(x[2], "m"), x[3],
				eng(x[5], "m"), eng(x[6], "m"), pile, eng(bi.F0, "Hz"), bi.Q, eng(bi.StoredJ, "J"), eng(bi.EnergyDensity, "J/m³"), ord)
		}
		if len(sr.Front) == 0 {
			p("_No valid design found._\n\n")
			continue
		}
		hiQ, hiE := sr.Front[0], sr.Front[0]
		for _, pt := range sr.Front {
			if pt.Bifilar.EnergyDensity > hiE.Bifilar.EnergyDensity {
				hiE = pt
			}
		}
		p("\n")
		for _, c := range []struct {
			name string
			pt   point
		}{{"Highest Q", hiQ}, {"Highest energy per volume", hiE}} {
			r := c.pt.Bifilar
			d := r.Design
			p("**%s** — %d turns, L = %s, coil self-capacitance %s bifilar (ordinary would be %s; ratio %.4g = N²/4 at this N), ",
				c.name, d.Coil.Turns, eng(r.L, "H"), eng(r.Cself, "F"),
				eng(withBifilar(d, false).Coil.SelfCapacitance(), "F"),
				r.Cself/withBifilar(d, false).Coil.SelfCapacitance())
			p("coil resistance at f0 %s. Element EMF %s. ", eng(r.Rcoil, "Ω"), eng(r.Element.EMF, "V"))
			if d.Variant.seebeck() {
				p("At the steady state no current flows (the plate capacitance blocks DC): the EMF only charges the plates, "+
					"with time constant (R_int + R_coil)·C_plate = %s; a switch-on transient shorter than that ends with the plates partly charged. ",
					eng((r.Element.Rint+r.Rcoil)*r.Element.Cplate, "s"))
				if d.Variant == Thermopile {
					p("Thermopile: %d pairs, legs %s × %s, internal resistance %s; holding ΔT = %g K costs %s of heat conducted through the legs (not part of the electrical ledger). ",
						d.Pairs, eng(d.LegDiameter, "m"), eng(d.LegLength, "m"), eng(r.Element.Rint, "Ω"), d.A.DeltaT, eng(r.Element.HeatLeak, "W"))
				}
			} else {
				p("At the steady state the cell drives %s through the coil; source power %s, all of it dissipated (%s). ",
					eng(r.SourcePower/r.Element.EMF, "A"), eng(r.SourcePower, "W"), eng(r.Heat, "W"))
			}
			if d.Variant == Wet {
				p("Aluminium plate capacity %s (%.3g Ah): at the steady-state current it lasts %s. ",
					eng(r.Element.CapacityC, "C"), r.Element.CapacityC/3600, lifetime(r))
			}
			p("Resonance %s is the n = %.6g subharmonic of %s, offset %.3g ppb (any frequency is within 1/(2n) = %.3g ppb of some subharmonic).\n\n",
				eng(r.F0, "Hz"), r.SubN, a.Line.Element, r.SubOffset*1e9, 0.5/r.SubN*1e9)
			if rd, err := simulateStep(d, r); err == nil {
				transients++
				worstResidual = math.Max(worstResidual, rd.Residual)
				maxSupplied = math.Max(maxSupplied, math.Abs(rd.Supplied))
				if rd.EnergyIn > 0 {
					worstMismatchRel = math.Max(worstMismatchRel, math.Abs(rd.Mismatch)/rd.EnergyIn)
				}
				p("  Switch-on transient over %s (%d steps): energy in from the %s port %s; dissipated %s; stored at the end %s; "+
					"integration mismatch %s; worst ledger residual %.2g; model-supplied %s. %s\n\n",
					eng(rd.Stop, "s"), rd.Steps, d.channelName(), eng(rd.EnergyIn, "J"), eng(rd.Dissipated, "J"), eng(rd.StoredEnd, "J"),
					eng(rd.Mismatch, "J"), rd.Residual, eng(rd.Supplied, "J"), ringText(rd, r))
			} else {
				p("  Switch-on transient: %v\n\n", err)
			}
		}
	}
	p("## Reading the results\n\n")
	for _, sr := range rs {
		if len(sr.Front) == 0 {
			continue
		}
		lowerF, higherQ := 0, 0
		for _, pt := range sr.Front {
			if pt.OrdErr != "" {
				continue
			}
			if pt.Bifilar.F0 < pt.Ordinary.F0*(1-1e-3) {
				lowerF++
			}
			if pt.Bifilar.Q > pt.Ordinary.Q*(1+1e-3) {
				higherQ++
			}
		}
		p("- **%s, bifilar against ordinary.** On %d of %d front designs the bifilar winding lowers the resonance by more than 0.1 %%; on %d it raises Q by more than 0.1 %%. ",
			sr.Variant, lowerF, len(sr.Front), higherQ)
		p("The stored energy at the steady state is the same for both windings: at DC no voltage stands across the coil's turn-to-turn capacitance, so the winding style changes the resonance, not the stored energy.\n")
		best := sr.Front[0]
		for _, pt := range sr.Front {
			if pt.Bifilar.EnergyDensity > best.Bifilar.EnergyDensity {
				best = pt
			}
		}
		atBound := []string{}
		sp := space(sr.Variant)
		for i, g := range sp {
			if g.Min == g.Max {
				continue
			}
			x := best.X[i]
			if x <= g.Min*(1+1e-3) || x >= g.Max*(1-1e-3) {
				atBound = append(atBound, g.Name)
			}
		}
		if len(atBound) > 0 {
			p("- **%s:** the highest energy per volume sits at the edge of the allowed range for %s; a wider range would move it further, so that optimum is set by the chosen limits.\n",
				sr.Variant, strings.Join(atBound, ", "))
		}
		if sr.Variant.seebeck() {
			r := best.Bifilar
			d := r.Design
			e := r.Element.EMF / d.Gap
			bound := 0.5 * material.Eps0 * d.DielectricEps * e * e * d.Gap / (d.Gap + 2*d.A.PlateThickness)
			p("- **%s, physical limit.** The only store is the field in the gap, ½·ε·(EMF/gap)², diluted by the plates: for the best design %s, against %s simulated (coil volume adds to the denominator). ",
				sr.Variant, eng(bound, "J/m³"), eng(r.EnergyDensity, "J/m³"))
			if sr.Variant == Dry {
				p("With a millivolt-level Seebeck EMF the energy per volume is necessarily microscopic; a thermopile of many junctions in series, or a larger ΔT, raises it as EMF².\n")
			} else {
				p("With %d pairs the EMF is %s, %.4g times one junction's, and the stored energy rises with its square; the legs' heat leak (%s) is the price of holding ΔT.\n",
					d.Pairs, eng(r.Element.EMF, "V"), float64(d.Pairs), eng(r.Element.HeatLeak, "W"))
			}
		} else {
			r := best.Bifilar
			jd := r.SourcePower / r.Element.EMF / (r.Design.PlateSide * r.Design.PlateSide)
			p("- **Wet, what limits it.** The stored energy is mostly ½·L·I² from the steady cell current; the current is set by the assumed charge-transfer resistance (%g Ω·cm² per electrode) and reaches %s (%.3g mA/cm²) in the best design. ",
				r.Design.A.ChargeTransfer*1e4, eng(r.SourcePower/r.Element.EMF, "A"), jd/10)
			p("A real aluminium electrode passivates well below such current densities, so these numbers are upper bounds until the cell is measured. The aluminium plate's charge lasts the same few hours for every design because both capacity and current scale with plate area.\n")
		}
	}
	ledgerFails := 0
	for _, sr := range rs {
		for k, n := range sr.Failures {
			if strings.HasPrefix(k, "ledger") {
				ledgerFails += n
			}
		}
	}
	p("- Energy accounting: %d evaluated designs failed the operating-point ledger check. Over the %d switch-on transients: worst per-point ledger residual %.2g (relative to the run's peak power), ", ledgerFails, transients, worstResidual)
	p("worst integration mismatch %.2g of the energy taken in, model-supplied energy %s. ", worstMismatchRel, eng(maxSupplied, "J"))
	if worstResidual <= 1e-9 && maxSupplied == 0 && ledgerFails == 0 {
		p("So all stored and dissipated energy is accounted for by the element's declared source — heat for the dry variant, the aluminium plate's chemistry for the wet one — and nothing else.\n")
	} else {
		p("**The ledger does not close within tolerance: this must be investigated as a model or numerical error before any result is used** (CLAUDE.md §1, rule 4).\n")
	}
	p("- The bifilar winding raises the coil's self-capacitance by about N²/4 compared with an ordinary winding of the same geometry (Tesla's 250,000 at N = 1,000), which lowers the self-resonance; the inductance is unchanged.\n")
	p("- The idealized turn-to-turn model is a first approximation (it underestimates an ordinary coil's measured self-capacitance; Medhurst 1947). Field-solved values arrive with `em/extract` (M9i).\n")
	p("- The neutrino channel is not yet computed (the environment model, M9k); it is listed as not modeled.\n")
	return b.String(), nil
}

func withBifilar(d Design, bif bool) Design {
	d.Coil.Bifilar = bif
	return d
}

func (d Design) channelName() string {
	if d.Variant == Wet {
		return "chemical"
	}
	return "thermal"
}

func lifetime(r Result) string {
	i := r.SourcePower / math.Max(r.Element.EMF, 1e-30)
	if !(i > 0) {
		return "indefinitely (no steady current)"
	}
	return (time.Duration(r.Element.CapacityC/i) * time.Second).Round(time.Second).String()
}

func ringText(rd ringDown, r Result) string {
	if rd.RingHz == 0 {
		return fmt.Sprintf("The ringing dies out within about two cycles (Q = %.3g).", r.Q)
	}
	return fmt.Sprintf("It rings at %s (resonance %s, difference %.2g %%).", eng(rd.RingHz, "Hz"), eng(r.F0, "Hz"),
		100*(rd.RingHz-r.F0)/r.F0)
}

// failureReason shortens an error message to its cause, for tallies.
func failureReason(msg string) string {
	switch {
	case strings.Contains(msg, "not quasi-static"):
		return "wire longer than λ/10 at resonance"
	case strings.Contains(msg, "Wheeler"):
		return "coil too short for Wheeler's formula"
	case strings.Contains(msg, "pitch"):
		return "pitch below wire diameter"
	case strings.Contains(msg, "1/5 of the coil diameter"):
		return "wire too thick for the coil"
	case strings.Contains(msg, "gap"):
		return "plate gap not small against the plates"
	case strings.Contains(msg, "no parallel resonance"):
		return "no resonance"
	}
	if i := strings.LastIndex(msg, ": "); i >= 0 && len(msg)-i < 80 {
		return msg[i+2:]
	}
	return msg
}
