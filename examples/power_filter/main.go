// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command power_filter evolves passive harmonic filters for a three-phase
// six-pulse diode rectifier, scored with geometric-algebra power theory
// (CLAUDE.md §3b item 3, milestone M9f).
//
// A 60 Hz, 120 V (phase) supply with 0.5 mH of line inductance feeds a
// diode bridge with a 2200 µF, 10 Ω DC load. At the point of common
// coupling (PCC) each phase carries two series L–C branches to a grounded
// star point, nominally for the 5th and 7th harmonics; the genome is their
// four values (E12 series). Fitness is 1 − PF at the PCC, with
// PF = ΣP/ΣS over the phases and each phase's P and S the scalar part and
// the magnitude of the geometric power M = V·I (analysis.GeometricPower).
// The report shows, before and after, the IEEE 1459 quantities and the
// pairs of harmonics that carry the most distortion power, and checks the
// quaternion instantaneous power of the three phases against ΣP.
//
//	go run ./examples/power_filter [-gens 10] [-pop 24] [-seed 1]
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/output"
)

const (
	f0      = 60.0
	vPhase  = 120.0 // RMS
	cycles  = 14    // simulated
	window  = 2     // analysed (the last cycles)
	maxHarm = 25
)

// Filter is one design: the two tuned branches per phase (henries,
// farads). A zero L or C omits that branch.
type Filter struct{ L5, C5, L7, C7 float64 }

// Circuit builds the rectifier with filter f.
func Circuit(f Filter) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("six-pulse rectifier with tuned filters")
	b.Model("dr", "D", map[string]float64{"is": 1e-12, "n": 1, "rs": 1e-3})
	for k, ph := range []string{"a", "b", "c"} {
		deg := -120 * float64(k)
		b.VoltageWave("V"+ph, "s"+ph, "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, math.Sqrt2 * vPhase, f0, 0, 0, deg}})
		b.Resistor("RS"+ph, "s"+ph, "m"+ph, 0.05)
		b.Inductor("LS"+ph, "m"+ph, "p"+ph, 0.5e-3)
		b.Diode("DU"+ph, "p"+ph, "pos", "dr")
		b.Diode("DD"+ph, "neg", "p"+ph, "dr")
		for _, br := range []struct {
			n    string
			l, c float64
		}{{"5", f.L5, f.C5}, {"7", f.L7, f.C7}} {
			if br.l > 0 && br.c > 0 {
				b.Inductor("LF"+br.n+ph, "p"+ph, "f"+br.n+ph, br.l)
				b.Capacitor("CF"+br.n+ph, "f"+br.n+ph, "r"+br.n+ph, br.c)
				b.Resistor("RF"+br.n+ph, "r"+br.n+ph, "0", 0.05)
			}
		}
	}
	b.Capacitor("CDC", "pos", "neg", 2200e-6)
	b.Resistor("RDC", "pos", "neg", 10)
	b.Resistor("RREF", "neg", "0", 1e5)
	return b.Build()
}

// Phase holds the analysis of one phase at the PCC.
type Phase struct {
	GA   *analysis.GAPower
	IEEE analysis.IEEE1459
}

// Result is the analysis of a design.
type Result struct {
	Phases   [3]Phase
	PF       float64 // ΣP / ΣS
	P        float64 // ΣP [W]
	QuatMean float64 // mean of the quaternion instantaneous active power
}

// Analyze simulates a design and evaluates its power at the PCC.
func Analyze(ctx context.Context, f Filter, opt analysis.Options) (*Result, error) {
	c, err := Circuit(f)
	if err != nil {
		return nil, err
	}
	return analyze(ctx, c, opt, nil)
}

func analyze(ctx context.Context, c *circuit.Circuit, opt analysis.Options, tran func(analysis.TranSpec) (*output.Waveforms, error)) (*Result, error) {
	h := 1 / (f0 * 600)
	spec := analysis.TranSpec{Step: h, Stop: cycles / f0, MaxStep: h}
	if tran == nil {
		opt.Context = ctx
		tran = func(s analysis.TranSpec) (*output.Waveforms, error) {
			return analysis.Tran(c, s, analysis.TranOptions{Options: opt, MaxSteps: 400000})
		}
	}
	wf, err := tran(spec)
	if err != nil {
		return nil, err
	}
	r := &Result{}
	var sumS float64
	var v3, i3 [3][]float64
	for k, ph := range []string{"a", "b", "c"} {
		v, ok := wf.NodeVoltage("p" + ph)
		iv, ok2 := wf.ElementCurrent("V" + ph)
		if !ok || !ok2 {
			return nil, fmt.Errorf("power_filter: missing waveforms for phase %s", ph)
		}
		i := make([]float64, len(iv))
		for n := range iv {
			i[n] = -iv[n] // into the PCC from the supply
		}
		v3[k], i3[k] = v, i
		hm, err := analysis.HarmonicsOf(wf.Time, v, i, f0, window, maxHarm)
		if err != nil {
			return nil, err
		}
		g, err := analysis.GeometricPower(hm)
		if err != nil {
			return nil, err
		}
		q, err := analysis.IEEE1459Of(wf.Time, v, i, f0, window)
		if err != nil {
			return nil, err
		}
		r.Phases[k] = Phase{GA: g, IEEE: q}
		r.P += g.P
		sumS += g.S
	}
	if sumS > 0 {
		r.PF = r.P / sumS
	}
	// Quaternion instantaneous power of the three phases, averaged over the
	// window (trapezoids on the solver's time points).
	t0 := wf.Time[len(wf.Time)-1] - window/f0
	var acc, span float64
	prev := math.NaN()
	for n := range wf.Time {
		if wf.Time[n] < t0 {
			continue
		}
		s := kron.InstantPower(kron.Phases(v3[0][n], v3[1][n], v3[2][n]), kron.Phases(i3[0][n], i3[1][n], i3[2][n])).W
		if !math.IsNaN(prev) {
			dt := wf.Time[n] - wf.Time[n-1]
			acc += 0.5 * (s + prev) * dt
			span += dt
		}
		prev = s
	}
	if span > 0 {
		r.QuatMean = acc / span
	}
	return r, nil
}

// Problem is the GA problem: four genes, fitness 1 − PF.
func Problem() evo.Problem {
	return evo.Problem{
		Space: evo.Space{
			{Name: "L5", Min: 0.1e-3, Max: 10e-3, Series: "E12"},
			{Name: "C5", Min: 10e-6, Max: 1000e-6, Series: "E12"},
			{Name: "L7", Min: 0.1e-3, Max: 10e-3, Series: "E12"},
			{Name: "C7", Min: 10e-6, Max: 1000e-6, Series: "E12"},
		},
		Build: func(v []float64) (*circuit.Circuit, error) {
			return Circuit(Filter{L5: v[0], C5: v[1], L7: v[2], C7: v[3]})
		},
		Objectives: []evo.Objective{{
			Name: "1-PF", Weight: 1,
			Eval: func(e *evo.Eval) (float64, error) {
				r, err := analyze(e.Ctx, e.Circuit, e.AnalysisOptions(), e.Tran)
				if err != nil {
					return 0, err
				}
				return 1 - r.PF, nil
			},
		}},
	}
}

// Report prints a design's power quantities.
func Report(name string, r *Result) {
	fmt.Printf("%s: PF %.4f, P %.0f W (quaternion mean p %.0f W)\n", name, r.PF, r.P, r.QuatMean)
	for k, ph := range []string{"a", "b", "c"} {
		q, g := r.Phases[k].IEEE, r.Phases[k].GA
		fmt.Printf("  phase %s: S %.0f VA, P %.0f W, Q1 %.0f var, D_I %.0f VA, THD_I %.1f %%, THD_V %.1f %%\n",
			ph, q.S, q.P, q.Q1, q.DI, 100*q.THDI, 100*q.THDV)
		fmt.Printf("           GA: Q %.0f var, D %.0f VA; largest pairs:", g.Q, g.D)
		for n := 0; n < 3 && n < len(g.Distortion); n++ {
			d := g.Distortion[n]
			fmt.Printf(" %d–%d %.0f VA;", d.H, d.K, d.Value)
		}
		fmt.Println()
	}
}

func main() {
	gens := flag.Int("gens", 10, "generations")
	pop := flag.Int("pop", 24, "population")
	seed := flag.Uint64("seed", 1, "seed")
	flag.Parse()
	ctx := context.Background()
	base, err := Analyze(ctx, Filter{}, analysis.Options{})
	if err != nil {
		log.Fatal(err)
	}
	Report("no filter", base)
	s, err := evo.Run(ctx, Problem(), evo.Options{Population: *pop, Generations: *gens, Seed: *seed}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	v := s.Best.Values
	f := Filter{L5: v[0], C5: v[1], L7: v[2], C7: v[3]}
	best, err := Analyze(ctx, f, analysis.Options{})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("\nevolved: L5 %.3g H, C5 %.3g F (%.0f Hz), L7 %.3g H, C7 %.3g F (%.0f Hz), %d evaluations\n",
		f.L5, f.C5, 1/(2*math.Pi*math.Sqrt(f.L5*f.C5)), f.L7, f.C7, 1/(2*math.Pi*math.Sqrt(f.L7*f.C7)), s.Evaluations)
	Report("evolved filter", best)
	if best.PF <= base.PF {
		fmt.Println("the evolved filter did not improve the power factor")
		os.Exit(1)
	}
}
