// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"sort"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
)

// The genes: the drive, the resonator, the ladder, the load, and the
// aluminium stacks (cells in series for voltage, plates in parallel for
// current). Ranges are practical: forming up to 700 V (capacitor foil),
// plates from 3 mm × 3 mm to 10 cm × 10 cm, up to 10 cells in series.
var pumpSpace = evo.Space{
	{Name: "drive_Hz", Min: 50, Max: 100e3},
	{Name: "drive_V", Min: 1, Max: 50},
	{Name: "tank_F", Min: 1e-9, Max: 10e-6},
	{Name: "stages", Min: 1, Max: 8, Integer: true},
	{Name: "ladder_F", Min: 10e-9, Max: 100e-6},
	{Name: "load_ohm", Min: 1e4, Max: 1e9},
	{Name: "cells_series", Min: 1, Max: 10, Integer: true},
	{Name: "plates_parallel", Min: 1, Max: 20, Integer: true},
	{Name: "plate_m2", Min: 1e-5, Max: 1e-2},
	{Name: "gap_m", Min: 1e-3, Max: 2e-2},
	{Name: "form_V", Min: 50, Max: maxForm},
}

// PumpFrom decodes gene values into a pump.
func PumpFrom(x []float64) Pump {
	p := DefaultPump()
	p.F, p.Vs, p.Ct = x[0], x[1], x[2]
	p.Stages = int(math.Round(x[3]))
	p.Cstage = x[4]
	p.Cbank = x[4]
	p.Rload = x[5]
	// The film is evolved as a plain capacitor (κ = 0): the varactor term
	// costs 16× the run time and changes the result by a few per cent; the
	// Pareto designs are re-measured with κ = 2 % (VaractorCheck).
	p.Cell = Stack{Ns: int(math.Round(x[6])), Np: int(math.Round(x[7])), Area: x[8], Gap: x[9], Vform: x[10]}
	p.Cycles, p.PerCycle, p.Window = 300, 80, 60
	return p
}

// pipeSpace: the same drive, resonator, ladder and load, and one tank of
// concentric pipes per rectifier: wired in parallel (alternating
// aluminium and iron on two bus bars) or in series (clad pipes, sealed
// ends — perfect seals while evolving), plain or graded to equal areas, the number of annuli, the rod's
// radius, the water gap, the wetted length and the forming voltage.
var pipeSpace = evo.Space{
	pumpSpace[0], pumpSpace[1], pumpSpace[2], pumpSpace[3], pumpSpace[4], pumpSpace[5],
	{Name: "wiring", Min: 1, Max: 3, Integer: true}, // 1 parallel, 2 series, 3 graded series
	{Name: "annuli", Min: 1, Max: 10, Integer: true},
	{Name: "rod_m", Min: 1e-3, Max: 30e-3},
	{Name: "gap_m", Min: 1e-3, Max: 2e-2},
	{Name: "length_m", Min: 5e-3, Max: 0.3},
	{Name: "form_V", Min: 50, Max: maxForm},
}

// PipesFrom decodes pipeSpace gene values into a pump.
func PipesFrom(x []float64) Pump {
	p := PumpFrom(append(append([]float64{}, x[:6]...), 1, 1, 1e-3, 1e-3, 100))
	w := math.Round(x[6])
	p.Pipes = &Pipes{Series: w >= 2, Graded: w == 3, N: int(math.Round(x[7])), R0: x[8], Gap: x[9], Wall: wallDef, Len: x[10], Vform: x[11]}
	return p
}

func measurePump(e *evo.Eval, from func([]float64) Pump) (Result, error) {
	r, err := e.Memo("pump", func() (any, error) {
		return from(e.Values).Run(e.AnalysisOptions())
	})
	if err != nil {
		return Result{}, err
	}
	return r.(Result), nil
}

// pumpProblem: three objectives, one constraint.
//
//   - 1 − η, η = power into the load and bank ÷ power from the source:
//     the "catalytic efficiency". Lifting an electron through |V| costs
//     e·|V| whatever the circuit; η = 1 means nothing more is spent.
//   - 1 kV / |Vout| (higher voltage);
//   - 1 mA / Iout (more current);
//   - constraints: no stack above 90 % of its holding voltage (an oxide
//     driven past its forming voltage sparks and is destroyed), and the
//     bank in steady state (its voltage drifting under 1 % over the
//     measured window), so energy still going into storage is not
//     counted as output.
func pumpProblem(space evo.Space, from func([]float64) Pump) evo.Problem {
	return evo.Problem{
		Space: space,
		Build: func(x []float64) (*circuit.Circuit, error) { return from(x).Build() },
		Objectives: []evo.Objective{
			{Name: "1 − efficiency", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measurePump(e, from)
				if err != nil {
					return 0, err
				}
				return math.Max(0, 1-r.Eff), nil
			}},
			{Name: "1 kV / |Vout|", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measurePump(e, from)
				if err != nil {
					return 0, err
				}
				return 1e3 / math.Abs(r.Out), nil
			}},
			{Name: "1 mA / Iout", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measurePump(e, from)
				if err != nil {
					return 0, err
				}
				return 1e-3 / r.Iout, nil
			}},
			{Name: "reverse voltage over 90 % of the oxide's", Weight: 1, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measurePump(e, from)
				if err != nil {
					return 0, err
				}
				return math.Max(0, r.RevRatio-0.9), nil
			}},
			{Name: "not in steady state (bank drift over 1 %)", Weight: 1, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := measurePump(e, from)
				if err != nil {
					return 0, err
				}
				return math.Max(0, r.Drift-0.01), nil
			}},
		},
	}
}

// Evolved is one Pareto design with its measurement.
type Evolved struct {
	Pump Pump
	R    Result
}

// EvolvePump runs NSGA-II and returns the feasible Pareto designs,
// re-measured.
func EvolvePump(pop, gens int, seed uint64, log func(string)) ([]Evolved, error) {
	return evolveWith(pumpSpace, PumpFrom, pop, gens, seed, log)
}

// EvolvePipes evolves the one-tank concentric-pipe rectifiers.
func EvolvePipes(pop, gens int, seed uint64, log func(string)) ([]Evolved, error) {
	return evolveWith(pipeSpace, PipesFrom, pop, gens, seed, log)
}

func evolveWith(space evo.Space, from func([]float64) Pump, pop, gens int, seed uint64, log func(string)) ([]Evolved, error) {
	opt := evo.Options{Population: pop, Generations: gens, Seed: seed, Algorithm: evo.NSGA2, Sigma: 0.2, Timeout: 120 * time.Second}
	st, err := evo.Run(context.Background(), pumpProblem(space, from), opt, nil, func(st *evo.State) error {
		if log != nil {
			b := st.History[len(st.History)-1]
			log(fmt.Sprintf("generation %d: front %d designs, failed %d, timed out %d", b.Generation, b.FrontSize, b.Failed, b.TimedOut))
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	var out []Evolved
	for _, ind := range evo.ParetoFront(st.Population) {
		p := from(ind.Values)
		r, err := p.Run(analysis.Options{})
		if err != nil {
			continue
		}
		out = append(out, Evolved{p, r})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].R.Eff > out[j].R.Eff })
	return out, nil
}
