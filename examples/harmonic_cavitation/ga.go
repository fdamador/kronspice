// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"sort"
	"time"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
)

// Stage 3: NSGA-II (evo) over the machine. Each candidate's Build returns
// its shaft circuit (ShaftCircuit), and the objectives read the operating
// point of that circuit, so the efficiency every design is ranked by is
// the power the engine books into RCOLLAPSE over the power VSHAFT gives.

// PShaftMax is the shaft-power limit for a home heater [W] (declared).
const PShaftMax = 3000.0

// dimpleSpace: band radius, rows, dimple Ø, depth/Ø, rotor pitch/Ø,
// housing dimples (1 none, 2 yes), housing pitch/Ø, gap, speed, outlet
// pressure.
var dimpleSpace = evo.Space{
	{Name: "R_m", Min: 0.03, Max: 0.10},
	{Name: "rows", Min: 1, Max: 10, Integer: true},
	{Name: "W_m", Min: 0.003, Max: 0.025},
	{Name: "depth_per_W", Min: 0.15, Max: 0.5},
	{Name: "rotor_pitch_per_W", Min: 1.2, Max: 6},
	{Name: "housing", Min: 1, Max: 2, Integer: true},
	{Name: "housing_pitch_per_W", Min: 1.2, Max: 12},
	{Name: "gap_m", Min: 0.8e-3, Max: 5e-3},
	{Name: "rpm", Min: 500, Max: 15000},
	{Name: "p_out_Pa", Min: 1e5, Max: 5e5},
}

// DimpleFrom decodes dimpleSpace values.
func DimpleFrom(x []float64) Machine {
	ph := 0.0
	if math.Round(x[5]) == 2 {
		ph = x[6]
	}
	return DimpleMachine(x[0], int(math.Round(x[1])), x[2], x[3], x[4], ph, x[7], x[8], x[9])
}

// bulbSpace: power n, band radius, rows, housing, housing pitch, gap,
// speed, outlet pressure. Latitude rings within ±0.5 rad, detail limit
// 0.3 mm (metal printing), MaxIt 10, polar convention.
var bulbSpace = evo.Space{
	{Name: "n", Min: 3, Max: 16, Integer: true},
	{Name: "R_m", Min: 0.03, Max: 0.10},
	{Name: "rows", Min: 1, Max: 6, Integer: true},
	{Name: "housing", Min: 1, Max: 2, Integer: true},
	{Name: "housing_pitch_per_W", Min: 1.2, Max: 12},
	{Name: "gap_m", Min: 0.8e-3, Max: 5e-3},
	{Name: "rpm", Min: 500, Max: 15000},
	{Name: "p_out_Pa", Min: 1e5, Max: 5e5},
}

// Detail limits for the bulb.
const (
	bulbLatMax = 0.5
	bulbRes    = 0.3e-3
	bulbMaxIt  = 10
)

// BulbFrom decodes bulbSpace values.
func BulbFrom(x []float64) Machine {
	ph := 0.0
	if math.Round(x[3]) == 2 {
		ph = x[4]
	}
	b := Bulb{N: math.Round(x[0]), MaxIt: bulbMaxIt}
	return BulbMachine(b, x[1], int(math.Round(x[2])), bulbLatMax, bulbRes, ph, x[5], x[6], x[7])
}

// Scored is one evaluated design.
type Scored struct {
	X []float64
	M Machine
	E Eval
}

func problem(p Params, space evo.Space, from func([]float64) Machine, landOf func([]float64) float64) evo.Problem {
	eval := func(e *evo.Eval) (Eval, error) {
		v, err := e.Memo("eval", func() (any, error) { return p.Evaluate(from(e.Values)) })
		if err != nil {
			return Eval{}, err
		}
		return v.(Eval), nil
	}
	return evo.Problem{
		Space: space,
		Build: func(x []float64) (*circuit.Circuit, error) {
			m := from(x)
			ev, err := p.Evaluate(m)
			if err != nil {
				return nil, err
			}
			return ShaftCircuit(ev, m.RPM)
		},
		Objectives: []evo.Objective{
			{Name: "1 − η (collapse energy per shaft joule, from the shaft circuit's OP)", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				r, err := e.OP()
				if err != nil {
					return 0, err
				}
				src, _ := r.Element("VSHAFT")
				col, ok := r.Element("RCOLLAPSE")
				if !ok || src.PowerAbsorbed >= 0 {
					return 1, nil
				}
				return 1 - col.PowerAbsorbed/(-src.PowerAbsorbed), nil
			}},
			{Name: "100 W / collapse power", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				ev, err := eval(e)
				if err != nil {
					return 0, err
				}
				return 100 / math.Max(ev.PCollapse, 1e-3), nil
			}},
			{Name: "shaft power over the limit", Weight: 1, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				ev, err := eval(e)
				if err != nil {
					return 0, err
				}
				return math.Max(0, ev.PShaft/PShaftMax-1), nil
			}},
			{Name: "collapse energy above the cavity drag (model inconsistency)", Weight: 1, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				ev, err := eval(e)
				if err != nil {
					return 0, err
				}
				if ev.Consistent {
					return 0, nil
				}
				return ev.PCollapse/ev.PDrag - 1, nil
			}},
			{Name: "land between cavities under 1 mm", Weight: 1, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				return math.Max(0, 1e-3-landOf(e.Values)) * 1e3, nil
			}},
		},
	}
}

func dimpleLand(x []float64) float64 {
	l := (x[4] - 1) * x[2]
	if math.Round(x[5]) == 2 {
		l = math.Min(l, (x[6]-1)*x[2])
	}
	return l
}

func bulbLand([]float64) float64 { return 1 } // the bulb's lobes carry their own lands

// Evolve runs NSGA-II and returns the feasible Pareto designs, sorted by η.
func Evolve(p Params, space evo.Space, from func([]float64) Machine, landOf func([]float64) float64,
	pop, gens int, seed uint64, log func(string)) ([]Scored, error) {
	opt := evo.Options{Population: pop, Generations: gens, Seed: seed, Algorithm: evo.NSGA2, Sigma: 0.15, Timeout: 120 * time.Second}
	st, err := evo.Run(context.Background(), problem(p, space, from, landOf), opt, nil, func(st *evo.State) error {
		if log != nil {
			h := st.History[len(st.History)-1]
			if h.Generation%20 == 0 {
				log(fmt.Sprintf("    generation %d: front %d, failed %d, timed out %d", h.Generation, h.FrontSize, h.Failed, h.TimedOut))
			}
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	var out []Scored
	for _, ind := range evo.ParetoFront(st.Population) {
		m := from(ind.Values)
		e, err := p.Evaluate(m)
		if err != nil || e.PShaft > PShaftMax*1.0001 || !e.Consistent || landOf(ind.Values) < 1e-3-1e-12 {
			continue
		}
		out = append(out, Scored{X: append([]float64(nil), ind.Values...), M: m, E: e})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].E.Eta > out[j].E.Eta })
	return out, nil
}

// SensRow is the best design of a search with one declared constant changed.
type SensRow struct {
	Name string
	Best Scored
	Err  error
}

// Sensitivity reruns the dimple search (60 × 60) with each declared
// constant moved to the ends of its range.
func Sensitivity(p Params, seed uint64) []SensRow {
	type alt struct {
		name string
		f    func(*Params)
	}
	var out []SensRow
	for _, a := range []alt{
		{"Cd 0.10", func(q *Params) { q.Cd = 0.10 }},
		{"Cd 0.40", func(q *Params) { q.Cd = 0.40 }},
		{"Cpulse 0.15", func(q *Params) { q.Cpulse = 0.15 }},
		{"Cpulse 0.50", func(q *Params) { q.Cpulse = 0.50 }},
		{"κ_v 0.10", func(q *Params) { q.KappaV = 0.10 }},
		{"St only 0.446 (Xue)", func(q *Params) { q.St = []float64{p.Cal.StLock} }},
	} {
		q := p
		q.St = append([]float64(nil), p.St...)
		a.f(&q)
		fr, err := Evolve(q, dimpleSpace, DimpleFrom, dimpleLand, 60, 60, seed, nil)
		r := SensRow{Name: a.name, Err: err}
		if err == nil && len(fr) == 0 {
			r.Err = fmt.Errorf("no feasible design")
		}
		if r.Err == nil {
			r.Best = fr[0]
		}
		out = append(out, r)
	}
	return out
}

// NRow is the best Mandelbulb band with the power n fixed.
type NRow struct {
	N        int
	Eta      float64
	Cavities int
}

// PerN runs the bulb search (50 × 50) for each n from 3 to 16.
func PerN(p Params, seed uint64) []NRow {
	var rows []NRow
	for n := 3; n <= 16; n++ {
		sp := append(evo.Space(nil), bulbSpace...)
		sp[0] = evo.Gene{Name: "n", Min: float64(n), Max: float64(n), Integer: true}
		fr, err := Evolve(p, sp, BulbFrom, bulbLand, 50, 50, seed, nil)
		if err != nil || len(fr) == 0 {
			rows = append(rows, NRow{n, math.NaN(), 0})
			continue
		}
		rows = append(rows, NRow{n, fr[0].E.Eta, fr[0].E.Cavities})
	}
	return rows
}

// Heater picks the home-heater design from a front sorted by η: the one
// with the most collapse power among those within 0.5 % of the best η
// (the drum front is flat in η, so it can be sized to the heating load).
func Heater(front []Scored) Scored {
	h := front[0]
	for _, s := range front {
		if s.E.Eta >= 0.995*front[0].E.Eta && s.E.PCollapse > h.E.PCollapse {
			h = s
		}
	}
	return h
}
