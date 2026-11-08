// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"math/cmplx"
	"sync"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
)

// Part 3: groups of output coils round the toroid, co-optimized with the
// sin/cos drive frequency by the genetic algorithm (evo, value genome,
// seeded). Genes are log-scaled by evo; quantities that should be
// searched evenly on a linear range (angles, positions) use a gene x in
// [1, 10] mapped through log10(x) ∈ [0, 1].
//
// Fitness: average power into the load(s) from the full AC solution of
// the sectional network with the groups in it (they load the drive), a
// resistive load R_L standing for each rectifier and bank, times
// BridgeFactor (the share left after two diode drops); the winner is
// then run in a kronspice transient with real Schottky bridges and the
// capacitor bank, which gives the bank's charging and the ledger.

// Space is the declared genome with its bounds.
func Part3Space() evo.Space { return Part3SpaceLimited(1, 0.06) }

// Part3SpaceLimited is the genome with at least minG groups and coil radii
// up to maxR (to search for distributed groups of local coils).
func Part3SpaceLimited(minG int, maxR float64) evo.Space {
	return evo.Space{
		{Name: "f [Hz]", Min: 1e3, Max: 3e7},
		{Name: "groups", Min: float64(minG), Max: 8, Integer: true},
		{Name: "θ0 (log10 x · 360°)", Min: 1, Max: 10},
		{Name: "Δθ (log10 x · 360°/groups)", Min: 1, Max: 10},
		{Name: "β (log10 x · 90°)", Min: 1, Max: 10},
		{Name: "ρ (log10 x · 9 cm)", Min: 1, Max: 10},
		{Name: "z (−7 cm + log10 x · 14 cm)", Min: 1, Max: 10},
		{Name: "coil radius [m]", Min: 3e-3, Max: maxR},
		{Name: "turns", Min: 1, Max: 200, Integer: true},
		{Name: "combine (1 separate, 2 series, 3 parallel)", Min: 1, Max: 3, Integer: true},
		{Name: "R_L [Ω]", Min: 0.1, Max: 1e4},
		{Name: "sense (1 toroidal, 2 dipole)", Min: 1, Max: 2, Integer: true},
	}
}

// Decode turns gene values into a layout.
func Decode(v []float64) Layout {
	lin := func(x float64) float64 { return math.Log10(x) }
	ng := int(v[1])
	l := Layout{F: v[0], NG: ng,
		Theta0: 2 * math.Pi * lin(v[2]),
		DTh:    2 * math.Pi / float64(ng) * lin(v[3]),
		Shape: GroupShape{Beta: math.Pi / 2 * lin(v[4]), Rho: 0.09 * lin(v[5]), Z: -0.07 + 0.14*lin(v[6]),
			Radius: v[7], Turns: int(v[8])},
		Comb: Combine(int(v[9])), RLoad: v[10], Sense: Toroidal}
	if int(v[11]) == 2 {
		l.Sense = Dipole
	}
	return l
}

// Feasible returns the worst clearance [m] of a layout: every coil clear
// of the bundles and outside the former, and the coils clear of each
// other (negative = violation).
func (m *Model) Feasible(l Layout) float64 {
	th := l.Thetas()
	worst := math.Inf(1)
	// a group coil must stay a lumped coil at f (wire ≤ λ/20) and its
	// tuning capacitor must be a real part (≥ MinTuneC), well above the
	// coil's own winding capacitance; violations count as metres short
	lw := float64(l.Shape.Turns) * 2 * math.Pi * l.Shape.Radius
	if d := em.SpeedOfLight/(20*l.F) - lw; d < 0 {
		worst = math.Min(worst, d*1e-2)
	}
	if ct := 1 / (math.Pow(2*math.Pi*l.F, 2) * l.Shape.SelfL()); ct < MinTuneC {
		worst = math.Min(worst, -1e-3*math.Log10(MinTuneC/ct)-1e-6)
	}
	loops := make([][]P3, len(th))
	for k, t := range th {
		worst = math.Min(worst, m.Clear(l.Shape, t))
		loops[k] = l.Shape.Loop(t)
	}
	for i := range loops {
		for j := i + 1; j < len(loops); j++ {
			d := math.Inf(1)
			for _, p := range loops[i] {
				for _, q := range loops[j] {
					d = math.Min(d, dist(p, q))
				}
			}
			worst = math.Min(worst, d-2*OutWireRadius*math.Sqrt(float64(l.Shape.Turns)/Fill)-Clearance)
		}
	}
	return worst
}

// MinTuneC is the smallest tuning capacitor allowed (declared): ten times
// a few pF of winding capacitance.
const MinTuneC = 30e-12

// Part3 holds the two driven nets (one per sense).
type Part3 struct {
	M    *Model
	Nets map[int]*Net
	mu   sync.Mutex
}

// NewPart3 builds the driven nets.
func NewPart3(m *Model) *Part3 {
	return &Part3{M: m, Nets: map[int]*Net{Toroidal: m.NewNet(DrivenPair, Toroidal, 1e5), Dipole: m.NewNet(DrivenPair, Dipole, 1e5)}}
}

// net returns a private copy of the driven net (SetFreq writes R).
func (p *Part3) net(sense int) *Net {
	p.mu.Lock()
	defer p.mu.Unlock()
	src := p.Nets[sense]
	c := *src
	c.R = append([]float64(nil), src.R...)
	return &c
}

// Score evaluates a layout: load power [W] and the AC result.
func (p *Part3) Score(l Layout) (float64, *ACResult, error) {
	if c := p.M.Feasible(l); c < 0 {
		return 0, nil, fmt.Errorf("infeasible by %.3g m", -c)
	}
	r, err := p.M.Evaluate(p.net(l.Sense), l)
	if err != nil {
		return 0, nil, err
	}
	return r.PLoad, r, nil
}

// BridgeDrop is the forward drop of one Schottky diode at the bank's
// charging current (declared, 1N5819-class).
const BridgeDrop = 0.3

// BridgeFactor estimates the share of the AC load power a bridge passes
// into the bank, smoothly: (V_pk/(V_pk + 2·V_d))², near 1 when the load
// voltage amplitude is far above two diode drops and falling as the
// square of V_pk below them (the smallest group voltage counts when there
// are several bridges). It keeps the search away from outputs too weak to
// rectify while leaving it a gradient; the kronspice transient then
// measures the real charging.
func BridgeFactor(r *ACResult) float64 {
	f := math.Inf(1)
	for _, v := range r.VLoad {
		a := cmplx.Abs(v)
		f = math.Min(f, math.Pow(a/(a+2*BridgeDrop), 2))
	}
	if math.IsInf(f, 1) {
		return 0
	}
	return f
}

// GAOptions size the search.
type GAOptions struct {
	Population, Generations int
	Seed                    uint64
	MinGroups               int     // 0 = 1
	MaxRadius               float64 // 0 = 6 cm
}

// GAResult is the outcome.
type GAResult struct {
	Best    Layout
	Power   float64
	AC      *ACResult
	History []evo.GenStats
	Evals   int
}

// Evolve runs the genetic algorithm. Fitness to minimize:
// 6 − log10(P_load·BridgeFactor / 1 W) for a feasible layout (≥ 0 for any P below 1 MW),
// 100 + 10³ × (clearance violation in m) for an infeasible one.
func (p *Part3) Evolve(o GAOptions) (*GAResult, error) {
	dummy := func([]float64) (*circuit.Circuit, error) {
		return circuit.NewBuilder("placeholder").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "0", 1).Build()
	}
	space := Part3Space()
	if o.MinGroups > 1 || o.MaxRadius > 0 {
		space = Part3SpaceLimited(max(1, o.MinGroups), map[bool]float64{true: o.MaxRadius, false: 0.06}[o.MaxRadius > 0])
	}
	prob := evo.Problem{Build: dummy, Space: space, Objectives: []evo.Objective{{Name: "6 − log10(P_load/W)", Weight: 1,
		Eval: func(e *evo.Eval) (float64, error) {
			l := Decode(e.Values)
			if c := p.M.Feasible(l); c < 0 {
				return 100 + 1e3*(-c), nil
			}
			pw, r, err := p.Score(l)
			if err != nil {
				return 0, err
			}
			return 6 - math.Log10(pw*BridgeFactor(r)+1e-30), nil
		}}}}
	st, err := evo.Run(context.Background(), prob, evo.Options{Population: o.Population, Generations: o.Generations, Seed: o.Seed, Timeout: 600e9}, nil, nil)
	if err != nil {
		return nil, err
	}
	if st.Best.Failure != "" {
		return nil, fmt.Errorf("best failed: %s %s", st.Best.Failure, st.Best.Error)
	}
	best := Decode(st.Best.Values)
	pw, ac, err := p.Score(best)
	if err != nil {
		return nil, err
	}
	return &GAResult{Best: best, Power: pw, AC: ac, History: st.History, Evals: st.Evaluations}, nil
}

// GroupEMFs returns each group's current phase [deg] in the AC solution:
// out-of-phase groups would cancel in a shared bridge.
func GroupPhases(r *ACResult) []float64 {
	out := make([]float64, len(r.IG))
	for k, i := range r.IG {
		out[k] = cmplx.Phase(i) * 180 / math.Pi
	}
	return out
}

func (l Layout) String() string {
	s := "toroidal"
	if l.Sense == Dipole {
		s = "dipole"
	}
	return fmt.Sprintf("f = %.4g Hz, sense %s, %d group(s) from θ0 = %.1f° every %.1f°, coil radius %.1f mm, %d turns, β = %.1f°, at ρ = %.1f mm, z = %.1f mm, %s, R_L = %.3g Ω",
		l.F, s, l.NG, l.Theta0*180/math.Pi, l.DTh*180/math.Pi, l.Shape.Radius*1e3, l.Shape.Turns, l.Shape.Beta*180/math.Pi,
		l.Shape.Rho*1e3, l.Shape.Z*1e3, l.Comb, l.RLoad)
}
