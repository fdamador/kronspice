// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 5: the multipactor pump driven by a standard microwave-oven
// magnetron, evolved by the genetic algorithm (evo).
//
// Chain: mains → magnetron (2.45 GHz CW) → circulator → coupling port →
// re-entrant copper cavity whose capacitive gap is the multipactor gap →
// grid → collector. The circulator is standard practice in industrial
// 2.45 GHz systems: the magnetron must see a matched load (the LG 2M214
// sheet rates VSWR ≤ 1.1 for its quoted output, 4 absolute maximum), so
// the power the cavity does not take goes to a water dummy load.
//
// Circuit model of the cavity (kronspice twin in Build): gap capacitance
// C = ε0·A/d, L resonating it at 2.45 GHz, wall loss R0 = Q0/(ωC), the
// port R_ext = R0/β (β = coupling coefficient), driven by a Norton source
// of available power P: amplitude √(8P/R_ext). The electrons are the sheet
// of pumps.go, loading the gap through Ramo's current.
//
// The genetic algorithm (evo.Run, value genome, seeded, deterministic)
// sets gap d, gap area A, coupling β and grid transparency τ, and picks
// the surface (copper as received or conditioned, the measured yields we
// have), to collect the most electrons per watt drawn from the mains.
// Constraints: the gap must stay lumped (side ≤ λ/10, gene range) and the
// peak gap field below Kilpatrick's vacuum-spark limit.

import (
	"context"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
)

// Magnetron is a CW oven magnetron.
type Magnetron struct {
	Name       string
	F          float64 // [Hz]
	Pout       float64 // RF output [W]
	Efficiency float64 // RF out / anode input (declared)
	Filament   float64 // heater power [W]
}

// LG2M214: 2450 MHz, 900 W (860–950 W), heater 3.3 V × 10.5 A, from the
// LG 2M214 specification sheet. The sheet's efficiency is only in a
// chart; 70 % is a declared typical value for oven magnetrons, not a
// datasheet number.
var LG2M214 = Magnetron{Name: "LG 2M214 (oven magnetron)", F: 2.45e9, Pout: 900, Efficiency: 0.70, Filament: 3.3 * 10.5}

// StellantCWM100: 915 MHz, 100 kW CW, efficiency 85 %, filament 12.6 V ×
// 69 A (Stellant Systems CWM-100kW sell sheet, May 2023): an industrial
// heating magnetron.
var StellantCWM100 = Magnetron{Name: "Stellant CWM-100kW (industrial magnetron)", F: 915e6, Pout: 100e3, Efficiency: 0.85, Filament: 12.6 * 69}

// Mains is the power drawn from the supply [W].
func (m Magnetron) Mains() float64 { return m.Pout/m.Efficiency + m.Filament }

// MagPump is one design.
type MagPump struct {
	Mag     Magnetron
	Gap     float64 // [m]
	Area    float64 // [m²]
	Beta    float64 // coupling coefficient R0/R_ext
	Grid    float64 // grid transparency
	Surface em.SEY
	Q0      float64 // unloaded cavity Q (copper re-entrant cavity, declared)
	// N identical gaps share the magnetron through an ideal, isolated power
	// divider (declared lossless): each gets Pout/N available. 0 means 1.
	N int
}

func (p MagPump) gaps() float64 {
	if p.N < 1 {
		return 1
	}
	return float64(p.N)
}

// Surfaces the genome can pick (measured: Baglin et al. 2002).
var Surfaces = []em.SEY{em.CopperAsReceived, em.CopperConditioned}

// Cavity returns C, L, R0, R_ext and the Norton amplitude.
func (p MagPump) Cavity() (C, L, R0, Rext, I float64) {
	w := 2 * math.Pi * p.Mag.F
	C = em.Eps0 * p.Area / p.Gap
	L = 1 / (w * w * C)
	R0 = p.Q0 / (w * C)
	Rext = R0 / p.Beta
	I = math.Sqrt(8 * p.Mag.Pout / p.gaps() / Rext)
	return
}

// Multipactor is the gap model of pumps.go for this design: loaded Q and
// the unloaded steady-state amplitude I·(R0 ∥ R_ext); strays only.
func (p MagPump) Multipactor(steps int) Multipactor {
	_, _, R0, Rext, I := p.Cavity()
	Rl := R0 * Rext / (R0 + Rext)
	return Multipactor{Gap: p.Gap, Area: p.Area, F: p.Mag.F, Q: p.Q0 / (1 + p.Beta), V0: I * Rl,
		Grid: p.Grid, Seed: 0, SeedEvery: 49, Surface: p.Surface, Rext: Rext, Steps: steps}
}

// Build is the cavity as a kronspice circuit (no electrons): Norton
// source with the port resistance, wall resistance, L and the gap C.
func (p MagPump) Build() (*circuit.Circuit, error) {
	C, L, R0, Rext, I := p.Cavity()
	b := circuit.NewBuilder("magnetron-driven multipactor cavity")
	b.CurrentWave("IMAG", "0", "g", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, I, p.Mag.F}})
	b.Resistor("REXT", "g", "0", Rext)
	b.Resistor("RWALL", "g", "0", R0)
	b.Inductor("LCAV", "g", "0", L)
	b.Capacitor("CGAP", "g", "0", C)
	return b.BuildTransient()
}

// MagResult is the outcome for one design.
type MagResult struct {
	Design                  MagPump
	V0, PeakV, PeakField    float64 // unloaded amplitude, largest gap voltage [V], field [V/m]
	Kilpatrick              float64 // [V/m]
	Current                 float64 // average collector current [A]
	ElectronsPerJoule       float64 // collected per joule from the mains
	Mains, Delivered, Dummy float64 // [W]: mains, into the cavity, to the dummy load
	CavityHeat, Electrons   float64 // [W]: wall heat, electron impacts + collector
	Residual                float64
	N                       int // gaps sharing the magnetron
	Run                     MPResult
}

// KilpatrickField [V/m]: f[MHz] = 1.64·E²·exp(−8.5/E), E in MV/m (W. D.
// Kilpatrick, Rev. Sci. Instrum. 28, 824, 1957; form of T. J. Boyd,
// LA-UR-82-1166, 1982), as in examples/double_hull.
func KilpatrickField(f float64) float64 {
	fm := f / 1e6
	lo, hi := 0.01, 1000.0
	for i := 0; i < 100; i++ {
		e := 0.5 * (lo + hi)
		if 1.64*e*e*math.Exp(-8.5/e) < fm {
			lo = e
		} else {
			hi = e
		}
	}
	return 0.5 * (lo + hi) * 1e6
}

// Evaluate runs the design for `periods` RF periods.
func (p MagPump) Evaluate(periods, steps int) MagResult {
	m := p.Multipactor(steps)
	r := m.Run(periods, true)
	T := float64(periods) / p.Mag.F
	res := MagResult{Design: p, V0: m.V0, PeakV: r.PeakGapV, PeakField: r.PeakGapV / p.Gap, Kilpatrick: KilpatrickField(p.Mag.F), Run: r}
	n := p.gaps()
	res.N = int(n)
	res.Current = n * E * r.Collected / T
	res.Mains = p.Mag.Mains()
	delivered := r.SourceWork - r.PortOut // into the cavity, net of what the port sends back
	res.Delivered = n * delivered / T
	res.Dummy = p.Mag.Pout - res.Delivered
	res.CavityHeat = n * (r.TankHeat - r.PortOut) / T
	res.Electrons = n * (r.WallHeat + r.Collector) / T
	res.ElectronsPerJoule = n * r.Collected / (res.Mains * T)
	res.Residual = r.Residual
	return res
}

// space: gap, area, β, grid, surface (1-based index) and, for arrays, the
// number of gaps. Gap range and largest area scale with the wavelength
// (area side ≤ λ/10: 1.2 cm at 2.45 GHz, 3.3 cm at 915 MHz).
func (o EvolveOptions) space() evo.Space {
	sc := 2.45e9 / o.mag().F
	s := evo.Space{
		{Name: "gap [m]", Min: 0.1e-3 * sc, Max: 3e-3 * sc},
		{Name: "area [m²]", Min: 0.05e-4 * sc * sc, Max: 1.5e-4 * sc * sc},
		{Name: "coupling β", Min: 1e-6, Max: 1},
		{Name: "grid", Min: 0.02, Max: 0.9},
		{Name: "surface", Min: 1, Max: 2, Integer: true},
	}
	if o.MaxGaps > 1 {
		s = append(s, evo.Gene{Name: "gaps", Min: 1, Max: float64(o.MaxGaps), Integer: true})
	}
	return s
}

func (o EvolveOptions) mag() Magnetron {
	if o.Mag.F == 0 {
		return LG2M214
	}
	return o.Mag
}

func (o EvolveOptions) design(v []float64) MagPump {
	p := MagPump{Mag: o.mag(), Gap: v[0], Area: v[1], Beta: v[2], Grid: v[3], Surface: Surfaces[int(v[4])-1], Q0: 2000}
	if len(v) > 5 {
		p.N = int(v[5])
	}
	return p
}

func magDesign(v []float64) MagPump { return EvolveOptions{}.design(v) }

// EvolveOptions sizes the search.
type EvolveOptions struct {
	Population, Generations, Periods, Steps int
	Seed                                    uint64
	Mag                                     Magnetron // zero = LG2M214
	MaxGaps                                 int       // > 1 adds the number of gaps as a gene
}

// DefaultEvolve: 40 designs × 30 generations, 1500 periods (612 ns) at
// 200 steps per period for the search; the winner is re-run at 800 steps
// for 2000 periods.
func DefaultEvolve() EvolveOptions {
	return EvolveOptions{Population: 40, Generations: 30, Periods: 1500, Steps: 200, Seed: 2450}
}

// Evolution is the outcome of Evolve.
type Evolution struct {
	First, Best MagResult // best of generation 0, best overall (re-run finer)
	History     []evo.GenStats
}

// MaxGapV keeps the electron sheet non-relativistic: 10 kV gives γ =
// 1.02, so the classical sheet errs by about 2 % at most; above it the
// design is outside the model's validity and is rejected, not scored.
const MaxGapV = 10e3

// Fitness to minimize: mains watts per collected microamp. A design that
// collects less than 1 nA per gap (strays only, no multiplication) scores
// 10⁶ plus a guide toward the multipactor zone, 10⁶·(max(0, 1.5 − δ̄) +
// |ln(V0/Vz)|), where δ̄ is em.Multipactor's best mean yield per impact at
// the cavity's unloaded amplitude V0 and Vz = m·ω²·d²/(2e) marks the
// first-order zone. Past Kilpatrick's limit the design scores 10⁹ × the excess,
// and past MaxGapV 10⁸ × the excess.
func magObjectives(o EvolveOptions) []evo.Objective {
	return []evo.Objective{
		{Name: "W per µA collected", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
			p := o.design(e.Values)
			r := p.Evaluate(o.Periods, o.Steps)
			if r.PeakField > r.Kilpatrick {
				return 1e9 * r.PeakField / r.Kilpatrick, nil
			}
			if r.V0 > MaxGapV {
				return 1e8 * r.V0 / MaxGapV, nil
			}
			if r.Current < 1e-9*p.gaps() {
				g, _ := em.Multipactor(p.Gap, p.Mag.F, r.V0, 0, p.Surface)
				// distance from the first-order zone, whose voltage scales
				// as m·ω²·d²/e (Hatch & Williams); a smooth guide when the
				// yield test sees nothing
				vz := 0.5 * mE * math.Pow(2*math.Pi*p.Mag.F*p.Gap, 2) / E
				return 1e6 * (1 + math.Max(0, 1.5-g) + math.Abs(math.Log(r.V0/vz))), nil
			}
			return r.Mains / (r.Current * 1e6), nil
		}},
	}
}

// Evolve runs the genetic algorithm.
func Evolve(o EvolveOptions) (Evolution, error) {
	prob := evo.Problem{
		Build: func(v []float64) (*circuit.Circuit, error) { return o.design(v).Build() },
		Space: o.space(), Objectives: magObjectives(o),
	}
	var first []float64
	st, err := evo.Run(context.Background(), prob, evo.Options{Population: o.Population, Generations: o.Generations, Seed: o.Seed, Timeout: 600e9},
		nil, func(st *evo.State) error {
			if st.Generation == 0 {
				first = append([]float64(nil), st.Best.Values...)
			}
			return nil
		})
	if err != nil {
		return Evolution{}, err
	}
	if st.Best.Failure != "" {
		return Evolution{}, fmt.Errorf("best design failed: %s %s", st.Best.Failure, st.Best.Error)
	}
	ev := Evolution{History: st.History}
	ev.First = o.design(first).Evaluate(2000, 800)
	ev.Best = o.design(st.Best.Values).Evaluate(2000, 800)
	return ev, nil
}
