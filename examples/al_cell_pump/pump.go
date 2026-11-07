// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// electron charge [C]
const qe = 1.602176634e-19

// Pump is the resonant-driven negative Cockcroft–Walton ground pump of
// examples/photon_electron (GroundPump) with its diodes replaced by
// aluminium-water stacks (or kept as silicon diodes, for comparison):
//
//	VS (Vs at F) – RS – L – tank node "ac" – Ct to earth;
//	n stages: push capacitors from "ac", hold capacitors from earth, a
//	rectifier from each hold node to the next push node and from each push
//	node to its hold node; the last hold node charges the bank, with the
//	load across it, returned to earth.
//
// Earth lends the electrons that pile up on the bank's negative plate; the
// energy comes only from the source. The series resonator lifts a small
// drive to a large tank swing: RS = ω·L/QL with the inductor's declared
// Q (QL = 50, an air-core coil).
type Pump struct {
	Vs, F    float64 // drive amplitude [V] and frequency [Hz]
	Ct       float64 // tank capacitor [F]
	QL       float64 // inductor Q
	Stages   int
	Cstage   float64 // ladder capacitors [F]
	Cbank    float64 // [F]
	Rload    float64 // [Ω]
	Diode    bool    // silicon diodes instead of aluminium stacks
	Cell     Stack
	Pipes    *Pipes // concentric pipes in one tank instead of Cell (nil: plates)
	Cycles   int    // simulated drive cycles
	PerCycle int    // time steps per cycle
	Window   int    // the last cycles measured
}

// DefaultPump is photon_electron's ground pump (10 V, 100 kHz, 1 µF tank,
// 3 stages of 100 nF, 100 nF bank, 10 MΩ load) with QL = 50 instead of the
// tank's Q of 20, and aluminium stacks of 2 cells × 1 plate of 10 cm²,
// 5 mm gap, formed at 600 V.
func DefaultPump() Pump {
	return Pump{Vs: 10, F: 100e3, Ct: 1e-6, QL: 50, Stages: 3, Cstage: 100e-9, Cbank: 100e-9, Rload: 10e6,
		Cell: Stack{Ns: 2, Np: 1, Area: 1e-3, Gap: 5e-3, Vform: 600, Kappa: kappaDef}, Cycles: 600, PerCycle: 100, Window: 100}
}

// Build makes the circuit.
func (p Pump) Build() (*circuit.Circuit, error) {
	w := 2 * math.Pi * p.F
	L := 1 / (w * w * p.Ct)
	rs := w * L / p.QL
	b := circuit.NewBuilder("aluminium-cell ground pump")
	b.Model("dhv", "D", map[string]float64{"is": 1e-9, "n": 1.8, "rs": 0.05})
	b.VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, p.Vs, p.F}})
	b.Resistor("RS", "s", "m", rs)
	b.Inductor("LT", "m", "ac", L)
	b.Capacitor("CT", "ac", "0", p.Ct)
	rect := func(name, anode, cathode string) {
		if p.Diode {
			b.Diode(name, anode, cathode, "dhv")
			return
		}
		p.rect().Add(b, name, anode, cathode)
	}
	prevPush, prevHold := "ac", "0"
	for k := 1; k <= p.Stages; k++ {
		push, hold := fmt.Sprintf("p%d", k), fmt.Sprintf("h%d", k)
		b.Capacitor(fmt.Sprintf("CP%d", k), prevPush, push, p.Cstage)
		b.Capacitor(fmt.Sprintf("CH%d", k), prevHold, hold, p.Cstage)
		rect(fmt.Sprintf("DA%d", k), push, prevHold)
		rect(fmt.Sprintf("DB%d", k), hold, push)
		prevPush, prevHold = push, hold
	}
	b.Capacitor("CBANK", prevHold, "0", p.Cbank)
	b.Resistor("RLOAD", prevHold, "0", p.Rload)
	return b.BuildTransient()
}

// Result is one run, measured over the last Window cycles.
type Result struct {
	Out  float64 // mean bank voltage [V] (negative: electrons on the bank)
	Iout float64 // mean load current [A] (magnitude)
	Pin  float64 // mean power from the source [W]
	Pout float64 // mean power into the load [W]
	Eff  float64 // Pout/Pin in steady state: the "catalytic efficiency" (1 = no loss at all)
	// Where the input goes [W]: the resonator coil's resistance, the
	// rectifiers conducting (electrolysis onset and electrolyte), the
	// rectifiers leaking (reverse), and the rest (energy still going into
	// the capacitors, numerical).
	LossCoil, LossForward, LossLeak, Rest float64
	CellRev                               []float64 // pipes: largest reverse on each annulus / its forming voltage, over all tanks
	LossSeal                              float64   // series pipes: through the leaks around the sealed pipe ends [W]
	ElPerJ                                float64   // electrons delivered per joule from the source
	Ideal                                 float64   // electrons per joule if nothing were lost: 1/(e·|Out|)
	TankPeak                              float64   // [V]
	RevRatio                              float64   // largest reverse voltage on any rectifier / its holding voltage
	Drift                                 float64   // relative change of the bank voltage over the window (steady state check)
	Ledger                                float64   // relative integration mismatch of the whole run
	Time                                  []float64
	Bank                                  []float64
}

// ErrNoPump: the bank did not charge.
var ErrNoPump = errors.New("the bank did not charge")

// Run simulates and measures the pump.
func (p Pump) Run(opts analysis.Options) (Result, error) {
	var r Result
	c, err := p.Build()
	if err != nil {
		return r, err
	}
	T := 1 / p.F
	dt := T / float64(p.PerCycle)
	stop := float64(p.Cycles) * T
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: stop, MaxStep: dt}, analysis.TranOptions{Options: opts, MaxSteps: 5_000_000})
	if err != nil {
		return r, err
	}
	l := wf.Ledger
	r.Ledger = math.Abs(l.IntegrationMismatch) / math.Max(l.EnergyIn, 1e-300)
	t := wf.Time
	out, _ := wf.NodeVoltage(fmt.Sprintf("h%d", p.Stages))
	ac, _ := wf.NodeVoltage("ac")
	vs, _ := wf.NodeVoltage("s")
	is, _ := wf.ElementCurrent("VS")
	r.Time, r.Bank = t, out
	t0 := stop - float64(p.Window)*T
	var ein, eload, span float64
	var vStart, vEnd float64
	first := true
	for i := 1; i < len(t); i++ {
		if t[i] <= t0 {
			continue
		}
		h := t[i] - t[i-1]
		if first {
			vStart, first = out[i-1], false
		}
		// The source delivers −v·i (its current is counted into its + node).
		ein += h * 0.5 * (-vs[i]*is[i] - vs[i-1]*is[i-1])
		eload += h * 0.5 * (out[i]*out[i] + out[i-1]*out[i-1]) / p.Rload
		r.Out += h * 0.5 * (out[i] + out[i-1])
		span += h
		r.TankPeak = math.Max(r.TankPeak, math.Abs(ac[i]))
		vEnd = out[i]
	}
	r.Out /= span
	r.Iout = math.Abs(r.Out) / p.Rload
	r.Pin = ein / span
	r.Pout = eload / span
	r.Drift = math.Abs(vEnd-vStart) / math.Max(math.Abs(r.Out), 1e-12)
	if r.Pin > 0 {
		r.Eff = r.Pout / r.Pin
		r.ElPerJ = r.Iout / qe / r.Pin
	}
	if math.Abs(r.Out) > 0 {
		r.Ideal = 1 / (qe * math.Abs(r.Out))
	}
	// Loss breakdown and the largest reverse voltage on any oxide over the
	// window (element v·i, trapezoidal).
	m := meter{
		node: func(n string) []float64 {
			if n == "0" {
				return make([]float64, len(t))
			}
			v, _ := wf.NodeVoltage(n)
			return v
		},
		cur: func(e string) []float64 {
			i, ok := wf.ElementCurrent(e)
			if !ok {
				return nil
			}
			return i
		},
		avg: func(f func(i int) float64) float64 {
			var acc float64
			for i := 1; i < len(t); i++ {
				if t[i] > t0 {
					acc += (t[i] - t[i-1]) * 0.5 * (f(i) + f(i-1))
				}
			}
			return acc / span
		},
		max: func(f func(i int) float64) float64 {
			x := math.Inf(-1)
			for i := range t {
				if t[i] > t0 {
					x = math.Max(x, f(i))
				}
			}
			return x
		},
	}
	if irs, ok := wf.ElementCurrent("RS"); ok {
		w := 2 * math.Pi * p.F
		rs := w * (1 / (w * w * p.Ct)) / p.QL
		r.LossCoil = m.avg(func(i int) float64 { return irs[i] * irs[i] * rs })
	}
	prevHold := "0"
	for k := 1; k <= p.Stages; k++ {
		push, hold2 := fmt.Sprintf("p%d", k), fmt.Sprintf("h%d", k)
		for _, d := range []struct{ name, a, c string }{{fmt.Sprintf("DA%d", k), push, prevHold}, {fmt.Sprintf("DB%d", k), hold2, push}} {
			var f, l, sl, rv float64
			if p.Diode {
				f, l = twoTerm(m, d.name, d.a, d.c)
			} else {
				f, l, sl, rv = p.rect().measure(m, d.name, d.a, d.c)
			}
			r.LossForward += f
			r.LossLeak += l
			r.LossSeal += sl
			r.RevRatio = math.Max(r.RevRatio, rv)
			if p.Pipes != nil && !p.Diode {
				cr := p.Pipes.cellRev(m, d.name, d.a, d.c)
				if r.CellRev == nil {
					r.CellRev = make([]float64, len(cr))
				}
				for j, x := range cr {
					r.CellRev[j] = math.Max(r.CellRev[j], x)
				}
			}
		}
		prevHold = hold2
	}
	r.Rest = r.Pin - r.Pout - r.LossCoil - r.LossForward - r.LossLeak - r.LossSeal
	if math.Abs(r.Out) < 1 {
		return r, ErrNoPump
	}
	return r, nil
}

// rectifier is a valve the pump can be built with.
type rectifier interface {
	Add(b *circuit.Builder, name, anode, cathode string)
	Hold() float64
	measure(m meter, name, anode, cathode string) (fwd, leak, seal, rev float64)
}

func (p Pump) rect() rectifier {
	if p.Pipes != nil {
		return *p.Pipes
	}
	return p.Cell
}

func b2f(b bool) float64 {
	if b {
		return 1
	}
	return 0
}
