// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Diode is anything that can sit in a rectifier arm between anode a and
// cathode k.
type Diode interface {
	Add(b *circuit.Builder, name, a, k string)
}

// Arm is n identical diodes in series.
type Arm struct {
	D Diode
	N int
}

func (r Arm) add(b *circuit.Builder, name, a, k string) {
	n := max(r.N, 1)
	prev := a
	for i := 0; i < n; i++ {
		next := k
		if i < n-1 {
			next = fmt.Sprintf("%s_%d", name, i)
		}
		r.D.Add(b, fmt.Sprintf("%s%d", name, i), prev, next)
		prev = next
	}
}

// Rectifier is an AC source (amplitude Vp, frequency F, source resistance
// Rs) feeding a half-wave rectifier or a four-arm bridge into a
// capacitor-filtered resistive load.
type Rectifier struct {
	Name         string
	Bridge       bool
	Arm          Arm
	Vp, F, Rs    float64
	Cload, Rload float64
	Cycles, Skip int // cycles simulated; leading cycles left out of the averages
}

// RectResult is a rectifier run.
type RectResult struct {
	Vdc, Ripple            float64 // mean load voltage, peak-to-peak ripple [V]
	Idc                    float64 // mean load current [A]
	Pload, Psource         float64 // [W] over the averaged cycles
	Efficiency             float64
	ReversePeak            float64 // largest reverse voltage across one arm [V]
	EnergyIn, Heat, Stored float64 // whole run [J]
	Mismatch               float64
	Time, Vout             []float64
	Forward, Backward      int // spark-gap strikes in the forward and the reverse direction
}

// Build makes the circuit; the load sits between nodes "p" and "n" (n is
// ground for the half-wave rectifier).
func (r Rectifier) Build() (*circuit.Circuit, error) {
	b := circuit.NewBuilder(r.Name)
	b.VoltageWave("VAC", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, r.Vp, r.F}})
	b.Resistor("RS", "s", "ac1", r.Rs)
	if r.Bridge {
		r.Arm.add(b, "DA", "ac1", "p")
		r.Arm.add(b, "DB", "0", "p")
		r.Arm.add(b, "DC", "n", "ac1")
		r.Arm.add(b, "DD", "n", "0")
		b.Resistor("RREF", "n", "0", 1e9) // keeps the floating output referenced
	} else {
		r.Arm.add(b, "DA", "ac1", "p")
	}
	out := "0"
	if r.Bridge {
		out = "n"
	}
	if r.Cload > 0 {
		b.Capacitor("CL", "p", out, r.Cload)
	}
	b.Resistor("RL", "p", out, r.Rload)
	return b.BuildTransient()
}

// Run simulates the rectifier.
func (r Rectifier) Run() (RectResult, error) {
	ckt, err := r.Build()
	if err != nil {
		return RectResult{}, err
	}
	per := 1 / r.F
	step := per / 400
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: float64(r.Cycles) * per, MaxStep: step}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return RectResult{}, err
	}
	vp, _ := wv.NodeVoltage("p")
	vn, _ := wv.NodeVoltage("n")
	if !r.Bridge {
		vn = make([]float64, len(vp))
	}
	vs, _ := wv.NodeVoltage("s")
	va, _ := wv.NodeVoltage("ac1")
	var res RectResult
	t0 := float64(r.Skip) * per
	lo, hi := math.Inf(1), math.Inf(-1)
	var T, sv, spl, sps float64
	for i := 1; i < len(wv.Time); i++ {
		t, dt := wv.Time[i], wv.Time[i]-wv.Time[i-1]
		vo := vp[i] - vn[i]
		res.Time = append(res.Time, t)
		res.Vout = append(res.Vout, vo)
		// reverse voltage on arm A: cathode p above anode ac1
		res.ReversePeak = math.Max(res.ReversePeak, vp[i]-va[i])
		if t < t0 {
			continue
		}
		lo, hi = math.Min(lo, vo), math.Max(hi, vo)
		is := (vs[i] - va[i]) / r.Rs
		T += dt
		sv += vo * dt
		spl += vo * vo / r.Rload * dt
		sps += vs[i] * is * dt
	}
	res.Vdc = sv / T
	res.Ripple = hi - lo
	res.Idc = res.Vdc / r.Rload
	res.Pload, res.Psource = spl/T, sps/T
	res.Efficiency = res.Pload / res.Psource
	l := wv.Ledger
	res.EnergyIn, res.Heat, res.Stored, res.Mismatch = l.EnergyIn, l.EnergyDissipated, l.StoredEnd-l.StoredStart, l.IntegrationMismatch
	return res, nil
}

// IVTwin sweeps a cell from −vmax to vmax in kronspice (.DC) and returns
// the largest difference from Cell.Current, relative to the largest current.
func IVTwin(c Cell, vmax float64) (float64, error) {
	b := circuit.NewBuilder("cell I–V")
	b.VoltageSource("VS", "a", "0", 0)
	b.NonlinearVC("X", "a", "0", c.curve(20*c.Vspark))
	ckt, err := b.Build()
	if err != nil {
		return 0, err
	}
	res, err := analysis.DC(ckt, []analysis.DCSweep{{Source: "VS", Start: -vmax, Stop: vmax, Step: vmax / 100}}, analysis.Options{})
	if err != nil {
		return 0, err
	}
	worst, top := 0.0, 0.0
	for k, v := range res.Values[0] {
		var i float64
		for e, n := range res.ElementNames {
			if n == "vs" || n == "VS" {
				i = -res.I[e][k] // source current flows out of its + terminal into the cell
			}
		}
		want := c.Current(v)
		top = math.Max(top, math.Abs(want))
		worst = math.Max(worst, math.Abs(i-want))
	}
	return worst / top, nil
}

// Leak is a plain resistor arm (a spark gap below breakdown).
type Leak struct{ R float64 }

// Add puts the resistor between a and k.
func (l Leak) Add(b *circuit.Builder, name, a, k string) { b.Resistor(name, a, k, l.R) }
