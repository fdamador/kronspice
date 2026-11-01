// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 3: tune the flash circuit to the coil, as in a spark-gap Tesla coil.
//
// In Part 2 the flash capacitor (10 nF) rings with the 6 µH primary at
// about 650 kHz, 16 times the coil's 41.6 kHz, so almost nothing reaches
// the coil. Tuning means Lp·Cp = Ls·Cs: the primary circuit rings at the
// coil's own frequency, and the energy swings from primary to coil and
// back at the beat rate. With coupling k the first complete transfer is
// after half a beat, t = 1/(2·k·f1) (two coupled tuned circuits; e.g.
// F. E. Terman, Radio Engineers' Handbook, 1943, §3, coupled circuits).
//
// The spark gap matters then:
//   - an ordinary arc keeps conducting through the current zeros until the
//     current has stayed small for its de-ionization time (declared
//     20 µs), so the energy swings back into the primary and burns in the
//     arc;
//   - a quenched gap (Tesla's quenched gap; Wien's quenched spark) opens at
//     the first current zero after the transfer and leaves the energy in
//     the coil.
//
// Both water gaps of Part 2 are run untuned, tuned with an arc and tuned
// and quenched; then the genetic algorithm tunes the mist gap's capacitor,
// primary turns and quench time for the largest share of each flash's
// energy delivered to the coil.

import (
	"context"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
)

// WithTurns returns the coil with an n-turn primary (radius 6 cm, pitch
// 2 mm, 1 mm wire, as in examples/pulsed_bifilar). Lp from the turns'
// self and mutual inductances; M scaled from the measured 5-turn value
// (28.1 µH) in proportion to the turns (declared: the turns sit close
// together at the coil's base).
func (c Coil) WithTurns(n int) Coil {
	const (
		radius = 0.06
		pitch  = 2e-3
		wire   = 1e-3
	)
	lp := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				lp += em.Mu0 * radius * (math.Log(8*radius/wire) - 2)
			} else {
				lp += em.MutualCircles(radius, radius, float64(i-j)*pitch)
			}
		}
	}
	c.Lp = lp
	c.M = 28.1e-6 * float64(n) / 5
	return c
}

// TunedC is the flash capacitor that rings with the primary at f1.
func (c Coil) TunedC() float64 {
	w := 2 * math.Pi * c.F1
	return 1 / (w * w * c.Lp)
}

// K is the coupling coefficient.
func (c Coil) K() float64 { return c.M / math.Sqrt(c.Lp*c.Ls) }

// Transfer is the time of the first full transfer, 1/(2·k·f1).
func (c Coil) Transfer() float64 { return 1 / (2 * c.K() * c.F1) }

// TunedRow is one Part 3 run.
type TunedRow struct {
	Label                     string
	C                         float64
	Turns                     int
	Run                       CoilRun
	PerFlash, Share, OfSupply float64 // energy to the coil per flash [J]; of flash energy; of supply energy
}

// tunedCase returns the gap case retuned: capacitor cp, bipolar spark
// with the given quench, and a charging resistor large enough that the
// supply cannot hold the arc by itself (Vs/R ≤ Iext/2; a smaller one gives
// a power-follow arc that never goes out), and at least 50 ms·cp.
func tunedCase(base Case, cp float64, quenchAt float64) Case {
	c := base
	c.C = cp
	c.R = math.Max(0.05/cp, 2*c.Vs/c.Lamp.Iext)
	c.Lamp.Bipolar = true
	c.Lamp.Quench = 20e-6
	c.Lamp.QuenchAt = quenchAt
	return c
}

func row(label string, c Case, coil Coil, turns int, seconds float64) (TunedRow, error) {
	r, err := c.RunCoil(coil, seconds, 7)
	if err != nil {
		return TunedRow{}, err
	}
	tr := TunedRow{Label: label, C: c.C, Turns: turns, Run: r}
	if r.Flashes > 0 {
		tr.PerFlash = r.CoilEnergy / float64(r.Flashes)
		tr.Share = r.CoilEnergy / r.FlashEnergy
		tr.OfSupply = r.CoilEnergy / r.Supply
	}
	return tr, nil
}

// TunedComparison runs the liquid and mist gaps untuned, tuned with an
// ordinary arc, and tuned with a quenched gap.
func TunedComparison() ([]TunedRow, error) {
	coil := DefaultCoil()
	var rows []TunedRow
	for _, tag := range []string{"liquid", "mist"} {
		var base Case
		for _, c := range Cases() {
			if c.Tag == tag {
				base = c
			}
		}
		r0, err := row(base.Name+": untuned (10 nF, Part 2)", base, coil, 5, 0.3)
		cpT := coil.TunedC()
		span := 6 * tunedCase(base, cpT, 0).R * cpT
		if err != nil {
			return nil, err
		}
		cp := coil.TunedC()
		r1, err := row(base.Name+": tuned, ordinary arc", tunedCase(base, cp, 0), coil, 5, span)
		if err != nil {
			return nil, err
		}
		r2, err := row(base.Name+": tuned, quenched gap", tunedCase(base, cp, coil.Transfer()), coil, 5, span)
		if err != nil {
			return nil, err
		}
		rows = append(rows, r0, r1, r2)
	}
	return rows, nil
}

// --- evolution of the mist gap's tuning -------------------------------------

var tuneSpace = evo.Space{
	{Name: "flash capacitor [F]", Min: 1e-9, Max: 20e-6},
	{Name: "primary turns", Min: 1, Max: 20, Integer: true},
	{Name: "quench after [beat halves]", Min: 0.2, Max: 3},
}

func tuneDesign(v []float64) (Case, Coil) {
	var base Case
	for _, c := range Cases() {
		if c.Tag == "mist" {
			base = c
		}
	}
	coil := DefaultCoil().WithTurns(int(v[1]))
	return tunedCase(base, v[0], v[2]*coil.Transfer()), coil
}

// TuneResult is the outcome of the evolution.
type TuneResult struct {
	Best    TunedRow
	Quench  float64 // in beat halves
	History []evo.GenStats
}

// EvolveTuning lets the genetic algorithm (evo.Run, seeded) pick the
// flash capacitor, primary turns and quench time for the mist gap,
// maximizing the share of each flash's energy that reaches the coil
// (fitness 1 − share; at least 3 flashes in six charging time constants).
func EvolveTuning(pop, gens int, seed uint64) (TuneResult, error) {
	prob := evo.Problem{
		Build: func(v []float64) (*circuit.Circuit, error) {
			_, coil := tuneDesign(v)
			b := circuit.NewBuilder("tuned primary and bifilar coil")
			b.CurrentWave("IP", "0", "p", circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 0.2e-6, 0.2e-6, 2e-6, 1}})
			b.Inductor("LP", "p", "0", coil.Lp)
			b.Inductor("LS", "s", "0", coil.Ls)
			b.Capacitor("CS", "s", "0", coil.Cs)
			b.Resistor("RS", "s", "0", coil.Rs)
			b.Mutual("K1", "LP", "LS", coil.K())
			return b.BuildTransient()
		},
		Space: tuneSpace,
		Objectives: []evo.Objective{{Name: "1 − share of flash energy to the coil", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
			c, coil := tuneDesign(e.Values)
			r, err := c.RunCoil(coil, 6*c.R*c.C, 7)
			if err != nil {
				return 0, err
			}
			if r.Flashes < 3 {
				return 2 - float64(r.Flashes)/3, nil
			}
			return 1 - r.CoilEnergy/r.FlashEnergy, nil
		}}},
	}
	st, err := evo.Run(context.Background(), prob, evo.Options{Population: pop, Generations: gens, Seed: seed, Timeout: 600e9}, nil, nil)
	if err != nil {
		return TuneResult{}, err
	}
	if st.Best.Failure != "" {
		return TuneResult{}, fmt.Errorf("best design failed: %s %s", st.Best.Failure, st.Best.Error)
	}
	c, coil := tuneDesign(st.Best.Values)
	r, err := row("evolved", c, coil, int(st.Best.Values[1]), 6*c.R*c.C)
	if err != nil {
		return TuneResult{}, err
	}
	return TuneResult{Best: r, Quench: st.Best.Values[2], History: st.History}, nil
}
