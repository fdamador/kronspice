// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 4: plates in parallel on the spark-gap side.
//
//  (a) Several spark gaps in parallel across one capacitor. The first gap to
//      break down collapses the voltage, so only one arc carries the flash;
//      what N gaps do change is the number of places a starting electron
//      can appear, N times the rate: shorter, steadier waits.
//  (b) More capacitor plates in parallel (n sections of the evolved
//      capacitor). Each flash holds n times the energy and the spark current
//      rises; to stay tuned to the coil the primary needs fewer turns
//      (Lp·Cp = Ls·Cs), chosen here from 1–20 turns.

import "math"

// ParallelRow is one Part 4 run.
type ParallelRow struct {
	N, Turns                  int
	C                         float64
	MeanStrike, SDStrike      float64
	MeanInterval, Jitter      float64
	PeakIp, PerFlash, ToCoil  float64 // A, J per flash, J per flash into the coil
	Share, Residual, Detuning float64
}

func parallelRun(c Case, coil Coil, n, turns int) (ParallelRow, error) {
	r, err := c.RunCoil(coil, 6*c.R*c.C, 7)
	if err != nil {
		return ParallelRow{}, err
	}
	row := ParallelRow{N: n, Turns: turns, C: c.C, MeanStrike: r.MeanStrike, SDStrike: r.SDStrike,
		MeanInterval: r.MeanInterval, Jitter: r.JitterInterval, PeakIp: r.PeakIp, Residual: r.Residual}
	if r.Flashes > 0 {
		row.PerFlash = r.FlashEnergy / float64(r.Flashes)
		row.ToCoil = r.CoilEnergy / float64(r.Flashes)
		row.Share = r.CoilEnergy / r.FlashEnergy
	}
	row.Detuning = 1/(2*math.Pi*math.Sqrt(coil.Lp*c.C))/coil.F1 - 1
	return row, nil
}

// ParallelGaps runs N = 1, 2, 4, 8 mist gaps across the evolved tuned
// capacitor (design: capacitor, turns, quench in beat halves).
func ParallelGaps(cp float64, turns int, quench float64) ([]ParallelRow, error) {
	var rows []ParallelRow
	for _, n := range []int{1, 2, 4, 8} {
		c, coil := tuneDesign([]float64{cp, float64(turns), quench})
		c.Lamp.Seeds *= float64(n)
		r, err := parallelRun(c, coil, n, turns)
		if err != nil {
			return nil, err
		}
		rows = append(rows, r)
	}
	return rows, nil
}

// ParallelPlates runs n = 1, 2, 4, 8 capacitor sections of cp in parallel,
// each with the primary turn count that best keeps the tuning.
func ParallelPlates(cp float64, quench float64) ([]ParallelRow, error) {
	var rows []ParallelRow
	for _, n := range []int{1, 2, 4, 8} {
		C := cp * float64(n)
		best, bt := math.Inf(1), 1
		for t := 1; t <= 20; t++ {
			k := DefaultCoil().WithTurns(t)
			if d := math.Abs(math.Log(k.Lp * C * math.Pow(2*math.Pi*k.F1, 2))); d < best {
				best, bt = d, t
			}
		}
		c, coil := tuneDesign([]float64{C, float64(bt), quench})
		r, err := parallelRun(c, coil, n, bt)
		if err != nil {
			return nil, err
		}
		rows = append(rows, r)
	}
	return rows, nil
}
