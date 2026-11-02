// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Twin is the electrical side of a stroke solved in kronspice: the leaf
// pair (1 m²) as a capacitor varying with the gap the stroke followed,
// C(t) = 4ε0/(L − a(t)) (VaryingCapacitor), fed by the servo voltage V(t)
// through the supply resistance.
type Twin struct {
	ToSupply float64 // energy into the supply [J]
	Charge   float64 // charge returned to the supply [C]
	PumpWork float64 // work done on the capacitor by the moving leaves, −∫½v²·dC/dt dt [J]
	Heat     float64
	Mismatch float64 // kronspice ledger mismatch ÷ largest flow
}

// KronTwin runs the twin for stroke s of leaf l.
func KronTwin(l Leaf, s Stroke) (Twin, error) {
	var tw Twin
	n := len(s.Time)
	c0 := l.Cap(s.Gap[0])
	scale := make([]float64, 0, 2*n)
	volt := make([]float64, 0, 2*n)
	for i := 0; i < n; i++ {
		scale = append(scale, s.Time[i], l.Cap(s.Gap[i])/c0)
		volt = append(volt, s.Time[i], s.Volt[i])
	}
	b := circuit.NewBuilder("Forward leaf pair").
		VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WavePWL, Args: volt}).
		Resistor("RS", "s", "p", l.R).
		VaryingCapacitor("CP", "p", "0", c0, circuit.Waveform{Kind: circuit.WavePWL, Args: scale})
	ckt, err := b.BuildTransient()
	if err != nil {
		return tw, err
	}
	stop := s.Time[n-1]
	step := stop / float64(4*n)
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return tw, fmt.Errorf("kronspice: %w", err)
	}
	vs, _ := wf.NodeVoltage("s")
	vp, _ := wf.NodeVoltage("p")
	isrc, _ := wf.ElementCurrent("VS")
	// interpolated C(t) and its rate for the pump work
	cAt := func(t float64) float64 {
		k := 0
		for k+1 < n && s.Time[k+1] < t {
			k++
		}
		if k+1 >= n {
			return l.Cap(s.Gap[n-1])
		}
		f := (t - s.Time[k]) / (s.Time[k+1] - s.Time[k])
		return l.Cap(s.Gap[k]) + f*(l.Cap(s.Gap[k+1])-l.Cap(s.Gap[k]))
	}
	for k := 1; k < len(wf.Time); k++ {
		dt := wf.Time[k] - wf.Time[k-1]
		// SPICE sign: the source current flows from + through the source to −;
		// power into the source is v·i.
		tw.ToSupply += 0.5 * dt * (vs[k]*isrc[k] + vs[k-1]*isrc[k-1])
		tw.Charge += 0.5 * dt * (isrc[k] + isrc[k-1])
		ir := (vs[k] - vp[k]) / l.R
		tw.Heat += l.R * ir * ir * dt
		dc := cAt(wf.Time[k]) - cAt(wf.Time[k-1])
		vm := 0.5 * (vp[k] + vp[k-1])
		tw.PumpWork -= 0.5 * vm * vm * dc
	}
	led := wf.Ledger
	tw.Mismatch = math.Abs(led.IntegrationMismatch) / math.Max(math.Abs(tw.ToSupply), 1e-300)
	return tw, nil
}
