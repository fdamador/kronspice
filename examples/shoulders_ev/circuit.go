// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Bench is the measuring circuit both runs share: the anode (witness
// plate or collector) has capacitance CA to ground and feeds a 50 Ω scope
// line. Declared values.
type Bench struct {
	CA    float64 // anode capacitance [F]
	Scope float64 // scope termination [Ω]
}

// DefaultBench is a 1 pF anode into 50 Ω.
func DefaultBench() Bench { return Bench{CA: 1e-12, Scope: 50} }

// Shot is one kronspice run seen at the scope.
type Shot struct {
	Time, Volt []float64 // scope voltage
	Charge     float64   // charge through the scope [C]
	ToScope    float64   // energy into the scope [J]
	EnergyIn   float64   // energy from independent sources [J]
	Stored0    float64   // energy stored at t = 0 [J]
	Heat       float64   // all heat [J]
	Mismatch   float64   // ledger mismatch ÷ largest energy
	PeakV      float64
}

// Source runs Shoulders' smallest EV source (2-7): a cathode of
// capacitance cs charged to v0, discharging across the gap (an arc of
// resistance rg in series with inductance lg; declared) into the bench.
// Everything that reaches the anode, EVs included, has to come from the
// charge cs·v0 and the energy ½cs·v0².
func (b Bench) Source(cs, v0, rg, lg, stop float64) (Shot, error) {
	ckt, err := circuit.NewBuilder("EV source").
		CapacitorIC("CS", "k", "0", cs, v0).
		Resistor("RG", "k", "g", rg).
		Inductor("LG", "g", "a", lg).
		Capacitor("CA", "a", "0", b.CA).
		Resistor("RSCOPE", "a", "0", b.Scope).
		BuildTransient()
	if err != nil {
		return Shot{}, err
	}
	return b.run(ckt, stop)
}

// Strike runs an EV of n electrons arriving at the anode in a triangular
// pulse lasting width: an ideal current source carries the charge in. The
// ledger then shows what energy the arriving charge must bring.
func (b Bench) Strike(n, width, stop float64) (Shot, error) {
	q := n * Qe
	ip := 2 * q / width
	w := circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, width / 2, ip, width, 0, stop, 0}}
	ckt, err := circuit.NewBuilder("EV strike").
		CurrentWave("IEV", "0", "a", w).
		Capacitor("CA", "a", "0", b.CA).
		Resistor("RSCOPE", "a", "0", b.Scope).
		BuildTransient()
	if err != nil {
		return Shot{}, err
	}
	return b.run(ckt, stop)
}

func (b Bench) run(ckt *circuit.Circuit, stop float64) (Shot, error) {
	step := stop / 4000
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step, UIC: true}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return Shot{}, fmt.Errorf("kronspice: %w", err)
	}
	va, _ := wf.NodeVoltage("a")
	s := Shot{Time: wf.Time, Volt: va}
	for k := 1; k < len(wf.Time); k++ {
		dt := wf.Time[k] - wf.Time[k-1]
		s.Charge += 0.5 * dt * (va[k] + va[k-1]) / b.Scope
		s.ToScope += 0.5 * dt * (va[k]*va[k] + va[k-1]*va[k-1]) / b.Scope
	}
	for _, v := range va {
		s.PeakV = math.Max(s.PeakV, math.Abs(v))
	}
	l := wf.Ledger
	s.EnergyIn, s.Stored0, s.Heat = l.EnergyIn, l.StoredStart, l.EnergyDissipated
	big := math.Max(math.Max(l.EnergyIn, l.StoredStart), 1e-300)
	s.Mismatch = math.Abs(l.IntegrationMismatch) / big
	return s, nil
}
