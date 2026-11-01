// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo_test

import (
	"context"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
)

// Evolve the resistor of an RC low-pass (capacitor fixed at 100 nF) for a
// 1 kHz cutoff, with resistors restricted to the E24 series.
func ExampleRun() {
	tpl, _ := circuit.NewBuilder("rc").
		Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "in", "0").
		Resistor("R1", "in", "out", 10e3).
		Capacitor("C1", "out", "0", 100e-9).
		Build()
	spec := analysis.ACSpec{Variation: "dec", Points: 200, FStart: 10, FStop: 1e5}
	p := evo.Problem{
		Template:   tpl,
		Space:      evo.Space{{Element: "R1", Min: 10, Max: 1e6, Series: "E24"}},
		Objectives: []evo.Objective{evo.Cutoff("out", "0", spec, 1000, true)},
	}
	st, err := evo.Run(context.Background(), p, evo.Options{Population: 20, Generations: 20, Seed: 1}, nil, nil)
	if err != nil {
		fmt.Println(err)
		return
	}
	r := st.Best.Values[0]
	fmt.Printf("R1 = %g Ω, fc = %.0f Hz\n", r, 1/(2*math.Pi*r*100e-9))
	// Output: R1 = 1600 Ω, fc = 995 Hz
}
