// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/kron"
)

// FromTorsion builds kronspice's machine model directly from Einstein's
// torsion of the stationary-coil 4-bein: the constant inductances L̄
// (Hoffmann eq. 29) and the speed matrix
//
//	G_λν = 2·(Ω_θλν + Ω_νλθ),   Ω_λμν = L̄_νρ·Λ^ρ_λμ
//
// (the speed-voltage part of Hoffmann's eq. 5′, written with Einstein's
// frame torsion Λ). No inductance derivative and no kron.RotatingFrame is
// used — only the frame and its torsion.
func FromTorsion(th float64) *circuit.MachineSpec {
	f := Frame{}
	lb := f.NewMetric(th)
	A := f.Anholonomy(th)
	var om [4][4][4]float64
	for l := 0; l < 4; l++ {
		for m := 0; m < 4; m++ {
			for n := 0; n < 4; n++ {
				for r := 0; r < 4; r++ {
					om[l][m][n] += lb[n][r] * A[r][l][m]
				}
			}
		}
	}
	s := circuit.NewMachineSpec(3)
	for l := 0; l < 3; l++ {
		for n := 0; n < 3; n++ {
			s.L0[l][n] = lb[l][n]
			s.G[l][n] = 2 * (om[3][l][n] + om[n][l][3])
			if math.Abs(s.G[l][n]) < 1e-15 {
				s.G[l][n] = 0
			}
		}
	}
	return s
}

// Run is a start-up of the motor (Hoffmann's single-phase motor, 120 V,
// 60 Hz, stator 2 Ω, rotor coils 3 Ω, J, friction 10⁻⁴ N·m·s), pushed to
// w0 [rad/s], for the given seconds.
type Run struct {
	Name           string
	Time, Speed    []float64
	Final          float64
	Supplied, Heat float64
	Mismatch       float64
}

// Start runs spec s in kronspice.
func Start(name string, s *circuit.MachineSpec, w0, seconds float64) (Run, error) {
	ckt, err := circuit.NewBuilder("Hoffmann motor").Domain(exchange.Mechanical, "w").
		MachineModel("mm", s).
		VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 120 * math.Sqrt2, 60}}).
		Resistor("RS", "s", "i1", 2).
		Resistor("R2", "r2", "0", 3).Resistor("R3", "r3", "0", 3).
		Machine("N1", "mm", [][2]string{{"i1", "0"}, {"r2", "0"}, {"r3", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", J, w0).Resistor("RF", "w", "0", 1e4).BuildTransient()
	if err != nil {
		return Run{}, err
	}
	step := 1 / (60.0 * 200)
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: seconds, MaxStep: step, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}, MaxSteps: 5_000_000})
	if err != nil {
		return Run{}, fmt.Errorf("%s: %w", name, err)
	}
	w, _ := wf.NodeVoltage("w")
	l := wf.Ledger
	return Run{Name: name, Time: wf.Time, Speed: w, Final: w[len(w)-1], Supplied: l.EnergyIn - l.EnergyOut, Heat: l.EnergyDissipated,
		Mismatch: l.IntegrationMismatch / (math.Abs(l.EnergyIn) + math.Abs(l.StoredEnd) + math.Abs(l.StoredStart))}, nil
}

// Specs returns the three machines: frame I, Kron's frame II(b)
// (kron.RotatingFrame), and the one built from Einstein's torsion; plus
// the same constant inductances with the torsion removed.
func Specs() (frameI, kronII, einstein, noTorsion *circuit.MachineSpec, err error) {
	frameI = Motor()
	kronII, err = kron.RotatingFrame(frameI, Frame{}.Moving())
	if err != nil {
		return
	}
	einstein = FromTorsion(0.4)
	noTorsion = circuit.NewMachineSpec(3)
	for i := range einstein.L0 {
		copy(noTorsion.L0[i], einstein.L0[i])
	}
	return
}
