// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/material"
)

// The λ/10 limit: a 1 m conductor is at exactly λ/10 at c/10 m ≈ 30 MHz;
// above that it is an error naming the part, below it passes; a
// dielectric shortens λ by √εr; WarnOnly reports without failing.
func TestQuasiStaticLimit(t *testing.T) {
	c, err := circuit.NewBuilder("line").VoltageSource("V1", "a", "0", 1).
		Resistor("R1", "a", "b", 50).Resistor("R2", "b", "0", 50).Build()
	if err != nil {
		t.Fatal(err)
	}
	f10 := SpeedOfLight / 10 // d/λ = 0.1 for d = 1 m
	if rep, err := CheckQuasiStatic(c, map[string]float64{"R1": 1}, 0.05, 0.99*f10, SizeOptions{}); err != nil || rep.Entries[0].Name != "R1" || math.Abs(rep.Entries[0].Ratio-0.099) > 1e-12 {
		t.Errorf("below the limit: %+v, %v", rep, err)
	}
	_, err = CheckQuasiStatic(c, map[string]float64{"R1": 1}, 0, 1.01*f10, SizeOptions{})
	if !errors.Is(err, ErrNotQuasiStatic) || !strings.Contains(err.Error(), "R1") {
		t.Errorf("above the limit: %v", err)
	}
	if _, err := CheckQuasiStatic(c, map[string]float64{"R1": 1}, 0, 0.6*f10, SizeOptions{EpsR: 4}); !errors.Is(err, ErrNotQuasiStatic) {
		t.Errorf("εr = 4 halves λ: 0.6·f10 should fail, got %v", err)
	}
	if rep, err := CheckQuasiStatic(c, nil, 2, f10, SizeOptions{WarnOnly: true}); err != nil || !rep.Entries[0].TooLong {
		t.Errorf("WarnOnly: %+v, %v", rep, err)
	}
	if _, err := CheckQuasiStatic(c, map[string]float64{"RX": 1}, 0, f10, SizeOptions{}); err == nil {
		t.Error("unknown element accepted")
	}
}

// Source bandwidth: SIN at its frequency, PULSE and PWL edges at the knee
// 0.5/t_r, EXP at the knee of 2.2τ, defaulted PULSE edges at the time step.
func TestSourceBandwidth(t *testing.T) {
	for _, tc := range []struct {
		w    circuit.Waveform
		edge float64
		want float64
	}{
		{circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 1, 1e3}}, 0, 1e3},
		{circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-9, 2e-9, 1e-6, 2e-6}}, 0, 0.5e9},
		{circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1}}, 1e-8, 0.5e8},
		{circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{0, 1, 0, 1e-6, 1e-3, 1e-5}}, 0, 0.5 / 2.2e-6},
		{circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-6, 1, 1e-3, 1, 1.0001e-3, 0}}, 0, 0.5 / 1e-7},
	} {
		c, err := circuit.NewBuilder("src").VoltageWave("V1", "a", "0", tc.w).Resistor("R1", "a", "0", 1).Build()
		if err != nil {
			t.Fatal(err)
		}
		if got := SourceBandwidth(c, tc.edge); math.Abs(got-tc.want) > 1e-9*tc.want {
			t.Errorf("%v: %g Hz, want %g", tc.w, got, tc.want)
		}
	}
}

// A coil's conductor length is its wire: a 100-turn, 5 cm coil has 15.7 m
// of wire, so as a lumped inductor it is valid only to about 1.9 MHz —
// beyond that it needs a distributed model (em/tline, Knight's helix, M9i).
func TestCoilElectricalLength(t *testing.T) {
	coil := Coil{Turns: 100, Diameter: 0.05, WireDiameter: 0.5e-3, Pitch: 1e-3, InsulationEps: 3, Conductor: material.Copper}
	c, _ := circuit.NewBuilder("coil").VoltageSource("V1", "a", "0", 1).Inductor("L1", "a", "0", coil.Inductance()).Build()
	f := 0.1 * SpeedOfLight / coil.WireLength() // the frequency at which the wire is λ/10
	if _, err := CheckQuasiStatic(c, map[string]float64{"L1": coil.WireLength()}, 0, 1.01*f, SizeOptions{}); !errors.Is(err, ErrNotQuasiStatic) {
		t.Errorf("coil wire at %g Hz: %v", f, err)
	}
	t.Logf("%.1f m of wire: lumped up to %.3g MHz at λ/10", coil.WireLength(), f/1e6)
}
