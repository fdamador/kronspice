// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"errors"
	"testing"

	"github.com/fdamador/kronspice/netlist"
)

func TestBuildDivider(t *testing.T) {
	c, err := NewBuilder("div").
		VoltageSource("V1", "in", "0", 10).
		Resistor("R1", "in", "out", 1e3).
		Resistor("R2", "out", "gnd", 1e3).
		Build()
	if err != nil {
		t.Fatal(err)
	}
	if c.NumNodes() != 3 {
		t.Fatalf("nodes = %v", c.NodeNames)
	}
	if c.Elements[2].Nodes[1] != 0 {
		t.Error("gnd should map to node 0")
	}
}

func TestValidationErrors(t *testing.T) {
	cases := []struct {
		name string
		b    *Builder
		want error
	}{
		{"empty", NewBuilder("x"), ErrEmptyCircuit},
		{"no ground", NewBuilder("x").Resistor("R1", "a", "b", 1), ErrNoGround},
		{"floating", NewBuilder("x").Resistor("R1", "a", "0", 1).Resistor("R2", "b", "c", 1), ErrFloatingNode},
		{"current source only path", NewBuilder("x").CurrentSource("I1", "a", "0", 1), ErrFloatingNode},
		{"duplicate", NewBuilder("x").Resistor("R1", "a", "0", 1).Resistor("r1", "a", "0", 2), ErrDuplicateName},
		{"zero resistor", NewBuilder("x").Resistor("R1", "a", "0", 0), ErrBadValue},
		{"self loop", NewBuilder("x").Resistor("R1", "a", "a", 1).Resistor("R2", "a", "0", 1), ErrSelfConnected},
	}
	for _, tc := range cases {
		if _, err := tc.b.Build(); !errors.Is(err, tc.want) {
			t.Errorf("%s: err = %v, want %v", tc.name, err, tc.want)
		}
	}
}

func TestFromDeck(t *testing.T) {
	d, err := netlist.ParseString("t\nV1 1 0 5\nR1 1 0 1k\n.op\n")
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	if len(c.Elements) != 2 || c.Elements[1].Value != 1000 {
		t.Fatalf("elements = %+v", c.Elements)
	}
}

func TestCapacitorPaths(t *testing.T) {
	// Node "b" reaches ground only through a capacitor: no DC path, but fine
	// for a transient analysis that skips the operating point.
	mk := func() *Builder {
		return NewBuilder("c").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "0", 1).Capacitor("C1", "b", "0", 1e-6)
	}
	if _, err := mk().Build(); !errors.Is(err, ErrFloatingNode) {
		t.Fatalf("DC build err = %v, want ErrFloatingNode", err)
	}
	if _, err := mk().BuildTransient(); err != nil {
		t.Fatalf("transient build err = %v", err)
	}
	// Inductors are DC paths.
	if _, err := NewBuilder("l").VoltageSource("V1", "a", "0", 1).Inductor("L1", "a", "b", 1e-3).Resistor("R1", "b", "0", 1).Build(); err != nil {
		t.Fatal(err)
	}
}

func TestWaveformValidation(t *testing.T) {
	bad := []Waveform{
		{Kind: WavePulse, Args: []float64{1}},
		{Kind: WavePulse, Args: []float64{0, 1, -1}},
		{Kind: WaveSin, Args: []float64{0, 1, -5}},
		{Kind: WavePWL, Args: []float64{0, 1, 2}},
		{Kind: WavePWL, Args: []float64{1, 0, 0.5, 1}},
		{Kind: WaveExp, Args: []float64{0, 1, 0, -1}},
	}
	for _, w := range bad {
		_, err := NewBuilder("w").VoltageWave("V1", "a", "0", w).Resistor("R1", "a", "0", 1).Build()
		if !errors.Is(err, ErrBadWaveform) {
			t.Errorf("%v %v: err = %v, want ErrBadWaveform", w.Kind, w.Args, err)
		}
	}
	good := Waveform{Kind: WavePulse, Args: []float64{0, 1, 0, 1e-9, 1e-9, 1e-6, 2e-6}}
	if _, err := NewBuilder("w").VoltageWave("V1", "a", "0", good).Resistor("R1", "a", "0", 1).Build(); err != nil {
		t.Fatal(err)
	}
	if _, err := NewBuilder("w").Capacitor("C1", "a", "0", 0).Resistor("R1", "a", "0", 1).Build(); !errors.Is(err, ErrBadValue) {
		t.Fatalf("zero capacitor err = %v", err)
	}
}

// K elements: two distinct inductors, 0 < k ≤ 1, one K per pair.
func TestMutualValidation(t *testing.T) {
	base := func() *Builder {
		return NewBuilder("k").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "b", 1).
			Inductor("L1", "b", "0", 1e-3).Inductor("L2", "c", "0", 1e-3).Resistor("R2", "c", "0", 1)
	}
	if _, err := base().Mutual("K1", "L1", "L2", 1).Build(); err != nil {
		t.Fatalf("k = 1: %v", err)
	}
	for _, bad := range []*Builder{
		base().Mutual("K1", "L1", "L2", 0),
		base().Mutual("K1", "L1", "L2", 1.5),
		base().Mutual("K1", "L1", "L1", 0.5),
		base().Mutual("K1", "L1", "R1", 0.5),
		base().Mutual("K1", "L1", "L9", 0.5),
		base().Mutual("K1", "L1", "L2", 0.5).Mutual("K2", "L2", "L1", 0.3),
	} {
		if _, err := bad.Build(); !errors.Is(err, ErrBadValue) {
			t.Errorf("want ErrBadValue, got %v", err)
		}
	}
}
