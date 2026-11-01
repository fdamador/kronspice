// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/material"
)

func coil(n int, bif bool) Coil {
	return Coil{Turns: n, Diameter: 0.05, WireDiameter: 0.5e-3, Pitch: 0.6e-3, InsulationEps: 3, Bifilar: bif,
		Conductor: material.Copper}
}

// Tesla, US 512,340: 1,000 turns at 100 V → 0.1 V between turns ordinary,
// 50 V bifilar; energy between turns (50/0.1)² = 250,000 times larger.
func TestTeslaRatio(t *testing.T) {
	o, b := coil(1000, false), coil(1000, true)
	if v := o.TurnVoltage(100); math.Abs(v-0.1) > 1e-15 {
		t.Errorf("ordinary turn voltage %g", v)
	}
	if v := b.TurnVoltage(100); v != 50 {
		t.Errorf("bifilar turn voltage %g", v)
	}
	if r := b.TurnToTurnEnergy(100) / o.TurnToTurnEnergy(100); math.Abs(r-250000) > 1e-6 {
		t.Errorf("energy ratio %g, want 250000", r)
	}
	if o.Inductance() != b.Inductance() {
		t.Error("series-aiding bifilar coil should have the ordinary coil's inductance")
	}
}

// Wheeler's own check value: r = 1 in, ℓ = 2 in, N = 100 ⇒
// L = 1·10⁴/(9 + 20) = 344.8 µH.
func TestWheeler(t *testing.T) {
	c := Coil{Turns: 100, Diameter: 2 * 0.0254, WireDiameter: 0.3e-3, Pitch: 2 * 0.0254 / 100, InsulationEps: 1,
		Conductor: material.Copper}
	if l := c.Inductance(); math.Abs(l-1e4/29*1e-6) > 1e-12 {
		t.Fatalf("L = %g H", l)
	}
	// Long-solenoid limit: µ0 N² A / ℓ within ~2 % for ℓ = 20 D.
	long := Coil{Turns: 2000, Diameter: 0.01, WireDiameter: 0.05e-3, Pitch: 0.1e-3, InsulationEps: 1,
		Conductor: material.Copper}
	ideal := material.Mu0 * 4e6 * math.Pi * 0.005 * 0.005 / long.Length()
	if r := long.Inductance() / ideal; r < 0.96 || r > 1.0 {
		t.Fatalf("long solenoid L/L_ideal = %g", r)
	}
}

func TestTurnCapacitanceAndResistance(t *testing.T) {
	c := coil(100, false)
	// Two-wire line per unit length: π ε0 εr / acosh(p/d).
	want := math.Pi * material.Eps0 * 3 / math.Acosh(1.2) * math.Pi * 0.05
	if got := c.TurnCapacitance(); math.Abs(got-want) > 1e-24 {
		t.Fatalf("C_tt = %g", got)
	}
	// DC: ρ·N·πD/(πd²/4).
	rdc := 1.68e-8 * 100 * math.Pi * 0.05 / (math.Pi * 0.25e-6 / 4)
	if got := c.Resistance(0); math.Abs(got-rdc) > 1e-12*rdc {
		t.Fatalf("R_dc = %g, want %g", got, rdc)
	}
	// Skin depth in copper at 1 MHz ≈ 65 µm ⇒ R ≈ R_dc·0.5 mm/(4·65 µm) ≈ 1.92·R_dc.
	if r := c.Resistance(1e6) / rdc; r < 1.85 || r > 2.0 {
		t.Fatalf("R(1 MHz)/R_dc = %g", r)
	}
	if c.Resistance(10) != rdc {
		t.Fatal("no skin effect expected at 10 Hz")
	}
}

func TestGeometryLimits(t *testing.T) {
	bad := []Coil{
		{Turns: 1, Diameter: 0.05, WireDiameter: 1e-3, Pitch: 2e-3, InsulationEps: 1},
		{Turns: 100, Diameter: 0.05, WireDiameter: 1e-3, Pitch: 0.9e-3, InsulationEps: 1},
		{Turns: 5, Diameter: 0.05, WireDiameter: 1e-3, Pitch: 2e-3, InsulationEps: 1}, // ℓ = 10 mm < 20 mm
	}
	for _, c := range bad {
		if err := c.Validate(); !errors.Is(err, ErrOutsideModelValidity) {
			t.Errorf("%+v: %v", c, err)
		}
	}
	if err := coil(100, true).Validate(); err != nil {
		t.Fatal(err)
	}
	if c, err := ParallelPlate(1e-2, 1e-4, 4); err != nil || math.Abs(c-material.Eps0*4*100) > 1e-20 {
		t.Fatalf("C = %g, %v", c, err)
	}
	if _, err := ParallelPlate(1e-2, 0.05, 4); !errors.Is(err, ErrOutsideModelValidity) {
		t.Fatal("wide gap accepted")
	}
}
