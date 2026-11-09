// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// IF97 region 4 verification value, and the tables' saturation points.
func TestSatPressure(t *testing.T) {
	if p := SatPressure(300); math.Abs(p/3536.58941-1) > 1e-8 {
		t.Fatalf("p_s(300 K) = %.6f Pa", p)
	}
	for _, s := range []Sat{sat20k, sat101k, sat150k, sat200k} {
		if p := SatPressure(s.T + 273.15); math.Abs(p/s.P-1) > 3e-3 {
			t.Fatalf("table %.0f Pa at %.2f °C, IF97 gives %.0f", s.P, s.T, p)
		}
		if math.Abs(s.Hf+s.Hfg-s.Hg) > 0.25 || math.Abs(s.Sf+s.Sfg-s.Sg) > 2e-3 {
			t.Fatalf("table row %+v inconsistent", s)
		}
	}
}

// The loop's books close, the turbine gives back less than Carnot allows,
// and the kronspice shaft circuit has the motor supply the rest.
func TestLoop(t *testing.T) {
	p := DefaultParams()
	m := DimpleMachine(0.0843, 2, 0.006, 0.5, 3.33, 1.2, 0.8e-3, 1828, 1.18e5)
	for _, r := range []LoopResult{
		LoopA(p, m, sat150k, sat20k, 0.3, 3000),
		LoopB(p, m, 3e5, sat200k, sat101k, sat20k, 0.3, 3000),
	} {
		if math.Abs(r.Residual) > 1e-9*r.PElec {
			t.Fatalf("%s: books residual %.3g W", r.Name, r.Residual)
		}
		carnot := 1 - (sat20k.T+273.15)/(sat150k.T+273.15)
		if r.PTurbine/r.PShaft > carnot {
			t.Fatalf("%s: turbine returns %.3f of the heat, above Carnot %.3f", r.Name, r.PTurbine/r.PShaft, carnot)
		}
		if !(r.PElec > 0.8*r.PShaft) {
			t.Fatalf("%s: the loop cannot run itself; electricity %.0f W for %.0f W shaft", r.Name, r.PElec, r.PShaft)
		}
		mot, res, err := LoopCircuit(r, r.RPM)
		if err != nil || math.Abs(mot/r.PMotorShaft-1) > 1e-9 || res > 1e-9 {
			t.Fatalf("%s: circuit motor %.3f vs %.3f, residual %.1e, err %v", r.Name, mot, r.PMotorShaft, res, err)
		}
	}
}
