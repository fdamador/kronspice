// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Ledger efficiency as a GA fitness (M9g): with a 50 Ω source the DC
// efficiency R_L/(R_L + 50) rises with R_L, so the GA drives R_L to its
// upper bound; and a Build that checks the lumped limit penalizes an
// electrically large candidate as invalid.
func TestEfficiencyObjective(t *testing.T) {
	build := func(v []float64) (*circuit.Circuit, error) {
		c, err := circuit.NewBuilder("load").VoltageSource("V1", "s", "0", 10).
			Resistor("RS", "s", "a", 50).Resistor("RL", "a", "0", v[0]).Build()
		if err != nil {
			return nil, err
		}
		// The load's leads, v[1] metres long, carry a 100 MHz signal.
		if _, err := em.CheckQuasiStatic(c, map[string]float64{"RL": v[1]}, 0, 100e6, em.SizeOptions{}); err != nil {
			return nil, err
		}
		return c, nil
	}
	p := Problem{
		Space:      Space{{Name: "RL", Min: 1, Max: 1000}, {Name: "lead", Min: 0.01, Max: 1}},
		Build:      build,
		Objectives: []Objective{OPEfficiency("RL")},
	}
	s, err := Run(context.Background(), p, Options{Population: 24, Generations: 15, Seed: 2}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	rl := s.Best.Values[0]
	if want := 1 - rl/(rl+50); math.Abs(s.Best.Fitness-want) > 1e-12 || rl < 500 {
		t.Errorf("best R_L %g with fitness %g (want 1 − R_L/(R_L+50) = %g, R_L near 1 kΩ)", rl, s.Best.Fitness, want)
	}
	// 100 MHz: λ = 3 m, so leads above 0.3 m fail.
	_, err = build([]float64{100, 0.5})
	if _, class := Penalty(err); class != "invalid" {
		t.Errorf("0.5 m lead at 100 MHz: class %q (%v)", class, err)
	}
	invalid := 0
	for _, ind := range s.Population {
		if ind.Failure == "invalid" {
			invalid++
		}
	}
	if s.Best.Values[1] > 0.3 {
		t.Errorf("best candidate has a %.2f m lead", s.Best.Values[1])
	}
	t.Logf("best R_L %.1f Ω, η = %.4f; %d of %d in the last population electrically large", rl, 1-s.Best.Fitness, invalid, len(s.Population))
}
