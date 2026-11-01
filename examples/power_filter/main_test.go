// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"math"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/evo"
)

// The unfiltered rectifier: GA and IEEE 1459 agree per phase, the largest
// distortion is the fundamental voltage with the 5th-harmonic current, and
// the quaternion instantaneous power averages to ΣP. A short GA run then
// finds a filter with a better power factor (M9f: GA power as a GA
// fitness).
func TestPowerFilter(t *testing.T) {
	ctx := context.Background()
	base, err := Analyze(ctx, Filter{}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	for k, ph := range base.Phases {
		g, q := ph.GA, ph.IEEE
		// The GA product stops at the 25th harmonic; the time-domain RMS
		// values include the commutation notches above it (≈0.3 % of S).
		if math.Abs(g.P-q.P) > 1e-3*q.S || math.Abs(g.S1-q.S1) > 1e-9*q.S1 || math.Abs(g.S-q.S) > 5e-3*q.S {
			t.Errorf("phase %d: GA P %g S %g S1 %g, IEEE 1459 P %g S %g S1 %g", k, g.P, g.S, g.S1, q.P, q.S, q.S1)
		}
		if d := g.Distortion[0]; d.H != 1 || d.K != 5 {
			t.Errorf("phase %d: largest distortion pair %d–%d, want 1–5", k, d.H, d.K)
		}
	}
	if math.Abs(base.QuatMean-base.P) > 1e-3*base.P {
		t.Errorf("quaternion mean power %g, ΣP %g", base.QuatMean, base.P)
	}
	if !(base.PF > 0.85 && base.PF < 0.95) {
		t.Errorf("unfiltered PF %g outside the expected 0.85–0.95", base.PF)
	}
	if testing.Short() {
		t.Skip("GA run skipped in -short mode")
	}
	// A generous per-evaluation deadline: under the race detector a
	// simulation runs ~20× slower, and a timed-out candidate would change
	// the run.
	s, err := evo.Run(ctx, Problem(), evo.Options{Population: 12, Generations: 4, Seed: 3, Timeout: 10 * time.Minute}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if pf := 1 - s.Best.Fitness; !(pf > base.PF+0.01) {
		t.Errorf("evolved PF %g, unfiltered %g", pf, base.PF)
	}
}
