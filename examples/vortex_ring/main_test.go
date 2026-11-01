// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// Kelvin's ring speed from the engine's inductance: the numerical speed
// approaches Γ/4πR(ln 8R/a − ¼) as the core thins, the error falling
// roughly as (a/R)², the order of the terms the thin-ring formula drops.
func TestKelvinSpeed(t *testing.T) {
	prev := math.Inf(1)
	for _, a := range []float64{0.2, 0.1, 0.05} {
		rel := math.Abs(Speed(1, 1, a, 24, 48)/KelvinSpeed(1, 1, a) - 1)
		lrel := math.Abs(SelfInductance(1, a, 24, 48)/ThinRingInductance(1, a) - 1)
		t.Logf("a = %.2f: speed %.3f %%, inductance %.3f %% from the thin-ring forms", a, 100*rel, 100*lrel)
		if rel > prev/2.5 {
			t.Errorf("speed error should shrink ~a²: %.4f after %.4f", rel, prev)
		}
		prev = rel
	}
	if prev > 2e-3 {
		t.Errorf("a = 0.05: speed %.3f %% from Kelvin's formula", 100*prev)
	}
}

// Kelvin p. 10: impulse ↔ magnetic moment; the flow at the centre of a
// ring of circulation Γ is Γ/2R (the field at the centre of a loop).
func TestAnalogy(t *testing.T) {
	loop := em.Loop{P: em.Circle(v3{}, v3{0, 0, 1}, 0.7, 720), I: 1.3}
	if u := Velocity(loop, v3{})[2]; math.Abs(u/(1.3/(2*0.7))-1) > 1e-4 {
		t.Errorf("centre flow %.6f, want Γ/2R = %.6f", u, 1.3/1.4)
	}
}

// Helmholtz's leapfrogging rings: they pass through each other again and
// again while impulse and energy stay constant.
func TestLeapfrog(t *testing.T) {
	lf := RunLeapfrog(1, [2]Ring{{R: 1, Z: 0, A: 0.1}, {R: 1, Z: 0.6, A: 0.1}}, 30, 0.02, 360)
	t.Logf("%d passes, period %.2f, impulse drift %.1e, energy drift %.1e", lf.Passes, lf.Period, lf.ImpulseDrift, lf.EnergyDrift)
	if lf.Passes < 8 {
		t.Errorf("expected repeated leapfrogging, got %d passes", lf.Passes)
	}
	if lf.ImpulseDrift > 1e-4 || lf.EnergyDrift > 1e-4 {
		t.Errorf("invariants drift: impulse %.1e, energy %.1e", lf.ImpulseDrift, lf.EnergyDrift)
	}
}
