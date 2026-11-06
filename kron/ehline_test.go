// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// Light-by-light scattering in its classical low-energy form on Kron's
// network with the Euler–Heisenberg vacuum as its nonlinear primitive
// (EHLine): the phase a weak probe gains crossing a strong field.
//
// Expected, from the Euler–Heisenberg Lagrangian (derived from
// em.EHLinearize; Toll 1952; Adler 1971; R. Battesti, C. Rizzo, Rep. Prog.
// Phys. 76, 016401, 2013), with a = (α/90π)(E/E_s)² or (α/90π)(B/B_c)²:
//
//	counter-propagating plane-wave pump, probe E along the pump's E:  n − 1 = 16a
//	                                     probe E across it:            n − 1 = 28a
//	(the (α/45π)(11 ∓ 3)(E/E_s)² of the literature)
//	co-propagating pump, probe riding inside it (F = G = 0):           n − 1 = 0
//	static magnetic background B₀ across the beam, probe E along B₀:  n − 1 = 7a
//	                                               probe E across B₀:  n − 1 = 4a
//
// The probe's field is the run with pump and probe less the run with the
// pump alone; its phase is compared with the probe run alone. The linear
// line (Maxwell's vacuum) superposes exactly: no phase at all.
//
// Accuracy. The lattice spoils the exact cancellations at O((kΔz)²): on
// it a wave's E/B̄ is not exactly c. At λ = 20 cells the counter-propagating
// shifts come out 0.4 % low and a probe riding inside a co-moving pump picks
// up 0.13a (∥) and 0.23a (⊥) instead of 0; at λ = 40 cells both errors fall
// fourfold (0.1 %; 0.035a, 0.06a), and halving Δt changes them by only
// 15 %. The tolerances (3 %; 2 % of the counter-propagating shift over the
// same path) leave room for that and nothing else.
//
// Not here: photon splitting (one photon into two in a magnetic
// background). Its collinear amplitude vanishes at the quartic order this
// Lagrangian has; the first non-zero term is the hexagon, sixth order in
// the fields (S. L. Adler, Ann. Phys. 67, 599, 1971), and it is
// suppressed further by (ħω/m_e c²)⁵, so it needs both the next term of the
// effective Lagrangian and quantum photon energies.
//
// Scale: Δz = 50 nm, the probe a 1 µm wave packet (ħω = 1.2 eV ≪ m_e c²,
// the classical low-frequency limit), the pump at 0.09 E_s: an intensity
// of about 2×10²⁷ W/cm², some 10⁴ beyond the strongest lasers (≈10²³).
func TestEHLineLightByLight(t *testing.T) {
	if testing.Short() {
		t.Skip("two dozen 6000-cell nonlinear runs (about a minute)")
	}
	const (
		n       = 6000
		dz      = 50e-9
		courant = 0.5
		lambda  = 20.0 // probe wavelength [cells]
		sigma   = 40.0 // probe envelope [cells]
		z0      = 300.0
		zs      = 600.0  // pump start [cells]
		lp      = 1600.0 // pump plateau [cells]
		ramp    = 100.0
	)
	c := em.SpeedOfLight
	ep := 0.09 * em.SchwingerField
	k0 := 2 * math.Pi / (lambda * dz)
	probeAmp := 1e-4 * em.SchwingerField
	steps := int((zs + lp + 2*ramp - z0 + 8*sigma) / 2 / courant)

	// Pump profile in cells from its trailing… leading edge: sin² ramps.
	pump := func(x float64) float64 { // x = position within the pump [cells]
		switch {
		case x < 0 || x > lp+2*ramp:
			return 0
		case x < ramp:
			s := math.Sin(math.Pi * x / (2 * ramp))
			return s * s
		case x > lp+ramp:
			s := math.Sin(math.Pi * (lp + 2*ramp - x) / (2 * ramp))
			return s * s
		}
		return 1
	}
	envelope := func(z float64) float64 { // probe, z in cells
		u := (z - z0) / sigma
		return math.Exp(-u*u/2) * math.Cos(2*math.Pi*(z-z0)/lambda)
	}
	half := courant / 2 // B is half a step ahead: c·Δt/2 = courant/2 cells

	type config struct {
		pumpDir  float64    // +1, −1, or 0 for none
		pumpAt   float64    // pump's start [cells] (0 means zs)
		bg       [3]float64 // static background B [T]
		probePol int        // 0 = x, 1 = y; −1 for no probe
		linear   bool
	}
	start := func(cf config) float64 {
		if cf.pumpAt != 0 {
			return cf.pumpAt
		}
		return zs
	}
	run := func(cf config) (ex, ey []float64) {
		l, err := NewEHLine(n, dz, courant)
		if err != nil {
			t.Fatal(err)
		}
		l.Linear = cf.linear
		e := func(z float64) [3]float64 {
			zc := z / dz
			var f [3]float64
			if cf.pumpDir != 0 {
				f[0] = ep * pump(zc-start(cf))
			}
			if cf.probePol >= 0 {
				f[cf.probePol] += probeAmp * envelope(zc)
			}
			return f
		}
		b := func(z float64) [3]float64 {
			zc := z / dz
			f := cf.bg
			if cf.pumpDir != 0 {
				// a wave along ±z with E along x has B_y = ±E_x/c
				f[1] += cf.pumpDir * ep * pump(zc-cf.pumpDir*half-start(cf)) / c
			}
			switch cf.probePol { // the probe runs along +z
			case 0:
				f[1] += probeAmp * envelope(zc-half) / c
			case 1:
				f[0] -= probeAmp * envelope(zc-half) / c
			}
			return f
		}
		if err := l.SetFields(e, b); err != nil {
			t.Fatal(err)
		}
		w0 := l.Energy()
		for s := 0; s < steps; s++ {
			if err := l.Step(); err != nil {
				t.Fatal(err)
			}
		}
		if cf.pumpDir == -1 && cf.probePol == 0 && !cf.linear {
			// The Euler–Heisenberg part of the energy is ≈ 6a ≈ 1e-6 of it.
			if r := l.Energy()/w0 - 1; math.Abs(r) > 2e-8 {
				t.Errorf("energy %.10g → %.10g J/m² (%.2g)", w0, l.Energy(), r)
			}
			t.Logf("energy %.10g → %.10g J/m² (%.2g; the Euler–Heisenberg part is ≈ %.2g of it)", w0, l.Energy(), l.Energy()/w0-1, 6*em.EHBirefringenceParameter(ep/c))
		}
		return append([]float64(nil), l.Ex...), append([]float64(nil), l.Ey...)
	}
	pick := func(cf config, pol int) []float64 {
		ex, ey := run(cf)
		if pol == 1 {
			return ey
		}
		return ex
	}
	// Complex amplitude of the probe's carrier round where it is.
	zEnd := z0 + float64(steps)*courant
	amp := func(f []float64) complex128 {
		var s complex128
		for k := int(zEnd - 6*sigma); k <= int(zEnd+6*sigma); k++ {
			s += complex(f[(k+n)%n], 0) * cmplx.Exp(complex(0, -k0*float64(k)*dz))
		}
		return s
	}
	phase := func(cf config) float64 {
		withPump := pick(cf, cf.probePol)
		alone := cf
		alone.pumpDir, alone.bg = 0, [3]float64{}
		ref := amp(pick(alone, cf.probePol))
		if cf.pumpDir != 0 || cf.bg != ([3]float64{}) {
			only := cf
			only.probePol = -1
			po := pick(only, cf.probePol)
			for k := range withPump {
				withPump[k] -= po[k]
			}
		}
		return cmplx.Phase(amp(withPump) / ref)
	}

	// Interaction lengths [m]: half the pump's ∫(E/E_p)² for a counter-
	// propagating pump (they close at 2c); the whole path in a background.
	lint := 0.0
	for x := 0.0; x <= lp+2*ramp; x++ {
		lint += pump(x) * pump(x)
	}
	lint *= dz / 2
	path := float64(steps) * courant * dz
	aE := em.EHBirefringenceParameter(ep / c)
	b0 := 0.09 * em.CriticalB
	aB := em.EHBirefringenceParameter(b0)

	cases := []struct {
		name string
		cf   config
		want float64
	}{
		{"counter-propagating pump, probe ∥ pump E (16a)", config{pumpDir: -1, probePol: 0}, k0 * 16 * aE * lint},
		{"counter-propagating pump, probe ⊥ pump E (28a)", config{pumpDir: -1, probePol: 1}, k0 * 28 * aE * lint},
		{"co-propagating pump, probe ∥ pump E inside it (0)", config{pumpDir: +1, pumpAt: -150, probePol: 0}, 0},
		{"co-propagating pump, probe ⊥ pump E inside it (0)", config{pumpDir: +1, pumpAt: -150, probePol: 1}, 0},
		{"magnetic background, probe E ∥ B₀ (7a)", config{bg: [3]float64{b0, 0, 0}, probePol: 0}, k0 * 7 * aB * path},
		{"magnetic background, probe E ⊥ B₀ (4a)", config{bg: [3]float64{b0, 0, 0}, probePol: 1}, k0 * 4 * aB * path},
		{"linear vacuum, counter-propagating (0)", config{pumpDir: -1, probePol: 0, linear: true}, 0},
	}
	// Zero is judged against the counter-propagating shift over the same
	// path: a co-moving probe sits in the pump for the whole run.
	scale := k0 * 16 * aE * path
	for _, cs := range cases {
		got := phase(cs.cf)
		t.Logf("%-50s Δφ = %.6g rad, Euler–Heisenberg %.6g rad", cs.name, got, cs.want)
		if cs.want == 0 {
			if math.Abs(got) > 2e-2*scale {
				t.Errorf("%s: Δφ %.3g rad, want 0 (against %.3g counter-propagating)", cs.name, got, scale)
			}
			continue
		}
		if math.Abs(got/cs.want-1) > 3e-2 {
			t.Errorf("%s: Δφ %.6g rad, Euler–Heisenberg %.6g (%.3g)", cs.name, got, cs.want, got/cs.want-1)
		}
	}
	t.Logf("a(pump %.2g E_s) = %.4g, interaction length %.4g µm; a(B₀ %.2g B_c) = %.4g over %.4g µm",
		ep/em.SchwingerField, aE, lint*1e6, b0/em.CriticalB, aB, path*1e6)
}

// The line's plumbing, fast: a plane wave on the linear line keeps its
// shape and energy and moves at the lattice's speed; the nonlinear line
// at a weak field (10⁻⁴ E_s, a ≈ 3×10⁻¹⁵) gives the same fields; a field
// past the weak-field limit is refused.
func TestEHLineBasics(t *testing.T) {
	const n, dz, courant = 400, 50e-9, 0.5
	c := em.SpeedOfLight
	amp := 1e-4 * em.SchwingerField
	pulse := func(z float64) float64 { u := (z/dz - 100) / 15; return amp * math.Exp(-u*u/2) }
	mk := func(linear bool) *EHLine {
		l, err := NewEHLine(n, dz, courant)
		if err != nil {
			t.Fatal(err)
		}
		l.Linear = linear
		if err := l.SetFields(func(z float64) [3]float64 { return [3]float64{pulse(z), 0, 0} },
			func(z float64) [3]float64 { return [3]float64{0, pulse(z-courant/2*dz) / c, 0} }); err != nil {
			t.Fatal(err)
		}
		return l
	}
	lin, nl := mk(true), mk(false)
	w0 := lin.Energy()
	for s := 0; s < 400; s++ {
		if err := lin.Step(); err != nil {
			t.Fatal(err)
		}
		if err := nl.Step(); err != nil {
			t.Fatal(err)
		}
	}
	var m0, m1, diff float64
	for k := 0; k < n; k++ {
		m0 += lin.Ex[k] * lin.Ex[k]
		m1 += float64(k) * lin.Ex[k] * lin.Ex[k]
		diff = math.Max(diff, math.Abs(nl.Ex[k]-lin.Ex[k])/amp)
	}
	if moved := m1/m0 - 100; math.Abs(moved-200) > 1 {
		t.Errorf("pulse moved %.3f cells in 400 steps of half a cell", moved)
	}
	if r := lin.Energy()/w0 - 1; math.Abs(r) > 1e-3 {
		t.Errorf("linear energy changed by %.3g", r)
	}
	if diff > 1e-10 {
		t.Errorf("weak-field nonlinear line differs from the linear one by %.3g of the field", diff)
	}
	t.Logf("pulse moved %.4f cells; energy %.2g; weak-field nonlinear − linear %.2g of the field", m1/m0-100, lin.Energy()/w0-1, diff)
	strong, _ := NewEHLine(n, dz, courant)
	if err := strong.SetFields(func(z float64) [3]float64 { return [3]float64{0.2 * em.SchwingerField} },
		func(z float64) [3]float64 { return [3]float64{} }); !errors.Is(err, em.ErrEHStrongField) {
		t.Errorf("0.2 E_s accepted: %v", err)
	}
}
