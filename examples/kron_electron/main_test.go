// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s = %.6g, want %.6g ± %.2g", what, got, want, tol)
	}
}

// Kron's frame reproduces Kelvin's ring: moments of momentum ∓eΦ/2π, the
// fluid at Larmor's rate, kinetic energy equal to the bench's, the ledger
// closed. Same ring and coil as examples/kelvin_ring's TestInduction.
func TestKelvinTwin(t *testing.T) {
	rg := Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.02}
	c := Coil{R: 0.1, Z: 0.05, NI: 2e4}
	for _, si := range []bool{false, true} {
		r, err := RunKelvin(rg, c, 1, si)
		if err != nil {
			t.Fatal(err)
		}
		near(t, "fluid L / (eΦ/2π)", r.LFluid/r.LKelvin, 1, 1e-6)
		near(t, "tube L / (−eΦ/2π)", -r.LTube/r.LKelvin, 1, 1e-6)
		near(t, "fluid ω / Larmor", r.OmegaFluid/r.Larmor, 1, 1e-6)
		// The bench's ledger: the source's work against the ring's back-EMF
		// is the kinetic energy plus the (negative) mutual-field energy.
		near(t, "ring ledger residual", r.RingResidual, 0, 1e-6)
		// The engine's ledger is dominated by the coil's own ½L_c I² (124 J
		// against 5e-16 J of kinetic energy); its mismatch is the
		// integration error over the ramp's breakpoints.
		near(t, "engine mismatch / source energy", r.Mismatch/r.EnergyIn, 0, 1e-3)
	}
	// Kelvin's values for the main program's ring (both masses 10 g), from
	// the Lorentz-force bench: 1.11207e-10 kg m²/s and 4.44829e-06 rad/s.
	// The bench drives its charges with the Biot–Savart flux, 5.7e-6 above
	// Maxwell's elliptic value used here, so its kinetic energy is 1.1e-5
	// higher.
	r, err := RunKelvin(Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.01}, c, 1, false)
	if err != nil {
		t.Fatal(err)
	}
	near(t, "fluid L vs bench", r.LFluid/1.11207e-10, 1, 1e-5)
	near(t, "fluid ω vs bench", r.OmegaFluid/4.44829e-06, 1, 1e-5)
	near(t, "KE vs bench", r.KE/4.94683e-16, 1, 2e-5)
}

// Natural frequencies: the Compton ring's circumferential frequency is the
// Compton frequency m_e c²/h by identity; precession is g·eB/4πm; the
// fitted massless ring precesses at the measured ESR frequency by
// construction; the tilt libration of Kelvin's ring is the part-2 period.
func TestFrequencies(t *testing.T) {
	c := RingModel{A: LambdaC, Core: LambdaC / 100, Mass: ME}
	f := c.Frequencies(HBAR/2, 1)
	near(t, "c/2πλ̄_C / (m_e c²/h)", f.Circumferential/ComptonHz, 1, 1e-9)
	n := c.AtSpin(HBAR / 2)
	near(t, "precession / (g eB/4πm)", f.Precession/(n.G*QE/(4*math.Pi*ME)), 1, 1e-12)
	fit := FitMassless(0.01).Frequencies(HBAR/2, 1)
	near(t, "fitted ring precession / measured ESR", fit.Precession/GammaElectronHz, 1, 1e-4)
	if fit.Libration != 0 {
		t.Errorf("massless ring has no libration, got %g", fit.Libration)
	}
	tl := Tilt{A: 0.05, E: 1e-6, Mass: 0.01, OmegaSpin: 1000, B: 1, Theta0: 0.05}
	near(t, "Kelvin ring libration period", 1/(math.Sqrt(tl.Moment()*tl.B/tl.Inertia())/(2*math.Pi)), 19.869, 1e-2)
}

// A massless thin ring with S = ħ/2 moves at v/c = π/(2α(ln(8a/r) − 7/4))
// whatever its radius.
func TestMasslessSpeed(t *testing.T) {
	for _, a := range []float64{1e-15, LambdaC, 1e-10} {
		n := RingModel{A: a, Core: a / 100, Mass: 0}.AtSpin(HBAR / 2)
		want := math.Pi / (2 * FineStructure * (math.Log(800) - 1.75))
		near(t, "v/c", n.V/C0, want, 1e-9)
	}
}

// The torque ½iᵀ(∂L/∂θ)i of frame I is −μB sin θ: the tilting ring's
// period equals the pendulum value with its amplitude correction.
func TestTiltPeriod(t *testing.T) {
	r, err := RunTilt(Tilt{A: 0.05, E: 1e-6, Mass: 0.01, OmegaSpin: 1000, B: 1, Theta0: 0.05})
	if err != nil {
		t.Fatal(err)
	}
	near(t, "period / analytic", r.Period/r.Analytic, 1, 2e-4)
	near(t, "ring current drift", r.CurrentDrift, 0, 1e-9)
}

// The tilt's harmonics are the pendulum's: odd only, with the Jacobi
// series amplitudes.
func TestTiltHarmonics(t *testing.T) {
	h, err := TiltHarmonics(Tilt{A: 0.05, E: 1e-6, Mass: 0.01, OmegaSpin: 1000, B: 1, Theta0: 1})
	if err != nil {
		t.Fatal(err)
	}
	near(t, "3rd harmonic ratio", h.Measured[0], h.Exact[0], 2e-4)
	near(t, "5th harmonic ratio", h.Measured[1], h.Exact[1], 2e-4)
	near(t, "2nd harmonic ratio", h.Even, 0, 2e-4)
	near(t, "exact 3rd at small angle / θ₀²/64", PendulumHarmonics(0.05)[0]/(0.05*0.05/64), 1, 1e-3)
}

// Ring electron arithmetic: a ring whose mass circulates with its charge
// has g = 1/(1 + L e²/4π² m a²) ≤ 1; the massless fit carries the measured
// moment at S = ħ/2 but only at femtometre radius and v > c; every ring
// lies above the size bound.
func TestElectronRing(t *testing.T) {
	m := RingModel{A: LambdaC, Core: LambdaC / 100, Mass: ME}
	n := m.AtSpin(HBAR / 2)
	near(t, "g·(1 + L e²/4π² m a²)", n.G*(1+n.L*QE*QE/(4*math.Pi*math.Pi*ME*LambdaC*LambdaC)), 1, 1e-12)
	near(t, "precession / Larmor (g = 1 limit)", n.Precession/(QE/(2*ME)), n.G, 1e-12)
	if n.V > C0 {
		t.Errorf("Compton ring with mass: v = %.3g c, want below c", n.V/C0)
	}
	fit := FitMassless(0.01).AtSpin(HBAR / 2)
	near(t, "fitted moment / measured", fit.Mu/MuMeasured, 1, 1e-9)
	if fit.V < C0 {
		t.Errorf("massless fit: v = %.3g c, expected superluminal", fit.V/C0)
	}
	for _, sc := range Scenarios() {
		if sc.A < SizeBound {
			t.Errorf("%s: a = %.3g m below the size bound", sc.Name, sc.A)
		}
	}
}

// Field lines of the ring current close on themselves.
func TestFieldLinesClose(t *testing.T) {
	for _, rho := range []float64{0.3, 0.6, 0.88} {
		ln := FieldLine(1, rho, 0, 40000)
		end := ln[len(ln)-1]
		near(t, "line closes", math.Hypot(end[0]-rho, end[1]), 0, 1e-9)
		if len(ln) >= 40000 {
			t.Errorf("rho %.2f: line did not close in 40000 steps", rho)
		}
	}
}
