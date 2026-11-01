// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6g, want %.6g ± %.2g", what, got, want, tol)
	}
}

// Kelvin §§1–2 with his editor's correction: each half of the ring takes
// moment of momentum ∓eΦ/2π; the fluid turns at Larmor's rate; energy and
// angular momentum close; switching off returns everything to rest.
func TestInduction(t *testing.T) {
	rg := Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.02}
	r := rg.Run(Coil{R: 0.1, Z: 0.05, NI: 2e4}, 1, 2000)
	near(t, "flux Biot–Savart / Maxwell", r.Flux/r.FluxMaxwell, 1, 1e-4)
	near(t, "fluid L / (eΦ/2π)", r.LFluid/r.LKelvin, 1, 1e-6)
	near(t, "tube L / (−eΦ/2π)", -r.LTube/r.LKelvin, 1, 1e-6)
	near(t, "fluid ω / Larmor", r.OmegaFluid/r.Larmor, 1, 1e-6)
	near(t, "net angular impulse / eΦ/2π", r.AngularImpulseNet/r.LKelvin, 0, 1e-9)
	near(t, "energy ledger residual", r.LedgerResidual, 0, 1e-9)
	near(t, "fluid L after switch-off / eΦ/2π", r.LFluidEnd/r.LKelvin, 0, 1e-9)
	near(t, "net source work / KE", r.WSourceEnd/r.KE, 0, 1e-9)
	if r.Moment >= 0 || r.AxialForce <= 0 {
		t.Errorf("ring should be diamagnetic and repelled: moment %g, force %g", r.Moment, r.AxialForce)
	}
}

// The Zeeman test. Lorentz's bound charge gives the sharp normal triplet,
// ±ωL, opposite circular polarisations along B, π/σ across B. Kelvin's
// gyrostat as he wrote it splits by ωL/2·cos θ, so the ensemble line is
// broad — his own reason for discarding it. The full ring splits by
// eB/2(w+w′) = ωL/2 at every orientation (|G| constant), so its lines are
// sharp, but its modes circle the reflection of B, so the polarisation is
// not Zeeman's.
func TestZeeman(t *testing.T) {
	if testing.Short() {
		t.Skip("ensemble simulation")
	}
	z := DefaultZeeman()
	res := z.Resolution()
	lor := z.Run(Lorentz)
	gyr := z.Run(Gyrostat)
	ful := z.Run(FullRing)

	for _, o := range lor.PerOrient {
		near(t, "Lorentz splitting", o.Split, z.OmegaL, res/2)
	}
	near(t, "Lorentz upper |V|", math.Abs(lor.RadPol.UpperV), 1, 0.02)
	near(t, "Lorentz V upper = −V lower", lor.RadPol.UpperV+lor.RadPol.LowerV, 0, 0.02)
	near(t, "Lorentz centre ∥B", lor.RadPol.CentrePar, 1, 0.02)
	near(t, "Lorentz sides ⊥B", lor.RadPol.SidePerp, 1, 0.02)

	for _, o := range gyr.PerOrient {
		if want := z.OmegaL / 2 * o.CosTheta; want > 2*res {
			near(t, "gyrostat splitting", o.Split, want, res/2)
		}
	}
	if gyr.UpperWidth < 4*res {
		t.Errorf("gyrostat ensemble line should be broad: FWHM %.4f, resolution %.4f", gyr.UpperWidth, res)
	}

	for _, o := range ful.PerOrient {
		near(t, "full-ring splitting", o.Split, z.OmegaL/2, res/2)
	}
	if ful.UpperWidth > 2*res || ful.SharpFrac < 0.9 {
		t.Errorf("full ring should give sharp lines: FWHM %.4f, main pair %.2f", ful.UpperWidth, ful.SharpFrac)
	}
	near(t, "full-ring motion |V| vs 5/11", math.Abs(ful.ModePol.UpperV), MotionPolarisationFull, 0.06)
	if math.Abs(ful.RadPol.UpperV) > 0.2 {
		t.Errorf("full ring radiation along B should be nearly unpolarised, V = %.2f", ful.RadPol.UpperV)
	}
	t.Logf("splitting: Lorentz %.4f, gyrostat %.4f…%.4f (FWHM %.4f), full ring %.4f…%.4f (FWHM %.4f)",
		z.OmegaL, gyr.SplitMin, gyr.SplitMax, gyr.UpperWidth, ful.SplitMin, ful.SplitMax, ful.UpperWidth)
}

// G = (ea²/2)(B − 2(n·B)n) has the same length for every n.
func TestGyroVector(t *testing.T) {
	B := v3{0, 0, 0.3}
	for _, c := range []float64{0, 0.3, 0.7, 1} {
		n := v3{math.Sqrt(1 - c*c), 0, c}
		near(t, "|G|", norm(GyroVector(2, 0.5, n, B)), 2*0.25*0.3/2, 1e-15)
	}
}
