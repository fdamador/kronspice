// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import (
	"math"
	"testing"
)

// TestGravityProbeB: a gyroscope in a 642 km polar orbit (circular, a =
// equatorial radius + 642 km) precesses −6,606.1 mas/yr (geodetic) and
// −39.2 mas/yr (frame dragging) in general relativity; Gravity Probe B
// measured −6,601.8 ± 18.3 and −37.2 ± 7.2 (C. W. F. Everitt et al.,
// Phys. Rev. Lett. 106, 221101, 2011; orbit from the mission parameters).
// The model here is a spherical Earth without J2 and a circular orbit, so
// it is held to 0.5 % (geodetic) and 3 % (frame dragging).
func TestGravityProbeB(t *testing.T) {
	a := REarthEq + 642e3
	geo, fd := OrbitAveragedPrecession(Earth, a, Vec{1, 0, 0}, Vec{0, 0, 1}, 3600)
	// GP-B read the precession against the guide star IM Pegasi, which
	// lies in the orbit plane at declination +16°50′28.54″ (Hipparcos,
	// via Bartel et al., arXiv:1509.07529, Table 1). Geodetic precession is
	// along the orbit normal, already perpendicular to the star; frame
	// dragging is along Earth's axis, so GP-B saw its component
	// perpendicular to the star line, smaller by cos δ.
	dec := (16 + 50.0/60 + 28.54/3600) * math.Pi / 180
	star := Vec{math.Cos(dec), 0, math.Sin(dec)}
	g := MasPerYear(geo.sub(star.scale(geo.dot(star))).norm())
	t.Logf("full frame dragging along Earth's axis %.2f mas/yr", MasPerYear(fd.norm()))
	f := MasPerYear(fd.sub(star.scale(fd.dot(star))).norm())
	t.Logf("geodetic %.1f mas/yr (GR 6606.1, GP-B 6601.8 ± 18.3); frame dragging %.2f mas/yr (GR 39.2, GP-B 37.2 ± 7.2)", g, f)
	if math.Abs(g-6606.1)/6606.1 > 0.005 {
		t.Errorf("geodetic %.1f mas/yr", g)
	}
	if math.Abs(f-39.2)/39.2 > 0.03 {
		t.Errorf("frame dragging %.2f mas/yr", f)
	}
	// The two are perpendicular in a polar orbit: geodetic in the orbit
	// plane's normal, frame dragging along the spin axis's projection.
	if c := geo.dot(fd) / (geo.norm() * fd.norm()); math.Abs(c) > 1e-6 {
		t.Errorf("geodetic and frame-dragging directions not perpendicular (cos %.3g)", c)
	}
	// Measured values within their errors (2σ, since the model has no J2).
	if math.Abs(g-6601.8) > 2*18.3 || math.Abs(f-37.2) > 2*7.2 {
		t.Errorf("outside GP-B measurement: %.1f, %.2f", g, f)
	}
}

// TestRimRateTwin: the exact rim rate (Schwarzschild hub × local special
// relativity) equals the weak-field metric's first-order rate to second
// order in the small quantities.
func TestRimRateTwin(t *testing.T) {
	for _, c := range []struct{ r, vHub, omega, rho float64 }{
		{REarthEq, 465, 1000, 0.3},
		{REarthEq + 642e3, 7500, 6000, 0.25},
		{1e9, 0, 1e5, 0.01},
	} {
		ex, err := RimRateInField(Earth, c.r, c.vHub, c.omega, c.rho)
		if err != nil {
			t.Fatal(err)
		}
		wk := RimRateWeak(Earth, c.r, c.vHub, c.omega, c.rho)
		small := Earth.GM/(c.r*C*C) + (c.vHub*c.vHub+c.omega*c.omega*c.rho*c.rho)/(C*C)
		if math.Abs(ex-wk) > 2*small*small+4e-16 { // plus float64 rounding of numbers near 1
			t.Errorf("r %g: exact %.18f weak %.18f", c.r, ex, wk)
		}
	}
}

// TestSagnacEarthGyro: a 1 m² horizontal loop on the ground at 45°
// latitude sees Earth's rotation component along the vertical, 4Aω sinφ/c²,
// and the local frame's precession changes it by about 10⁻⁹ of that
// (geodetic and Lense–Thirring terms for an Earth-bound gyro are ~10⁻⁹ of
// Earth's rate).
func TestSagnacEarthGyro(t *testing.T) {
	lat := math.Pi / 4
	rep, err := AnalyzeMachine(Earth, EarthSurface(lat), Machine{Name: "ring", Radius: 0.5, Omega: 0, Loop: 1})
	if err != nil {
		t.Fatal(err)
	}
	want := 4 * 1 * OmegaEarth * math.Sin(lat) / (C * C)
	if math.Abs(rep.SagnacPlatform-want)/want > 1e-12 {
		t.Errorf("Sagnac %.6g s, want %.6g", rep.SagnacPlatform, want)
	}
	ratio := math.Abs(rep.SagnacFrameCorrection / rep.SagnacPlatform)
	t.Logf("Earth-surface ring: Sagnac %.4g s; local-frame correction %.3g of it", rep.SagnacPlatform, ratio)
	if ratio > 1e-8 || ratio < 1e-11 {
		t.Errorf("frame correction ratio %.3g outside the expected 1e-11–1e-8", ratio)
	}
}
