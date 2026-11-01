// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import (
	"errors"
	"math"
	"testing"
)

const day = 86400.0

// GPS (N. Ashby, "Relativity in the Global Positioning System," Living
// Rev. Relativ. 6, 1, 2003, §5): a satellite clock on a 26,561.75 km
// orbit gains 45.7 µs/day from the weaker potential, loses 7.1 µs/day
// from its speed, net 38.6 µs/day against a clock on the rotating geoid.
// The equator clock here (spherical Earth, with rotation) stands in for
// the geoid; Earth's oblateness (J2) is left out, which shifts the result
// by about 0.03 µs/day.
func TestGPSClock(t *testing.T) {
	a := 26561.75e3
	sat, err := OrbitRate(GMEarth, a)
	if err != nil {
		t.Fatal(err)
	}
	ground, _ := MovingRate(GMEarth, REarthEq, OmegaEarth*REarthEq)
	net := (sat/ground - 1) * day * 1e6
	if math.Abs(net-38.6) > 0.15 {
		t.Errorf("net %.3f µs/day, Ashby 38.6", net)
	}
	gs, _ := StaticRate(GMEarth, a)
	g0, _ := StaticRate(GMEarth, REarthEq)
	if grav := (gs/g0 - 1) * day * 1e6; math.Abs(grav-45.7) > 0.1 {
		t.Errorf("gravitational %.3f µs/day, Ashby 45.7", grav)
	}
	t.Logf("GPS satellite clock: %+.2f µs/day net (Ashby 2003: +38.6)", net)
}

// Pound–Rebka–Snider: 22.5 m tower, predicted fractional shift gh/c²
// (each way) = 2.46×10⁻¹⁵; measured/predicted = 0.9990 ± 0.0076
// (R. V. Pound, J. L. Snider, Phys. Rev. 140, B788, 1965).
func TestPoundRebka(t *testing.T) {
	h := 22.5
	r := REarthEq
	shift, err := StaticShift(GMEarth, r+h, r) // emitted at the top, received at the bottom: blue shift
	if err != nil {
		t.Fatal(err)
	}
	g := GMEarth / (r * r)
	want := g * h / (C * C) * r / (r + h) // GM(1/r − 1/(r+h))/c² to first order
	if math.Abs(shift/want-1) > 1e-5 {
		t.Errorf("shift %g, gh/c² %g", shift, want)
	}
	if math.Abs(want-2.45e-15) > 0.02e-15 {
		t.Errorf("gh/c² = %g, expected ≈ 2.45e-15", want)
	}
}

// Shapiro delay: the closed form equals the path integral of the weak-
// field index; the exact isotropic index agrees in the weak field.
func TestShapiroDelay(t *testing.T) {
	b, x1, x2 := 2*RSun, 1.5e11, 2e11 // Earth-ish to Mars-ish geometry past the Sun
	want := ShapiroStraight(GMSun, b, x1, x2, 1)
	got := PathDelay(func(r float64) float64 { return IndexWeak(GMSun, r, 1) }, b, x1, x2, 4000)
	if math.Abs(got/want-1) > 1e-8 {
		t.Errorf("path integral %g s, closed form %g s", got, want)
	}
	iso := PathDelay(func(r float64) float64 { return IndexIsotropic(GMSun, r) }, b, x1, x2, 4000)
	if math.Abs(iso/want-1) > 1e-5 {
		t.Errorf("isotropic index %g s vs weak %g s", iso, want)
	}
	t.Logf("one-way Shapiro delay past the Sun at 2 R☉: %.1f µs", want*1e6)
}

// Light bending at the solar limb: 4GM/(R c²) = 1.751″; VLBI measures
// (1+γ)/2 = 0.99992 ± 0.00023 (S. S. Shapiro et al., Phys. Rev. Lett. 92,
// 121101, 2004).
func TestSolarDeflection(t *testing.T) {
	d := Deflection(func(r float64) float64 { return IndexWeak(GMSun, r, 1) }, RSun, 4000)
	want := 4 * GMSun / (RSun * C * C)
	if math.Abs(d/want-1) > 1e-4 {
		t.Errorf("deflection %g rad, 4GM/bc² %g", d, want)
	}
	if as := d * 180 / math.Pi * 3600; math.Abs(as-1.751) > 0.002 {
		t.Errorf("%.4f arcsec, expected 1.751", as)
	}
}

func TestHorizonAndOrbitErrors(t *testing.T) {
	rs := SchwarzschildRadius(GMSun)
	if _, err := StaticRate(GMSun, rs); !errors.Is(err, ErrInsideHorizon) {
		t.Errorf("at r_s: %v", err)
	}
	if _, err := OrbitRate(GMSun, 1.4*rs); !errors.Is(err, ErrNoOrbit) {
		t.Errorf("inside photon sphere: %v", err)
	}
	// r = 3 r_s static: √(2/3).
	if r, _ := StaticRate(GMSun, 3*rs); math.Abs(r-math.Sqrt(2.0/3)) > 1e-15 {
		t.Errorf("rate at 3 r_s = %g", r)
	}
	if _, err := MovingRate(GMEarth, REarthEq, C); !errors.Is(err, ErrSuperluminal) {
		t.Errorf("v = c: %v", err)
	}
}

// Tidal stretching at Earth's surface: 2GM/R³ ≈ 3.08×10⁻⁶ (m/s²) per metre.
func TestTidal(t *testing.T) {
	if a := Tidal(GMEarth, REarthEq, 1); math.Abs(a-3.08e-6) > 0.02e-6 {
		t.Errorf("tidal %g", a)
	}
}

// Rotor clocks: for source at the centre and absorber at the rim, the
// standard k = 1/2 shift equals the exact rim time dilation to order
// (v/c)⁴ (Hay, Schiffer, Cranshaw, Egelstaff, Phys. Rev. Lett. 4, 165,
// 1960). The reported anomaly k ≈ 0.68 is 36 % larger.
func TestRotorShift(t *testing.T) {
	omega, r := 2*math.Pi*500, 0.093 // ≈ Kündig's rotor scale
	b2 := omega * omega * r * r / (C * C)
	std := RotorShift(omega, 0, r, 0.5)
	if exact := math.Expm1(-0.5 * math.Log1p(-b2)); math.Abs(std/exact-1) > 1e-10 { // 1/rate − 1
		t.Errorf("k=1/2 shift %g vs exact %g", std, exact)
	}
	anom := RotorShift(omega, 0, r, 0.68)
	if math.Abs(anom/std-1.36) > 1e-12 {
		t.Errorf("anomaly ratio %g", anom/std)
	}
}

// Sagnac: the exact counter-propagation difference reduces to 4Aω/c².
func TestSagnac(t *testing.T) {
	r, omega := 5.0, 100.0
	ex, err := SagnacRingExact(r, omega)
	if err != nil {
		t.Fatal(err)
	}
	if ap := SagnacDelay(math.Pi*r*r, omega); math.Abs(ex/ap-1) > 1e-9 {
		t.Errorf("exact %g vs 4Aω/c² %g", ex, ap)
	}
	// Ring splitting = Sagnac delay × f / (round-trip time): Δf = 4Aω/(λP).
	lambda := 8.0
	p := 2 * math.Pi * r
	if df, want := SagnacSplitting(math.Pi*r*r, p, lambda, omega), SagnacDelay(math.Pi*r*r, omega)*(C/lambda)/(p/C); math.Abs(df/want-1) > 1e-12 {
		t.Errorf("splitting %g vs %g", df, want)
	}
}

// Burst limit: T1000G fibre hoop bursts at √(6.37 GPa / 1800 kg/m³) ≈
// 1,881 m/s, where (v/c)² ≈ 3.9×10⁻¹¹ — rotating machines are deep in the
// weak-field regime.
func TestBurstAndRotor(t *testing.T) {
	v := BurstSpeed(T1000G)
	if math.Abs(v-1881.2) > 1 {
		t.Errorf("burst speed %g", v)
	}
	rot, err := AnalyzeRotor(v/0.5, 0.5, T1000G)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(rot.BurstMargin-1) > 1e-12 || math.Abs(rot.Beta2-3.94e-11) > 0.01e-11 {
		t.Errorf("%+v", rot)
	}
	if _, err := RimRate(C, 1); !errors.Is(err, ErrSuperluminal) {
		t.Error("rim at c should be rejected")
	}
	t.Logf("T1000G hoop at burst: v = %.0f m/s, (v/c)² = %.3g, rim clock lags %.2f µs/day", rot.RimSpeed, rot.Beta2, rot.LossPerDay*1e6)
}
