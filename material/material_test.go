// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package material

import (
	"math"
	"testing"
)

func TestKnownValues(t *testing.T) {
	// Bi–Al thermocouple: 75.5 µV/K, so 7.55 mV at ΔT = 100 K.
	if e := ThermocoupleEMF(Bismuth, Aluminium, 100); math.Abs(e-7.55e-3) > 1e-12 {
		t.Errorf("Bi–Al EMF at 100 K = %g V", e)
	}
	// Al anode, Bi cathode: 0.308 − (−1.662) = 1.970 V.
	if e := CellEMF(Aluminium, Bismuth); math.Abs(e-1.970) > 1e-12 {
		t.Errorf("Al|Bi cell EMF = %g V", e)
	}
	// 1 g of aluminium: 3·96485/26.98 = 10 728 C ≈ 2.98 Ah.
	if q := ChargeCapacity(Aluminium, 1e-3); math.Abs(q/3600-2.98) > 0.005 {
		t.Errorf("1 g Al = %g Ah", q/3600)
	}
	if math.Abs(Mu0*Eps0*299792458*299792458-1) > 1e-9 {
		t.Error("µ0·ε0·c² ≠ 1")
	}
}

func TestSpectralLines(t *testing.T) {
	// Air index near 500 nm is ≈ 1.000279 (Birch & Downs).
	if n := AirIndex(2); math.Abs(n-1.000279) > 2e-6 {
		t.Errorf("n_air(500 nm) = %.7f", n)
	}
	// Cu 324.754 nm (air) → 324.848 nm (vacuum), f ≈ 922.9 THz.
	v := LineCopper.VacuumWavelength()
	if math.Abs(v-324.848e-9) > 0.005e-9 {
		t.Errorf("Cu vacuum wavelength %g", v)
	}
	if f := LineCopper.Frequency(); math.Abs(f/1e12-922.87) > 0.1 {
		t.Errorf("Cu line %g THz", f/1e12)
	}
	n, off := LineAluminium.Subharmonic(1e6)
	if math.Abs(off) > 0.5/n+1e-18 || n < 7e8 {
		t.Errorf("n = %g, offset %g", n, off)
	}
}
