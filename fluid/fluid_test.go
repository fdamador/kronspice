// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package fluid

import (
	"math"
	"math/cmplx"
	"testing"
)

func TestMedia(t *testing.T) {
	if math.Abs(Water.Impedance()-1.4963e6) > 1e3 {
		t.Errorf("ρc of water %g rayl", Water.Impedance())
	}
	// D₂O: denser (≈ 1104 vs 997 kg/m³), slower (1404 vs 1502 m/s).
	if !(HeavyWater.Density > 1.1*Water.Density && HeavyWater.SoundSpeed < 0.94*Water.SoundSpeed) {
		t.Errorf("D2O %+v against H2O %+v", HeavyWater, Water)
	}
}

// The exact piston impedance has the known limits: (ka)²/2 and 8ka/(3π) for
// ka ≪ 1, and ρc/S with vanishing reactance for ka ≫ 1; the lumped network
// matches the exact impedance to O((ka)²).
func TestPistonRadiation(t *testing.T) {
	a := 0.01
	for _, m := range []Medium{Water, HeavyWater} {
		S := math.Pi * a * a
		zc := m.Impedance() / S
		Ma, Rpar, fmax, err := PistonRadiation(a, m)
		if err != nil {
			t.Fatal(err)
		}
		for _, ka := range []float64{0.01, 0.05, 0.2} {
			f := ka * m.SoundSpeed / (2 * math.Pi * a)
			z := PistonImpedance(a, m, f)
			if d := math.Abs(real(z)/zc - ka*ka/2); d > 0.2*ka*ka*ka*ka+1e-15 {
				t.Errorf("%s ka = %g: R1 = %g, (ka)²/2 = %g", m.Name, ka, real(z)/zc, ka*ka/2)
			}
			if d := math.Abs(imag(z)/zc - 8*ka/(3*math.Pi)); d > 0.25*ka*ka*ka { // next series term: (32/45π)·(ka)³
				t.Errorf("%s ka = %g: X1 = %g, 8ka/3π = %g", m.Name, ka, imag(z)/zc, 8*ka/(3*math.Pi))
			}
			w := 2 * math.Pi * f
			zn := complex(0, w*Ma) * complex(Rpar, 0) / complex(Rpar, w*Ma)
			if e := cmplx.Abs(zn-z) / cmplx.Abs(z); e > 2*ka*ka {
				t.Errorf("%s ka = %g: network %v, exact %v (%.2g)", m.Name, ka, zn, z, e)
			}
		}
		if hi := PistonImpedance(a, m, 3*m.SoundSpeed/(2*math.Pi*a)); math.Abs(real(hi)/zc-1) > 0.15 {
			t.Errorf("%s ka = 3: R1 = %g, want ≈ 1", m.Name, real(hi)/zc)
		}
		if math.Abs(fmax*2*math.Pi*a/m.SoundSpeed-0.5) > 1e-12 {
			t.Errorf("fmax %g", fmax)
		}
	}
	if _, _, _, err := PistonRadiation(0, Water); err == nil {
		t.Error("zero radius accepted")
	}
}

// The breathing mode's lumped compliance equals the mode's potential
// energy integrated numerically, and the frequency is c/(2R): 25.0 kHz for
// a 3 cm water sphere, lower in D₂O.
func TestSphericalResonator(t *testing.T) {
	R := 0.03
	f0, Ca, Ma, Rpar, err := SphericalResonator(R, Water, 1000)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(f0-Water.SoundSpeed/(2*R)) > 1e-9 || math.Abs(f0-25025) > 10 {
		t.Errorf("f0 = %g Hz", f0)
	}
	k := math.Pi / R
	n := 20000
	sum := 0.0
	for i := 0; i < n; i++ {
		r := (float64(i) + 0.5) * R / float64(n)
		p := math.Sin(k*r) / (k * r)
		sum += p * p * 4 * math.Pi * r * r * R / float64(n)
	}
	W := sum / (2 * Water.Density * Water.SoundSpeed * Water.SoundSpeed) // for P = 1 Pa
	if math.Abs(W-0.5*Ca) > 1e-6*W {
		t.Errorf("∫p²/(2ρc²)dV = %g J, ½·Ca = %g J", W, 0.5*Ca)
	}
	w0 := 2 * math.Pi * f0
	if math.Abs(w0*w0*Ma*Ca-1) > 1e-12 || math.Abs(Rpar*w0*Ca-1000) > 1e-9 {
		t.Errorf("Ma, Rpar inconsistent")
	}
	fd, _, _, _, _ := SphericalResonator(R, HeavyWater, 1000)
	if !(fd < 0.94*f0) {
		t.Errorf("D2O f0 = %g Hz, H2O %g Hz", fd, f0)
	}
}
