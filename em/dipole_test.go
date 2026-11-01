// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"testing"
)

// Triple twin (§3c rule 2) for the radiating dipole: the Poynting flux of
// Jackson's full fields (near and far terms) through a sphere, computed
// with the GA product ½FF† and with the vector cross product, equals
// Larmor's formula. The flux is independent of the radius (no field energy
// is lost between spheres in vacuum), so it is checked in the near zone
// (kR = 0.5) and the far zone (kR = 50).
func TestDipoleRadiation(t *testing.T) {
	p := [3]complex128{complex(1e-30, 0), complex(0, 2e-30), complex(0.5e-30, 0.3e-30)} // elliptical, C·m
	omega := 2 * math.Pi * SpeedOfLight / 6e-6                                          // 6 µm light (~1667 cm⁻¹)
	k := omega / SpeedOfLight
	want := LarmorPower(p, omega)
	for _, kr := range []float64{0.5, 50} {
		R := kr / k
		ga := RadiatedPower(p, omega, R, 24, 48, true)
		vec := RadiatedPower(p, omega, R, 24, 48, false)
		if math.Abs(ga-want) > 1e-9*want || math.Abs(vec-want) > 1e-9*want {
			t.Errorf("kR=%g: GA %.12g, vector %.12g, Larmor %.12g W", kr, ga, vec, want)
		}
	}
}

// Far from the dipole the field is a transverse wave with |E| = c|B| and
// both perpendicular to the direction of travel (lumped radiator, §3g
// standard prediction).
func TestDipoleFarField(t *testing.T) {
	p := [3]complex128{0, 0, complex(1e-30, 0)}
	omega := 2 * math.Pi * SpeedOfLight / 6e-6
	r := [3]float64{0.3, 0.4, 0.2}
	E, B := DipoleFields(p, r, omega)
	e2, b2, en := 0.0, 0.0, complex(0, 0)
	rr := math.Sqrt(0.29)
	for d := 0; d < 3; d++ {
		e2 += real(E[d] * complex(real(E[d]), -imag(E[d])))
		b2 += real(B[d] * complex(real(B[d]), -imag(B[d])))
		en += E[d] * complex(r[d]/rr, 0)
	}
	if math.Abs(math.Sqrt(e2/b2)/SpeedOfLight-1) > 1e-5 {
		t.Errorf("|E|/(c|B|) = %g", math.Sqrt(e2/b2)/SpeedOfLight)
	}
	if math.Hypot(real(en), imag(en)) > 1e-5*math.Sqrt(e2) {
		t.Errorf("longitudinal E %g of %g", math.Hypot(real(en), imag(en)), math.Sqrt(e2))
	}
}
