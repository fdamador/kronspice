// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// A machine whose inductances do not vary with θ is flat.
func TestCurvatureFlat(t *testing.T) {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0], s.L0[1][1], s.L0[0][1], s.L0[1][0] = 3e-3, 4e-3, 1e-3, 1e-3
	r, err := Curvature(s, 1e-3, 0.7)
	if err != nil {
		t.Fatal(err)
	}
	for a := range r {
		for b := range r[a] {
			for c := range r[a][b] {
				for d := range r[a][b][c] {
					if r[a][b][c][d] != 0 {
						t.Fatalf("R[%d][%d][%d][%d] = %g", a, b, c, d, r[a][b][c][d])
					}
				}
			}
		}
	}
}

// One salient-pole winding, L(θ) = L₀ + a·cos 2θ, with inertia J: the
// metric L(θ)·dq² + J·dθ² is a surface of Gaussian curvature
// K = −f″/(f·J), f = √L, so R_qθqθ = K·L·J; the tensor has the
// antisymmetries R_abcd = −R_bacd = −R_abdc and pair symmetry.
func TestCurvatureSalientWinding(t *testing.T) {
	const l0, a, J = 10e-3, 3e-3, 2e-4
	s := circuit.NewMachineSpec(1)
	s.L0[0][0] = l0
	s.Harmonic(2).Cos[0][0] = a
	for _, th := range []float64{0, 0.3, 1.0, 2.2} {
		r, err := Curvature(s, J, th)
		if err != nil {
			t.Fatal(err)
		}
		E := l0 + a*math.Cos(2*th)
		E1 := -2 * a * math.Sin(2*th)
		E2 := -4 * a * math.Cos(2*th)
		f := math.Sqrt(E)
		f2 := E2/(2*f) - E1*E1/(4*f*f*f)
		K := -f2 / (f * J)
		want := K * E * J
		if got := r[0][1][0][1]; math.Abs(got-want) > 1e-12*math.Abs(want)+1e-18 {
			t.Errorf("θ=%g: R_qθqθ = %g, want K·det g = %g", th, got, want)
		}
		if r[1][0][0][1] != -r[0][1][0][1] || r[0][1][1][0] != -r[0][1][0][1] || r[1][0][1][0] != r[0][1][0][1] {
			t.Errorf("θ=%g: symmetries broken", th)
		}
	}
}
