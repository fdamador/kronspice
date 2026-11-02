// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// hoffmannMotor is Hoffmann's single-phase squirrel-cage motor in frame I
// (1949, eq. 42): stator coil 1, rotor coils 2, 3 turning with the rotor.
func hoffmannMotor(l11, m, ls, ld float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(3)
	s.L0[0][0], s.L0[1][1], s.L0[2][2] = l11, ls, ls
	h1 := s.Harmonic(1)
	h1.Cos[0][1], h1.Cos[1][0] = m, m
	h1.Sin[0][2], h1.Sin[2][0] = -m, -m
	h2 := s.Harmonic(2)
	h2.Cos[1][1], h2.Cos[2][2] = -ld, ld
	h2.Sin[1][2], h2.Sin[2][1] = ld, ld
	return s
}

// TestHoffmannTwin: eq. 5′ in the stationary-coil frame, with the torsion
// of eq. 19, gives the frame-I Lagrange–Maxwell forces transformed by
// C_λ^α (voltages Cᵀ·e, the same torque), at random states; without the
// torsion term it does not ("wrong results", p. 540).
func TestHoffmannTwin(t *testing.T) {
	rng := rand.New(rand.NewSource(1949))
	for _, ld := range []float64{0, 0.01} {
		s := hoffmannMotor(0.25, 0.23, 0.25, ld)
		f := RotorFrame(3, [][2]int{{1, 2}}, 1)
		const J = 1e-3
		worst, wrong := 0.0, 0.0
		for k := 0; k < 50; k++ {
			th := 2 * math.Pi * rng.Float64()
			ip := []float64{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
			dip := []float64{rng.NormFloat64() * 100, rng.NormFloat64() * 100, rng.NormFloat64() * 100}
			om, dom := rng.NormFloat64()*300, rng.NormFloat64()*1000
			c, d := f.C(th), f.D(th)
			i, di := make([]float64, 3), make([]float64, 3)
			for a := 0; a < 3; a++ {
				for b := 0; b < 3; b++ {
					i[a] += c[a][b] * ip[b]
					di[a] += c[a][b]*dip[b] + d[a][b]*ip[b]*om
				}
			}
			old := Geodesic(s, J, th, i, di, om, dom)
			want := make([]float64, 4)
			for la := 0; la < 3; la++ {
				for a := 0; a < 3; a++ {
					want[la] += c[a][la] * old[a]
				}
			}
			want[3] = old[3]
			got, err := NonRiemannian(s, J, f, th, ip, dip, om, dom, true)
			if err != nil {
				t.Fatal(err)
			}
			bad, _ := NonRiemannian(s, J, f, th, ip, dip, om, dom, false)
			scale := 0.0
			for _, v := range want {
				scale = math.Max(scale, math.Abs(v))
			}
			for la := range want {
				worst = math.Max(worst, math.Abs(got[la]-want[la])/scale)
				wrong = math.Max(wrong, math.Abs(bad[la]-want[la])/scale)
			}
		}
		if worst > 1e-12 {
			t.Errorf("L_D = %g: eq. 5′ differs from frame I by %g", ld, worst)
		}
		if wrong < 0.1 {
			t.Errorf("L_D = %g: leaving out the torsion should be wrong, differs only %g", ld, wrong)
		}
	}
}

// TestHoffmannTorsion: Ω is antisymmetric in its first two indices and the
// new metric is constant in the stationary-coil frame (eq. 29).
func TestHoffmannTorsion(t *testing.T) {
	s := hoffmannMotor(0.25, 0.23, 0.25, 0.01)
	f := RotorFrame(3, [][2]int{{1, 2}}, 1)
	var first [][]float64
	for _, th := range []float64{0, 0.3, 1.7, 4} {
		om, err := Torsion(s, 1e-3, f, th)
		if err != nil {
			t.Fatal(err)
		}
		for a := 0; a < 4; a++ {
			for b := 0; b < 4; b++ {
				for c := 0; c < 4; c++ {
					if math.Abs(om[a][b][c]+om[b][a][c]) > 1e-15 {
						t.Fatalf("Ω not antisymmetric at %d%d%d", a, b, c)
					}
				}
			}
		}
		lb, dlb, _ := NewMetric(s, 1e-3, f, th)
		for a := range lb {
			for b := range lb {
				if math.Abs(dlb[a][b]) > 1e-12 || (first != nil && math.Abs(lb[a][b]-first[a][b]) > 1e-12) {
					t.Fatalf("new metric varies with θ at %d,%d", a, b)
				}
			}
		}
		first = lb
	}
}
