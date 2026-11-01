// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/rand/v2"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/hyper/quat"
)

func apply3(m [3][3]float64, x [3]float64) [3]float64 {
	var y [3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			y[i] += m[i][j] * x[j]
		}
	}
	return y
}

// The Clarke and Park rotors equal the standard matrices, and r·x·r*
// equals the matrix transform of random three-phase values (CLAUDE.md §3b
// rule 1: quaternion Clarke/Park vs the matrix transforms, 1e−9).
func TestClarkeParkRotors(t *testing.T) {
	cm := ClarkeMatrix()
	if d := maxDiff3(ClarkeRotor().Matrix(), cm); d > 1e-15 {
		t.Errorf("Clarke rotor matrix differs by %g", d)
	}
	r := rand.New(rand.NewPCG(21, 22))
	for n := 0; n < 500; n++ {
		th := 4 * math.Pi * (r.Float64() - 0.5)
		pm := ParkMatrix(th)
		if d := maxDiff3(ParkRotor(th).Matrix(), pm); d > 1e-14 {
			t.Fatalf("Park rotor matrix differs by %g at θ = %g", d, th)
		}
		x := [3]float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		want := apply3(pm, x)
		got := Rotate3(ParkRotor(th), Phases(x[0], x[1], x[2]))
		if math.Abs(got.X-want[0])+math.Abs(got.Y-want[1])+math.Abs(got.Z-want[2]) > 1e-14*(1+math.Abs(want[0])+math.Abs(want[1])+math.Abs(want[2])) || math.Abs(got.W) > 1e-14 {
			t.Fatalf("r·x·r* = %v, matrix gives %v", got, want)
		}
	}
	// A balanced set √2·V·cos(ωt + φ − 2πk/3) is constant in dq0 at the
	// synchronous angle ωt: d + jq = √3·V·e^(jφ), 0 = 0.
	V, phi := 120.0, 0.3
	for n := 0; n < 50; n++ {
		wt := 2 * math.Pi * r.Float64()
		x := Phases(math.Sqrt2*V*math.Cos(wt+phi), math.Sqrt2*V*math.Cos(wt+phi-2*math.Pi/3), math.Sqrt2*V*math.Cos(wt+phi+2*math.Pi/3))
		dq := Rotate3(ParkRotor(wt), x)
		if math.Abs(dq.X-math.Sqrt(3)*V*math.Cos(phi)) > 1e-12*V || math.Abs(dq.Y-math.Sqrt(3)*V*math.Sin(phi)) > 1e-12*V || math.Abs(dq.Z) > 1e-12*V {
			t.Fatalf("balanced set in dq0 = (%g, %g, %g)", dq.X, dq.Y, dq.Z)
		}
	}
	// Zero sequence lies on the 0 axis alone: (1, 1, 1) → (0, 0, √3).
	z := Rotate3(ClarkeRotor(), Phases(1, 1, 1))
	if math.Abs(z.X)+math.Abs(z.Y) > 1e-15 || math.Abs(z.Z-math.Sqrt(3)) > 1e-15 {
		t.Errorf("zero sequence → %v", z)
	}
}

func maxDiff3(a, b [3][3]float64) float64 {
	m := 0.0
	for i := range a {
		for j := range a[i] {
			m = math.Max(m, math.Abs(a[i][j]-b[i][j]))
		}
	}
	return m
}

// Quaternion power: S.s = Σ v·i, V.s = −v×i, and both are frame
// invariants (p equal and |V.s| equal in abc, αβ0 and dq0); for a balanced
// sinusoidal set p = 3·V·I·cos φ and |V.s| = 3·V·I·sin φ, both constant
// (Akagi's p and q).
func TestQuaternionPower(t *testing.T) {
	r := rand.New(rand.NewPCG(23, 24))
	for n := 0; n < 500; n++ {
		va, vb, vc := r.NormFloat64(), r.NormFloat64(), r.NormFloat64()
		ia, ib, ic := r.NormFloat64(), r.NormFloat64(), r.NormFloat64()
		v, i := Phases(va, vb, vc), Phases(ia, ib, ic)
		s := InstantPower(v, i)
		p := va*ia + vb*ib + vc*ic
		cx, cy, cz := -(vb*ic - vc*ib), -(vc*ia - va*ic), -(va*ib - vb*ia)
		if math.Abs(s.W-p) > 1e-14*(1+math.Abs(p)) || math.Abs(s.X-cx)+math.Abs(s.Y-cy)+math.Abs(s.Z-cz) > 1e-14*(1+math.Abs(p)) {
			t.Fatalf("s = %v, want p %g and −v×i (%g, %g, %g)", s, p, cx, cy, cz)
		}
		th := 2 * math.Pi * r.Float64()
		for _, rot := range []quat.Q{ClarkeRotor(), ParkRotor(th)} {
			s2 := InstantPower(Rotate3(rot, v), Rotate3(rot, i))
			if math.Abs(s2.W-s.W) > 1e-13*(1+math.Abs(p)) {
				t.Fatalf("p not invariant: %g vs %g", s2.W, s.W)
			}
			want := Rotate3(rot, quat.Vec(s.X, s.Y, s.Z))
			if math.Abs(s2.X-want.X)+math.Abs(s2.Y-want.Y)+math.Abs(s2.Z-want.Z) > 1e-13*(1+s.Abs()) {
				t.Fatalf("non-active vector does not rotate with the frame")
			}
		}
	}
	V, I, phi := 120.0, 10.0, 0.4
	for n := 0; n < 50; n++ {
		wt := 2 * math.Pi * r.Float64()
		ph := func(a, k float64) float64 { return math.Sqrt2 * a * math.Cos(wt-2*math.Pi*k/3) }
		v := Phases(ph(V, 0), ph(V, 1), ph(V, 2))
		ii := func(k float64) float64 { return math.Sqrt2 * I * math.Cos(wt-phi-2*math.Pi*k/3) }
		s := InstantPower(v, Phases(ii(0), ii(1), ii(2)))
		if math.Abs(s.W-3*V*I*math.Cos(phi)) > 1e-10*V*I || math.Abs(math.Sqrt(s.X*s.X+s.Y*s.Y+s.Z*s.Z)-3*V*I*math.Sin(phi)) > 1e-10*V*I {
			t.Fatalf("balanced p = %g, |q| = %g; want %g, %g", s.W, math.Sqrt(s.X*s.X+s.Y*s.Y+s.Z*s.Z), 3*V*I*math.Cos(phi), 3*V*I*math.Sin(phi))
		}
	}
}

// salientStator is the stator of a salient-pole machine (P. C. Krause,
// O. Wasynczuk, S. D. Sudhoff, Analysis of Electric Machinery, 2nd ed.,
// IEEE Press 2002, Ch. 5, the stator self- and mutual inductances of a
// salient-pole synchronous machine; θ the electrical rotor angle): with
// phase angles φ = 0, 2π/3, −2π/3, L_xx = Ls + Lm·cos 2(θ − φx) and
// L_xy = −Ms + Lm·cos 2(θ − (φx + φy)/2), so L_xy = −Ms − Lm·cos 2(θ + π/6)
// for a, b as in Krause.
func salientStator(ls, ms, lm float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(3)
	h := s.Harmonic(2)
	ang := []float64{0, 2 * math.Pi / 3, -2 * math.Pi / 3}
	for x := 0; x < 3; x++ {
		for y := 0; y < 3; y++ {
			// cos 2(θ − (φx + φy)/2) = cos2θ·cos(φx+φy) + sin2θ·sin(φx+φy)
			a := ang[x] + ang[y]
			if x == y {
				s.L0[x][y] = ls
				h.Cos[x][y] = lm * math.Cos(a)
				h.Sin[x][y] = lm * math.Sin(a)
			} else {
				s.L0[x][y] = -ms
				h.Cos[x][y] = lm * math.Cos(a)
				h.Sin[x][y] = lm * math.Sin(a)
			}
		}
	}
	return s
}

// Kron's rotating frame with the quaternion Park frame: a salient-pole
// stator's θ-dependent abc inductances become the constant
// Ld = Ls + Ms + (3/2)Lm, Lq = Ls + Ms − (3/2)Lm, L0 = Ls − 2Ms (Krause
// et al. Ch. 5; the same in power-invariant axes), and the speed matrix
// carries the speed voltages e_d = −ω·Lq·i_q, e_q = ω·Ld·i_d — checked
// against flux linkage differentiated numerically in the abc frame.
func TestParkFrameSalientMachine(t *testing.T) {
	ls, ms, lm := 10e-3, 4e-3, 2e-3
	s := salientStator(ls, ms, lm)
	if err := s.CheckPositive(); err != nil {
		t.Fatal(err)
	}
	f := ParkFrame(3, [][3]int{{0, 1, 2}}, 1)
	// D is dC/dθ.
	for _, th := range []float64{0.1, 1.3, 2.9} {
		h := 1e-6
		cp, cm, d := f.C(th+h), f.C(th-h), f.D(th)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				if fd := (cp[i][j] - cm[i][j]) / (2 * h); math.Abs(fd-d[i][j]) > 1e-8 {
					t.Fatalf("D differs from dC/dθ at (%d,%d): %g vs %g", i, j, d[i][j], fd)
				}
			}
		}
	}
	out, err := RotatingFrame(s, f)
	if err != nil {
		t.Fatal(err)
	}
	ld, lq, l0 := ls+ms+1.5*lm, ls+ms-1.5*lm, ls-2*ms
	want := [][]float64{{ld, 0, 0}, {0, lq, 0}, {0, 0, l0}}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if math.Abs(out.L0[i][j]-want[i][j]) > 1e-15 {
				t.Errorf("L′[%d][%d] = %g, want %g", i, j, out.L0[i][j], want[i][j])
			}
		}
	}
	// Speed voltages: with steady i_dq0 and ω, e_abc = d/dt(L(θ)·Pᵀ·i′),
	// e′ = P·e_abc must equal ω·G′·i′.
	ip := []float64{3, -2, 0.5}
	omega, th := 377.0, 0.7
	psi := func(th float64) []float64 {
		l := square(3)
		s.Inductance(th, l, nil, nil)
		c := f.C(th)
		iabc := make([]float64, 3)
		for r := 0; r < 3; r++ {
			for k := 0; k < 3; k++ {
				iabc[r] += c[r][k] * ip[k]
			}
		}
		out := make([]float64, 3)
		for r := 0; r < 3; r++ {
			for k := 0; k < 3; k++ {
				out[r] += l[r][k] * iabc[k]
			}
		}
		return out
	}
	h := 1e-6
	a, b := psi(th+h), psi(th-h)
	c := f.C(th)
	for k := 0; k < 3; k++ {
		ek := 0.0
		for r := 0; r < 3; r++ {
			ek += c[r][k] * omega * (a[r] - b[r]) / (2 * h)
		}
		gk := 0.0
		for j := 0; j < 3; j++ {
			gk += omega * out.G[k][j] * ip[j]
		}
		if math.Abs(ek-gk) > 1e-6*math.Abs(omega*ld*3) {
			t.Errorf("speed voltage on axis %d: %g from flux, %g from G′", k, ek, gk)
		}
	}
	if math.Abs(out.G[0][1]+lq) > 1e-15 || math.Abs(out.G[1][0]-ld) > 1e-15 {
		t.Errorf("G′ = %v, want e_d = −ωLq·i_q and e_q = ωLd·i_d", out.G)
	}
}
