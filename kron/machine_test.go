// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/rand/v2"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// primitiveMachine is Kron's two-phase primitive machine: stator windings
// α, β (1, 2) fixed, rotor windings a, b (3, 4) turning with the rotor at
// electrical angle p·θ; the stator–rotor mutuals form M·R(p·θ).
func primitiveMachine(ls, lr, m float64, p int) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(4)
	s.L0[0][0], s.L0[1][1], s.L0[2][2], s.L0[3][3] = ls, ls, lr, lr
	h := s.Harmonic(p)
	set := func(mat [][]float64, i, j int, v float64) { mat[i][j], mat[j][i] = v, v }
	// L_αa = M cos, L_αb = −M sin, L_βa = M sin, L_βb = M cos.
	set(h.Cos, 0, 2, m)
	set(h.Sin, 0, 3, -m)
	set(h.Sin, 1, 2, m)
	set(h.Cos, 1, 3, m)
	return s
}

// Kron's eq. 3 with the Christoffel symbols of eq. 4 reproduces the N
// element's frame-I equations: winding voltages L·di/dt + ω·L′·i and the
// shaft torque J·dω/dt − ½·iᵀ·L′·i, at random states.
func TestChristoffelEquations(t *testing.T) {
	s := primitiveMachine(10e-3, 12e-3, 8e-3, 2)
	r := rand.New(rand.NewPCG(3, 4))
	const J = 2e-3
	for trial := 0; trial < 50; trial++ {
		th, om, dom := r.Float64()*7, r.NormFloat64()*100, r.NormFloat64()*1e3
		i, di := make([]float64, 4), make([]float64, 4)
		for k := range i {
			i[k], di[k] = r.NormFloat64(), r.NormFloat64()*1e3
		}
		e := Geodesic(s, J, th, i, di, om, dom)
		l, d1 := square(4), square(4)
		s.Inductance(th, l, d1, nil)
		for k := 0; k < 4; k++ {
			want := 0.0
			for b := 0; b < 4; b++ {
				want += l[k][b]*di[b] + om*d1[k][b]*i[b]
			}
			if math.Abs(e[k]-want) > 1e-12*(1+math.Abs(want)) {
				t.Fatalf("winding %d: Kron %g, element %g", k+1, e[k], want)
			}
		}
		T := 0.0
		for a := 0; a < 4; a++ {
			for b := 0; b < 4; b++ {
				T += 0.5 * i[a] * d1[a][b] * i[b]
			}
		}
		if want := J*dom - T; math.Abs(e[4]-want) > 1e-12*(1+math.Abs(want)) {
			t.Fatalf("shaft: Kron %g, element %g", e[4], want)
		}
	}
}

// The rotor frame makes the primitive machine's inductances constant
// (L′ = [[Ls·I, M·I], [M·I, Lr·I]]), puts speed voltages in the rotor
// rows only, and keeps the torque: i′ᵀ·G′·i′ = ½·iᵀ·L′(θ)·i with i = C·i′.
func TestRotatingFrame(t *testing.T) {
	const ls, lr, m = 10e-3, 12e-3, 8e-3
	for _, p := range []int{1, 2} {
		s := primitiveMachine(ls, lr, m, p)
		f := RotorFrame(4, [][2]int{{2, 3}}, float64(p))
		d, err := RotatingFrame(s, f)
		if err != nil {
			t.Fatal(err)
		}
		want := [][]float64{{ls, 0, m, 0}, {0, ls, 0, m}, {m, 0, lr, 0}, {0, m, 0, lr}}
		for i := range want {
			for j := range want[i] {
				if math.Abs(d.L0[i][j]-want[i][j]) > 1e-15 {
					t.Errorf("p=%d L′[%d][%d] = %g, want %g", p, i, j, d.L0[i][j], want[i][j])
				}
			}
			if i < 2 {
				for j := 0; j < 4; j++ {
					if d.G[i][j] != 0 {
						t.Errorf("p=%d: stator row %d has speed term %g", p, i+1, d.G[i][j])
					}
				}
			}
		}
		r := rand.New(rand.NewPCG(5, uint64(p)))
		for trial := 0; trial < 20; trial++ {
			th := r.Float64() * 7
			ip := []float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
			c := f.C(th)
			i := make([]float64, 4)
			for a := 0; a < 4; a++ {
				for b := 0; b < 4; b++ {
					i[a] += c[a][b] * ip[b]
				}
			}
			d1 := square(4)
			s.Inductance(th, nil, d1, nil)
			var t1, t2 float64
			for a := 0; a < 4; a++ {
				for b := 0; b < 4; b++ {
					t1 += 0.5 * i[a] * d1[a][b] * i[b]
					t2 += ip[a] * d.G[a][b] * ip[b]
				}
			}
			if math.Abs(t1-t2) > 1e-12*(1+math.Abs(t1)) {
				t.Fatalf("p=%d torque %g in the original frame, %g in the rotor frame", p, t1, t2)
			}
		}
	}
	// A frame that does not remove θ is refused.
	s := primitiveMachine(ls, lr, m, 1)
	if _, err := RotatingFrame(s, RotorFrame(4, nil, 1)); err == nil || !strings.Contains(err.Error(), "still depend on θ") {
		t.Errorf("identity frame: %v", err)
	}
}

// inductionMotor is a two-phase induction motor starting from rest:
// stator fed with a balanced two-phase supply through Rs, rotor windings
// shorted through Rr, a rotor of inertia J with a friction load.
func inductionMotor(t *testing.T, s *circuit.MachineSpec) *circuit.Circuit {
	t.Helper()
	const (
		vm, f  = 100.0, 50.0
		rs, rr = 0.5, 0.8
		J, fr  = 0.01, 50.0 // inertia [kg·m²], mobility of the friction [rad/s per N·m]
	)
	w := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vm, f, 0, 0, 90}}
	c, err := circuit.NewBuilder("two-phase induction motor").Domain(exchange.Mechanical, "w").
		MachineModel("mm", s).
		VoltageWave("VA", "sa", "0", w).
		VoltageWave("VB", "sb", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vm, f, 0, 0, 0}}).
		Resistor("RSA", "sa", "ia", rs).Resistor("RSB", "sb", "ib", rs).
		Resistor("RRA", "ra", "0", rr).Resistor("RRB", "rb", "0", rr).
		Machine("N1", "mm", [][2]string{{"ia", "0"}, {"ib", "0"}, {"ra", "0"}, {"rb", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", J, 0).Resistor("RF", "w", "0", fr).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// Kron's transformation twin: the same induction motor solved in the
// holonomic frame I (rotor windings fixed to the rotor, L(θ) trigonometric)
// and in the rotor frame II(b) produced by RotatingFrame (constant L,
// speed voltages) gives the same speed and stator currents through the
// whole start-up, and both ledgers close.
func TestMachineFrameTwin(t *testing.T) {
	const ls, lr, m = 0.12, 0.12, 0.11
	s1 := primitiveMachine(ls, lr, m, 1)
	s2, err := RotatingFrame(s1, RotorFrame(4, [][2]int{{2, 3}}, 1))
	if err != nil {
		t.Fatal(err)
	}
	spec := analysis.TranSpec{Step: 2e-4, Stop: 0.6, UIC: true}
	var res [2]struct {
		tm, w, ia []float64
		mm        float64
	}
	for k, s := range []*circuit.MachineSpec{s1, s2} {
		c := inductionMotor(t, s)
		wf, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}})
		if err != nil {
			t.Fatalf("frame %d: %v", k+1, err)
		}
		res[k].tm = wf.Time
		res[k].w, _ = wf.NodeVoltage("w")
		res[k].ia, _ = wf.ElementCurrent("RSA")
		l := wf.Ledger
		res[k].mm = math.Abs(l.IntegrationMismatch) / (math.Abs(l.EnergyIn) + math.Abs(l.StoredEnd))
	}
	interp := func(tm, v []float64, t0 float64) float64 {
		j := 1
		for j < len(tm)-1 && tm[j] < t0 {
			j++
		}
		a := (t0 - tm[j-1]) / (tm[j] - tm[j-1])
		return v[j-1] + a*(v[j]-v[j-1])
	}
	wmax, imax, dw, di := 0.0, 0.0, 0.0, 0.0
	for i, t0 := range res[0].tm {
		w2, i2 := interp(res[1].tm, res[1].w, t0), interp(res[1].tm, res[1].ia, t0)
		wmax, imax = math.Max(wmax, math.Abs(res[0].w[i])), math.Max(imax, math.Abs(res[0].ia[i]))
		dw, di = math.Max(dw, math.Abs(res[0].w[i]-w2)), math.Max(di, math.Abs(res[0].ia[i]-i2))
	}
	if dw > 2e-3*wmax || di > 1e-2*imax {
		t.Errorf("frames differ: speed by %.3g of %.4g rad/s, stator current by %.3g of %.4g A", dw/wmax, wmax, di/imax, imax)
	}
	for k := range res {
		if res[k].mm > 2e-3 {
			t.Errorf("frame %d: energy mismatch %.3g", k+1, res[k].mm)
		}
	}
	t.Logf("final speed %.4g rad/s (frame I) %.4g rad/s (frame II); worst Δω %.2g, worst Δi %.2g of peak; %d and %d points; energy mismatch %.2g, %.2g",
		res[0].w[len(res[0].w)-1], res[1].w[len(res[1].w)-1], dw/wmax, di/imax, len(res[0].tm), len(res[1].tm), res[0].mm, res[1].mm)
}
