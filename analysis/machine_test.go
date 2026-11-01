// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/exchange"
)

// rotatingMutual is a two-winding machine whose mutual inductance turns
// with the rotor: L = [[L1, M·cos θ], [M·cos θ, L2]].
func rotatingMutual(l1, l2, m, theta0 float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0], s.L0[1][1] = l1, l2
	h := s.Harmonic(1)
	h.Cos[0][1], h.Cos[1][0] = m, m
	s.Theta0 = theta0
	return s
}

// withIC sets a spec's initial winding currents (UIC).
func withIC(s *circuit.MachineSpec, ic ...float64) *circuit.MachineSpec {
	copy(s.IC, ic)
	return s
}

// At standstill with DC winding currents the torque is
// T = ½·iᵀ·L′·i = −M·sin θ0·I1·I2; through a damper of mobility Rd the
// shaft settles at ω = Rd·T, and each winding shows the speed voltage
// v1 = ω·(L′·i)_1 = −ω·M·sin θ0·I2.
func TestMachineTorqueOP(t *testing.T) {
	const (
		l1, l2, m = 10e-3, 20e-3, 5e-3
		i1, i2    = 2.0, 3.0
		th0       = 0.7
		rd        = 4.0 // rad/s per N·m
	)
	c, err := circuit.NewBuilder("torque").Domain(exchange.Mechanical, "w").
		MachineModel("mm", rotatingMutual(l1, l2, m, th0)).
		CurrentSource("I1", "0", "a", i1).CurrentSource("I2", "0", "b", i2).
		Resistor("RA", "a", "0", 1e6).Resistor("RB", "b", "0", 1e6).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}}, "w", "0").
		Resistor("RD", "w", "0", rd).Build()
	if err != nil {
		t.Fatal(err)
	}
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	T := -m * math.Sin(th0) * i1 * i2
	w := r.NodeVoltages[c.NodeIndex("w")]
	if math.Abs(w-rd*T) > 1e-9*math.Abs(rd*T) {
		t.Errorf("ω = %g, want %g", w, rd*T)
	}
	va := r.NodeVoltages[c.NodeIndex("a")]
	want := -w * m * math.Sin(th0) * i2
	if math.Abs(va-want) > 1e-6*math.Abs(want) {
		t.Errorf("speed voltage %g, want %g", va, want)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-9 {
		t.Errorf("ledger residual %g", rr)
	}
	_ = device.Machine{}
}

// A rotor of inertia J between two DC-fed windings is a pendulum:
// J·θ̈ = −M·I1·I2·sin θ. The transient ω(t) matches the pendulum equation
// integrated independently (RK4), and the ledger closes: the sources'
// work equals the rotor's kinetic energy plus the change of ½·iᵀ·L(θ)·i.
func TestMachinePendulum(t *testing.T) {
	const (
		l1, l2, m = 10e-3, 20e-3, 5e-3
		i1, i2    = 2.0, 3.0
		th0       = 1.0
		J         = 1e-4
	)
	c, err := circuit.NewBuilder("pendulum").Domain(exchange.Mechanical, "w").
		MachineModel("mm", withIC(rotatingMutual(l1, l2, m, th0), i1, i2)).
		CurrentSource("I1", "0", "a", i1).CurrentSource("I2", "0", "b", i2).
		// 1 GΩ across each winding fixes its node voltage at the UIC
		// start, where the winding current is imposed (as for an inductor).
		Resistor("RA", "a", "0", 1e9).Resistor("RB", "b", "0", 1e9).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", J, 0).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	k := m * i1 * i2 / J
	period := 2 * math.Pi / math.Sqrt(k) // small-angle period
	stop := 2 * period
	wf, err := Tran(c, TranSpec{Step: period / 200, Stop: stop, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	w, _ := wf.NodeVoltage("w")
	// Reference: RK4 on θ̈ = −k·sin θ with a fine fixed step.
	ref := func(tt float64) float64 {
		n := int(math.Ceil(tt / (period / 20000)))
		if n == 0 {
			return 0
		}
		h := tt / float64(n)
		th, om := th0, 0.0
		f := func(th float64) float64 { return -k * math.Sin(th) }
		for i := 0; i < n; i++ {
			k1t, k1w := om, f(th)
			k2t, k2w := om+h/2*k1w, f(th+h/2*k1t)
			k3t, k3w := om+h/2*k2w, f(th+h/2*k2t)
			k4t, k4w := om+h*k3w, f(th+h*k3t)
			th += h / 6 * (k1t + 2*k2t + 2*k3t + k4t)
			om += h / 6 * (k1w + 2*k2w + 2*k3w + k4w)
		}
		return om
	}
	wmax, worst := 0.0, 0.0
	for i, tt := range wf.Time {
		r := ref(tt)
		wmax = math.Max(wmax, math.Abs(r))
		worst = math.Max(worst, math.Abs(w[i]-r))
	}
	if worst > 2e-3*wmax {
		t.Errorf("ω differs from the pendulum by %.3g of its peak %.3g rad/s", worst/wmax, wmax)
	}
	l := wf.Ledger
	if l.MaxRelativeResidual > 1e-9 {
		t.Errorf("ledger residual %g", l.MaxRelativeResidual)
	}
	// Energy over the run: the sources' work equals the change of stored
	// energy (½·iᵀ·L(θ)·i in the machine, ½·J·ω² in the inertia), to the
	// accuracy of integrating the sampled power.
	scale := math.Abs(l.EnergyIn) + math.Abs(l.StoredStart) + math.Abs(l.StoredEnd)
	if mm := math.Abs(l.IntegrationMismatch) / scale; mm > 1e-3 {
		t.Errorf("energy mismatch %.3g of %.3g J", mm, scale)
	}
	t.Logf("%d points, ω peak %.4g rad/s, worst error %.2g of peak; sources' net work %.4g J, stored %.4g → %.4g J, mismatch %.2g",
		len(wf.Time), wmax, worst/wmax, l.EnergyIn-l.EnergyOut, l.StoredStart, l.StoredEnd, l.IntegrationMismatch/scale)
}

// A temporal (parametric) inductance: one winding with
// L(θ) = L0·(1 + m·cos 2θ), its coordinate driven at constant speed ω0 by a
// mechanical velocity source, so L(t) = L0·(1 + m·cos 2ω0t) pumps an LC
// tank tuned to ω0. The flux obeys Mathieu's equation
// φ̈ + ω0²·(1 − m·cos 2ω0t)·φ = 0 (to first order), whose principal
// instability grows as e^(m·ω0/4·t); a series resistance takes R/(2L0)
// off that rate. The run checks the growth and decay rates and that the
// energy gained by the tank comes from the pump, through the mechanical
// port.
func TestMachineParametric(t *testing.T) {
	const (
		l0 = 1e-3
		m  = 0.05
		f0 = 1e3
	)
	w0 := 2 * math.Pi * f0
	cc := 1 / (w0 * w0 * l0)
	for _, tc := range []struct {
		name string
		r    float64 // series resistance
	}{
		{"growing", m * w0 / 8 * 2 * l0},  // damping half the pump rate
		{"decaying", m * w0 / 2 * 2 * l0}, // damping twice the pump rate
	} {
		s := circuit.NewMachineSpec(1)
		s.L0[0][0] = l0
		s.Harmonic(2).Cos[0][0] = m * l0
		c, err := circuit.NewBuilder("parametric").Domain(exchange.Mechanical, "w").
			MachineModel("pump", s).
			Machine("N1", "pump", [][2]string{{"a", "b"}}, "w", "0").
			Resistor("R1", "b", "0", tc.r).
			CapacitorIC("C1", "a", "0", cc, 1).
			VoltageSource("VW", "w", "0", w0).BuildTransient()
		if err != nil {
			t.Fatal(err)
		}
		period := 1 / f0
		stop := 60 * period
		wf, err := Tran(c, TranSpec{Step: period / 100, Stop: stop, UIC: true}, TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		v, _ := wf.NodeVoltage("a")
		// Peak |v| per period over the last 30 periods, then the slope of
		// its logarithm.
		var tp, lp []float64
		for p := 30; p < 60; p++ {
			pk, at := 0.0, 0.0
			for i, tt := range wf.Time {
				if tt >= float64(p)*period && tt < float64(p+1)*period && math.Abs(v[i]) > pk {
					pk, at = math.Abs(v[i]), tt
				}
			}
			tp, lp = append(tp, at), append(lp, math.Log(pk))
		}
		n := float64(len(tp))
		var st, sl, stt, stl float64
		for i := range tp {
			st, sl, stt, stl = st+tp[i], sl+lp[i], stt+tp[i]*tp[i], stl+tp[i]*lp[i]
		}
		rate := (n*stl - st*sl) / (n*stt - st*st)
		want := m*w0/4 - tc.r/(2*l0)
		if math.Abs(rate-want) > 0.05*math.Abs(want) {
			t.Errorf("%s: growth rate %.4g /s, want %.4g /s", tc.name, rate, want)
		}
		// Energy: the pump's mechanical work is what the tank gained plus
		// what R dissipated.
		l := wf.Ledger
		scale := math.Abs(l.EnergyIn) + math.Abs(l.EnergyOut) + math.Abs(l.StoredStart) + math.Abs(l.StoredEnd)
		if mm := math.Abs(l.IntegrationMismatch) / scale; mm > 2e-3 {
			t.Errorf("%s: energy mismatch %.3g", tc.name, mm)
		}
		t.Logf("%s: rate %.4g /s (Mathieu %.4g); pump net work %.4g J, R heat %.4g J, tank %.4g → %.4g J, mismatch %.2g",
			tc.name, rate, want, l.EnergyIn-l.EnergyOut, l.EnergyDissipated, l.StoredStart, l.StoredEnd, l.IntegrationMismatch/scale)
	}
}

// BenchmarkMachinePendulum is the pendulum transient (two windings, the
// rotor, about 400 steps).
func BenchmarkMachinePendulum(b *testing.B) {
	s := withIC(rotatingMutual(10e-3, 20e-3, 5e-3, 1), 2, 3)
	c, err := circuit.NewBuilder("pendulum").Domain(exchange.Mechanical, "w").MachineModel("mm", s).
		CurrentSource("I1", "0", "a", 2).CurrentSource("I2", "0", "b", 3).
		Resistor("RA", "a", "0", 1e9).Resistor("RB", "b", "0", 1e9).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", 1e-4, 0).BuildTransient()
	if err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := Tran(c, TranSpec{Step: 1e-4, Stop: 0.05, UIC: true}, TranOptions{}); err != nil {
			b.Fatal(err)
		}
	}
}
