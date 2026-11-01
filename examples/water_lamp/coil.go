// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Coil is the Tesla bifilar solenoid of examples/pulsed_bifilar reduced to
// its first mode: inductance Ls with the capacitance Cs that resonates it
// at f1 and the parallel resistance Rs that gives its Q, coupled to the
// 5-turn primary (Lp) by the primary mutual inductance M. Values from
// examples/pulsed_bifilar/run_output.txt (L = 4.853 mH, f1 = 41.61 kHz,
// Q ≈ 83) and its primary-to-coil mutual (28.1 µH); Lp is the 5-turn
// primary of that example (radius 6 cm, pitch 2 mm, 1 mm wire), computed
// the same way. A one-mode model: the coil's higher modes are left out.
type Coil struct {
	Ls, Cs, Rs, Lp, M, F1, Q float64
}

// DefaultCoil builds the reduced bifilar solenoid.
func DefaultCoil() Coil {
	const (
		turns  = 5
		radius = 0.06
		pitch  = 2e-3
		wire   = 1e-3
	)
	lp := 0.0
	for i := 0; i < turns; i++ {
		for j := 0; j < turns; j++ {
			if i == j {
				lp += em.Mu0 * radius * (math.Log(8*radius/wire) - 2)
			} else {
				lp += em.MutualCircles(radius, radius, float64(i-j)*pitch)
			}
		}
	}
	c := Coil{Ls: 4.853e-3, F1: 41.61e3, Q: 83, Lp: lp, M: 28.1e-6}
	w := 2 * math.Pi * c.F1
	c.Cs = 1 / (w * w * c.Ls)
	c.Rs = c.Q * w * c.Ls
	return c
}

// CoilRun is a Part 2 result.
type CoilRun struct {
	Run
	CoilPeak    float64 // largest coil voltage [V]
	CoilEnergy  float64 // heat in the coil + energy left in it [J]
	FlashEnergy float64 // energy the capacitor gave up in flashes [J]
	PeakIp      float64 // largest spark (primary) current [A]
}

// RunCoil drives the coil from this lamp circuit.
func (c Case) RunCoil(coil Coil, seconds float64, seed uint64) (CoilRun, error) {
	cc := coil
	return simulate(c, &cc, seconds, seed)
}

// Twin drives the reduced coil's primary with a 1 A, 2 µs current pulse in
// kronspice and in the model's own RK4 equations, and returns the largest
// difference in coil voltage over 10 periods, relative to its peak.
func (c Coil) Twin() (float64, error) {
	pulse := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 0.2e-6, 0.2e-6, 2e-6, 1}}
	b := circuit.NewBuilder("bifilar coil, first mode")
	b.CurrentWave("IP", "0", "p", pulse)
	b.Inductor("LP", "p", "0", c.Lp)
	b.Inductor("LS", "s", "0", c.Ls)
	b.Capacitor("CS", "s", "0", c.Cs)
	b.Resistor("RS", "s", "0", c.Rs)
	b.Mutual("K1", "LP", "LS", c.M/math.Sqrt(c.Lp*c.Ls))
	ckt, err := b.BuildTransient()
	if err != nil {
		return 0, err
	}
	stop := 10 / c.F1
	step := 1 / c.F1 / 2000
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return 0, err
	}
	vs, _ := wv.NodeVoltage("s")
	// own equations with the primary current imposed: Ls dIs/dt = Vc − M dIp/dt
	ip := func(t float64) (i, di float64) {
		switch {
		case t < 0.2e-6:
			return t / 0.2e-6, 1 / 0.2e-6
		case t < 2.2e-6:
			return 1, 0
		case t < 2.4e-6:
			return (2.4e-6 - t) / 0.2e-6, -1 / 0.2e-6
		}
		return 0, 0
	}
	// Ls·dIs/dt + M·dIp/dt = V, Cs·dV/dt = −Is − V/Rs (the same sign convention as simulate)
	f := func(t, is, v float64) (float64, float64) {
		_, di := ip(t)
		return (v - c.M*di) / c.Ls, (-is - v/c.Rs) / c.Cs
	}
	dt := step / 20
	is, v, t := 0.0, 0.0, 0.0
	peak, worst := 0.0, 0.0
	k := 1
	for k < len(wv.Time) {
		a1, b1 := f(t, is, v)
		a2, b2 := f(t+dt/2, is+dt/2*a1, v+dt/2*b1)
		a3, b3 := f(t+dt/2, is+dt/2*a2, v+dt/2*b2)
		a4, b4 := f(t+dt, is+dt*a3, v+dt*b3)
		is += dt / 6 * (a1 + 2*a2 + 2*a3 + a4)
		v += dt / 6 * (b1 + 2*b2 + 2*b3 + b4)
		t += dt
		for k < len(wv.Time) && wv.Time[k] <= t+dt/2 {
			if math.Abs(wv.Time[k]-t) < dt {
				peak = math.Max(peak, math.Abs(v))
				worst = math.Max(worst, math.Abs(vs[k]-v))
			}
			k++
		}
	}
	return worst / peak, nil
}
