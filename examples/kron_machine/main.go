// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kron_machine shows Kron's dynamical view of moving inductors
// (G. Kron, Equivalent Circuits of Electric Machinery, Wiley 1951,
// Epilogue): one two-phase induction motor solved in the holonomic frame I
// (rotor windings on the rotor, inductances trigonometric in θ) and in the
// non-holonomic rotor frame II(b) produced by kron.RotatingFrame (constant
// inductances, speed voltages G·i·ω), starting from rest; then a
// parametric (time-varying) inductance pumping an LC tank, against
// Mathieu's growth rate.
//
//	go run ./examples/kron_machine
package main

import (
	"fmt"
	"log"
	"math"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/kron"
)

// Primitive is Kron's two-phase primitive machine: stator α, β (1, 2),
// rotor a, b (3, 4) at angle θ, stator–rotor mutuals M·R(θ).
func Primitive(ls, lr, m float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(4)
	s.L0[0][0], s.L0[1][1], s.L0[2][2], s.L0[3][3] = ls, ls, lr, lr
	h := s.Harmonic(1)
	set := func(mat [][]float64, i, j int, v float64) { mat[i][j], mat[j][i] = v, v }
	set(h.Cos, 0, 2, m)
	set(h.Sin, 0, 3, -m)
	set(h.Sin, 1, 2, m)
	set(h.Cos, 1, 3, m)
	return s
}

// Motor is the motor circuit: a balanced two-phase 100 V, 50 Hz supply
// through 0.5 Ω, rotor windings shorted through 0.8 Ω, rotor inertia
// 0.002 kg·m² and a friction load of 0.2 mN·m per rad/s.
func Motor(s *circuit.MachineSpec) (*circuit.Circuit, error) {
	cosw := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 100, 50, 0, 0, 90}}
	sinw := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 100, 50, 0, 0, 0}}
	return circuit.NewBuilder("two-phase induction motor").Domain(exchange.Mechanical, "w").
		MachineModel("mm", s).
		VoltageWave("VA", "sa", "0", cosw).VoltageWave("VB", "sb", "0", sinw).
		Resistor("RSA", "sa", "ia", 0.5).Resistor("RSB", "sb", "ib", 0.5).
		Resistor("RRA", "ra", "0", 0.8).Resistor("RRB", "rb", "0", 0.8).
		Machine("N1", "mm", [][2]string{{"ia", "0"}, {"ib", "0"}, {"ra", "0"}, {"rb", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", 0.002, 0).Resistor("RF", "w", "0", 5000).BuildTransient()
}

func main() {
	s1 := Primitive(0.12, 0.12, 0.11)
	s2, err := kron.RotatingFrame(s1, kron.RotorFrame(4, [][2]int{{2, 3}}, 1))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("Kron's rotor frame (constant inductances, speed matrix G):")
	for i := range s2.L0 {
		fmt.Printf("  L′ %v   G′ %v\n", fmtRow(s2.L0[i]), fmtRow(s2.G[i]))
	}
	fmt.Println("\nStart-up from rest, 1.5 s (synchronous speed 314.16 rad/s):")
	spec := analysis.TranSpec{Step: 5e-4, Stop: 1.5, UIC: true}
	for _, v := range []struct {
		name string
		s    *circuit.MachineSpec
	}{{"frame I  (holonomic, L(θ))", s1}, {"frame II (rotor frame, G·i·ω)", s2}} {
		c, err := Motor(v.s)
		if err != nil {
			log.Fatal(err)
		}
		t0 := time.Now()
		wf, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}})
		if err != nil {
			log.Fatal(err)
		}
		w, _ := wf.NodeVoltage("w")
		l := wf.Ledger
		fmt.Printf("  %-30s final speed %8.3f rad/s, %5d points, %5.2f s; supply %.2f J, heat %.2f J, stored %.2f J, mismatch %.1e\n",
			v.name, w[len(w)-1], len(wf.Time), time.Since(t0).Seconds(), l.EnergyIn-l.EnergyOut, l.EnergyDissipated, l.StoredEnd-l.StoredStart,
			l.IntegrationMismatch/(math.Abs(l.EnergyIn)+math.Abs(l.StoredEnd)))
	}

	fmt.Println("\nParametric inductance L(t) = L0·(1 + m·cos 2ω0t) pumping an LC tank at ω0 (1 kHz, m = 0.05):")
	for _, damp := range []float64{0, 0.5, 2} {
		rate, pump, heat := parametric(1e-3, 0.05, 1e3, damp)
		w0 := 2 * math.Pi * 1e3
		want := 0.05 * w0 / 4 * (1 - damp)
		fmt.Printf("  damping %.1f × pump rate: growth %8.3f /s (Mathieu %8.3f /s); pump work %.4g J, heat %.4g J\n", damp, rate, want, pump, heat)
	}
}

// parametric runs the pumped tank and returns the envelope growth rate,
// the pump's net work and the resistor's heat. damp is R/(2L0) as a
// multiple of the pump rate m·ω0/4.
func parametric(l0, m, f0, damp float64) (rate, pump, heat float64) {
	w0 := 2 * math.Pi * f0
	s := circuit.NewMachineSpec(1)
	s.L0[0][0] = l0
	s.Harmonic(2).Cos[0][0] = m * l0
	r := damp * m * w0 / 4 * 2 * l0
	b := circuit.NewBuilder("parametric").Domain(exchange.Mechanical, "w").MachineModel("pump", s).
		Machine("N1", "pump", [][2]string{{"a", "b"}}, "w", "0").
		CapacitorIC("C1", "a", "0", 1/(w0*w0*l0), 1).VoltageSource("VW", "w", "0", w0)
	b.Resistor("R1", "b", "0", math.Max(r, 1e-9))
	c, err := b.BuildTransient()
	if err != nil {
		log.Fatal(err)
	}
	period := 1 / f0
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: period / 100, Stop: 60 * period, UIC: true}, analysis.TranOptions{})
	if err != nil {
		log.Fatal(err)
	}
	v, _ := wf.NodeVoltage("a")
	var ts, ls []float64
	for p := 30; p < 60; p++ {
		pk, at := 0.0, 0.0
		for i, tt := range wf.Time {
			if tt >= float64(p)*period && tt < float64(p+1)*period && math.Abs(v[i]) > pk {
				pk, at = math.Abs(v[i]), tt
			}
		}
		ts, ls = append(ts, at), append(ls, math.Log(pk))
	}
	n := float64(len(ts))
	var st, sl, stt, stl float64
	for i := range ts {
		st, sl, stt, stl = st+ts[i], sl+ls[i], stt+ts[i]*ts[i], stl+ts[i]*ls[i]
	}
	return (n*stl - st*sl) / (n*stt - st*st), wf.Ledger.EnergyIn - wf.Ledger.EnergyOut, wf.Ledger.EnergyDissipated
}

func fmtRow(r []float64) string {
	out := "["
	for i, v := range r {
		if i > 0 {
			out += " "
		}
		out += fmt.Sprintf("%7.3g", v)
	}
	return out + "]"
}
