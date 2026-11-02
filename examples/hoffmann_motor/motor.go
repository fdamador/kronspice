// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/kron"
)

// Params is Hoffmann's single-phase squirrel-cage motor (Rev. Mod. Phys.
// 21, 535, 1949, §3, Figs. 1–3): stator coil 1, the cage represented by two
// perpendicular rotor coils 2, 3, rotor angle θ (two-pole: mechanical and
// electrical angle are the same). Inductances follow his eq. 42: mutuals
// M·cos θ, −M·sin θ; rotor self and mutual L_S ∓ L_D·cos 2θ, L_D·sin 2θ,
// where L_D ≠ 0 means a salient pole (an uneven air gap). All values are
// declared: a fractional-horsepower motor on 120 V, 60 Hz.
type Params struct {
	L11, M, LS, LD float64 // [H]
	Rs, Rr         float64 // stator, rotor-coil resistance [Ω]
	J, Friction    float64 // inertia [kg·m²], friction [N·m per rad/s]
	Vpeak, F       float64 // supply [V], [Hz]
}

// Default motor (declared).
func Default() Params {
	return Params{L11: 0.25, M: 0.23, LS: 0.25, Rs: 2, Rr: 3, J: 1e-3, Friction: 1e-4, Vpeak: 120 * math.Sqrt2, F: 60}
}

// Sync is the synchronous speed [rad/s].
func (p Params) Sync() float64 { return 2 * math.Pi * p.F }

// FrameI is the machine with the rotor coils turning with the rotor
// (Hoffmann eqs. 24, 42: the holonomic, Riemannian frame).
func (p Params) FrameI() *circuit.MachineSpec {
	s := circuit.NewMachineSpec(3)
	s.L0[0][0], s.L0[1][1], s.L0[2][2] = p.L11, p.LS, p.LS
	h1 := s.Harmonic(1)
	h1.Cos[0][1], h1.Cos[1][0] = p.M, p.M
	h1.Sin[0][2], h1.Sin[2][0] = -p.M, -p.M
	if p.LD != 0 {
		h2 := s.Harmonic(2)
		h2.Cos[1][1], h2.Cos[2][2] = -p.LD, p.LD
		h2.Sin[1][2], h2.Sin[2][1] = p.LD, p.LD
	}
	return s
}

// Frame is Hoffmann's change to stationary rotor coils 2̄, 3̄ (eq. 35).
func Frame() kron.MovingFrame { return kron.RotorFrame(3, [][2]int{{1, 2}}, 1) }

// FrameIIb is the machine on stationary rotor coils: constant inductances
// (eq. 29) and the speed matrix G that carries Hoffmann's torsion (kron.RotatingFrame).
func (p Params) FrameIIb() (*circuit.MachineSpec, error) {
	return kron.RotatingFrame(p.FrameI(), Frame())
}

// Naive is Hoffmann's warning case: the constant inductances of eq. 29 used
// in the ordinary Lagrange–Maxwell equations, with no speed matrix (no
// torsion). It is not the motor.
func (p Params) Naive() (*circuit.MachineSpec, error) {
	s, err := p.FrameIIb()
	if err != nil {
		return nil, err
	}
	n := circuit.NewMachineSpec(3)
	for i := range s.L0 {
		copy(n.L0[i], s.L0[i])
	}
	return n, nil
}

// build makes the motor circuit. With fixedSpeed true the shaft is held at
// omega by an ideal speed source (VW, whose current is the motor torque);
// otherwise it starts at omega with inertia J and friction.
func (p Params) build(s *circuit.MachineSpec, omega float64, fixedSpeed bool) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("Hoffmann single-phase squirrel-cage motor").Domain(exchange.Mechanical, "w").
		MachineModel("mm", s).
		VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, p.Vpeak, p.F}}).
		Resistor("RS", "s", "i1", p.Rs).
		Resistor("R2", "r2", "0", p.Rr).Resistor("R3", "r3", "0", p.Rr).
		Machine("N1", "mm", [][2]string{{"i1", "0"}, {"r2", "0"}, {"r3", "0"}}, "w", "0")
	if fixedSpeed {
		b.VoltageSource("VW", "w", "0", omega)
	} else {
		b.CapacitorIC("CJ", "w", "0", p.J, omega).Resistor("RF", "w", "0", 1/p.Friction)
	}
	return b.BuildTransient()
}

// Run is a start-up run.
type Run struct {
	Name                   string
	Time, Speed            []float64
	Final                  float64 // speed at the end [rad/s]
	Supplied, Heat, Stored float64 // [J]
	Mismatch               float64 // ledger mismatch ÷ energy scale
}

// Start runs the motor for seconds from initial speed omega0.
func (p Params) Start(name string, s *circuit.MachineSpec, omega0, seconds float64) (Run, error) {
	c, err := p.build(s, omega0, false)
	if err != nil {
		return Run{}, err
	}
	step := 1 / (p.F * 200)
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: step, Stop: seconds, MaxStep: step, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}, MaxSteps: 5_000_000})
	if err != nil {
		return Run{}, fmt.Errorf("%s: %w", name, err)
	}
	w, _ := wf.NodeVoltage("w")
	l := wf.Ledger
	r := Run{Name: name, Time: wf.Time, Speed: w, Final: w[len(w)-1],
		Supplied: l.EnergyIn - l.EnergyOut, Heat: l.EnergyDissipated, Stored: l.StoredEnd - l.StoredStart}
	r.Mismatch = l.IntegrationMismatch / (math.Abs(l.EnergyIn) + math.Abs(l.StoredEnd) + math.Abs(l.StoredStart))
	return r, nil
}

// Torque is the mean motor torque [N·m] with the shaft held at omega,
// averaged over the last 10 supply cycles of a 1 s run.
func (p Params) Torque(s *circuit.MachineSpec, omega float64) (float64, error) {
	c, err := p.build(s, omega, true)
	if err != nil {
		return 0, err
	}
	step := 1 / (p.F * 200)
	stop := 1.0
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}, MaxSteps: 5_000_000})
	if err != nil {
		return 0, err
	}
	iw, ok := wf.ElementCurrent("VW")
	if !ok {
		return 0, fmt.Errorf("no current for VW")
	}
	t0 := stop - 10/p.F
	var T, sum float64
	for k := 1; k < len(wf.Time); k++ {
		if wf.Time[k] <= t0 {
			continue
		}
		dt := wf.Time[k] - wf.Time[k-1]
		T += dt
		sum += iw[k] * dt // current out of the speed source's + node = torque the motor applies
	}
	return sum / T, nil
}
