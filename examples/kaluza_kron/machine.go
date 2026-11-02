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

// Kaluza's construction inside a Kron machine. Kron's coordinates are the
// winding charges q1, q2 and the rotor angle θ; the inductances depend on θ
// only, never on the charges — Kaluza's cylinder condition, as Hoffmann
// (Rev. Mod. Phys. 21, 535, 1949, p. 536) calls it. With the windings
// shorted and lossless (superconducting) the charges are cyclic
// coordinates, so their momenta — the flux linkages ψ = L(θ)·i — are
// conserved, as Kaluza's p_w is. Eliminating them (Routh, as Kaluza–Klein
// reduction eliminates w) leaves the rotor alone in an effective potential
// set by the conserved fluxes:
//
//	J·θ̈ = −∂W/∂θ |ψ,   W(θ) = ½·ψᵀ·L(θ)⁻¹·ψ
//
// the machine's counterpart of the Lorentz force the fifth dimension leaves
// behind. Three ways: Kron's geodesic in (q1, q2, θ) (kron.Geodesic), the
// reduced rotor equation, and kronspice's machine element with the
// windings shorted through a tiny resistance.

// Machine is a stator coil and a rotor coil, mutual M·cos θ (declared).
type Machine struct {
	L1, L2, M, J float64
	I1, I2       float64 // initial currents [A]
	Theta0, W0   float64 // initial angle [rad], speed [rad/s]
}

// DefaultMachine (declared values).
func DefaultMachine() Machine {
	return Machine{L1: 1, L2: 1, M: 0.5, J: 0.01, I1: 1, I2: 0.5, Theta0: 0.3}
}

// Spec is the kronspice machine model (frame I).
func (m Machine) Spec() *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0], s.L0[1][1] = m.L1, m.L2
	h := s.Harmonic(1)
	h.Cos[0][1], h.Cos[1][0] = m.M, m.M
	s.Theta0 = m.Theta0
	s.IC = []float64{m.I1, m.I2}
	return s
}

func (m Machine) lmat(th float64) (l11, l12, l22 float64) {
	return m.L1, m.M * math.Cos(th), m.L2
}

// Fluxes are the conserved momenta ψ = L(θ)·i at the start.
func (m Machine) Fluxes() (p1, p2 float64) {
	a, b, d := m.lmat(m.Theta0)
	return a*m.I1 + b*m.I2, b*m.I1 + d*m.I2
}

// MTrace is a machine run.
type MTrace struct {
	Time, Theta, Omega []float64
	FluxDrift          float64 // largest relative change of ψ1, ψ2 (geodesic run)
	EnergyDrift        float64
}

// Geodesic integrates Kron's eq. 3 with no resistance and no applied
// voltage in all three coordinates (RK4).
func (m Machine) Geodesic(stop, dt float64) MTrace {
	s := m.Spec()
	p1, p2 := m.Fluxes()
	y := []float64{m.I1, m.I2, m.W0, m.Theta0}
	f := func(y []float64) []float64 {
		// Γ·v·v from Geodesic with zero accelerations, then g·ẍ = −Γvv
		e := kron.Geodesic(s, m.J, y[3], []float64{y[0], y[1]}, []float64{0, 0}, y[2], 0)
		a, b, d := m.lmat(y[3])
		det := a*d - b*b
		di1 := (-d*e[0] + b*e[1]) / det
		di2 := (b*e[0] - a*e[1]) / det
		return []float64{di1, di2, -e[2] / m.J, y[2]}
	}
	var tr MTrace
	e0 := m.energy(y)
	steps := int(math.Round(stop / dt))
	for k := 0; k <= steps; k++ {
		if k%10 == 0 {
			tr.Time, tr.Theta, tr.Omega = append(tr.Time, float64(k)*dt), append(tr.Theta, y[3]), append(tr.Omega, y[2])
		}
		a, b, d := m.lmat(y[3])
		q1, q2 := a*y[0]+b*y[1], b*y[0]+d*y[1]
		tr.FluxDrift = math.Max(tr.FluxDrift, math.Max(math.Abs(q1-p1)/math.Abs(p1), math.Abs(q2-p2)/math.Max(math.Abs(p1), 1e-300)))
		tr.EnergyDrift = math.Max(tr.EnergyDrift, math.Abs(m.energy(y)-e0)/e0)
		rk4(f, y, dt)
	}
	return tr
}

func (m Machine) energy(y []float64) float64 {
	a, b, d := m.lmat(y[3])
	return 0.5*(a*y[0]*y[0]+2*b*y[0]*y[1]+d*y[1]*y[1]) + 0.5*m.J*y[2]*y[2]
}

// Reduced integrates the rotor alone with the fluxes held at their initial
// values: J·θ̈ = −∂W/∂θ|ψ (Routh).
func (m Machine) Reduced(stop, dt float64) MTrace {
	p1, p2 := m.Fluxes()
	W := func(th float64) float64 {
		a, b, d := m.lmat(th)
		det := a*d - b*b
		return 0.5 * (d*p1*p1 - 2*b*p1*p2 + a*p2*p2) / det
	}
	f := func(y []float64) []float64 {
		h := 1e-6
		return []float64{-(W(y[1]+h) - W(y[1]-h)) / (2 * h) / m.J, y[0]}
	}
	y := []float64{m.W0, m.Theta0}
	var tr MTrace
	steps := int(math.Round(stop / dt))
	for k := 0; k <= steps; k++ {
		if k%10 == 0 {
			tr.Time, tr.Theta, tr.Omega = append(tr.Time, float64(k)*dt), append(tr.Theta, y[1]), append(tr.Omega, y[0])
		}
		rk4(f, y, dt)
	}
	return tr
}

func rk4(f func([]float64) []float64, y []float64, h float64) {
	n := len(y)
	t := make([]float64, n)
	k1 := f(y)
	for i := range y {
		t[i] = y[i] + h/2*k1[i]
	}
	k2 := f(t)
	for i := range y {
		t[i] = y[i] + h/2*k2[i]
	}
	k3 := f(t)
	for i := range y {
		t[i] = y[i] + h*k3[i]
	}
	k4 := f(t)
	for i := range y {
		y[i] += h / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
	}
}

// KronSpice runs the machine in kronspice: windings shorted through r (a
// near-superconducting short; L/r far longer than the run), rotor inertia J.
func (m Machine) KronSpice(stop, step, r float64) (MTrace, float64, error) {
	ckt, err := circuit.NewBuilder("cyclic windings").Domain(exchange.Mechanical, "w").
		MachineModel("mm", m.Spec()).
		Resistor("R1", "a", "0", r).Resistor("R2", "b", "0", r).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", m.J, m.W0).BuildTransient()
	if err != nil {
		return MTrace{}, 0, err
	}
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-6}, MaxSteps: 5_000_000})
	if err != nil {
		return MTrace{}, 0, fmt.Errorf("kronspice: %w", err)
	}
	w, _ := wf.NodeVoltage("w")
	var tr MTrace
	th := m.Theta0
	for k := range wf.Time {
		if k > 0 {
			th += 0.5 * (w[k] + w[k-1]) * (wf.Time[k] - wf.Time[k-1])
		}
		tr.Time, tr.Theta, tr.Omega = append(tr.Time, wf.Time[k]), append(tr.Theta, th), append(tr.Omega, w[k])
	}
	l := wf.Ledger
	mis := l.IntegrationMismatch / math.Max(l.StoredStart, 1e-300)
	return tr, mis, nil
}
