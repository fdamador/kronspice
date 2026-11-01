// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command tether models an electrodynamic tether in low Earth orbit — a
// conductor pushing and pulling against a planetary magnetic field
// (doc/em-push-pull.md; CLAUDE.md §3f; milestone M9n).
//
// A vertical aluminium wire of length L moves east at orbital speed v
// through a northward field B. The Lorentz force on its charges, by grade
// selection (em.MotionalEMF), gives an EMF (v×B)·L along it. The circuit
// closes through the ionosphere by plasma contactors at the two ends:
//
//	EMF — tether resistance — contactors — load (generator)
//	EMF — tether resistance — contactors — supply above the EMF (thruster)
//
// The current feels I·L×B. As a generator the force is drag: the power the
// load receives comes from the orbit's kinetic energy (F·v = −EMF·I). With
// a supply driving the current against the EMF the force turns into
// thrust. The momentum ledger books the force on the tether with its
// reaction on the geomagnetic field and ionosphere (the exterior source).
//
// Inputs are declared, not looked up: v = 7.7 km/s, B = 30 µT (order of
// magnitude in LEO; a real orbit takes B from a cited field model such as
// IGRF), contactor resistances 100 Ω each (an illustrative linear stand-in
// for the contactors' nonlinear current–voltage characteristics).
//
//	go run ./examples/tether
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/material"
)

// Tether is one configuration.
type Tether struct {
	Length, WireDiameter float64 // [m]
	V, B                 float64 // orbital speed east [m/s], field north [T]
	RContactor           float64 // each end [Ω]
	Load                 float64 // generator load [Ω] (0: thruster mode)
	Supply               float64 // thruster supply [V], opposing the EMF
}

// Defaults is a 5 km, 1 mm aluminium tether.
func Defaults() Tether {
	return Tether{Length: 5000, WireDiameter: 1e-3, V: 7700, B: 30e-6, RContactor: 100, Load: 500}
}

// EMF is ∫(v×B)·dl along the tether (upward), from the Lorentz force.
func (t Tether) EMF() float64 {
	var path [][3]float64
	for k := 0; k <= 20; k++ {
		path = append(path, [3]float64{0, 0, t.Length * float64(k) / 20})
	}
	return em.MotionalEMF(path,
		func([3]float64) [3]float64 { return [3]float64{t.V, 0, 0} },
		func([3]float64) [3]float64 { return [3]float64{0, t.B, 0} })
}

// Resistance of the wire.
func (t Tether) Resistance() float64 {
	a := math.Pi * t.WireDiameter * t.WireDiameter / 4
	return material.Aluminium.Resistivity * t.Length / a
}

// Result of one configuration.
type Result struct {
	EMF, Current  float64 // current upward along the tether [A]
	Force         [3]float64
	PElectrical   float64 // power into the load (generator) or from the supply (thruster) [W]
	PMechanical   float64 // F·v: < 0 drag (orbit loses energy), > 0 thrust [W]
	Momentum      em.MomentumLedger
	LedgerResidue float64
}

// Solve builds the circuit, solves it and books the force.
func (t Tether) Solve() (*Result, error) {
	emf := t.EMF()
	b := circuit.NewBuilder("electrodynamic tether")
	// Bottom end "b", top end "top": the EMF raises the top.
	b.VoltageSource("VEMF", "top", "m", emf).
		Resistor("RT", "m", "bot", t.Resistance()).
		Resistor("RCT", "top", "x", t.RContactor).
		Resistor("RCB", "bot", "y", t.RContactor)
	if t.Load > 0 {
		b.Resistor("RL", "x", "y", t.Load)
	} else {
		// The supply drives current down the tether, against the EMF.
		b.VoltageSource("VS", "x", "y", t.Supply)
	}
	b.Resistor("RREF", "y", "0", 1e9) // reference the floating loop
	c, err := b.Build()
	if err != nil {
		return nil, err
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return nil, err
	}
	r := &Result{EMF: emf}
	vemf, _ := op.Element("VEMF")
	// SPICE reports a source's current flowing n+ → n− through it (top →
	// m, downward), so the current up the tether is its negative.
	r.Current = -vemf.Current
	// Force on the tether ∫ I dl × B, dl = L·ẑ, B = B·ŷ: I·L·B·(ẑ×ŷ) = −I·L·B·x̂.
	f := [3]float64{-r.Current * t.Length * t.B, 0, 0}
	r.Force = f
	r.Momentum.AddExterior("tether", "geomagnetic field and ionosphere", f, false)
	r.PMechanical = f[0] * t.V
	if t.Load > 0 {
		l, _ := op.Element("RL")
		r.PElectrical = l.PowerAbsorbed
	} else {
		s, _ := op.Element("VS")
		r.PElectrical = -s.PowerAbsorbed
	}
	r.LedgerResidue = op.Ledger.RelativeResidual()
	return r, nil
}

func main() {
	t := Defaults()
	fmt.Printf("tether %.0f m, wire %.1f mm (%.0f Ω), v %.0f m/s, B %.0f µT\n", t.Length, t.WireDiameter*1e3, t.Resistance(), t.V, t.B*1e6)
	fmt.Printf("EMF from the Lorentz force %.1f V (vBL = %.1f V)\n\n", t.EMF(), t.V*t.B*t.Length)
	r, err := t.Solve()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("generator, load %.0f Ω: I = %.3f A, %.0f W to the load; force %.3f N along x (drag); F·v = %.0f W — the orbit pays\n",
		t.Load, r.Current, r.PElectrical, r.Force[0], r.PMechanical)
	th := t
	th.Load, th.Supply = 0, 3000
	r, err = th.Solve()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("thruster, supply %.0f V: I = %.3f A, %.0f W from the supply; thrust %.3f N; F·v = %.0f W into the orbit; %.2g N/kW\n",
		th.Supply, r.Current, r.PElectrical, r.Force[0], r.PMechanical, r.Force[0]/(r.PElectrical/1e3))
	for src, f := range r.Momentum.ByExterior() {
		fmt.Printf("momentum ledger: reaction %.3f N on %s\n", -f[0], src)
	}
}
