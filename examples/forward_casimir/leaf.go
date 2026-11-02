// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// One leaf pair of Forward's battery as a Kron electromechanical system,
// per square metre of leaf.
//
// Two parallel aluminium leaves (thickness t, mass m = ρ·t per area each)
// face each other across the gap a. Both are held at the supply voltage V
// (same polarity, Forward §IV), so there is no field in the gap; each leaf
// faces a grounded neighbour at distance D = (L − a)/2, so the pair is a
// capacitor to ground of C(a) = 2ε0/D = 4ε0/(L − a) per area. The field
// outside pulls each leaf outward — Forward's electrostatic repulsion.
//
// Kron's coordinates are the charge q (per area) and the gap a. The
// Lagrange–Maxwell equations (Kron 1934; Hoffmann 1949, eqs. 1–8) with
// kinetic energy ¼m·ȧ², potential q²/2C(a) + E_cas(a), dissipation
// ½r·q̇², and the supply voltage V(t) impressed on q:
//
//	r·q̇ + q/C(a) = V
//	½m·ä = ½q²·C′(a)/C² − E′_cas(a)          (outward push − Casimir pull)
//
// The supply is servo-controlled, as Forward describes: at each instant it
// sets V so the outward push equals the Casimir pull plus what is needed
// to move the leaves at a steady speed v (closing, then reopening).

// Leaf is the leaf pair and its servo.
type Leaf struct {
	Mat      Material
	T        float64 // leaf thickness [m]
	L        float64 // a + 2D: spacing between the grounded neighbours [m]
	R        float64 // supply resistance × area [Ω·m²]
	A0, AMin float64 // starting and closest gap [m]
	V        float64 // servo speed [m/s]
	Gain     float64 // servo speed gain [1/s]
	tab      Table
}

// DefaultLeaf (declared): 100 nm aluminium leaves, 2 µm to the grounded
// neighbours, from 1 µm down to 20 nm at 1 mm/s.
func DefaultLeaf(m Material) Leaf {
	return Leaf{Mat: m, T: 100e-9, L: 3e-6, R: 1e-6, A0: 1e-6, AMin: 20e-9, V: 1e-3, Gain: 1e5}
}

func (l *Leaf) table() Table {
	if l.tab.la == nil {
		l.tab = NewTable(l.Mat, 5e-9, 10e-6, 400)
	}
	return l.tab
}

// Cap is C(a) per area [F/m²].
func (l Leaf) Cap(a float64) float64 { return 4 * Eps0 / (l.L - a) }

// Stroke is one run of the leaf pair.
type Stroke struct {
	Time, Gap, Volt, Charge []float64
	Casimir                 float64 // work done by the Casimir force [J/m²]
	ToSupply                float64 // energy delivered into the supply [J/m²]
	Heat                    float64 // in the supply resistance [J/m²]
	Stored                  float64 // change of electric energy ½q²/C [J/m²]
	Kinetic                 float64 // change of kinetic energy [J/m²]
	Residual                float64 // Casimir − (supply + heat + stored + kinetic)
	Net                     float64 // to the supply over charge-up, stroke and discharge: ToSupply + Stored [J/m²]
	ChargeIn                float64 // charge the supply sends into the leaves during the stroke [C/m²]
	Electrons               float64 // electrons on the leaves at the end (all returned at discharge) per m²
	PerElectron             float64 // Net ÷ (end charge): energy gained per electron cycled [eV]
	PeakVolt                float64
}

// Run moves the gap from a0 to a1 (closing if a1 < a0) under the servo,
// with the electrical side solved exactly for each step (backward Euler on
// q; the RC time r·C is far shorter than the motion).
func (l *Leaf) Run(a0, a1 float64) Stroke {
	tab := l.table()
	m := l.Mat.Density * l.T
	dir := math.Copysign(1, a1-a0)
	a, v := a0, 0.0
	cv := func(a float64) float64 { // servo voltage for gap a and speed v
		want := dir * l.V
		acc := l.Gain * (want - v) // desired ä
		push := -tab.Pressure(a) + 0.5*m*acc
		if push < 0 {
			push = 0
		}
		// push = ½V²·C′/C·... per area: ½q²C′/C² = ½V²C′
		cp := 4 * Eps0 / ((l.L - a) * (l.L - a))
		return math.Sqrt(2 * push / cp)
	}
	V := cv(a)
	q := l.Cap(a) * V
	var s Stroke
	e0 := 0.5 * q * q / l.Cap(a)
	s.Time, s.Gap, s.Volt, s.Charge = []float64{0}, []float64{a}, []float64{V}, []float64{q} // the start point
	t := 0.0
	dt := 1e-9
	for (a1-a)*dir > 0 {
		// adapt the step: move at most a/200000 per step (forces vary as a⁻⁴)
		dt = math.Min(1e-7, math.Max(1e-13, a/(200000*math.Max(math.Abs(v), 1e-6))))
		V = cv(a)
		cA := l.Cap(a)
		// mechanics: semi-implicit Euler with the current charge
		cp := 4 * Eps0 / ((l.L - a) * (l.L - a))
		push := 0.5 * q * q * cp / (cA * cA)
		pull := tab.Pressure(a) // negative
		acc := 2 * (push + pull) / m
		vn := v + acc*dt
		an := a + vn*dt
		ec := tab.Energy(a) - tab.Energy(an) // Casimir work = −ΔE
		s.Casimir += ec
		// electrical: r·(qn − q)/dt + qn/C(an) = Vn
		cn := l.Cap(an)
		Vn := cv(an)
		qn := (Vn + l.R*q/dt) / (l.R/dt + 1/cn)
		i := (qn - q) / dt
		s.Heat += l.R * i * i * dt
		s.ToSupply -= Vn * i * dt // the supply receives −V·i
		// mechanical work by the electric push on the gap
		s.Kinetic += 0.25 * m * (vn*vn - v*v)
		s.ChargeIn += qn - q
		a, v, q = an, vn, qn
		t += dt
		s.PeakVolt = math.Max(s.PeakVolt, Vn)
		if t-s.Time[len(s.Time)-1] > 2e-7 || math.Abs(an-s.Gap[len(s.Gap)-1]) > 0.002*an {
			s.Time, s.Gap, s.Volt, s.Charge = append(s.Time, t), append(s.Gap, a), append(s.Volt, Vn), append(s.Charge, q)
		}
	}
	s.Time, s.Gap, s.Volt, s.Charge = append(s.Time, t), append(s.Gap, a), append(s.Volt, q/l.Cap(a)), append(s.Charge, q) // the end point
	s.Stored = 0.5*q*q/l.Cap(a) - e0
	s.Residual = s.Casimir - (s.ToSupply + s.Heat + s.Stored + s.Kinetic)
	// The supply first charges the pair to its starting voltage (paying the
	// starting field energy) and, after the stroke, with the leaves held at
	// a mechanical stop, takes the final field energy back by lowering V
	// slowly (lossless in the limit). Over the whole operation it gains
	// ToSupply − e_start + e_end = ToSupply + Stored.
	s.Net = s.ToSupply + s.Stored
	s.Electrons = q / Qe
	s.PerElectron = s.Net / q
	return s
}
