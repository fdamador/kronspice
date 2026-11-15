// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Helix is the core field of one pass decomposed in the frame of the
// spring's axis. The flux runs along the wire, whose tangent turns once per
// pitch: in the plane across the axis the field's direction rotates with
// height (a spatially circular "polarization", wavenumber κ = 2π/p), and a
// small part points along the axis.
type Helix struct {
	B        float64 // flux density along the wire, Φ/A [T]
	Alpha    float64 // pitch angle, atan(p/2πR) [rad]
	Trans    float64 // transverse amplitude B·cos α (B(0) of the B Cyclic theorem) [T]
	Axial    float64 // axial component in the wire, B·sin α [T]
	Kappa    float64 // rotation wavenumber 2π/p [1/m]
	B3Cyclic float64 // |B(3)| from B(1)×B(2) = i·B(0)·B(3)* : equals B(0) [T]
	Product  float64 // standard reading: |B(1)×B(2)| = B(0)² [T²], an axial vector, not a field
}

// HelixField decomposes the field of pass flux phi.
func (s Spring) HelixField(phi float64) Helix {
	b := phi / s.Area()
	a := math.Atan(s.P / (2 * math.Pi * s.R))
	h := Helix{B: b, Alpha: a, Trans: b * math.Cos(a), Axial: b * math.Sin(a), Kappa: 2 * math.Pi / s.P}
	h.B3Cyclic = h.Trans
	h.Product = h.Trans * h.Trans
	return h
}

// PickupPerDPhi is the EMF per unit dΦ/dt of the three ways a coil can meet
// the helical flux: a turn encircling the whole spring links the net axial
// flux once; the straight wire up the bore links the circulation N times;
// a wall winding turn links the n passes it encloses.
func (s Spring) PickupPerDPhi(n int) (encircling, bore, wall float64) {
	return 1, float64(s.N), float64(n)
}

// MagAmp is a magnetic-amplifier estimate for the spring core: a control
// winding and a gate (load) winding, each round some of the passes. The
// core is linear in this example, so the numbers are first-order: the
// ampere-turns that bring the core to saturation, the gate's volt-second
// capacity, and the gains from the ampere-turn balance of a
// self-saturating magnetic amplifier (Mali 1960, pp. 34–35, 50–55):
// I_load·N_g·n_g ≈ I_c·N_c·n_c.
type MagAmp struct {
	Fsat    float64 // effective ampere-turns round the core to reach Bs [A]
	IcSat   float64 // control current for Fsat [A]
	Vgate   float64 // largest AC voltage (rms) the gate blocks before saturating, at f [V]
	Ki      float64 // current gain (N_c·n_c)/(N_g·n_g)
	Kp      float64 // power gain Ki²·R_load/R_control
	Rc, Rg  float64 // control and gate winding resistance [Ω]
	KiPlain float64 // the same windings on an ordinary core: N_c/N_g
}

// MagAmpEstimate evaluates control winding c and gate winding g at f into
// a load rl.
func (s Spring) MagAmpEstimate(c, g Winding, f, rl, cu float64) MagAmp {
	phiSat := s.Mat.Bsat * s.Area()
	rtot := float64(s.N*s.S)*s.SegReluctance() + s.ReturnReluctance()
	var m MagAmp
	m.Fsat = phiSat * rtot
	m.IcSat = m.Fsat / float64(c.Turns*c.Passes)
	// reset to −Bs, set to +Bs in one half cycle: 2·λsat = ∫v dt over T/2
	m.Vgate = 4.44 * f * float64(g.Turns*g.Passes) * phiSat
	m.Ki = float64(c.Turns*c.Passes) / float64(g.Turns*g.Passes)
	m.KiPlain = float64(c.Turns) / float64(g.Turns)
	m.Rc, m.Rg = s.WindingResistance(c, cu), s.WindingResistance(g, cu)
	m.Kp = m.Ki * m.Ki * rl / m.Rc
	return m
}
