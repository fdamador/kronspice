// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/material"
)

// Stack is an aluminium-water rectifier built from plates: Ns cells in
// series (to hold more voltage), each cell Np plates in parallel (to carry
// more current). One cell is a formed aluminium plate facing an iron (or
// lead) plate across a gap of electrolyte — the valve cell of
// examples/sparkgap_diode (Pollak, Nodon; EB 1911), with the oxide's
// leakage and capacitance taken from capacitor-foil data:
//
//   - Forward (aluminium negative): electrolysis, onset Vd per cell, then
//     the electrolyte's resistance ρ·gap/(area·Np) per cell.
//   - Reverse (aluminium positive): the anodic barrier oxide, formed at
//     Vform, blocks. Its leakage is ionic conduction by the
//     Güntherschulze–Betz law j = j0·exp(β·E), E = v/d, d = k·Vform: it
//     climbs steeply as the voltage nears Vform, as capacitor makers'
//     curves show. β = 3·10⁻⁸ m/V and k = 1.4 nm/V are in the published
//     ranges for anodic alumina (β 10⁻⁹…5·10⁻⁸ m/V; 1.2–1.5 nm/V); j0 is
//     set by the leakage density at Vform, JForm (declared 0.1 A/m²,
//     10 µA/cm², a well-formed foil). Above Vform the oxide sparks: the
//     cell conducts through RSpark (and the design is rejected).
//   - The oxide is a capacitor across each cell, ε_r·ε0·area·Np/d with
//     ε_r = 9.5 (EPCOS data). A barrier oxide has a fixed thickness, so its
//     capacitance hardly changes with voltage; Kappa, the fractional change
//     of capacitance from 0 to Vform, is declared small (2 %) and makes the
//     film a weak varactor: q(v) = C0·(v + Kappa·v·|v|/(2·Vform)).
type Stack struct {
	Ns, Np int
	Area   float64 // one plate [m²]
	Gap    float64 // electrolyte gap [m]
	Vform  float64 // forming voltage of the oxide [V]
	Kappa  float64 // varactor coefficient of the oxide film
}

// Declared electrolyte and oxide constants.
const (
	rhoWater = 0.1    // Ω·m, a borate/phosphate forming electrolyte (declared; as sparkgap_diode)
	vdCell   = 2.0    // V, electrolysis onset per cell (as sparkgap_diode)
	betaGB   = 3e-8   // m/V, Güntherschulze–Betz field coefficient
	kOx      = 1.4e-9 // m/V, oxide thickness per forming volt
	jForm    = 0.1    // A/m², leakage density at the forming voltage (declared)
	epsOx    = 9.5    // relative permittivity of Al2O3 (EPCOS)
	rSpark   = 2.0    // Ω per cell while sparking (as sparkgap_diode)
	kappaDef = 0.02   // varactor coefficient (declared)
	maxForm  = 700.0  // V, highest forming voltage in capacitor-foil practice (≈ 700–760 V)
)

// Rw is the forward (electrolyte) resistance of the stack [Ω].
func (s Stack) Rw() float64 { return float64(s.Ns) * rhoWater * s.Gap / (s.Area * float64(s.Np)) }

// C0 is the oxide capacitance of the stack (Ns cells in series) [F].
func (s Stack) C0() float64 {
	return epsOx * material.Eps0 * s.Area * float64(s.Np) / (kOx * s.Vform) / float64(s.Ns)
}

// Hold is the largest reverse voltage the stack holds before sparking [V].
func (s Stack) Hold() float64 { return float64(s.Ns) * s.Vform }

// Current is the static current [A] at stack voltage v, positive in the
// conducting direction (v > 0: aluminium negative).
func (s Stack) Current(v float64) float64 {
	ns := float64(s.Ns)
	if v >= 0 {
		vc := v/ns - vdCell
		if vc <= 0 {
			return 0
		}
		return vc * ns / s.Rw()
	}
	vr := -v / ns // reverse voltage per cell
	area := s.Area * float64(s.Np)
	// β·E at the forming voltage is β/k (the field there is 1/k), so
	// j = jForm·exp((β/k)·(vr/Vform − 1)), less its value at zero volts.
	bef := betaGB / kOx
	if vr <= s.Vform {
		return -area * jForm * (math.Exp(bef*(vr/s.Vform-1)) - math.Exp(-bef))
	}
	// Sparking: the excess voltage drives current through the sparks.
	excess := (vr - s.Vform) * ns
	return -(area*jForm*(1-math.Exp(-bef)) + excess/(ns*rSpark+s.Rw()))
}

// curve samples the stack's i(v) densely (the leakage is exponential).
func (s Stack) curve() []float64 {
	var out []float64
	hold := s.Hold()
	nf := float64(s.Ns) * vdCell
	vmax := 1.5 * hold
	// reverse: 300 points to −1.5·hold, denser near the forming voltage
	for j := 300; j >= 1; j-- {
		u := float64(j) / 300
		v := -vmax * u * u
		out = append(out, v, s.Current(v))
	}
	out = append(out, 0, 0)
	for j := 1; j <= 200; j++ {
		v := nf + (2*hold-nf)*math.Pow(float64(j)/200, 2)
		if j == 1 {
			out = append(out, nf, 0)
		}
		out = append(out, v, s.Current(v))
	}
	return out
}

// Add puts the stack in like a diode from anode to cathode: it conducts
// when V(anode) > V(cathode) (the aluminium plates face the cathode side:
// aluminium negative is the conducting direction). The oxide film is
// across it.
func (s Stack) Add(b *circuit.Builder, name, anode, cathode string) {
	b.NonlinearVC(name, anode, cathode, s.curve())
	c0 := s.C0()
	if s.Kappa == 0 {
		b.Capacitor("C"+name, anode, cathode, c0)
		return
	}
	var q []float64
	hold := s.Hold()
	for j := -300; j <= 300; j++ {
		v := 1.5 * hold * float64(j) / 300
		q = append(q, v, c0*(v+s.Kappa*v*math.Abs(v)/(2*hold)))
	}
	b.NonlinearCapacitor("C"+name, anode, cathode, q, circuit.Waveform{})
}

// Describe prints the stack.
func (s Stack) Describe() string {
	return fmt.Sprintf("%s, formed at %.0f V: holds %.0f V, forward %.1f V + %.3g Ω, film %.3g µF",
		s.Short(), s.Vform, s.Hold(), float64(s.Ns)*vdCell, s.Rw(), s.C0()*1e6)
}

// Short names the stack's geometry: cells in series, plates in parallel,
// plate area and gap.
func (s Stack) Short() string {
	return fmt.Sprintf("%d cells × %d plates of %.3g mm², gap %.1f mm", s.Ns, s.Np, s.Area*1e6, s.Gap*1e3)
}
