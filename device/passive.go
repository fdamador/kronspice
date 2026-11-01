// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Resistor is a linear resistor.
//
// Maxwell derivation (CLAUDE.md §3c table, row "Resistor"): Ohm's law
// 𝔎 = C𝔈 — conduction current proportional to electromotive intensity —
// integrated along a uniform conductor in the quasi-static limit gives v = R·i.
//
// Ports: thermal (out) is intrinsic — the Joule heat i²R leaves the circuit to
// the thermal exterior. Electrical terminals connect inside the network and
// exchange nothing with the exterior. All other channels: not modeled.
//
// Negative resistance (R < 0). Under sign inversion (time reversal of the
// current, or a negative-impedance converter) v = R·i becomes v = −|R|·i:
// dissipation turns into generation, and the element delivers v²/|R| to the
// circuit. It is an idealized active element, like a controlled source, so
// that power is reported as model-supplied (device.ModelSupplier), not as
// heat flowing in through the thermal port: a real negative-resistance
// device (tunnel or Gunn diode, op-amp converter) takes it from a bias
// supply that the ideal element leaves out. Its thermal port is therefore
// not modeled. (The network dual of a resistor, swapping the roles of v and
// i, is the conductance 1/R — a positive resistor; capacitor and inductor
// are each other's duals.)
type Resistor struct {
	Base
	R  float64
	g  float64
	ch exchange.ChannelID // where the absorbed power goes (Thermal by default)
}

// NewResistor builds a resistor from a circuit element. Its absorbed
// power leaves through its thermal port, or through e.Channel when one is
// set: an acoustic radiation resistance sends sound out through the
// acoustic port, a mechanical load does work outside the circuit through
// the mechanical port (CLAUDE.md §3a.1 item 3).
func NewResistor(e circuit.Element) *Resistor {
	r := &Resistor{Base: newBase(e), R: e.Value, ch: exchange.Thermal}
	if e.Channel != exchange.Electrical {
		r.ch = e.Channel
	}
	if r.R < 0 {
		_ = r.SetIntrinsic(exchange.Electrical, exchange.Both,
			"idealized negative resistance: delivers v²/|R| as model-supplied power")
		return r
	}
	_ = r.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network")
	if r.ch == exchange.Thermal {
		_ = r.SetIntrinsic(exchange.Thermal, exchange.Out, "Joule heating i²R to surroundings")
	} else {
		_ = r.SetIntrinsic(r.ch, exchange.Out, "absorbed power leaves through the "+r.ch.String()+" exterior")
	}
	return r
}

// Setup maps the resistor's nodes.
func (r *Resistor) Setup(ctx *SetupContext) error {
	r.setupNodes(ctx)
	r.g = 1 / r.R
	return nil
}

// Load stamps the conductance (Ho, Ruehli & Brennan 1975, MNA stamp for a conductance).
func (r *Resistor) Load(ctx *LoadContext) error {
	a, b := r.u[0], r.u[1]
	ctx.M.Add(a, a, r.g)
	ctx.M.Add(b, b, r.g)
	ctx.M.Add(a, b, -r.g)
	ctx.M.Add(b, a, -r.g)
	return nil
}

// Operating returns voltage, current and absorbed power.
func (r *Resistor) Operating(x []float64) OperatingPoint {
	v := r.voltage(x)
	i := v * r.g
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// ExteriorPower: heat leaves through the thermal port (R > 0 only).
func (r *Resistor) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == r.ch && r.R > 0 {
		return -op.PowerAbsorbed
	}
	return 0
}

// ModelSuppliedPower is the power a negative resistance delivers (0 for R > 0).
func (r *Resistor) ModelSuppliedPower(op OperatingPoint) float64 {
	if r.R < 0 {
		return -op.PowerAbsorbed
	}
	return 0
}

// LoadAC stamps the conductance.
func (r *Resistor) LoadAC(ctx *ACContext) error {
	stampY(ctx.M, r.u[0], r.u[1], complex(r.g, 0))
	return nil
}

// OperatingAC returns the resistor's phasors.
func (r *Resistor) OperatingAC(x []complex128) ACPoint {
	v := r.voltageAC(x)
	return ACPoint{Voltage: v, Current: v * complex(r.g, 0)}
}
