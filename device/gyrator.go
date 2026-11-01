// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Gyrator is an ideal, lossless two-port transducer between two energy
// domains (or within one):
//
//	i1 = g·V2,   i2 = −g·V1     (each current p+ → p− through its port)
//
// so the power absorbed at one port is delivered at the other,
// V1·i1 + V2·i2 = 0, and the element neither stores nor dissipates energy.
//
// In the mobility analogy (across = velocity, through = force for the
// mechanical domain; across = pressure, through = volume velocity for the
// acoustic one; exchange channel table) the classic transducers are
// gyrators (F. A. Firestone, "A new analogy between mechanical and
// electrical systems," J. Acoust. Soc. Am. 4, 1933; L. L. Beranek,
// T. J. Mellow, Acoustics: Sound Fields and Transducers, 2012, Ch. 3):
//
//   - a piezoelectric or electrostatic transducer with force factor φ
//     [N/V]: motional current i = φ·v, blocked force F = φ·V (g = φ);
//   - a piston of area S between a mechanical and an acoustic port:
//     volume velocity U = S·v, force F = S·p (g = S).
//
// Maxwell derivation (CLAUDE.md §3c): the electromechanical coupling is
// the force of the field on the charges of the transducer's electrodes and
// the current their motion drives — the Lorentz-force row; in the lumped
// limit it is linear, and lossless because the field's energy is booked in
// the separate capacitor (C0) that the user puts across the electrical
// port.
//
// MNA: both ports are voltage-controlled current sources, stamped into the
// node rows (Ho, Ruehli & Brennan 1975); no branch unknown.
type Gyrator struct {
	Base
	g  float64
	n2 [2]int // port-2 nodes
	u2 [2]int // port-2 unknowns
}

// NewGyrator builds a transducer from a circuit element.
func NewGyrator(e circuit.Element) *Gyrator {
	y := &Gyrator{Base: newBase(e), g: e.Value, n2: [2]int{e.Extra[0], e.Extra[1]}}
	_ = y.SetIntrinsic(exchange.Electrical, exchange.Both,
		"lossless transducer: power in at one port leaves at the other, inside the network")
	return y
}

// Setup maps both ports' nodes.
func (y *Gyrator) Setup(ctx *SetupContext) error {
	y.setupNodes(ctx)
	y.u2 = [2]int{ctx.Unknown(y.n2[0]), ctx.Unknown(y.n2[1])}
	return nil
}

// stamp writes the two controlled currents.
func (y *Gyrator) stamp(add func(i, j int, v float64)) {
	a, b, c, d := y.u[0], y.u[1], y.u2[0], y.u2[1]
	g := y.g
	// i1 = g·(V(c) − V(d)) leaves a, enters b.
	add(a, c, g)
	add(a, d, -g)
	add(b, c, -g)
	add(b, d, g)
	// i2 = −g·(V(a) − V(b)) leaves c, enters d.
	add(c, a, -g)
	add(c, b, g)
	add(d, a, g)
	add(d, b, -g)
}

// Load stamps the linear two-port.
func (y *Gyrator) Load(ctx *LoadContext) error {
	y.stamp(ctx.M.Add)
	return nil
}

func (y *Gyrator) v2(x []float64) float64 { return at(x, y.u2[0]) - at(x, y.u2[1]) }

// Operating reports port 1: its voltage and current, and the element's
// total absorbed power V1·i1 + V2·i2 (zero up to round-off).
func (y *Gyrator) Operating(x []float64) OperatingPoint {
	v1, v2 := y.voltage(x), y.v2(x)
	i1, i2 := y.g*v2, -y.g*v1
	return OperatingPoint{Voltage: v1, Current: i1, PowerAbsorbed: v1*i1 + v2*i2}
}

// PortPower returns the power absorbed at port 1 (and so delivered at port
// 2) at solution x: the power the transducer converts.
func (y *Gyrator) PortPower(x []float64) float64 {
	return y.voltage(x) * y.g * y.v2(x)
}

// ExteriorPower: the transducer exchanges nothing with the exterior.
func (y *Gyrator) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// LoadAC stamps the same two-port.
func (y *Gyrator) LoadAC(ctx *ACContext) error {
	y.stamp(func(i, j int, v float64) { ctx.M.Add(i, j, complex(v, 0)) })
	return nil
}

func (y *Gyrator) v2AC(x []complex128) complex128 {
	var v complex128
	if y.u2[0] >= 0 {
		v += x[y.u2[0]]
	}
	if y.u2[1] >= 0 {
		v -= x[y.u2[1]]
	}
	return v
}

// OperatingAC returns port 1's phasors.
func (y *Gyrator) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: y.voltageAC(x), Current: complex(y.g, 0) * y.v2AC(x)}
}

// ACPower is ½·(V1·i1* + V2·i2*), zero up to round-off.
func (y *Gyrator) ACPower(x []complex128) complex128 {
	v1, v2 := y.voltageAC(x), y.v2AC(x)
	g := complex(y.g, 0)
	return 0.5 * (v1*cmplx.Conj(g*v2) + v2*cmplx.Conj(-g*v1))
}
