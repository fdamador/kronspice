// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Controlled is a linear controlled source (SPICE3 E, F, G, H; Quarles et
// al., SPICE3f3 User's Manual §3.2; stamps as in Vlach & Singhal, Computer
// Methods for Circuit Analysis and Design, ch. 4, Table 4.3):
//
//	E (VCVS):  V(n+) − V(n−) = k·(V(nc+) − V(nc−))   branch unknown i
//	G (VCCS):  i(n+ → n−)    = k·(V(nc+) − V(nc−))
//	F (CCCS):  i(n+ → n−)    = k·I(Vc)
//	H (CCVS):  V(n+) − V(n−) = k·I(Vc)                 branch unknown i
//
// I(Vc) is the branch current of the controlling voltage source, positive
// into its n+ terminal (SPICE convention). The control port draws no
// current and no power.
//
// Maxwell derivation (CLAUDE.md §3c): none on its own — a controlled source
// is an idealization of an active device whose supply is left out of the
// circuit. Its power, of either sign, is therefore booked as model-supplied
// in the ledger, never as inflow from the exterior.
type Controlled struct {
	Base
	gain     float64
	cu       [2]int // control node unknowns (E, G)
	ctrlName string // control source (F, H)
	ck       int    // control branch unknown (F, H)
	branch   int    // own branch unknown (E, H)
	ac       acStamper
}

// acStamper adapts a complex stamper to real stamps, without allocating.
type acStamper struct{ m CStamper }

func (a *acStamper) Add(i, j int, v float64) { a.m.Add(i, j, complex(v, 0)) }

// NewControlled builds an E, F, G or H element.
func NewControlled(e circuit.Element) *Controlled {
	c := &Controlled{Base: newBase(e), gain: e.Value, cu: e.CtrlNodes, ctrlName: e.CtrlSource, ck: -1, branch: -1}
	_ = c.SetIntrinsic(exchange.Electrical, exchange.Both,
		"idealized controlled source: its power, of either sign, is model-supplied")
	return c
}

// Setup maps nodes and allocates the branch of a voltage output.
func (c *Controlled) Setup(ctx *SetupContext) error {
	c.setupNodes(ctx)
	if c.kind.VoltageControlled() {
		c.cu = [2]int{ctx.Unknown(c.cu[0]), ctx.Unknown(c.cu[1])}
	}
	if c.kind.VoltageOutput() {
		c.branch = ctx.AllocBranch()
		// A CCCS or CCVS may sense this source's current, as in ngspice.
		ctx.RegisterBranch(c.name, c.branch)
	}
	return nil
}

// Link finds the controlling source's branch (F, H).
func (c *Controlled) Link(ctx *SetupContext) error {
	if c.kind.VoltageControlled() {
		return nil
	}
	k, ok := ctx.BranchOf(c.ctrlName)
	if !ok {
		return fmt.Errorf("%w: %s: controlling source %s not found", ErrUnsupported, c.name, c.ctrlName)
	}
	c.ck = k
	return nil
}

// stamp writes the element's matrix entries (it has no right-hand side).
func (c *Controlled) stamp(m Stamper) {
	add := m.Add
	p, n := c.u[0], c.u[1]
	g := c.gain
	if c.kind.VoltageOutput() {
		k := c.branch
		add(p, k, 1)
		add(n, k, -1)
		add(k, p, 1)
		add(k, n, -1)
		if c.kind == circuit.VCVS {
			add(k, c.cu[0], -g)
			add(k, c.cu[1], g)
		} else {
			add(k, c.ck, -g)
		}
		return
	}
	if c.kind == circuit.VCCS {
		add(p, c.cu[0], g)
		add(p, c.cu[1], -g)
		add(n, c.cu[0], -g)
		add(n, c.cu[1], g)
	} else {
		add(p, c.ck, g)
		add(n, c.ck, -g)
	}
}

// Load stamps the linear element; it has no state and no right-hand side.
func (c *Controlled) Load(ctx *LoadContext) error {
	c.stamp(ctx.M)
	return nil
}

func (c *Controlled) control(x []float64) float64 {
	if c.kind.VoltageControlled() {
		return at(x, c.cu[0]) - at(x, c.cu[1])
	}
	return at(x, c.ck)
}

// Operating reports the output port's voltage and current (n+ → n−
// through the element).
func (c *Controlled) Operating(x []float64) OperatingPoint {
	v := c.voltage(x)
	var i float64
	if c.kind.VoltageOutput() {
		i = at(x, c.branch)
	} else {
		i = c.gain * c.control(x)
	}
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// ExteriorPower: a controlled source has no exterior port of its own.
func (c *Controlled) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// ModelSuppliedPower is all of the element's delivered power.
func (c *Controlled) ModelSuppliedPower(op OperatingPoint) float64 { return -op.PowerAbsorbed }

// LoadAC stamps the same linear relation with a real gain.
func (c *Controlled) LoadAC(ctx *ACContext) error {
	c.ac.m = ctx.M
	c.stamp(&c.ac)
	return nil
}

// OperatingAC returns the output port's phasors.
func (c *Controlled) OperatingAC(x []complex128) ACPoint {
	v := c.voltageAC(x)
	if c.kind.VoltageOutput() {
		return ACPoint{Voltage: v, Current: catC(x, c.branch)}
	}
	var ctl complex128
	if c.kind == circuit.VCCS {
		ctl = catC(x, c.cu[0]) - catC(x, c.cu[1])
	} else {
		ctl = catC(x, c.ck)
	}
	return ACPoint{Voltage: v, Current: complex(c.gain, 0) * ctl}
}
