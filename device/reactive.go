// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Capacitor is a linear capacitor.
//
// Maxwell derivation (CLAUDE.md §3c table, row "Capacitor"): the displacement
// current ∂𝔇/∂t with 𝔇 = K𝔈/4π, confined between the plates, integrates to
// i = C dv/dt in the quasi-static limit.
//
// Transient companion models (Nagel 1975; Vlach & Singhal ch. 10):
//
//	backward Euler: i_n = (C/h)(v_n − v_{n−1})
//	trapezoidal:    i_n = (2C/h)(v_n − v_{n−1}) − i_{n−1}
//
// each stamped as a conductance geq in parallel with a history current.
// The capacitor stores energy ½Cv² inside the circuit; its terminals connect
// inside the network and no exterior port is modeled (ideal, lossless).
type Capacitor struct {
	Base
	C     float64
	ic    float64
	hasIC bool

	geq, ihist float64 // companion model of the most recent Load
	tran       bool
	omega      float64 // angular frequency of the most recent LoadAC

	vPrev, iPrev float64 // accepted state
	hist         history // accepted charges for LTE
}

// NewCapacitor builds a capacitor from a circuit element.
func NewCapacitor(e circuit.Element) *Capacitor {
	c := &Capacitor{Base: newBase(e), C: e.Value, ic: e.IC, hasIC: e.HasIC}
	_ = c.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores ½Cv²")
	return c
}

// Setup maps the capacitor's nodes.
func (c *Capacitor) Setup(ctx *SetupContext) error {
	c.setupNodes(ctx)
	return nil
}

// uicConductance pins a capacitor to its initial voltage when reporting the
// t = 0 point of a UIC transient (a Norton equivalent of a voltage source).
const uicConductance = 1e9

// Load stamps an open circuit (DC), the UIC pin, or the companion model.
func (c *Capacitor) Load(ctx *LoadContext) error {
	a, b := c.u[0], c.u[1]
	switch ctx.Mode {
	case ModeOP, ModeTranOP:
		c.tran, c.geq, c.ihist = false, 0, 0
		return nil
	case ModeUICInit:
		c.tran = false
		g := uicConductance
		ctx.M.Add(a, a, g)
		ctx.M.Add(b, b, g)
		ctx.M.Add(a, b, -g)
		ctx.M.Add(b, a, -g)
		ctx.AddRHS(a, g*c.vPrev)
		ctx.AddRHS(b, -g*c.vPrev)
		c.geq, c.ihist = 0, 0
		return nil
	}
	c.tran = true
	if ctx.Method == Trapezoidal {
		c.geq = 2 * c.C / ctx.H
		c.ihist = -(c.geq*c.vPrev + c.iPrev)
	} else {
		c.geq = c.C / ctx.H
		c.ihist = -c.geq * c.vPrev
	}
	g := c.geq
	ctx.M.Add(a, a, g)
	ctx.M.Add(b, b, g)
	ctx.M.Add(a, b, -g)
	ctx.M.Add(b, a, -g)
	// Branch current i = g·v + ihist flows from n+ to n−.
	ctx.AddRHS(a, -c.ihist)
	ctx.AddRHS(b, c.ihist)
	return nil
}

// Operating returns voltage, current and absorbed power.
func (c *Capacitor) Operating(x []float64) OperatingPoint {
	v := c.voltage(x)
	i := 0.0
	if c.tran {
		i = c.geq*v + c.ihist
	}
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i, StoragePower: v * i}
}

// Accept commits the state. For UIC the initial voltage is the IC value.
func (c *Capacitor) Accept(x []float64, t float64, mode Mode) {
	switch mode {
	case ModeUICInit:
		// vPrev already holds the IC; keep it exactly.
		c.iPrev = 0
		c.hist.reset()
	case ModeOP, ModeTranOP:
		c.vPrev, c.iPrev = c.voltage(x), 0
		c.hist.reset()
	default:
		op := c.Operating(x)
		c.vPrev, c.iPrev = op.Voltage, op.Current
	}
	c.hist.push(t, c.C*c.vPrev)
}

// SetInitialCondition loads the UIC state (IC value, or 0 if none was given).
func (c *Capacitor) SetInitialCondition() {
	c.vPrev, c.iPrev = 0, 0
	if c.hasIC {
		c.vPrev = c.ic
	}
}

// StoredEnergy is ½Cv².
func (c *Capacitor) StoredEnergy(op OperatingPoint) float64 {
	return 0.5 * c.C * op.Voltage * op.Voltage
}

// Truncate estimates the step from the charge history.
func (c *Capacitor) Truncate(x []float64, tc *TruncContext) float64 {
	op := c.Operating(x)
	return lteStep(&c.hist, c.C*op.Voltage, op.Current, c.iPrev, tc.AbsTol, tc)
}

// ExteriorPower: an ideal capacitor exchanges nothing with the exterior.
func (c *Capacitor) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// Inductor is a linear inductor with a branch-current unknown (current from
// n+ to n− through the inductor).
//
// Maxwell derivation (CLAUDE.md §3c table, row "Inductor"): 𝔅 = V.∇𝔄 and the
// flux linkage of 𝔄 around the coil give v = L di/dt in the quasi-static limit.
//
// Branch equations:
//
//	DC:             v = 0
//	backward Euler: v_n = (L/h)(i_n − i_{n−1})
//	trapezoidal:    v_n = (2L/h)(i_n − i_{n−1}) − v_{n−1}
//
// The inductor stores ½Li² inside the circuit; no exterior port is modeled.
type Inductor struct {
	Base
	L      float64
	ic     float64
	hasIC  bool
	branch int

	vPrev, iPrev float64
	hist         history
}

// NewInductor builds an inductor from a circuit element.
func NewInductor(e circuit.Element) *Inductor {
	l := &Inductor{Base: newBase(e), L: e.Value, ic: e.IC, hasIC: e.HasIC}
	_ = l.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores ½Li²")
	return l
}

// Setup maps nodes and allocates the branch-current unknown.
func (l *Inductor) Setup(ctx *SetupContext) error {
	l.setupNodes(ctx)
	l.branch = ctx.AllocBranch()
	// A K element finds the inductor's branch by name (Mutual.Link).
	ctx.RegisterBranch(l.name, l.branch)
	return nil
}

// Load stamps the branch equation for the current mode.
func (l *Inductor) Load(ctx *LoadContext) error {
	a, b, k := l.u[0], l.u[1], l.branch
	ctx.M.Add(a, k, 1)
	ctx.M.Add(b, k, -1)
	switch ctx.Mode {
	case ModeOP, ModeTranOP:
		ctx.M.Add(k, a, 1)
		ctx.M.Add(k, b, -1)
		return nil
	case ModeUICInit:
		// Pin the branch current to the initial condition: i = IC.
		ctx.M.Add(k, k, 1)
		ctx.AddRHS(k, l.iPrev)
		return nil
	}
	req := l.L / ctx.H
	rhs := -req * l.iPrev
	if ctx.Method == Trapezoidal {
		req = 2 * l.L / ctx.H
		rhs = -req*l.iPrev - l.vPrev
	}
	ctx.M.Add(k, a, 1)
	ctx.M.Add(k, b, -1)
	ctx.M.Add(k, k, -req)
	ctx.AddRHS(k, rhs)
	return nil
}

// Operating returns voltage, current and absorbed power.
func (l *Inductor) Operating(x []float64) OperatingPoint {
	v := l.voltage(x)
	i := x[l.branch]
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i, StoragePower: v * i}
}

// Accept commits the state.
func (l *Inductor) Accept(x []float64, t float64, mode Mode) {
	switch mode {
	case ModeUICInit:
		// iPrev already holds the IC.
		l.vPrev = 0
		l.hist.reset()
	case ModeOP, ModeTranOP:
		l.vPrev, l.iPrev = 0, x[l.branch]
		l.hist.reset()
	default:
		op := l.Operating(x)
		l.vPrev, l.iPrev = op.Voltage, op.Current
	}
	l.hist.push(t, l.L*l.iPrev)
}

// SetInitialCondition loads the UIC state (IC value, or 0 if none was given).
func (l *Inductor) SetInitialCondition() {
	l.vPrev, l.iPrev = 0, 0
	if l.hasIC {
		l.iPrev = l.ic
	}
}

// StoredEnergy is ½Li².
func (l *Inductor) StoredEnergy(op OperatingPoint) float64 {
	return 0.5 * l.L * op.Current * op.Current
}

// Truncate estimates the step from the flux history; the "current" in the
// tolerance is the inductor voltage (the time derivative of flux).
func (l *Inductor) Truncate(x []float64, tc *TruncContext) float64 {
	op := l.Operating(x)
	return lteStep(&l.hist, l.L*op.Current, op.Voltage, l.vPrev, tc.VnTol, tc)
}

// ExteriorPower: an ideal inductor exchanges nothing with the exterior.
func (l *Inductor) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// LoadAC stamps the admittance jωC.
func (c *Capacitor) LoadAC(ctx *ACContext) error {
	stampY(ctx.M, c.u[0], c.u[1], complex(0, ctx.Omega*c.C))
	c.omega = ctx.Omega
	return nil
}

// OperatingAC returns the capacitor's phasors.
func (c *Capacitor) OperatingAC(x []complex128) ACPoint {
	v := c.voltageAC(x)
	return ACPoint{Voltage: v, Current: v * complex(0, c.omega*c.C)}
}

// LoadAC stamps the branch equation V(n+) − V(n−) − jωL·I = 0.
func (l *Inductor) LoadAC(ctx *ACContext) error {
	a, b, k := l.u[0], l.u[1], l.branch
	ctx.M.Add(a, k, 1)
	ctx.M.Add(b, k, -1)
	ctx.M.Add(k, a, 1)
	ctx.M.Add(k, b, -1)
	ctx.M.Add(k, k, complex(0, -ctx.Omega*l.L))
	return nil
}

// OperatingAC returns the inductor's phasors.
func (l *Inductor) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: l.voltageAC(x), Current: x[l.branch]}
}

// InitialConditioner devices accept UIC initial states.
type InitialConditioner interface {
	SetInitialCondition()
}
