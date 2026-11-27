// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Arc is the kronspice A element: a dynamic arc or spark gap whose
// conductance g follows the Mayr, Cassie or Habedank law (circuit.Arc).
//
// Unknowns: the two nodes and the arc current I (n+ → n−). The branch
// equation is g(I)·v − I = 0, written in conductance form so that a cold gap
// (g = GOFF) stays well conditioned.
//
// Time discretization. With the current-driven laws the state is linear in
// g (Mayr) or g² (Cassie) for a fixed current, so one step of length h from
// the accepted state is integrated exactly for a current held at the
// trapezoidal mean of i² over the step, q = (i₀² + i₁²)/2:
//
//	Mayr:    g₁  = q/P0  + (g₀  − q/P0)·e^(−h/τ)
//	Cassie:  g₁² = q/U0² + (g₀² − q/U0²)·e^(−2h/τ)
//
// (Mayr 1943; Cassie 1939; current-driven form as in Habedank 1993). This
// is unconditionally stable and exact for constant current. The step is
// implicit in i₁, so a fast discharge of a capacitor through the arc stays
// stable at steps far above the RC time; the Newton system folds over only
// when h exceeds the arc's growth time τ·P0/(v·i), which Truncate prevents. At an
// operating point the conductance is frozen at its state (G0 initially):
// the arc is a resistor there, which makes the DC solution unique.
//
// Maxwell derivation (CLAUDE.md §3c): Ohm's law 𝔎 = C𝔈 for the arc column,
// whose conductivity is set by the energy balance of the column (Mayr:
// constant cooling power; Cassie: constant current density, so a constant
// field). The element is passive, v·i = g·v² ≥ 0; the power leaves as heat
// and light through the thermal port (not split between them).
type Arc struct {
	Base
	p  circuit.ArcParams
	br int

	gm, gc float64 // accepted Mayr and Cassie conductances [S]
	i0, t0 float64 // accepted current and time

	broke int  // breakdowns so far
	fresh bool // a breakdown fired and no short step has followed yet

	tran              bool
	lastH             float64
	linG, linDG, linI float64
}

// NewArc builds an arc from a circuit element.
func NewArc(e circuit.Element) *Arc {
	a := &Arc{Base: newBase(e), p: circuit.ArcParamsOf(e)}
	a.setState(a.p.G0)
	_ = a.SetIntrinsic(exchange.Electrical, exchange.Both, "arc column between the terminals")
	_ = a.SetIntrinsic(exchange.Thermal, exchange.Out, "arc power g·v² as heat and light (not split)")
	return a
}

// setState sets the total conductance to g (each series part 2g).
func (a *Arc) setState(g float64) {
	switch a.p.Model {
	case circuit.ArcMayr:
		a.gm = g
	case circuit.ArcCassie:
		a.gc = g
	default:
		a.gm, a.gc = 2*g, 2*g
	}
}

// total returns the series conductance of the parts and its derivatives.
func (a *Arc) total(gm, dgm, gc, dgc float64) (float64, float64) {
	switch a.p.Model {
	case circuit.ArcMayr:
		return gm, dgm
	case circuit.ArcCassie:
		return gc, dgc
	}
	s := gc + gm
	return gc * gm / s, (gm*gm*dgc + gc*gc*dgm) / (s * s)
}

func (a *Arc) clamp(g, dg float64) (float64, float64) {
	if g < a.p.Goff {
		return a.p.Goff, 0
	}
	if g > a.p.Gmax {
		return a.p.Gmax, 0
	}
	return g, dg
}

// step integrates the state from the accepted point over h with the
// driving term q = i² taken as given. It returns the new part conductances
// and the total.
func (a *Arc) step(h, q float64) (gm, gc, g float64) {
	gm, gc, g, _ = a.stepD(h, q)
	return
}

// stepD is step with dg/dq.
func (a *Arc) stepD(h, q float64) (gm, gc, g, dg float64) {
	var dgm, dgc float64
	if a.p.Model != circuit.ArcCassie {
		tau := a.p.Tau
		if a.p.Model == circuit.ArcHabedank {
			tau = a.p.TauM
		}
		e := math.Exp(-h / tau)
		gm, dgm = a.clamp(q/a.p.P0+(a.gm-q/a.p.P0)*e, (1-e)/a.p.P0)
	}
	if a.p.Model != circuit.ArcMayr {
		tau := a.p.Tau
		if a.p.Model == circuit.ArcHabedank {
			tau = a.p.TauC
		}
		e := math.Exp(-2 * h / tau)
		u2 := a.p.U0 * a.p.U0
		gc = math.Sqrt(math.Max(0, q/u2+(a.gc*a.gc-q/u2)*e))
		if gc > 0 {
			dgc = (1 - e) / u2 / (2 * gc)
		}
		gc, dgc = a.clamp(gc, dgc)
	}
	g, dg = a.total(gm, dgm, gc, dgc)
	return gm, gc, g, dg
}

// Conductance returns the accepted total conductance [S].
func (a *Arc) Conductance() float64 {
	g, _ := a.total(a.gm, 0, a.gc, 0)
	return g
}

// Breakdowns returns the number of gap breakdowns so far.
func (a *Arc) Breakdowns() int { return a.broke }

// Setup maps the nodes and allocates the arc current.
func (a *Arc) Setup(ctx *SetupContext) error {
	a.setupNodes(ctx)
	a.br = ctx.AllocBranch()
	return nil
}

// Load stamps g·v − I = 0 with g integrated over the step, and the arc
// current into KCL.
func (a *Arc) Load(ctx *LoadContext) error {
	ua, ub, br := a.u[0], a.u[1], a.br
	g, dgdi := a.Conductance(), 0.0
	v0, i := at(ctx.X, ua)-at(ctx.X, ub), at(ctx.X, br)
	if ctx.Mode == ModeTran {
		a.tran, a.lastH = true, ctx.H
		var dq float64
		_, _, g, dq = a.stepD(ctx.H, a.q(i))
		dgdi = dq * i // dq/di = i
	} else {
		a.tran = false
	}
	a.linG, a.linDG, a.linI = g, dgdi, i
	m := ctx.M
	m.Add(ua, br, 1)
	m.Add(ub, br, -1)
	// g(I)·v − I ≈ g·v + g′·v0·(I − I0) − I = 0
	m.Add(br, ua, g)
	m.Add(br, ub, -g)
	m.Add(br, br, dgdi*v0-1)
	ctx.AddRHS(br, dgdi*v0*i)
	return nil
}

// Limited reports no junction limiting.
func (a *Arc) Limited() bool { return false }

// q is the driving term i² averaged over the step (trapezoidal).
func (a *Arc) q(i1 float64) float64 { return 0.5 * (a.i0*a.i0 + i1*i1) }

// ConvTest checks the branch equation g(I)·v = I at x.
func (a *Arc) ConvTest(x []float64, reltol, abstol float64) bool {
	if !a.tran {
		return true
	}
	v, i := a.voltage(x), at(x, a.br)
	_, _, g := a.step(a.lastH, a.q(i))
	return math.Abs(g*v-i) <= reltol*math.Max(math.Abs(g*v), math.Abs(i))+abstol
}

// Operating reports the arc voltage, current and absorbed power.
func (a *Arc) Operating(x []float64) OperatingPoint {
	v, i := a.voltage(x), at(x, a.br)
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// Accept commits the state at an accepted point and fires a breakdown.
func (a *Arc) Accept(x []float64, t float64, mode Mode) {
	v, i := a.voltage(x), at(x, a.br)
	if a.fresh && mode == ModeTran && t-a.t0 <= 1.01*a.restart() {
		a.fresh = false
	}
	if mode == ModeTran && t > a.t0 {
		a.gm, a.gc, _ = a.step(t-a.t0, a.q(i))
	}
	a.i0, a.t0 = i, t
	if (mode == ModeTran || mode == ModeTranOP) && a.p.Vbd > 0 && math.Abs(v) >= a.p.Vbd && a.Conductance() < a.p.Gon {
		a.setState(a.p.Gon)
		a.broke++
		a.fresh = true
	}
}

// ExteriorPower: the arc's power leaves through the thermal port.
func (a *Arc) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -op.PowerAbsorbed
	}
	return 0
}

// StoredEnergy is zero: the arc column's heat content is not modeled as
// circuit storage (its thermal inertia is the time constant τ).
func (a *Arc) StoredEnergy(OperatingPoint) float64 { return 0 }

// Truncate limits the step so that ln g changes by at most arcLogStep,
// which resolves the arc's own time constant, and restarts with a short
// step after a breakdown.
func (a *Arc) Truncate(x []float64, tc *TruncContext) float64 {
	h := inf
	if a.fresh {
		h = a.restart()
	}
	_, _, g := a.step(tc.H, a.q(at(x, a.br)))
	if d := math.Abs(math.Log(g / a.Conductance())); d > 1e-9 {
		h = math.Min(h, tc.H*arcLogStep/d)
	}
	return h
}

// restart is the step taken after a breakdown, which changes the
// circuit's state discontinuously: a thousandth of the shortest time
// constant.
func (a *Arc) restart() float64 {
	tau := a.p.Tau
	if a.p.Model == circuit.ArcHabedank {
		tau = math.Min(a.p.TauC, a.p.TauM)
	}
	return 1e-3 * tau
}

const arcLogStep = 0.25 // largest change of ln g in one step

// LoadAC stamps the arc as its accepted conductance (the state is frozen
// for small signals; the Mayr and Cassie small-signal impedances, with
// their τ-dependent phase, are not modeled).
func (a *Arc) LoadAC(ctx *ACContext) error {
	ua, ub, br := a.u[0], a.u[1], a.br
	g := complex(a.Conductance(), 0)
	ctx.M.Add(ua, br, 1)
	ctx.M.Add(ub, br, -1)
	ctx.M.Add(br, ua, g)
	ctx.M.Add(br, ub, -g)
	ctx.M.Add(br, br, -1)
	return nil
}

// OperatingAC returns the arc's phasors.
func (a *Arc) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: a.voltageAC(x), Current: x[a.br]}
}
