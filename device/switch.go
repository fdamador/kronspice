// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// StepLimiter is implemented by devices without stored energy that still
// need the transient step cut, such as a switch that should change state
// close to the instant its control crosses the threshold. The analysis
// treats a proposal below the step like a local-truncation-error
// rejection.
type StepLimiter interface {
	LimitStep(x []float64, tc *TruncContext) float64
	// Restart reports whether the solution just accepted changed the
	// device discontinuously (a switch changing state), so integration
	// restarts with backward Euler and a small step, as after a source
	// breakpoint.
	Restart() bool
}

// Switch is a voltage- (SPICE S) or current-controlled (SPICE W) switch,
// following the ngspice switch model (ngspice manual, "Switches"; Berkeley
// SPICE3 sw/csw devices).
//
// With hysteresis VH ≥ 0 the switch is a resistance RON or ROFF. It turns on
// when the control c exceeds VT + VH, off when c falls below VT − VH, and
// otherwise keeps its state. The state is part of the solution history:
// the accepted state of the previous solution (or the ON/OFF initial state
// for the first operating point) is the reference for the next.
//
//   - Operating point and DC sweep: the state is the one the iterate's
//     control selects from that reference. A Newton iterate that would
//     change the state is not converged, so the solve repeats with the
//     other resistance (SPICE3 SWload); a circuit with no consistent state
//     (a relaxation oscillator at DC) returns no convergence.
//   - Transient step: the state is held at its accepted value for the whole
//     step, so the step is a linear solve. If the control crossed a
//     threshold during the step, LimitStep has the step repeated, ending
//     just after the crossing found by linear interpolation of the
//     control, and Accept then changes the state for the next step. The
//     switching instant is resolved to a hundredth of the step that
//     crossed it, and a switch that feeds its own control (a relaxation
//     oscillator) cannot flip back and forth inside one step.
//
// With VH < 0 the switch is smooth: between VT − |VH| and VT + |VH|, ln R
// follows the cubic in x = (c − VT)/|VH| that runs from ln ROFF to ln RON
// with zero slope at both ends,
//
//	ln R = ln √(RON·ROFF) + ln(RON/ROFF)·(3x − x³)/4,
//
// which is ngspice's smooth switch (the only cubic with those end values and
// slopes). It is a nonlinear resistor in c and is linearized for Newton.
//
// Maxwell derivation (CLAUDE.md §3c table, row "Resistor"): an ohmic
// conductor whose conductivity is set externally. Ports: the dissipated
// v·i leaves through the thermal port, as for a resistor. The control is
// ideal (no current into its nodes, no power drawn), so the switch is
// passive: it supplies no model power. Changing state does no work, since
// the switch stores no energy.
type Switch struct {
	Base
	current        bool   // W: controlled by a source current
	cu             [2]int // control node unknowns (S)
	ctrlName       string // control source (W)
	ck             int    // control branch unknown (W)
	vt, vh         float64
	ron, roff      float64
	smooth         bool
	lmean, lratio  float64 // ln √(RON·ROFF), ln(RON/ROFF) (smooth)
	initOn         bool
	on, accOn      bool    // state stamped by the last Load; accepted state
	accepted       bool    // a solution has been accepted
	frozen         bool    // the last Load held the accepted state (transient step)
	flipped        bool    // the last Accept changed the state
	g, k, c0, cAcc float64 // linearization: i ≈ g·v + k·(c − c0)
	v0             float64
}

// NewSwitch builds a switch from its element and "sw" or "csw" model.
func NewSwitch(e circuit.Element, m circuit.Model) *Switch {
	s := &Switch{Base: newBase(e), current: e.Kind == circuit.CurrentSwitch, ctrlName: e.CtrlSource, ck: -1}
	s.vt, s.vh, s.ron, s.roff = circuit.SwitchParams(m)
	if s.vh < 0 {
		s.smooth = true
		s.vh = -s.vh
		s.lmean = math.Log(math.Sqrt(s.ron * s.roff))
		s.lratio = math.Log(s.ron / s.roff)
	}
	s.initOn = e.HasIC && e.IC == 1
	s.on, s.accOn = s.initOn, s.initOn
	s.cu = e.CtrlNodes // node indices until Setup
	_ = s.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network")
	_ = s.SetIntrinsic(exchange.Thermal, exchange.Out, "Joule heating i²R to surroundings")
	return s
}

// Setup maps nodes.
func (s *Switch) Setup(ctx *SetupContext) error {
	s.setupNodes(ctx)
	if !s.current {
		s.cu = [2]int{ctx.Unknown(s.cu[0]), ctx.Unknown(s.cu[1])}
	}
	return nil
}

// Link finds the controlling source's branch (W).
func (s *Switch) Link(ctx *SetupContext) error {
	if !s.current {
		return nil
	}
	k, ok := ctx.BranchOf(s.ctrlName)
	if !ok {
		return fmt.Errorf("%w: %s: controlling source %s not found", ErrUnsupported, s.name, s.ctrlName)
	}
	s.ck = k
	return nil
}

func (s *Switch) control(x []float64) float64 {
	if s.current {
		return at(x, s.ck)
	}
	return at(x, s.cu[0]) - at(x, s.cu[1])
}

// next is the hysteretic state after control c, from state prev.
func (s *Switch) next(c float64, prev bool) bool {
	switch {
	case c > s.vt+s.vh:
		return true
	case c < s.vt-s.vh:
		return false
	}
	return prev
}

// smoothG returns the smooth switch's conductance and dg/dc at control c.
func (s *Switch) smoothG(c float64) (g, dg float64) {
	x := (c - s.vt) / s.vh
	switch {
	case x >= 1:
		return 1 / s.ron, 0
	case x <= -1:
		return 1 / s.roff, 0
	}
	g = math.Exp(-(s.lmean + s.lratio*(3*x-x*x*x)/4))
	return g, -g * s.lratio * 3 * (1 - x*x) / (4 * s.vh)
}

func (s *Switch) cond(on bool) float64 {
	if on {
		return 1 / s.ron
	}
	return 1 / s.roff
}

// Load stamps the conductance of the state the control selects (abrupt
// switch) or the linearized smooth conductance.
func (s *Switch) Load(ctx *LoadContext) error {
	a, b := s.u[0], s.u[1]
	if s.smooth {
		c, v := s.control(ctx.X), s.vab(ctx.X)
		if ctx.InitJct {
			c, v = s.vt-s.vh, 0
			if s.initOn {
				c = s.vt + s.vh
			}
		}
		g, dg := s.smoothG(c)
		s.g, s.k, s.c0, s.v0 = g, dg*v, c, v
		ctx.M.Add(a, a, g)
		ctx.M.Add(b, b, g)
		ctx.M.Add(a, b, -g)
		ctx.M.Add(b, a, -g)
		if s.k != 0 {
			s.stampCtrl(ctx.M, s.k)
			j0 := -s.k * c // i = g·v + k·c + j0
			ctx.AddRHS(a, -j0)
			ctx.AddRHS(b, j0)
		}
		return nil
	}
	s.frozen = ctx.Mode == ModeTran
	if ctx.InitJct || s.frozen {
		s.on = s.accOn
	} else {
		s.on = s.next(s.control(ctx.X), s.accOn)
	}
	s.g, s.k = s.cond(s.on), 0
	ctx.M.Add(a, a, s.g)
	ctx.M.Add(b, b, s.g)
	ctx.M.Add(a, b, -s.g)
	ctx.M.Add(b, a, -s.g)
	return nil
}

// stampCtrl adds the transconductance k from the control to the switch
// current (a → b).
func (s *Switch) stampCtrl(m Stamper, k float64) {
	a, b := s.u[0], s.u[1]
	if s.current {
		m.Add(a, s.ck, k)
		m.Add(b, s.ck, -k)
		return
	}
	m.Add(a, s.cu[0], k)
	m.Add(a, s.cu[1], -k)
	m.Add(b, s.cu[0], -k)
	m.Add(b, s.cu[1], k)
}

// Limited is always false: the switch does not limit its control.
func (s *Switch) Limited() bool { return false }

// ConvTest: an abrupt switch has converged when the iterate keeps the
// state it was stamped with (always, in a transient step, where the state
// is held); a smooth one when the linearized current matches the curve.
func (s *Switch) ConvTest(x []float64, reltol, abstol float64) bool {
	c := s.control(x)
	if !s.smooth {
		return s.frozen || s.next(c, s.accOn) == s.on
	}
	v := s.vab(x)
	g, _ := s.smoothG(c)
	ihat := s.g*v + s.k*(c-s.c0)
	i := g * v
	return math.Abs(ihat-i) <= reltol*math.Max(math.Abs(i), math.Abs(ihat))+abstol
}

// Operating reports voltage, the current the solve used, and power.
func (s *Switch) Operating(x []float64) OperatingPoint {
	v := s.vab(x)
	i := s.g*v + s.k*(s.control(x)-s.c0)
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// Accept commits the state of the accepted solution.
func (s *Switch) Accept(x []float64, _ float64, _ Mode) {
	c := s.control(x)
	s.cAcc, s.accepted = c, true
	if !s.smooth {
		on := s.next(c, s.accOn)
		s.flipped = on != s.accOn
		s.accOn, s.on = on, on
	}
}

// Restart reports whether the last accepted solution changed the state.
func (s *Switch) Restart() bool { return s.flipped }

// LimitStep cuts a transient step in which an abrupt switch changed state,
// unless the crossing lies in the last hundredth of the step.
func (s *Switch) LimitStep(x []float64, tc *TruncContext) float64 {
	if s.smooth || !s.accepted {
		return inf
	}
	c := s.control(x)
	if s.next(c, s.accOn) == s.accOn || c == s.cAcc {
		return inf
	}
	th := s.vt - s.vh
	if !s.accOn {
		th = s.vt + s.vh
	}
	frac := (th - s.cAcc) / (c - s.cAcc)
	if !(frac > 0) || frac >= 0.99 {
		return inf
	}
	return tc.H * math.Min(frac+0.005, 1)
}

// ExteriorPower: the dissipated power leaves as heat.
func (s *Switch) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -op.PowerAbsorbed
	}
	return 0
}

// LoadAC stamps the small-signal conductance at the operating point (and,
// for a smooth switch, the transconductance from its control).
func (s *Switch) LoadAC(ctx *ACContext) error {
	stampY(ctx.M, s.u[0], s.u[1], complex(s.g, 0))
	if s.smooth && s.k != 0 {
		a, b, k := s.u[0], s.u[1], complex(s.k, 0)
		if s.current {
			ctx.M.Add(a, s.ck, k)
			ctx.M.Add(b, s.ck, -k)
		} else {
			ctx.M.Add(a, s.cu[0], k)
			ctx.M.Add(a, s.cu[1], -k)
			ctx.M.Add(b, s.cu[0], -k)
			ctx.M.Add(b, s.cu[1], k)
		}
	}
	return nil
}

// OperatingAC returns the switch's phasors.
func (s *Switch) OperatingAC(x []complex128) ACPoint {
	v := s.voltageAC(x)
	i := v * complex(s.g, 0)
	if s.smooth && s.k != 0 {
		var c complex128
		if s.current {
			c = catC(x, s.ck)
		} else {
			c = catC(x, s.cu[0]) - catC(x, s.cu[1])
		}
		i += complex(s.k, 0) * c
	}
	return ACPoint{Voltage: v, Current: i}
}

// On reports the accepted state of an abrupt switch (for a smooth one,
// whether the accepted control is past the threshold).
func (s *Switch) On() bool {
	if s.smooth {
		return s.cAcc > s.vt
	}
	return s.accOn
}

// vab is V(n+) − V(n−), reading 0 for unknowns outside x.
func (s *Switch) vab(x []float64) float64 { return at(x, s.u[0]) - at(x, s.u[1]) }
