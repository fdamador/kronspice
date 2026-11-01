// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"
	"sort"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Behavioral is a nonlinear element whose output is a piecewise-linear
// function y = f(c) of one controlling quantity c (ngspice B source with
// pwl). Output 'I' drives a current from n+ to n−; output 'V' forces
// V(n+) − V(n−) through a branch-current unknown. Control 'V' is a node
// voltage difference, control 'I' the branch current of a voltage source.
//
// Two of the four combinations are two-terminal nonlinear resistors
// (circuit.SelfControlled):
//
//   - Voltage-controlled, i = f(v) (N-type: tunnel diode). It lives entirely
//     in the node-voltage (Kron open-path) frame of MNA.
//   - Current-controlled, v = g(i) (S-type: arc, glow lamp). Its current is a
//     branch unknown, sensed by a 0 V source in series.
//
// Either may have negative-slope (negative differential resistance)
// regions; a circuit containing one can then have several DC operating
// points, and which one Newton finds depends on the starting point
// (Options.NodeSet).
//
// Maxwell derivation (CLAUDE.md §3c): a nonlinear constitutive relation
// 𝔎 = C(𝔈)·𝔈 for the conduction current in the quasi-static limit, as for
// the resistor row, with C now depending on the field.
//
// Ledger. A nonlinear resistor's absorbed power v·i leaves as heat through
// its thermal port when positive; when the curve enters the second or
// fourth quadrant (v·i < 0) the element delivers power, which is booked as
// model-supplied, as for a negative resistance. Any other combination is a
// controlled source: all of its power, of either sign, is model-supplied
// (CLAUDE.md §3c ledger table).
//
// Newton: the curve is linearized at the control value of the previous
// iterate; a change of the control value larger than four mean segment
// widths in one iteration is limited (the device then reports Limited), so
// Newton crosses a few segments at a time instead of jumping across the
// curve.
type Behavioral struct {
	Base
	out, ctrl byte
	self      bool
	cu        [2]int // control node unknowns (control 'V')
	ctrlName  string // control source (control 'I')
	ck        int    // control branch unknown (control 'I')
	branch    int    // own branch unknown (output 'V')
	xs, ys    []float64
	maxStep   float64

	// Linearization of the most recent Load.
	c0, y0, g float64
	limited   bool
	cPrev     float64 // control value used at the previous iteration
	cAcc      float64 // accepted control value (AC linearization)
}

// NewBehavioral builds a behavioral element; self marks a nonlinear
// resistor (see circuit.SelfControlled).
func NewBehavioral(e circuit.Element, self bool) *Behavioral {
	b := &Behavioral{Base: newBase(e), out: e.BOut, ctrl: e.BCtrl, self: self, ctrlName: e.CtrlSource,
		ck: -1, branch: -1}
	n := len(e.Curve) / 2
	b.xs, b.ys = make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		b.xs[i], b.ys[i] = e.Curve[2*i], e.Curve[2*i+1]
	}
	if n >= 2 {
		b.maxStep = 4 * (b.xs[n-1] - b.xs[0]) / float64(n-1)
	}
	b.cu = [2]int{e.CtrlNodes[0], e.CtrlNodes[1]} // node indices until Setup
	if self {
		_ = b.SetIntrinsic(exchange.Electrical, exchange.Both,
			"nonlinear resistor: power delivered in the 2nd/4th quadrant is model-supplied")
		_ = b.SetIntrinsic(exchange.Thermal, exchange.Out, "dissipation v·i > 0 to surroundings")
	} else {
		_ = b.SetIntrinsic(exchange.Electrical, exchange.Both,
			"idealized controlled source: its power, of either sign, is model-supplied")
	}
	return b
}

// Setup maps nodes and allocates the output branch.
func (b *Behavioral) Setup(ctx *SetupContext) error {
	b.setupNodes(ctx)
	if b.ctrl == 'V' {
		b.cu = [2]int{ctx.Unknown(b.cu[0]), ctx.Unknown(b.cu[1])}
	}
	if b.out == 'V' {
		b.branch = ctx.AllocBranch()
	}
	return nil
}

// Link finds the controlling source's branch.
func (b *Behavioral) Link(ctx *SetupContext) error {
	if b.ctrl != 'I' {
		return nil
	}
	k, ok := ctx.BranchOf(b.ctrlName)
	if !ok {
		return fmt.Errorf("%w: %s: controlling source %s not found", ErrUnsupported, b.name, b.ctrlName)
	}
	b.ck = k
	return nil
}

// eval returns f(c) and f'(c); outside the curve the end segments extend.
func (b *Behavioral) eval(c float64) (y, g float64) {
	n := len(b.xs)
	k := sort.SearchFloat64s(b.xs, c) // first index with xs[k] >= c
	switch {
	case k <= 0:
		k = 1
	case k >= n:
		k = n - 1
	}
	x0, x1, y0, y1 := b.xs[k-1], b.xs[k], b.ys[k-1], b.ys[k]
	g = (y1 - y0) / (x1 - x0)
	return y0 + g*(c-x0), g
}

func (b *Behavioral) control(x []float64) float64 {
	if b.ctrl == 'V' {
		return at(x, b.cu[0]) - at(x, b.cu[1])
	}
	return at(x, b.ck)
}

// Load linearizes y ≈ y0 + g·(c − c0) and stamps it.
func (b *Behavioral) Load(ctx *LoadContext) error {
	c := b.control(ctx.X)
	b.limited = false
	if ctx.InitJct {
		c = 0
	} else if d := c - b.cPrev; b.maxStep > 0 && math.Abs(d) > b.maxStep {
		c = b.cPrev + math.Copysign(b.maxStep, d)
		b.limited = true
	}
	b.cPrev = c
	b.c0 = c
	b.y0, b.g = b.eval(c)
	b.stamp(ctx, b.g, b.y0-b.g*c)
	return nil
}

// stamp adds y = g·c + j0 to the matrix.
func (b *Behavioral) stamp(ctx *LoadContext, g, j0 float64) {
	p, n := b.u[0], b.u[1]
	if b.out == 'I' {
		if b.ctrl == 'V' {
			ctx.M.Add(p, b.cu[0], g)
			ctx.M.Add(p, b.cu[1], -g)
			ctx.M.Add(n, b.cu[0], -g)
			ctx.M.Add(n, b.cu[1], g)
		} else {
			ctx.M.Add(p, b.ck, g)
			ctx.M.Add(n, b.ck, -g)
		}
		ctx.AddRHS(p, -j0)
		ctx.AddRHS(n, j0)
		return
	}
	k := b.branch
	ctx.M.Add(p, k, 1)
	ctx.M.Add(n, k, -1)
	ctx.M.Add(k, p, 1)
	ctx.M.Add(k, n, -1)
	if b.ctrl == 'V' {
		ctx.M.Add(k, b.cu[0], -g)
		ctx.M.Add(k, b.cu[1], g)
	} else {
		ctx.M.Add(k, b.ck, -g)
	}
	ctx.AddRHS(k, j0)
}

// Limited reports whether the last Load limited the control step.
func (b *Behavioral) Limited() bool { return b.limited }

// ConvTest checks that the linearized output at x matches the curve.
func (b *Behavioral) ConvTest(x []float64, reltol, abstol float64) bool {
	c := b.control(x)
	yhat := b.y0 + b.g*(c-b.c0)
	y, _ := b.eval(c)
	tol := abstol
	if b.out == 'V' {
		tol = 1e-6 // vntol for a voltage output
	}
	return math.Abs(yhat-y) <= reltol*math.Max(math.Abs(y), math.Abs(yhat))+tol
}

// Operating reports the element's voltage and current at the solution.
func (b *Behavioral) Operating(x []float64) OperatingPoint {
	v := b.voltage(x)
	var i float64
	if b.out == 'I' {
		i = b.y0 + b.g*(b.control(x)-b.c0) // the linearized current the solve used
	} else {
		i = at(x, b.branch)
	}
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// Accept records the converged control value.
func (b *Behavioral) Accept(x []float64, _ float64, _ Mode) {
	b.cAcc = b.control(x)
	b.cPrev = b.cAcc
}

// ExteriorPower: a nonlinear resistor's dissipation leaves as heat.
func (b *Behavioral) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal && b.self && op.PowerAbsorbed > 0 {
		return -op.PowerAbsorbed
	}
	return 0
}

// ModelSuppliedPower is the power delivered without a declared port: the
// 2nd/4th-quadrant power of a nonlinear resistor, or all power of a
// controlled source (negative when it absorbs).
func (b *Behavioral) ModelSuppliedPower(op OperatingPoint) float64 {
	if b.self {
		if op.PowerAbsorbed < 0 {
			return -op.PowerAbsorbed
		}
		return 0
	}
	return -op.PowerAbsorbed
}

// LoadAC stamps the small-signal slope at the operating point.
func (b *Behavioral) LoadAC(ctx *ACContext) error {
	_, g := b.eval(b.cAcc)
	p, n := b.u[0], b.u[1]
	cg := complex(g, 0)
	if b.out == 'I' {
		if b.ctrl == 'V' {
			ctx.M.Add(p, b.cu[0], cg)
			ctx.M.Add(p, b.cu[1], -cg)
			ctx.M.Add(n, b.cu[0], -cg)
			ctx.M.Add(n, b.cu[1], cg)
		} else {
			ctx.M.Add(p, b.ck, cg)
			ctx.M.Add(n, b.ck, -cg)
		}
		return nil
	}
	k := b.branch
	ctx.M.Add(p, k, 1)
	ctx.M.Add(n, k, -1)
	ctx.M.Add(k, p, 1)
	ctx.M.Add(k, n, -1)
	if b.ctrl == 'V' {
		ctx.M.Add(k, b.cu[0], -cg)
		ctx.M.Add(k, b.cu[1], cg)
	} else {
		ctx.M.Add(k, b.ck, -cg)
	}
	return nil
}

// OperatingAC returns the element's phasors.
func (b *Behavioral) OperatingAC(x []complex128) ACPoint {
	v := b.voltageAC(x)
	if b.out == 'V' {
		return ACPoint{Voltage: v, Current: catC(x, b.branch)}
	}
	_, g := b.eval(b.cAcc)
	var c complex128
	if b.ctrl == 'V' {
		c = catC(x, b.cu[0]) - catC(x, b.cu[1])
	} else {
		c = catC(x, b.ck)
	}
	return ACPoint{Voltage: v, Current: complex(g, 0) * c}
}

func catC(x []complex128, i int) complex128 {
	if i < 0 || i >= len(x) {
		return 0
	}
	return x[i]
}
