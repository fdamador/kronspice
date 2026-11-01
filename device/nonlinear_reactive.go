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

// stateLaw is the constitutive law of a charge- or flux-defined element,
//
//	y(x, t) = s(t)·f(x)
//
// with x the voltage and y the charge (capacitor), or x the current and y
// the flux linkage (inductor). f is a piecewise-linear curve through the
// origin with positive slope, or lin·x; s is a dimensionless waveform > 0
// (1 when the element does not vary).
//
// Energy (Woodson & Melcher, Electromechanical Dynamics, Part I, ch. 3:
// energy and co-energy of a lossless coupling field). At fixed s the energy
// stored is W = ∫₀^y x dy′ = s·(x·f(x) − F(x)), with F(x) = ∫₀^x f(u) du the
// co-energy at s = 1. The electrical power into the element splits as
//
//	x·dy/dt = dW/dt + ṡ·F(x)
//
// so ṡ·F(x) leaves through whatever varies s — the mechanical force on
// moving plates or a moving core does that work. For a linear capacitor
// this is the familiar v·i = d(½Cv²)/dt + ½v²·dC/dt.
type stateLaw struct {
	lin     float64   // f(x) = lin·x when there is no curve
	xs, ys  []float64 // curve knots
	pre     []float64 // ∫ from xs[0] to xs[k] of f
	g0      float64   // ∫ from xs[0] to 0 of f
	maxStep float64   // Newton step limit inside the curve's knot range
	knotT   float64   // time at which a cut step lands on a knot (0: none)

	spec   circuit.Waveform
	scale  wave
	varies bool
}

func newStateLaw(e circuit.Element) stateLaw {
	l := stateLaw{lin: e.Value, spec: e.Scale, varies: e.Scale.Kind != circuit.WaveNone}
	n := len(e.Curve) / 2
	if n < 2 {
		return l
	}
	l.xs, l.ys, l.pre = make([]float64, n), make([]float64, n), make([]float64, n)
	for i := 0; i < n; i++ {
		l.xs[i], l.ys[i] = e.Curve[2*i], e.Curve[2*i+1]
	}
	for k := 1; k < n; k++ {
		l.pre[k] = l.pre[k-1] + 0.5*(l.ys[k-1]+l.ys[k])*(l.xs[k]-l.xs[k-1])
	}
	l.g0 = l.cum(0)
	l.maxStep = 4 * (l.xs[n-1] - l.xs[0]) / float64(n-1)
	return l
}

// resolve fills the scale waveform's SPICE defaults for one analysis and
// checks that s stays > 0 with them (an EXP bound depends on its defaulted
// time constants, which circuit validation does not know).
func (l *stateLaw) resolve(name string, step, stop float64) error {
	if !l.varies {
		return nil
	}
	l.scale = resolveWave(l.spec, step, stop)
	if l.scale.kind == circuit.WaveExp {
		// v = v1 + (v2 − v1)(e2 − e1) after td2, with e1 − e2 below
		// e^(−(td2 − td1)/τ1); before td2, between v1 and v2.
		a := &l.scale.a
		lo := math.Min(a[0], a[1])
		if a[1] > a[0] {
			lo = math.Min(lo, a[0]-(a[1]-a[0])*math.Exp(-(a[4]-a[2])/a[3]))
		}
		if !(lo > 0) {
			return fmt.Errorf("%w: %s scale EXP can reach %g; it must stay > 0", circuit.ErrBadWaveform, name, lo)
		}
	}
	return nil
}

// segment returns k such that knots k−1, k bound x (end segments extend).
func (l *stateLaw) segment(x float64) int {
	n := len(l.xs)
	k := sort.SearchFloat64s(l.xs, x)
	switch {
	case k <= 0:
		return 1
	case k >= n:
		return n - 1
	}
	return k
}

// f returns f(x) and f′(x).
func (l *stateLaw) f(x float64) (y, dy float64) {
	if l.xs == nil {
		return l.lin * x, l.lin
	}
	k := l.segment(x)
	g := (l.ys[k] - l.ys[k-1]) / (l.xs[k] - l.xs[k-1])
	return l.ys[k-1] + g*(x-l.xs[k-1]), g
}

// cum is ∫ from xs[0] to x of f, exact for the piecewise-linear curve.
func (l *stateLaw) cum(x float64) float64 {
	k := l.segment(x)
	g := (l.ys[k] - l.ys[k-1]) / (l.xs[k] - l.xs[k-1])
	d := x - l.xs[k-1]
	return l.pre[k-1] + l.ys[k-1]*d + 0.5*g*d*d
}

// coenergy is F(x) = ∫₀^x f(u) du.
func (l *stateLaw) coenergy(x float64) float64 {
	if l.xs == nil {
		return 0.5 * l.lin * x * x
	}
	return l.cum(x) - l.g0
}

// energy is x·f(x) − F(x) = ∫₀^f(x) x dy, the stored energy at s = 1.
func (l *stateLaw) energy(x float64) float64 {
	y, _ := l.f(x)
	return x*y - l.coenergy(x)
}

func (l *stateLaw) s(t float64) float64 {
	if !l.varies {
		return 1
	}
	return l.scale.value(t)
}

func (l *stateLaw) sdot(t float64) float64 {
	if !l.varies {
		return 0
	}
	return l.scale.slope(t)
}

// limit bounds a Newton update of x inside the curve's knot range to a few
// segment widths, as for the behavioral element; updates entirely beyond
// one end of the curve follow a straight line and are not limited.
func (l *stateLaw) limit(x, prev float64) (float64, bool) {
	if l.maxStep == 0 || math.Abs(x-prev) <= l.maxStep {
		return x, false
	}
	lo, hi := l.xs[0], l.xs[len(l.xs)-1]
	if (x <= lo && prev <= lo) || (x >= hi && prev >= hi) {
		return x, false
	}
	return prev + math.Copysign(l.maxStep, x-prev), true
}

// knotStep proposes a step for a candidate that moved x from the accepted
// x0 to x1 across a knot of the curve. The state y = s·f(x) is smooth there,
// so its truncation error does not see the corner in x (a saturating core's
// current turns sharply while its flux does not, and the new segment's time
// constant can be far shorter). A step that crosses a knot before 90% of the
// way is cut to end at the knot, by linear interpolation, so the corner
// falls on a time point. Interpolation
// undershoots when x moves faster beyond the knot, so a value within 10⁻³ of
// the neighbouring segment's width counts as on the knot; otherwise the
// steps would shrink geometrically towards it without ever crossing.
func (l *stateLaw) knotStep(x0, x1, h float64) float64 {
	if l.xs == nil || x0 == x1 {
		return inf
	}
	n := len(l.xs)
	k := sort.SearchFloat64s(l.xs, x0) // first knot ≥ x0
	var xk float64
	if x1 > x0 {
		for k < n && l.xs[k] <= x0 {
			k++
		}
		if k == n || l.xs[k] >= x1 {
			return inf
		}
		xk = l.xs[k]
	} else {
		if k--; k < 0 || l.xs[k] <= x1 {
			return inf
		}
		xk = l.xs[k]
	}
	w := inf
	if k > 0 {
		w = l.xs[k] - l.xs[k-1]
	}
	if k < n-1 {
		w = math.Min(w, l.xs[k+1]-l.xs[k])
	}
	if math.Abs(xk-x0) <= 1e-3*w {
		return inf
	}
	f := (xk - x0) / (x1 - x0)
	if f > 0.9 {
		return inf // nearly there: the step is accepted as it is
	}
	return f * h
}

// truncate applies knotStep to a candidate step. A cut also registers the
// landing time as a breakpoint, so the engine lands on it exactly and
// restarts with backward Euler: across the corner the element's voltage
// (capacitor current) jumps, and the trapezoidal rule would carry an
// alternating error past it, as after a source corner.
func (l *stateLaw) truncate(x0, x1 float64, tc *TruncContext) float64 {
	p := l.knotStep(x0, x1, tc.H)
	if p < inf {
		l.knotT = tc.Time - tc.H + p
	}
	return p
}

// accepted clears a knot breakpoint once the run has reached it.
func (l *stateLaw) accepted(t float64) {
	if l.knotT != 0 && t >= l.knotT {
		l.knotT = 0
	}
}

func (l *stateLaw) nextBreakpoint(t float64) float64 {
	b := inf
	if l.varies {
		b = l.scale.nextBreakpoint(t)
	}
	if l.knotT > t && l.knotT < b {
		b = l.knotT
	}
	return b
}

// scaleChannel is where the power of a varying law goes: the element's
// declared channel, or the mechanical exterior by default.
func scaleChannel(e circuit.Element) exchange.ChannelID {
	if e.Channel == exchange.Electrical {
		return exchange.Mechanical
	}
	return e.Channel
}

// NonlinearCapacitor is a capacitor defined by its charge,
// q(v, t) = s(t)·Q(v) (circuit.NonlinearCapacitor). Q is a piecewise-linear
// curve or C·v; s is the scale waveform.
//
// Maxwell derivation (CLAUDE.md §3c table, row "Capacitor"): the plate
// charge is the flux of 𝔇, and charge continuity (the scalar part of
// ∇ℭ = 0 with ℭ ⊃ ∂𝔇/∂t) gives the terminal current i = dq/dt in the
// quasi-static limit. With q = s(t)·Q(v),
//
//	i = s·Q′(v)·dv/dt + ṡ·Q(v)
//
// which for Q = C·v is i = C dv/dt + v dC/dt. Writing the element in terms
// of charge, not i = C(v)·dv/dt, conserves charge exactly; the companion
// models integrate the charge, as SPICE3 does for junction charges (Nagel
// 1975 ch. 4; Kundert, The Designer's Guide to SPICE and Spectre, ch. 2):
//
//	backward Euler: i_n = (q_n − q_{n−1})/h
//	trapezoidal:    i_n = 2(q_n − q_{n−1})/h − i_{n−1}
//
// with q_n = s(t_n)·Q(v_n), linearized in v_n for Newton when Q is a curve.
//
// Ledger. Stored energy is s·∫₀^v u·Q′(u) du (stateLaw). When s varies, the
// power ṡ·∫₀^v Q(u) du is exchanged through the scale's port (mechanical by
// default: positive out of the circuit when the capacitance grows at fixed
// voltage, as the plates' attraction does work while they close), so energy
// pumped in by moving plates enters as declared inflow and never appears
// from nowhere.
type NonlinearCapacitor struct {
	Base
	law   stateLaw
	ch    exchange.ChannelID
	ic    float64
	hasIC bool

	// Linearization of the most recent Load: i ≈ iq + geq·(v − v0), and the
	// exact companion i(v) = a·sn·Q(v) + b.
	v0, iq, geq float64
	a, b, sn    float64
	sdot        float64
	t           float64
	tran, ac    bool
	limited     bool
	vIter       float64 // voltage used at the previous iteration (limiting)
	omega       float64

	vPrev, qPrev, iPrev float64 // accepted state
	vAcc, sAcc          float64 // accepted voltage and scale (AC linearization)
	hist                history
}

// curvedCapacitor is a NonlinearCapacitor with a charge curve, which needs
// Newton iteration; with Q = C·v each step is linear and needs one solve.
type curvedCapacitor struct{ NonlinearCapacitor }

// NewNonlinearCapacitor builds a nonlinear or time-varying capacitor.
func NewNonlinearCapacitor(e circuit.Element) Device {
	c := NonlinearCapacitor{Base: newBase(e), law: newStateLaw(e), ch: scaleChannel(e), ic: e.IC, hasIC: e.HasIC}
	_ = c.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores ∫v dq")
	if c.law.varies {
		_ = c.SetIntrinsic(c.ch, exchange.Both, "work of the varying capacitance, ṡ·∫₀ᵛ Q(u) du")
	}
	if c.law.xs != nil {
		return &curvedCapacitor{c}
	}
	return &c
}

// Setup maps the nodes and resolves the scale waveform.
func (c *NonlinearCapacitor) Setup(ctx *SetupContext) error {
	c.setupNodes(ctx)
	return c.law.resolve(c.name, ctx.Step, ctx.Stop)
}

// Load stamps an open circuit (DC), the UIC pin, or the charge companion.
func (c *NonlinearCapacitor) Load(ctx *LoadContext) error {
	p, n := c.u[0], c.u[1]
	c.ac, c.limited, c.t = false, false, ctx.Time
	c.tran, c.geq, c.iq = false, 0, 0
	switch ctx.Mode {
	case ModeOP, ModeTranOP:
		return nil
	case ModeUICInit:
		g := uicConductance
		ctx.M.Add(p, p, g)
		ctx.M.Add(n, n, g)
		ctx.M.Add(p, n, -g)
		ctx.M.Add(n, p, -g)
		ctx.AddRHS(p, g*c.vPrev)
		ctx.AddRHS(n, -g*c.vPrev)
		return nil
	}
	c.tran = true
	v, lim := c.law.limit(c.voltage(ctx.X), c.vIter)
	c.vIter, c.limited, c.v0 = v, lim, v
	c.sn, c.sdot = c.law.s(ctx.Time), c.law.sdot(ctx.Time)
	if ctx.Method == Trapezoidal {
		c.a, c.b = 2/ctx.H, -2*c.qPrev/ctx.H-c.iPrev
	} else {
		c.a, c.b = 1/ctx.H, -c.qPrev/ctx.H
	}
	q, dq := c.law.f(v)
	c.iq = c.a*c.sn*q + c.b
	c.geq = c.a * c.sn * dq
	g := c.geq
	ctx.M.Add(p, p, g)
	ctx.M.Add(n, n, g)
	ctx.M.Add(p, n, -g)
	ctx.M.Add(n, p, -g)
	// Branch current i = g·v + ieq flows from n+ to n−.
	ieq := c.iq - g*v
	ctx.AddRHS(p, -ieq)
	ctx.AddRHS(n, ieq)
	return nil
}

// scalePower is the power ṡ·F(v) leaving through the scale's port.
func (c *NonlinearCapacitor) scalePower(v float64) float64 {
	if !c.tran || c.ac || c.sdot == 0 {
		return 0
	}
	return c.sdot * c.law.coenergy(v)
}

// Operating returns voltage, current, absorbed power and the part of it
// that is stored (the rest leaves through the scale's port).
func (c *NonlinearCapacitor) Operating(x []float64) OperatingPoint {
	v := c.voltage(x)
	i := 0.0
	if c.tran {
		i = c.iq + c.geq*(v-c.v0)
	}
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i, StoragePower: v*i - c.scalePower(v)}
}

// Accept commits the state. For UIC the initial voltage is the IC value.
func (c *NonlinearCapacitor) Accept(x []float64, t float64, mode Mode) {
	v := c.voltage(x)
	switch mode {
	case ModeUICInit:
		v = c.vPrev
		c.iPrev = 0
		c.hist.reset()
	case ModeOP, ModeTranOP:
		c.iPrev = 0
		c.hist.reset()
	default:
		c.iPrev = c.iq + c.geq*(v-c.v0)
	}
	c.vPrev, c.vIter, c.vAcc = v, v, v
	c.sAcc = c.law.s(t)
	q, _ := c.law.f(v)
	c.qPrev = c.sAcc * q
	c.hist.push(t, c.qPrev)
	c.law.accepted(t)
}

// SetInitialCondition loads the UIC state (IC value, or 0 if none was given).
func (c *NonlinearCapacitor) SetInitialCondition() {
	c.vPrev, c.iPrev = 0, 0
	if c.hasIC {
		c.vPrev = c.ic
	}
}

// StoredEnergy is s·∫₀^v u·Q′(u) du at the time of the most recent Load.
func (c *NonlinearCapacitor) StoredEnergy(op OperatingPoint) float64 {
	return c.law.s(c.t) * c.law.energy(op.Voltage)
}

// Truncate estimates the step from the charge history, and lands a crossed
// knot of the charge curve on a time point.
func (c *NonlinearCapacitor) Truncate(x []float64, tc *TruncContext) float64 {
	op := c.Operating(x)
	q, _ := c.law.f(op.Voltage)
	h := lteStep(&c.hist, c.law.s(tc.Time)*q, op.Current, c.iPrev, tc.AbsTol, tc)
	return math.Min(h, c.law.truncate(c.vAcc, op.Voltage, tc))
}

// NextBreakpoint lands on the corners of the scale waveform and on a knot
// crossing found by Truncate.
func (c *NonlinearCapacitor) NextBreakpoint(t float64) float64 { return c.law.nextBreakpoint(t) }

// ExteriorPower is the work of the varying capacitance on its channel
// (+ into the circuit); 0 on every other channel and outside a transient.
func (c *NonlinearCapacitor) ExteriorPower(ch exchange.ChannelID, op OperatingPoint) float64 {
	if ch != c.ch || !c.law.varies {
		return 0
	}
	return -c.scalePower(op.Voltage)
}

// LoadAC stamps jω·s·Q′(v) at the accepted operating point.
func (c *NonlinearCapacitor) LoadAC(ctx *ACContext) error {
	c.ac, c.omega = true, ctx.Omega
	stampY(ctx.M, c.u[0], c.u[1], complex(0, ctx.Omega*c.incremental()))
	return nil
}

func (c *NonlinearCapacitor) incremental() float64 {
	_, dq := c.law.f(c.vAcc)
	return c.sAcc * dq
}

// OperatingAC returns the capacitor's phasors.
func (c *NonlinearCapacitor) OperatingAC(x []complex128) ACPoint {
	v := c.voltageAC(x)
	return ACPoint{Voltage: v, Current: v * complex(0, c.omega*c.incremental())}
}

// Limited reports whether the last Load limited the voltage update.
func (c *curvedCapacitor) Limited() bool { return c.limited }

// ConvTest checks that the linearized current at x matches the charge
// companion evaluated on the curve.
func (c *curvedCapacitor) ConvTest(x []float64, reltol, abstol float64) bool {
	if !c.tran {
		return true
	}
	v := c.voltage(x)
	ihat := c.iq + c.geq*(v-c.v0)
	q, _ := c.law.f(v)
	i := c.a*c.sn*q + c.b
	return math.Abs(ihat-i) <= reltol*math.Max(math.Abs(i), math.Abs(ihat))+abstol
}

// NonlinearInductor is an inductor defined by its flux linkage,
// λ(i, t) = s(t)·Λ(i) (circuit.NonlinearInductor), with a branch-current
// unknown (current from n+ to n− through the inductor).
//
// Maxwell derivation (CLAUDE.md §3c table, row "Inductor"): 𝔅 = V.∇𝔄, and
// the line integral of 𝔄 around the coil is its flux linkage λ, so the
// vector part of ∇𝔈 = −∂𝔅/∂t integrated around the loop gives v = dλ/dt
// in the quasi-static limit. With λ = s(t)·Λ(i),
//
//	v = s·Λ′(i)·di/dt + ṡ·Λ(i)
//
// which for Λ = L·i is v = L di/dt + i dL/dt. The branch equations
// integrate the flux:
//
//	DC:             v = 0
//	backward Euler: v_n = (λ_n − λ_{n−1})/h
//	trapezoidal:    v_n = 2(λ_n − λ_{n−1})/h − v_{n−1}
//
// with λ_n = s(t_n)·Λ(i_n), linearized in i_n when Λ is a curve. A
// saturating core with a single-valued B–H curve fits this element;
// hysteresis needs the core model of CLAUDE.md §3j.
//
// Ledger. Stored energy is s·∫₀^i u·Λ′(u) du; when s varies, ṡ·∫₀^i Λ(u) du
// is exchanged through the scale's port (mechanical by default: a moving
// core or a changing gap), as for NonlinearCapacitor.
type NonlinearInductor struct {
	Base
	law    stateLaw
	ch     exchange.ChannelID
	ic     float64
	hasIC  bool
	branch int

	// Linearization of the most recent Load: λ ≈ lam0 + ld·(i − i0); the
	// branch voltage is v = a·λ + b.
	i0, lam0, ld float64
	a, b, sn     float64
	sdot         float64
	t            float64
	tran, ac     bool
	limited      bool
	iIter        float64
	omega        float64

	vPrev, iPrev, lamPrev float64
	iAcc, sAcc            float64
	hist                  history
}

// curvedInductor is a NonlinearInductor with a flux curve (Newton needed).
type curvedInductor struct{ NonlinearInductor }

// NewNonlinearInductor builds a nonlinear or time-varying inductor.
func NewNonlinearInductor(e circuit.Element) Device {
	l := NonlinearInductor{Base: newBase(e), law: newStateLaw(e), ch: scaleChannel(e), ic: e.IC, hasIC: e.HasIC}
	_ = l.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores ∫i dλ")
	if l.law.varies {
		_ = l.SetIntrinsic(l.ch, exchange.Both, "work of the varying inductance, ṡ·∫₀ⁱ Λ(u) du")
	}
	if l.law.xs != nil {
		return &curvedInductor{l}
	}
	return &l
}

// Setup maps nodes, allocates the branch current and resolves the scale.
func (l *NonlinearInductor) Setup(ctx *SetupContext) error {
	l.setupNodes(ctx)
	l.branch = ctx.AllocBranch()
	return l.law.resolve(l.name, ctx.Step, ctx.Stop)
}

// Load stamps the branch equation for the current mode.
func (l *NonlinearInductor) Load(ctx *LoadContext) error {
	p, n, k := l.u[0], l.u[1], l.branch
	l.ac, l.limited, l.t, l.tran = false, false, ctx.Time, false
	ctx.M.Add(p, k, 1)
	ctx.M.Add(n, k, -1)
	switch ctx.Mode {
	case ModeOP, ModeTranOP:
		ctx.M.Add(k, p, 1)
		ctx.M.Add(k, n, -1)
		return nil
	case ModeUICInit:
		ctx.M.Add(k, k, 1)
		ctx.AddRHS(k, l.iPrev)
		return nil
	}
	l.tran = true
	i, lim := l.law.limit(at(ctx.X, k), l.iIter)
	l.iIter, l.limited, l.i0 = i, lim, i
	l.sn, l.sdot = l.law.s(ctx.Time), l.law.sdot(ctx.Time)
	if ctx.Method == Trapezoidal {
		l.a, l.b = 2/ctx.H, -2*l.lamPrev/ctx.H-l.vPrev
	} else {
		l.a, l.b = 1/ctx.H, -l.lamPrev/ctx.H
	}
	lam, dl := l.law.f(i)
	l.lam0, l.ld = l.sn*lam, l.sn*dl
	// V(n+) − V(n−) − a·ld·i = a·(lam0 − ld·i0) + b.
	ctx.M.Add(k, p, 1)
	ctx.M.Add(k, n, -1)
	ctx.M.Add(k, k, -l.a*l.ld)
	ctx.AddRHS(k, l.a*(l.lam0-l.ld*l.i0)+l.b)
	return nil
}

// scalePower is the power ṡ·F(i) leaving through the scale's port.
func (l *NonlinearInductor) scalePower(i float64) float64 {
	if !l.tran || l.ac || l.sdot == 0 {
		return 0
	}
	return l.sdot * l.law.coenergy(i)
}

// Operating returns voltage, current, absorbed power and the stored part.
func (l *NonlinearInductor) Operating(x []float64) OperatingPoint {
	v := l.voltage(x)
	i := x[l.branch]
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i, StoragePower: v*i - l.scalePower(i)}
}

// Accept commits the state. For UIC the initial current is the IC value.
func (l *NonlinearInductor) Accept(x []float64, t float64, mode Mode) {
	i, v := x[l.branch], l.voltage(x)
	switch mode {
	case ModeUICInit:
		i, v = l.iPrev, 0
		l.hist.reset()
	case ModeOP, ModeTranOP:
		v = 0
		l.hist.reset()
	}
	l.iPrev, l.vPrev, l.iIter, l.iAcc = i, v, i, i
	l.sAcc = l.law.s(t)
	lam, _ := l.law.f(i)
	l.lamPrev = l.sAcc * lam
	l.hist.push(t, l.lamPrev)
	l.law.accepted(t)
}

// SetInitialCondition loads the UIC state (IC value, or 0 if none was given).
func (l *NonlinearInductor) SetInitialCondition() {
	l.vPrev, l.iPrev = 0, 0
	if l.hasIC {
		l.iPrev = l.ic
	}
}

// StoredEnergy is s·∫₀^i u·Λ′(u) du at the time of the most recent Load.
func (l *NonlinearInductor) StoredEnergy(op OperatingPoint) float64 {
	return l.law.s(l.t) * l.law.energy(op.Current)
}

// Truncate estimates the step from the flux history; the "current" in the
// tolerance is the inductor voltage (the time derivative of flux). A crossed
// knot of the flux curve is landed on a time point.
func (l *NonlinearInductor) Truncate(x []float64, tc *TruncContext) float64 {
	op := l.Operating(x)
	lam, _ := l.law.f(op.Current)
	h := lteStep(&l.hist, l.law.s(tc.Time)*lam, op.Voltage, l.vPrev, tc.VnTol, tc)
	return math.Min(h, l.law.truncate(l.iAcc, op.Current, tc))
}

// NextBreakpoint lands on the corners of the scale waveform and on a knot
// crossing found by Truncate.
func (l *NonlinearInductor) NextBreakpoint(t float64) float64 { return l.law.nextBreakpoint(t) }

// ExteriorPower is the work of the varying inductance on its channel.
func (l *NonlinearInductor) ExteriorPower(ch exchange.ChannelID, op OperatingPoint) float64 {
	if ch != l.ch || !l.law.varies {
		return 0
	}
	return -l.scalePower(op.Current)
}

// LoadAC stamps V(n+) − V(n−) − jω·s·Λ′(i)·I = 0 at the operating point.
func (l *NonlinearInductor) LoadAC(ctx *ACContext) error {
	p, n, k := l.u[0], l.u[1], l.branch
	l.ac, l.omega = true, ctx.Omega
	_, dl := l.law.f(l.iAcc)
	ctx.M.Add(p, k, 1)
	ctx.M.Add(n, k, -1)
	ctx.M.Add(k, p, 1)
	ctx.M.Add(k, n, -1)
	ctx.M.Add(k, k, complex(0, -ctx.Omega*l.sAcc*dl))
	return nil
}

// OperatingAC returns the inductor's phasors.
func (l *NonlinearInductor) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: l.voltageAC(x), Current: x[l.branch]}
}

// Limited reports whether the last Load limited the current update.
func (l *curvedInductor) Limited() bool { return l.limited }

// ConvTest checks that the linearized branch voltage at x matches the flux
// companion evaluated on the curve.
func (l *curvedInductor) ConvTest(x []float64, reltol, _ float64) bool {
	if !l.tran {
		return true
	}
	i := x[l.branch]
	vhat := l.a*(l.lam0+l.ld*(i-l.i0)) + l.b
	lam, _ := l.law.f(i)
	v := l.a*l.sn*lam + l.b
	return math.Abs(vhat-v) <= reltol*math.Max(math.Abs(v), math.Abs(vhat))+1e-6 // vntol
}
