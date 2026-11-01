// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Bearing is the kronspice P element, the ball-bearing (Huber-effect) motor:
// a resistive current path a → b through a rotating bearing, and a shaft
// port (m+, m−; across = ω, through = torque) driven with
//
//	T = κ·I²·ψ(ω)
//
// See circuit.BearingMotor for ψ in the two modes and the sources.
//
// Unknowns: the nodes and one branch current I (a → b). Equations:
//
//	branch:  V_a − V_b = VF·tanh(I/1 mA) + R·I + g(ω)·I,
//	         g = κ·tanh(ω/ω0)·ω (EM), 0 (thermal)
//	port:    the element drives torque T into m+ (current −T from m+ to m−)
//
// Both are linearized at each Newton iterate. V·I − ω·T is the power that
// leaves as heat through the thermal port: H = VF·|I| + R·I² in EM mode,
// H − ω·T (≥ (1 − η)·H) in thermal mode. VF is the measured, almost
// constant DC voltage of a current-carrying ball bearing (0.5–1.5 V,
// Safdarzadeh et al., Lubricants 12, 148, 2024). So every joule of shaft work is drawn
// from the electrical circuit, and the ledger closes.
//
// At ω = 0 the torque is zero for any current (no start from rest, as in
// every report that a hand spin is needed: Huber 1959; M. Harrison,
// electricstuff.co.uk; McDonald 2020 fn. 6), but dT/dω = κ·I²/ω0 > 0, so
// the rest state is unstable once κ·I²/ω0 exceeds the shaft's damping:
// a small spin in either sense grows, the direction set by the spin
// (Marinov 1989: "left or right").
//
// Maxwell derivation (CLAUDE.md §3c): in EM mode the torque is the Lorentz
// force of the bearing's own current and field on the current crowded into
// the moving ball contacts (Gruenberg 1978; Hatzikonstantinou & Moyssides,
// J. Phys. A 23, 3183, 1990); its work is the back-EMF term (Poynting). The
// thermal mode is not electromagnetic: it is a heat engine on the contact
// heating (Watson, Williams & Crimp, IEE Proc. A 140, 1993; Marinov 1989),
// bounded by its measured efficiency.
type Bearing struct {
	Base
	r, vf, kappa, w0, eta float64
	thermal               bool
	extra                 [2]int
	p                     [2]int // shaft unknowns
	br                    int

	// linearization at the last Load
	linOK              bool
	linT, linTi, linTw float64
	linI0, linW0       float64
	linG, linDG        float64
	linD, linDD        float64 // contact drop and its slope
	iOp, wOp           float64
}

// dropAmps is the current over which the contact drop VF·tanh(I/dropAmps)
// changes sign: a numerical rounding of the measured, almost constant DC
// bearing voltage, not a physical parameter.
const dropAmps = 1e-3

// NewBearing builds a bearing motor from a circuit element.
func NewBearing(e circuit.Element) *Bearing {
	kappa, w0, thermal, eta := circuit.BearingParams(e)
	b := &Bearing{Base: newBase(e), r: e.Value, vf: e.Params["vf"], kappa: kappa, w0: w0, eta: eta, thermal: thermal}
	if len(e.Extra) == 2 {
		b.extra = [2]int{e.Extra[0], e.Extra[1]}
	}
	_ = b.SetIntrinsic(exchange.Electrical, exchange.Both, "bearing current path inside the network")
	_ = b.SetIntrinsic(exchange.Mechanical, exchange.Both,
		"Huber torque κI²ψ(ω) on the shaft: energy converted inside the network")
	_ = b.SetIntrinsic(exchange.Thermal, exchange.Out, "contact and race heating V·I − ω·T to surroundings")
	return b
}

// Setup maps the current path, the shaft port and the branch current.
func (b *Bearing) Setup(ctx *SetupContext) error {
	b.setupNodes(ctx)
	b.p = [2]int{ctx.Unknown(b.extra[0]), ctx.Unknown(b.extra[1])}
	b.br = ctx.AllocBranch()
	return nil
}

func (b *Bearing) omega(x []float64) float64 { return at(x, b.p[0]) - at(x, b.p[1]) }

// drop returns the contact drop VF·tanh(I/I0) and its slope.
func (b *Bearing) drop(i float64) (float64, float64) {
	if b.vf == 0 {
		return 0, 0
	}
	t := math.Tanh(i / dropAmps)
	return b.vf * t, b.vf * (1 - t*t) / dropAmps
}

// heat is the electrical loss before any shaft work: drop·I + R·I², and its
// derivative in I.
func (b *Bearing) heat(i float64) (float64, float64) {
	d, dd := b.drop(i)
	return d*i + b.r*i*i, dd*i + d + 2*b.r*i
}

// torque returns T(I, ω) and ∂T/∂I, ∂T/∂ω. EM mode: T = κ·I²·tanh(ω/ω0).
// Thermal mode: T = a/(1 + a·ω/(η·H)) with a the EM form and H the heat,
// so ω·T < η·H.
func (b *Bearing) torque(i, w float64) (t, ti, tw float64) {
	th := math.Tanh(w / b.w0)
	dth := (1 - th*th) / b.w0
	a, ai, aw := b.kappa*th*i*i, 2*b.kappa*th*i, b.kappa*dth*i*i
	if !b.thermal {
		return a, ai, aw
	}
	h, hi := b.heat(i)
	if h <= 0 || a == 0 {
		return a, ai, aw
	}
	q := a * w / (b.eta * h)
	qi := w * (ai*h - a*hi) / (b.eta * h * h)
	qw := (aw*w + a) / (b.eta * h)
	d := (1 + q) * (1 + q)
	return a / (1 + q), (ai*(1+q) - a*qi) / d, (aw*(1+q) - a*qw) / d
}

// emf returns g(ω) = κ·tanh(ω/ω0)·ω and dg/dω (0 in thermal mode).
func (b *Bearing) emf(w float64) (float64, float64) {
	if b.thermal {
		return 0, 0
	}
	th := math.Tanh(w / b.w0)
	return b.kappa * th * w, b.kappa * ((1-th*th)/b.w0*w + th)
}

// Torque is T(I, ω) at solution x [N·m].
func (b *Bearing) Torque(x []float64) float64 {
	t, _, _ := b.torque(at(x, b.br), b.omega(x))
	return t
}

// BackEMF is g(ω)·I at solution x [V].
func (b *Bearing) BackEMF(x []float64) float64 {
	g, _ := b.emf(b.omega(x))
	return g * at(x, b.br)
}

// Load stamps the linearized branch equation and the torque source.
func (b *Bearing) Load(ctx *LoadContext) error {
	x := ctx.X
	i0, w0 := at(x, b.br), b.omega(x)
	g, dg := b.emf(w0)
	dv, ddv := b.drop(i0)
	ua, ub, pp, pn, br := b.u[0], b.u[1], b.p[0], b.p[1], b.br
	m := ctx.M
	// KCL: I leaves a, enters b.
	m.Add(ua, br, 1)
	m.Add(ub, br, -1)
	// V_a − V_b = (R + g0)·I + D0 + D0'·(I − I0) + g0'·I0·(ω − ω0)
	m.Add(br, ua, 1)
	m.Add(br, ub, -1)
	m.Add(br, br, -(b.r + g + ddv))
	rhs := dv - ddv*i0
	if dg != 0 {
		m.Add(br, pp, -dg*i0)
		m.Add(br, pn, dg*i0)
		rhs -= dg * i0 * w0
	}
	ctx.AddRHS(br, rhs)
	// Torque T ≈ T0 + Ti·(I − I0) + Tw·(ω − ω0) driven into m+.
	t0, ti, tw := b.torque(i0, w0)
	m.Add(pp, br, -ti)
	m.Add(pp, pp, -tw)
	m.Add(pp, pn, tw)
	m.Add(pn, br, ti)
	m.Add(pn, pp, tw)
	m.Add(pn, pn, -tw)
	c := t0 - ti*i0 - tw*w0
	ctx.AddRHS(pp, c)
	ctx.AddRHS(pn, -c)
	b.linOK, b.linT, b.linTi, b.linTw, b.linI0, b.linW0 = true, t0, ti, tw, i0, w0
	b.linG, b.linDG, b.linD, b.linDD = g, dg, dv, ddv
	return nil
}

// portTorque is the torque the last linear solve used (so the port's KCL
// closes in the ledger), or the exact torque before any Load.
func (b *Bearing) portTorque(x []float64) float64 {
	if !b.linOK {
		return b.Torque(x)
	}
	return b.linT + b.linTi*(at(x, b.br)-b.linI0) + b.linTw*(b.omega(x)-b.linW0)
}

// Limited: no junction limiting (tanh is bounded).
func (b *Bearing) Limited() bool { return false }

// ConvTest checks that the linearized torque and back-EMF agree with the
// element's equations at x.
func (b *Bearing) ConvTest(x []float64, reltol, abstol float64) bool {
	t, tl := b.Torque(x), b.portTorque(x)
	if math.Abs(t-tl) > reltol*math.Max(math.Abs(t), math.Abs(tl))+abstol {
		return false
	}
	if b.linOK && b.vf != 0 {
		d, _ := b.drop(at(x, b.br))
		dl := b.linD + b.linDD*(at(x, b.br)-b.linI0)
		if math.Abs(d-dl) > reltol*math.Max(math.Abs(d), math.Abs(dl))+abstol {
			return false
		}
	}
	if b.linOK && !b.thermal {
		i, w := at(x, b.br), b.omega(x)
		e, el := b.BackEMF(x), b.linG*i+b.linDG*b.linI0*(w-b.linW0)
		if math.Abs(e-el) > reltol*math.Max(math.Abs(e), math.Abs(el))+abstol {
			return false
		}
	}
	return true
}

// Operating reports the current path (V, I) and the element's net absorbed
// power V·I − ω·T: the heat it sends out through its thermal port.
func (b *Bearing) Operating(x []float64) OperatingPoint {
	v, i, w := b.voltage(x), at(x, b.br), b.omega(x)
	b.iOp, b.wOp = i, w
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v*i - w*b.portTorque(x)}
}

// Accept records the operating point for the small-signal model.
func (b *Bearing) Accept(x []float64, _ float64, _ Mode) {
	b.iOp, b.wOp = at(x, b.br), b.omega(x)
}

// ExteriorPower: the heat V·I − ω·T leaves through the thermal port.
func (b *Bearing) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -op.PowerAbsorbed
	}
	return 0
}

// LoadAC stamps the model linearized about the operating point (I0, ω0).
func (b *Bearing) LoadAC(ctx *ACContext) error {
	i0, w0 := b.iOp, b.wOp
	g, dg := b.emf(w0)
	_, ddv := b.drop(i0)
	_, ti, tw := b.torque(i0, w0)
	ua, ub, pp, pn, br := b.u[0], b.u[1], b.p[0], b.p[1], b.br
	m := ctx.M
	m.Add(ua, br, 1)
	m.Add(ub, br, -1)
	m.Add(br, ua, 1)
	m.Add(br, ub, -1)
	m.Add(br, br, complex(-(b.r+g+ddv), 0))
	m.Add(br, pp, complex(-dg*i0, 0))
	m.Add(br, pn, complex(dg*i0, 0))
	m.Add(pp, br, complex(-ti, 0))
	m.Add(pp, pp, complex(-tw, 0))
	m.Add(pp, pn, complex(tw, 0))
	m.Add(pn, br, complex(ti, 0))
	m.Add(pn, pp, complex(tw, 0))
	m.Add(pn, pn, complex(-tw, 0))
	return nil
}

// OperatingAC returns the current path's phasors.
func (b *Bearing) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: b.voltageAC(x), Current: x[b.br]}
}

// ACPower is ½·(V·I* + Ω·(−T)*) over both ports.
func (b *Bearing) ACPower(x []complex128) complex128 {
	_, ti, tw := b.torque(b.iOp, b.wOp)
	i, w := x[b.br], cat(x, b.p[0])-cat(x, b.p[1])
	t := complex(ti, 0)*i + complex(tw, 0)*w
	return 0.5 * (b.voltageAC(x)*cmplx.Conj(i) - w*cmplx.Conj(t))
}
