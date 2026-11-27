// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// mu0 is the vacuum permeability [H/m] (CODATA 2018).
const mu0 = 1.25663706212e-6

// Core is the kronspice U element: windings on a saturable core whose
// magnetization follows the Jiles–Atherton model (circuit.Core).
//
// Unknowns: the nodes and one current per winding. With H = Σ N_j·i_j/le
// and B(H) the magnetization history integrated from the accepted state,
// each winding obeys Faraday's law v_k = N_k·Ae·dB/dt, discretized as an
// inductor's flux (Nagel 1975):
//
//	backward Euler: v_k = N_k·Ae·(B − B₀)/h
//	trapezoidal:    v_k = 2·N_k·Ae·(B − B₀)/h − v_k₀
//
// and linearized at the Newton iterate with dB/dH = μ0·(1 + dM/dH). At an
// operating point every winding is a short (v_k = 0), and the core follows
// the initial magnetization curve to the operating point's H.
//
// The Jiles–Atherton equations are integrated in H with the effective
// field He = H + αM (Jiles & Atherton 1986, eqs. 1–8; Jiles, Thoelke &
// Devine, IEEE Trans. Magn. 28, 27, 1992):
//
//	dMirr/dHe = (Man − Mirr)/(δk)   (0 when δ(Man − Mirr) < 0, Lederer 1999)
//	χe = (1 − c)·dMirr/dHe + c·dMan/dHe
//	dM/dH = χe/(1 − αχe),  dMirr/dH = (dMirr/dHe)/(1 − αχe)
//
// by the midpoint rule in sub-steps no larger than a/32. The Jacobian uses
// the end point's dM/dH, the derivative of the exact path (within the
// sub-step error).
//
// Energy (CLAUDE.md §3c ledger). The winding power Σ v_k·i_k splits into
// hysteresis heat, which leaves through the thermal port, and the rest,
// which is booked as stored field energy. The heat is the Jiles–Atherton
// pinning loss (Jiles & Atherton 1986, eq. 12: the energy dissipated
// against pinning sites is μ0·k·∫|dMirr| per volume, with the reversible
// fraction c removed),
//
//	P_loss = Ae·le·μ0·k·(1 − c)·|dMirr/dt| ≥ 0,
//
// which over a closed major loop equals the loop area ∮H dB times the core
// volume (0.1 % for the 1986 parameters). The stored energy is the running
// integral of the remainder, so the ledger closes at every point.
type Core struct {
	Base
	p     circuit.CoreParams
	nw    int
	pairs [][2]int // winding nodes
	w     [][2]int // winding node unknowns
	br    []int    // winding current unknowns

	// accepted state
	h0, m0, mirr0, chi0 float64
	v0                  []float64 // accepted winding voltages (trapezoidal history)
	wAcc, spAcc         float64   // stored energy and its rate at the accepted point
	lossAcc             float64   // heat rate over the accepted step
	justAccepted        bool      // no Load since the last Accept

	// linearization of the last Load
	tran                bool
	lastH               float64
	method              Method
	linH, linB, linDBdH float64

	spOp float64 // storage rate at the last Operating call
}

// NewCore builds a core from a circuit element.
func NewCore(e circuit.Element) *Core {
	c := &Core{Base: newBase(e), p: circuit.CoreParamsOf(e)}
	c.nw = len(c.p.Turns)
	c.pairs = circuit.CoreWindings(e)
	c.w = make([][2]int, c.nw)
	c.br = make([]int, c.nw)
	c.v0 = make([]float64, c.nw)
	c.chi0 = c.initialChi()
	_ = c.SetIntrinsic(exchange.Electrical, exchange.Both, "windings connect inside the network; core stores field energy")
	_ = c.SetIntrinsic(exchange.Magnetic, exchange.Both, "flux confined to the core: energy stays inside the circuit")
	_ = c.SetIntrinsic(exchange.Thermal, exchange.Out, "hysteresis loss to surroundings")
	return c
}

func (c *Core) initialChi() float64 {
	_, _, chi := c.deriv(0, 0, 0, 1)
	return chi
}

// Setup maps every winding and allocates its current.
func (c *Core) Setup(ctx *SetupContext) error {
	c.setupNodes(ctx)
	for k, p := range c.pairs {
		c.w[k] = [2]int{ctx.Unknown(p[0]), ctx.Unknown(p[1])}
		c.br[k] = ctx.AllocBranch()
	}
	return nil
}

// langevin returns L(x) = coth x − 1/x and L'(x).
func langevin(x float64) (float64, float64) {
	ax := math.Abs(x)
	if ax < 1e-2 {
		x2 := x * x
		return x * (1.0/3 - x2/45 + 2*x2*x2/945), 1.0/3 - x2/15 + 2*x2*x2/189
	}
	if ax > 30 {
		return math.Copysign(1, x) - 1/x, 1 / (x * x)
	}
	s := math.Sinh(x)
	return 1/math.Tanh(x) - 1/x, 1/(x*x) - 1/(s*s)
}

// deriv returns dM/dH, dMirr/dH and dM/dH again (for the Jacobian) at
// (H, M, Mirr) for the sweep direction delta (±1).
func (c *Core) deriv(h, m, mirr, delta float64) (dm, dmirr, chi float64) {
	p := &c.p
	he := h + p.Alpha*m
	l, dl := langevin(he / p.A)
	man, dman := p.Ms*l, p.Ms/p.A*dl
	d := man - mirr
	dmirrE := 0.0
	if delta*d > 0 {
		dmirrE = d / (delta * p.K)
	}
	chiE := (1-p.C)*dmirrE + p.C*dman
	den := 1 - p.Alpha*chiE
	if den < 0.05 { // keep the mean-field feedback from running away
		den = 0.05
	}
	return chiE / den, dmirrE / den, chiE / den
}

// integrate follows the Jiles–Atherton equations from (h0, m0, mirr0) to h1.
func (c *Core) integrate(h0, m0, mirr0, h1 float64) (m, mirr, chi float64) {
	dh := h1 - h0
	delta := 1.0
	if dh < 0 {
		delta = -1
	}
	if dh == 0 {
		_, _, chi = c.deriv(h0, m0, mirr0, delta)
		return m0, mirr0, chi
	}
	n := int(math.Ceil(math.Abs(dh) / (c.p.A / 32)))
	if n > 4000 {
		n = 4000
	}
	step := dh / float64(n)
	h, m, mirr := h0, m0, mirr0
	for s := 0; s < n; s++ {
		k1m, k1i, _ := c.deriv(h, m, mirr, delta)
		k2m, k2i, _ := c.deriv(h+step/2, m+step/2*k1m, mirr+step/2*k1i, delta)
		m += step * k2m
		mirr += step * k2i
		h += step
	}
	if lim := c.p.Ms; math.Abs(m) > lim {
		m = math.Copysign(lim, m)
	}
	_, _, chi = c.deriv(h1, m, mirr, delta)
	return m, mirr, chi
}

// field returns H from the winding currents in x.
func (c *Core) field(x []float64) float64 {
	h := 0.0
	for k, n := range c.p.Turns {
		h += n * at(x, c.br[k])
	}
	return h / c.p.Len
}

// Load stamps KCL and each winding's Faraday equation.
func (c *Core) Load(ctx *LoadContext) error {
	mtx := ctx.M
	for k := range c.p.Turns {
		ua, ub, br := c.w[k][0], c.w[k][1], c.br[k]
		mtx.Add(ua, br, 1)
		mtx.Add(ub, br, -1)
		mtx.Add(br, ua, 1)
		mtx.Add(br, ub, -1)
	}
	c.justAccepted = false
	if ctx.Mode != ModeTran {
		c.tran = false
		return nil
	}
	c.tran, c.lastH, c.method = true, ctx.H, ctx.Method
	h := c.field(ctx.X)
	m, _, chi := c.integrate(c.h0, c.m0, c.mirr0, h)
	b, b0 := mu0*(h+m), mu0*(c.h0+c.m0)
	dbdh := mu0 * (1 + chi)
	c.linH, c.linB, c.linDBdH = h, b, dbdh
	kap := 1 / ctx.H
	if ctx.Method == Trapezoidal {
		kap = 2 / ctx.H
	}
	for k, nk := range c.p.Turns {
		br := c.br[k]
		g := kap * nk * c.p.Area * dbdh / c.p.Len
		for j, nj := range c.p.Turns {
			mtx.Add(br, c.br[j], -g*nj)
		}
		rhs := kap * nk * c.p.Area * (b - dbdh*h - b0)
		if ctx.Method == Trapezoidal {
			rhs -= c.v0[k]
		}
		ctx.AddRHS(br, rhs)
	}
	return nil
}

// Limited reports no junction limiting.
func (c *Core) Limited() bool { return false }

// ConvTest checks that B(H) at x agrees with the linearization of the last
// Load.
func (c *Core) ConvTest(x []float64, reltol, _ float64) bool {
	if !c.tran {
		return true
	}
	h := c.field(x)
	m, _, _ := c.integrate(c.h0, c.m0, c.mirr0, h)
	b, bl := mu0*(h+m), c.linB+c.linDBdH*(h-c.linH)
	return math.Abs(b-bl) <= reltol*math.Max(math.Abs(b), math.Abs(bl))+1e-9
}

// power returns Σ v_k·i_k at x.
func (c *Core) power(x []float64) float64 {
	p := 0.0
	for k := range c.p.Turns {
		p += (at(x, c.w[k][0]) - at(x, c.w[k][1])) * at(x, c.br[k])
	}
	return p
}

// loss is the pinning heat rate of the step from the accepted state to x.
func (c *Core) loss(x []float64) float64 {
	if !c.tran {
		return 0
	}
	_, mirr, _ := c.integrate(c.h0, c.m0, c.mirr0, c.field(x))
	return c.p.Area * c.p.Len * mu0 * c.p.K * (1 - c.p.C) * math.Abs(mirr-c.mirr0) / c.lastH
}

// Operating reports winding 1, the core's total absorbed power and the part
// of it going into storage (the rest is the hysteresis heat).
func (c *Core) Operating(x []float64) OperatingPoint {
	p := c.power(x)
	l := c.lossAcc
	if !c.justAccepted {
		l = c.loss(x)
	}
	c.spOp = p - l
	return OperatingPoint{Voltage: at(x, c.w[0][0]) - at(x, c.w[0][1]), Current: at(x, c.br[0]),
		PowerAbsorbed: p, StoragePower: p - l}
}

// Accept commits the magnetization and the energy split of the accepted
// step.
func (c *Core) Accept(x []float64, _ float64, mode Mode) {
	p := c.power(x)
	l := 0.0
	if mode == ModeTran && c.tran {
		l = c.loss(x)
		sp := p - l
		c.wAcc += 0.5 * c.lastH * (c.spAcc + sp)
		c.spAcc = sp
	} else {
		c.spAcc = 0
	}
	c.lossAcc = l
	h := c.field(x)
	c.m0, c.mirr0, c.chi0 = c.integrate(c.h0, c.m0, c.mirr0, h)
	c.h0 = h
	for k := range c.p.Turns {
		c.v0[k] = at(x, c.w[k][0]) - at(x, c.w[k][1])
	}
	c.justAccepted = true
}

// State returns the accepted H [A/m], M [A/m] and B [T].
func (c *Core) State() (h, m, b float64) { return c.h0, c.m0, mu0 * (c.h0 + c.m0) }

// StoredEnergy returns the stored field energy booked so far [J] (at the
// last Operating call).
func (c *Core) StoredEnergy(OperatingPoint) float64 {
	if c.justAccepted {
		return c.wAcc
	}
	return c.wAcc + 0.5*c.lastH*(c.spAcc+c.spOp)
}

// ExteriorPower: the hysteresis heat (absorbed power less the storage
// rate) leaves through the thermal port.
func (c *Core) ExteriorPower(ch exchange.ChannelID, op OperatingPoint) float64 {
	if ch == exchange.Thermal {
		return -(op.PowerAbsorbed - op.StoragePower)
	}
	return 0
}

// Truncate keeps each step's change of H within a and of B within 2 % of
// μ0·Ms, so the B–H path is followed on the curve's own scale.
func (c *Core) Truncate(x []float64, tc *TruncContext) float64 {
	h := c.field(x)
	m, _, _ := c.integrate(c.h0, c.m0, c.mirr0, h)
	dh := math.Abs(h - c.h0)
	db := math.Abs(mu0 * (h + m - c.h0 - c.m0))
	r := math.Max(dh/c.p.A, db/(0.02*mu0*c.p.Ms))
	if r < 1e-9 {
		return inf
	}
	return tc.H / r
}

// LoadAC stamps the windings as coupled inductors with the core's
// differential permeability at the operating point.
func (c *Core) LoadAC(ctx *ACContext) error {
	dbdh := mu0 * (1 + c.chi0)
	for k, nk := range c.p.Turns {
		ua, ub, br := c.w[k][0], c.w[k][1], c.br[k]
		ctx.M.Add(ua, br, 1)
		ctx.M.Add(ub, br, -1)
		ctx.M.Add(br, ua, 1)
		ctx.M.Add(br, ub, -1)
		for j, nj := range c.p.Turns {
			l := nk * nj * c.p.Area * dbdh / c.p.Len
			ctx.M.Add(br, c.br[j], complex(0, -ctx.Omega*l))
		}
	}
	return nil
}

// OperatingAC returns winding 1's phasors.
func (c *Core) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: cat(x, c.w[0][0]) - cat(x, c.w[0][1]), Current: x[c.br[0]]}
}

// ACPower is ½·Σ V_k·I_k* over the windings.
func (c *Core) ACPower(x []complex128) complex128 {
	var s complex128
	for k := range c.p.Turns {
		s += (cat(x, c.w[k][0]) - cat(x, c.w[k][1])) * cmplx.Conj(x[c.br[k]])
	}
	return 0.5 * s
}
