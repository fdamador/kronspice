// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Physical constants (CODATA 2018, exact in SI since 2019).
const (
	ChargeQ   = 1.602176634e-19 // elementary charge [C]
	BoltzmanK = 1.380649e-23    // Boltzmann constant [J/K]
	// TNomK is the SPICE nominal temperature, 27 °C.
	TNomK = 300.15
)

// DiodeParams are the SPICE3 level-1 diode model parameters, already scaled
// by the element's area factor.
type DiodeParams struct {
	IS, N, RS, CJO, VJ, M, TT, BV, IBV, FC float64
	HasBV                                  bool
}

// DefaultDiodeParams returns the SPICE3 defaults (SPICE3f3 User's Manual §4.2).
func DefaultDiodeParams() DiodeParams {
	return DiodeParams{IS: 1e-14, N: 1, RS: 0, CJO: 0, VJ: 1, M: 0.5, TT: 0, IBV: 1e-3, FC: 0.5}
}

// Diode is the SPICE3 level-1 junction diode.
//
// Maxwell derivation (CLAUDE.md §3c): the junction current is a
// semiconductor constitutive relation (Shockley's law) for the conduction
// current 𝔎 in the quasi-static limit; the depletion and diffusion charges
// give a displacement-current term ∂q/∂t, as for a capacitor.
//
// Equations (ported from the SPICE3 equations, not its code; Nagel 1975
// ch. 4; SPICE3f5 dioload/DEVpnjlim; Massobrio & Antognetti, Semiconductor
// Device Modeling with SPICE, ch. 1), with vte = N·kT/q:
//
//	vd ≥ −3·vte:            id = IS(e^(vd/vte) − 1) + gmin·vd
//	−XBV ≤ vd < −3·vte:      id = −IS[1 + (3·vte/(e·vd))³] + gmin·vd
//	vd < −XBV (if BV given): id = −IS·e^(−(XBV + vd)/vte) + gmin·vd
//
// XBV is the breakdown voltage adjusted so that the reverse current at BV is
// IBV: IS·(e^((BV−XBV)/vt) − 1 + XBV/vt) = IBV with vt = kT/q, or XBV = BV
// when IBV < IS·BV/vt (SPICE3f5 temperature setup of the diode).
//
//	depletion charge, vd < FC·VJ:  qj = CJO·VJ·(1 − (1 − vd/VJ)^(1−M))/(1−M)
//	beyond FC·VJ the capacitance is continued linearly (SPICE F1, F2, F3)
//	diffusion charge:              qd = TT·id
//
// Newton limiting uses pnjlim with vcrit = vte·ln(vte/(√2·IS)).
//
// Series resistance RS sits between the anode terminal and an internal node.
// Ports: thermal (out) is intrinsic and carries the dissipated power (junction
// and RS); stored junction energy is internal.
type Diode struct {
	Base
	P     DiodeParams
	vte   float64
	vcrit float64
	xbv   float64 // adjusted breakdown voltage (see the type comment)
	grs   float64 // 1/RS, 0 if RS = 0
	ai    int     // unknown index of the internal anode node (= anode if RS = 0)

	// Linearization of the most recent Load.
	vd, cd, gd    float64 // limited junction voltage, junction current and conductance
	iq, geq       float64 // charge current and its companion conductance (transient)
	tran          bool
	limited       bool
	hasCharge     bool
	gmin          float64
	vdSol         float64 // junction voltage of the most recent Operating call
	qPrev, iqPrev float64 // accepted charge and charge current
	vdAcc         float64 // accepted junction voltage (AC linearization, limiting)
	hist          history
	omega         float64
}

// NewDiode builds a diode from a circuit element and its model.
func NewDiode(e circuit.Element, m circuit.Model) *Diode {
	p := DefaultDiodeParams()
	get := func(k string, dst *float64) {
		if v, ok := m.Params[k]; ok {
			*dst = v
		}
	}
	get("is", &p.IS)
	get("n", &p.N)
	get("rs", &p.RS)
	get("cjo", &p.CJO)
	get("cj0", &p.CJO)
	get("vj", &p.VJ)
	get("m", &p.M)
	get("tt", &p.TT)
	get("ibv", &p.IBV)
	get("fc", &p.FC)
	if v, ok := m.Params["bv"]; ok {
		p.BV, p.HasBV = v, true
	}
	area := e.Area
	if area == 0 {
		area = 1
	}
	p.IS *= area
	p.CJO *= area
	p.IBV *= area
	if p.RS > 0 {
		p.RS /= area
	}
	d := &Diode{Base: newBase(e), P: p}
	d.vte = p.N * BoltzmanK * TNomK / ChargeQ
	d.vcrit = d.vte * math.Log(d.vte/(math.Sqrt2*p.IS))
	if p.RS > 0 {
		d.grs = 1 / p.RS
	}
	d.hasCharge = p.CJO > 0 || p.TT > 0
	if p.HasBV {
		d.xbv = breakdownVoltage(p.IS, p.BV, p.IBV, BoltzmanK*TNomK/ChargeQ)
	}
	_ = d.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores junction charge")
	_ = d.SetIntrinsic(exchange.Thermal, exchange.Out, "junction and series-resistance heating to surroundings")
	return d
}

// Setup maps nodes and allocates the internal node when RS > 0.
func (d *Diode) Setup(ctx *SetupContext) error {
	d.setupNodes(ctx)
	d.ai = d.u[0]
	if d.grs > 0 {
		d.ai = ctx.AllocUnknown()
	}
	return nil
}

func at(x []float64, i int) float64 {
	if i < 0 || i >= len(x) {
		return 0
	}
	return x[i]
}

// current returns the junction current and conductance at vd (without gmin).
func (d *Diode) current(vd float64) (id, gd float64) {
	p, vte := &d.P, d.vte
	switch {
	case vd >= -3*vte:
		e := math.Exp(math.Min(vd/vte, 700))
		return p.IS * (e - 1), p.IS * e / vte
	case !p.HasBV || vd >= -d.xbv:
		a := 3 * vte / (vd * math.E)
		a = a * a * a
		return -p.IS * (1 + a), p.IS * 3 * a / vd
	default:
		e := math.Exp(math.Min(-(d.xbv+vd)/vte, 700))
		return -p.IS * e, p.IS * e / vte
	}
}

// breakdownVoltage solves IS·(e^((bv−x)/vt) − 1 + x/vt) = ibv for x by
// bisection; the left side decreases in x below bv, so the root is unique.
func breakdownVoltage(is, bv, ibv, vt float64) float64 {
	if ibv < is*bv/vt {
		return bv
	}
	g := func(x float64) float64 { return is*(math.Exp(math.Min((bv-x)/vt, 700))-1+x/vt) - ibv }
	lo, hi := math.Min(0, bv-vt*math.Log1p(ibv/is))-10*vt, bv
	for i := 0; i < 200 && g(lo) <= 0; i++ {
		lo -= 10 * vt
	}
	for i := 0; i < 200; i++ {
		mid := 0.5 * (lo + hi)
		if mid == lo || mid == hi {
			break
		}
		if g(mid) > 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// charge returns the junction charge and capacitance at vd for diffusion
// current id and conductance gd.
func (d *Diode) charge(vd, id, gd float64) (q, c float64) {
	p := &d.P
	q, c = p.TT*id, p.TT*gd
	if p.CJO == 0 {
		return q, c
	}
	fcv := p.FC * p.VJ
	if vd < fcv {
		arg := 1 - vd/p.VJ
		sarg := math.Pow(arg, -p.M)
		q += p.CJO * p.VJ * (1 - arg*sarg) / (1 - p.M)
		c += p.CJO * sarg
		return q, c
	}
	f1 := p.VJ * (1 - math.Pow(1-p.FC, 1-p.M)) / (1 - p.M)
	f2 := math.Pow(1-p.FC, 1+p.M)
	f3 := 1 - p.FC*(1+p.M)
	q += p.CJO*f1 + p.CJO/f2*(f3*(vd-fcv)+p.M/(2*p.VJ)*(vd*vd-fcv*fcv))
	c += p.CJO / f2 * (f3 + p.M*vd/p.VJ)
	return q, c
}

// pnjlim limits a junction-voltage update (Nagel 1975; SPICE3 DEVpnjlim).
func pnjlim(vnew, vold, vt, vcrit float64) (float64, bool) {
	if vnew > vcrit && math.Abs(vnew-vold) > 2*vt {
		if vold > 0 {
			arg := 1 + (vnew-vold)/vt
			if arg > 0 {
				return vold + vt*math.Log(arg), true
			}
			return vcrit, true
		}
		return vt * math.Log(vnew/vt), true
	}
	return vnew, false
}

// Load linearizes the diode at the current iterate and stamps it.
func (d *Diode) Load(ctx *LoadContext) error {
	a, k, ai := d.u[0], d.u[1], d.ai
	d.gmin = ctx.Gmin
	vraw := at(ctx.X, ai) - at(ctx.X, k)
	switch {
	case ctx.InitJct:
		d.vd, d.limited = d.vcrit, true
	default:
		d.vd, d.limited = pnjlim(vraw, d.vdAcc, d.vte, d.vcrit)
		if bv := d.xbv; d.P.HasBV && vraw < math.Min(0, -bv+10*d.vte) {
			// Limit in breakdown by reflecting about −BV (SPICE3 dioload).
			vt, lim := pnjlim(-(vraw + bv), -(d.vdAcc + bv), d.vte, d.vcrit)
			d.vd, d.limited = -(vt + bv), lim
		}
	}
	d.vdAcc = d.vd // the next iteration limits relative to this one
	id, gd := d.current(d.vd)
	d.cd, d.gd = id+d.gmin*d.vd, gd+d.gmin

	d.tran, d.iq, d.geq = false, 0, 0
	if ctx.Mode == ModeTran && d.hasCharge {
		q, c := d.charge(d.vd, id, gd)
		d.tran = true
		if ctx.Method == Trapezoidal {
			d.iq = 2*(q-d.qPrev)/ctx.H - d.iqPrev
			d.geq = 2 * c / ctx.H
		} else {
			d.iq = (q - d.qPrev) / ctx.H
			d.geq = c / ctx.H
		}
	}
	g := d.gd + d.geq
	ieq := d.cd + d.iq - g*d.vd
	ctx.M.Add(ai, ai, g)
	ctx.M.Add(k, k, g)
	ctx.M.Add(ai, k, -g)
	ctx.M.Add(k, ai, -g)
	ctx.AddRHS(ai, -ieq)
	ctx.AddRHS(k, ieq)
	if d.grs > 0 {
		ctx.M.Add(a, a, d.grs)
		ctx.M.Add(ai, ai, d.grs)
		ctx.M.Add(a, ai, -d.grs)
		ctx.M.Add(ai, a, -d.grs)
	}
	return nil
}

// Limited reports whether the last Load limited the junction voltage.
func (d *Diode) Limited() bool { return d.limited }

// ConvTest checks that the linearized current at x agrees with the current
// at the linearization point (SPICE3 DIOconvTest).
func (d *Diode) ConvTest(x []float64, reltol, abstol float64) bool {
	vd := at(x, d.ai) - at(x, d.u[1])
	cdhat := d.cd + d.gd*(vd-d.vd)
	id, _ := d.current(vd)
	cd := id + d.gmin*vd
	tol := reltol*math.Max(math.Abs(cdhat), math.Abs(cd)) + abstol
	return math.Abs(cdhat-cd) <= tol
}

// Operating evaluates the linearized model at the solution, so the reported
// currents satisfy KCL exactly and the energy ledger balances.
func (d *Diode) Operating(x []float64) OperatingPoint {
	vd := at(x, d.ai) - at(x, d.u[1])
	d.vdSol = vd
	ir := d.cd + d.gd*(vd-d.vd)  // resistive (conduction) part
	ic := d.iq + d.geq*(vd-d.vd) // charge part
	i := ir + ic
	v := d.voltage(x)
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i, StoragePower: vd * ic}
}

// Accept commits the converged junction state.
func (d *Diode) Accept(x []float64, t float64, mode Mode) {
	vd := at(x, d.ai) - at(x, d.u[1])
	d.vdAcc = vd
	id, gd := d.current(vd)
	q, _ := d.charge(vd, id, gd)
	switch mode {
	case ModeTran:
		if d.tran {
			d.iqPrev = d.iq + d.geq*(vd-d.vd)
		} else {
			d.iqPrev = 0
		}
	default:
		d.iqPrev = 0
		d.hist.reset()
	}
	d.qPrev = q
	d.hist.push(t, q)
}

// StoredEnergy is ∫₀^vd v·C(v) dv = vd·q(vd) − ∫₀^vd q(v) dv (integration by
// parts, with q(0) = 0), at the junction voltage of the most recent
// Operating call. Both terms are in closed form for every region of the
// charge and current equations, so it is exact and cheap enough to report
// at every time point.
func (d *Diode) StoredEnergy(OperatingPoint) float64 {
	if !d.hasCharge {
		return 0
	}
	vd := d.vdSol
	id, gd := d.current(vd)
	q, _ := d.charge(vd, id, gd)
	return vd*q - d.intDepletion(vd) - d.P.TT*d.intCurrent(vd)
}

// intDepletion is ∫₀^v qj(u) du for the depletion charge.
func (d *Diode) intDepletion(v float64) float64 {
	p := &d.P
	if p.CJO == 0 {
		return 0
	}
	fcv := p.FC * p.VJ
	// Below FC·VJ: qj = CJO·VJ·(1 − (1 − u/VJ)^(1−M))/(1 − M).
	below := func(u float64) float64 {
		return p.CJO * p.VJ / (1 - p.M) * (u + p.VJ/(2-p.M)*(math.Pow(1-u/p.VJ, 2-p.M)-1))
	}
	if v < fcv {
		return below(v)
	}
	f1 := p.VJ * (1 - math.Pow(1-p.FC, 1-p.M)) / (1 - p.M)
	f2 := math.Pow(1-p.FC, 1+p.M)
	f3 := 1 - p.FC*(1+p.M)
	dv := v - fcv
	return below(fcv) + p.CJO*f1*dv +
		p.CJO/f2*(f3*dv*dv/2+p.M/(2*p.VJ)*((v*v*v-fcv*fcv*fcv)/3-fcv*fcv*dv))
}

// intCurrent is ∫₀^v id(u) du for the junction current (without gmin),
// region by region.
func (d *Diode) intCurrent(v float64) float64 {
	p, vte := &d.P, d.vte
	fwd := func(u float64) float64 { return p.IS * (vte*math.Expm1(math.Min(u/vte, 700)) - u) }
	v1 := -3 * vte
	if v >= v1 {
		return fwd(v)
	}
	total := fwd(v1)
	// Reverse: id = −IS(1 + a³/u³), a = 3·vte/e; ∫ = −IS(u − a³/(2u²)).
	a3 := math.Pow(3*vte/math.E, 3)
	rev := func(u float64) float64 { return -p.IS * (u - a3/(2*u*u)) }
	v2 := v
	if p.HasBV && v < -d.xbv {
		v2 = -d.xbv
	}
	total += rev(v2) - rev(v1)
	if v2 == v {
		return total
	}
	// Breakdown: id = −IS·e^(−(XBV+u)/vte); ∫ = IS·vte·e^(−(XBV+u)/vte).
	brk := func(u float64) float64 { return p.IS * vte * math.Exp(math.Min(-(d.xbv+u)/vte, 700)) }
	return total + brk(v) - brk(v2)
}

// Truncate estimates the step from the junction charge history.
func (d *Diode) Truncate(x []float64, tc *TruncContext) float64 {
	if !d.hasCharge {
		return inf
	}
	vd := at(x, d.ai) - at(x, d.u[1])
	id, gd := d.current(vd)
	q, _ := d.charge(vd, id, gd)
	ic := d.iq + d.geq*(vd-d.vd)
	return lteStep(&d.hist, q, ic, d.iqPrev, tc.AbsTol, tc)
}

// ExteriorPower: dissipated power (absorbed minus stored) leaves as heat.
func (d *Diode) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -(op.PowerAbsorbed - op.StoragePower)
	}
	return 0
}

// LoadAC stamps the small-signal model at the accepted operating point:
// gd + jωC across the junction, plus RS.
func (d *Diode) LoadAC(ctx *ACContext) error {
	a, k, ai := d.u[0], d.u[1], d.ai
	id, gd := d.current(d.vdAcc)
	_, c := d.charge(d.vdAcc, id, gd)
	y := complex(gd+d.gmin, ctx.Omega*c)
	stampY(ctx.M, ai, k, y)
	if d.grs > 0 {
		stampY(ctx.M, a, ai, complex(d.grs, 0))
	}
	d.omega = ctx.Omega
	return nil
}

// OperatingAC returns the diode's phasors.
func (d *Diode) OperatingAC(x []complex128) ACPoint {
	v := d.voltageAC(x)
	var vj complex128
	if d.ai >= 0 && d.ai < len(x) {
		vj = x[d.ai]
	}
	if d.u[1] >= 0 {
		vj -= x[d.u[1]]
	}
	id, gd := d.current(d.vdAcc)
	_, c := d.charge(d.vdAcc, id, gd)
	return ACPoint{Voltage: v, Current: vj * complex(gd+d.gmin, d.omega*c)}
}
