// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// BJTParams are the SPICE3 Gummel–Poon parameters at the nominal
// temperature, already scaled by the area factor. Zero for VAF, VAR, IKF,
// IKR, IRB and VTF means infinite (the effect is off), as in SPICE.
type BJTParams struct {
	IS, BF, NF, VAF, IKF, ISE, NE float64
	BR, NR, VAR, IKR, ISC, NC     float64
	RB, IRB, RBM, RE, RC          float64
	CJE, VJE, MJE, TF, XTF, VTF   float64
	ITF, CJC, VJC, MJC, XCJC, TR  float64
	CJS, VJS, MJS, FC             float64
	ISS, NS                       float64 // substrate junction saturation current, emission coefficient
	// SUBS: +1 vertical (substrate junction to the collector), −1 lateral
	// (to the base). 0 means the ngspice default: vertical for NPN,
	// lateral for PNP.
	SUBS float64
	Type float64 // +1 NPN, −1 PNP
}

// DefaultBJTParams returns the SPICE3 defaults (SPICE3f3 User's Manual §4.3).
func DefaultBJTParams() BJTParams {
	return BJTParams{IS: 1e-16, BF: 100, NF: 1, NE: 1.5, BR: 1, NR: 1, NC: 2,
		VJE: 0.75, MJE: 0.33, VJC: 0.75, MJC: 0.33, XCJC: 1, VJS: 0.75, FC: 0.5, NS: 1, Type: 1}
}

// BJT is the SPICE3 Gummel–Poon bipolar junction transistor.
//
// Maxwell derivation (CLAUDE.md §3c): as for the diode, the junction and
// transport currents are semiconductor constitutive relations for the
// conduction current 𝔎 in the quasi-static limit, and the depletion and
// diffusion charges add displacement-current terms ∂q/∂t.
//
// Equations (from the published model, not ported code: Gummel & Poon,
// Bell Syst. Tech. J. 49, 1970; Nagel 1975 ch. 4; Massobrio & Antognetti,
// Semiconductor Device Modeling with SPICE, 2nd ed., ch. 2; the SPICE3 /
// ngspice bjtload equations for the region formulas and limiting). With
// vt = kT/q, vbe and vbc the internal junction voltages times the polarity
// (+1 NPN, −1 PNP):
//
//	ibf = IS·(e^(vbe/NF·vt) − 1),  ibr = IS·(e^(vbc/NR·vt) − 1)
//	  (below −3·N·vt the cubic reverse tail of the diode is used)
//	ile = ISE·(e^(vbe/NE·vt) − 1) + gmin·vbe,  ilc = ISC·(…NC…) + gmin·vbc
//	q1  = 1 / (1 − vbc/VAF − vbe/VAR)
//	q2  = ibf/IKF + ibr/IKR
//	qb  = q1·(1 + √(1 + 4·q2))/2
//	transport  it = (ibf − ibr)/qb               (collector' → emitter')
//	base       ib = ibf/BF + ile + ibr/BR + ilc
//	collector  ic = it − ibr/BR − ilc
//
// Base resistance falls from RB to RBM with current: rbb = RBM + (RB −
// RBM)/qb, or with IRB the crowding formula 3(RB − RBM)(tan z − z)/(z·tan²z).
// Charges: qbe = TF·ibf·(1 + XTF·(ibf/(ibf+ITF))²·e^(vbc/1.44·VTF))/qb plus
// the CJE depletion charge; qbc = TR·ibr plus XCJC of the CJC depletion
// charge on the internal base, the rest (1 − XCJC) from the external base;
// the substrate charge on CJS. Depletion charges are continued linearly
// above FC·VJ, as for the diode.
//
// Newton limiting: pnjlim on vbe and vbc (vcrit = vt·ln(vt/(√2·IS))).
//
// Terminals: Nodes = {collector, emitter} (the reported port: V_CE and the
// collector current), Extra = {base[, substrate]}; a missing substrate is
// ground. Ports: thermal out carries the dissipated power (junctions,
// transport current and the series resistances); stored junction energy is
// internal. Excess phase (PTF) is not implemented and is rejected at
// validation.
type BJT struct {
	Base
	P                            BJTParams
	vt                           float64
	vcrit                        float64
	rbpr                         float64 // RBM
	rbpi                         float64 // RB − RBM
	gcpr                         float64 // 1/RC
	gepr                         float64 // 1/RE
	ovtf                         float64 // 1/(1.44·VTF)
	f1e, f2e, f3e, f1c, f2c, f3c float64

	// Unknowns: external c, e (Base.u), b, s; internal cp, bp, ep.
	b, s, cp, bp, ep int
	// The substrate junction runs from sa to sb: (s, cp) vertical,
	// (bp, s) lateral, so vcs = T·(V(sa) − V(sb)) in both cases.
	sa, sb     int
	isub, gsub float64 // substrate junction DC current sa→sb and slope at vcs

	// Linearization of the most recent Load (or of AC/Accept evaluations).
	vbe, vbc, vbx, vcs float64
	dc                 bjtDC
	ch                 bjtCharge
	ag                 float64 // integration coefficient (2/h or 1/h), 0 outside a transient step
	iq                 [4]float64
	limited, tran      bool
	gmin               float64

	// Accepted state.
	vbeAcc, vbcAcc float64
	vbxAcc, vcsAcc float64
	qPrev, iqPrev  [4]float64
	hist           [4]history
	hasCharge      bool
	vSol           [4]float64 // vbe, vbc, vbx, vcs of the most recent Operating call
	omega          float64
	ac             bjtAC
}

// Charge indices.
const (
	qBE = iota
	qBC
	qBX
	qCS
)

// bjtDC holds the DC currents and their derivatives at (vbe, vbc).
type bjtDC struct {
	cbe, gbe, cbc, gbc float64 // ideal diode currents ibf, ibr and slopes
	it, itBE, itBC     float64 // transport current and its derivatives
	ibe, gpi           float64 // base–emitter base current and slope
	ibc, gmu           float64 // base–collector base current and slope
	qb, dqbdve, dqbdvc float64
	gx                 float64 // base-resistance conductance
}

// bjtCharge holds the junction charges and capacitances.
type bjtCharge struct {
	q        [4]float64
	c        [4]float64 // dq/dv of each charge on its own voltage
	dqbeDvbc float64    // cross term of qbe on vbc
}

// NewBJT builds a transistor from its element and model.
func NewBJT(e circuit.Element, m circuit.Model) *BJT {
	p := DefaultBJTParams()
	if m.Type == "pnp" {
		p.Type = -1
	}
	get := func(dst *float64, keys ...string) {
		for _, k := range keys {
			if v, ok := m.Params[k]; ok {
				*dst = v
				return
			}
		}
	}
	get(&p.IS, "is")
	get(&p.BF, "bf")
	get(&p.NF, "nf")
	get(&p.VAF, "vaf", "va")
	get(&p.IKF, "ikf", "ik")
	get(&p.ISE, "ise")
	get(&p.NE, "ne")
	get(&p.BR, "br")
	get(&p.NR, "nr")
	get(&p.VAR, "var", "vb")
	get(&p.IKR, "ikr")
	get(&p.ISC, "isc")
	get(&p.NC, "nc")
	get(&p.RB, "rb")
	get(&p.IRB, "irb")
	p.RBM = p.RB
	get(&p.RBM, "rbm")
	get(&p.RE, "re")
	get(&p.RC, "rc")
	get(&p.CJE, "cje")
	get(&p.VJE, "vje", "pe")
	get(&p.MJE, "mje", "me")
	get(&p.TF, "tf")
	get(&p.XTF, "xtf")
	get(&p.VTF, "vtf")
	get(&p.ITF, "itf")
	get(&p.CJC, "cjc")
	get(&p.VJC, "vjc", "pc")
	get(&p.MJC, "mjc", "mc")
	get(&p.XCJC, "xcjc")
	get(&p.TR, "tr")
	get(&p.CJS, "cjs", "ccs")
	get(&p.VJS, "vjs", "ps")
	get(&p.MJS, "mjs", "ms")
	get(&p.FC, "fc")
	get(&p.ISS, "iss")
	get(&p.NS, "ns")
	get(&p.SUBS, "subs")
	if p.SUBS == 0 {
		// ngspice's default (bjtsetup.c): vertical NPN, lateral PNP. Checked
		// against ngspice-39: a PNP at cutoff draws gmin·V_B through the
		// substrate junction to its base (testdata/netlists/bjt_pnp_dc.cir,
		// i(vin) at 9.75 V: −1.7394e-11 A, of which 9.75e-12 is this term).
		p.SUBS = 1
		if p.Type < 0 {
			p.SUBS = -1
		}
	}
	if p.RBM > p.RB {
		p.RBM = p.RB // SPICE: RBM ≤ RB
	}
	area := e.Area
	if area == 0 {
		area = 1
	}
	for _, x := range []*float64{&p.IS, &p.ISE, &p.ISC, &p.IKF, &p.IKR, &p.IRB, &p.ITF, &p.CJE, &p.CJC, &p.CJS, &p.ISS} {
		*x *= area
	}
	for _, x := range []*float64{&p.RB, &p.RBM, &p.RE, &p.RC} {
		*x /= area
	}
	q := &BJT{Base: newBase(e), P: p}
	q.b = e.Extra[0]
	q.s = 0
	if len(e.Extra) > 1 {
		q.s = e.Extra[1]
	}
	q.vt = BoltzmanK * TNomK / ChargeQ
	q.vcrit = q.vt * math.Log(q.vt/(math.Sqrt2*p.IS))
	q.rbpr, q.rbpi = p.RBM, p.RB-p.RBM
	if p.RC > 0 {
		q.gcpr = 1 / p.RC
	}
	if p.RE > 0 {
		q.gepr = 1 / p.RE
	}
	if p.VTF > 0 {
		q.ovtf = 1 / (1.44 * p.VTF)
	}
	q.f1e, q.f2e, q.f3e = deplF(p.VJE, p.MJE, p.FC)
	q.f1c, q.f2c, q.f3c = deplF(p.VJC, p.MJC, p.FC)
	q.hasCharge = p.CJE > 0 || p.CJC > 0 || p.CJS > 0 || p.TF > 0 || p.TR > 0
	_ = q.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores junction charge")
	_ = q.SetIntrinsic(exchange.Thermal, exchange.Out, "junction, transport and series-resistance heating to surroundings")
	return q
}

// deplF returns SPICE's F1, F2, F3 for the linear continuation of a
// depletion charge above FC·VJ.
func deplF(vj, m, fc float64) (f1, f2, f3 float64) {
	return vj * (1 - math.Pow(1-fc, 1-m)) / (1 - m), math.Pow(1-fc, 1+m), 1 - fc*(1+m)
}

// Setup maps terminals and allocates internal nodes behind RC, RB, RE.
func (q *BJT) Setup(ctx *SetupContext) error {
	q.setupNodes(ctx)
	q.b, q.s = ctx.Unknown(q.b), ctx.Unknown(q.s)
	q.cp, q.bp, q.ep = q.u[0], q.b, q.u[1]
	if q.P.RC > 0 {
		q.cp = ctx.AllocUnknown()
	}
	if q.P.RB > 0 {
		q.bp = ctx.AllocUnknown()
	}
	if q.P.RE > 0 {
		q.ep = ctx.AllocUnknown()
	}
	q.sa, q.sb = q.s, q.cp
	if q.P.SUBS < 0 {
		q.sa, q.sb = q.bp, q.s
	}
	return nil
}

// junction returns an ideal diode current IS(e^(v/nvt) − 1) and its slope,
// with SPICE's cubic reverse tail below −3·nvt.
func junction(is, v, nvt float64) (i, g float64) {
	if is == 0 {
		return 0, 0
	}
	if v >= -3*nvt {
		e := math.Exp(math.Min(v/nvt, 700))
		return is * (e - 1), is * e / nvt
	}
	a := 3 * nvt / (v * math.E)
	a = a * a * a
	return -is * (1 + a), is * 3 * a / v
}

// substrate returns the substrate junction current and slope at vcs:
// ISS·(e^(vcs/NS·vt) − 1) + gmin·vcs, with the cubic reverse tail. ngspice
// (bjtload.c) has this branch; SPICE3f5 had only the substrate charge.
func (q *BJT) substrate(vcs float64) (i, g float64) {
	i, g = junction(q.P.ISS, vcs, q.P.NS*q.vt)
	return i + q.gmin*vcs, g + q.gmin
}

// evalDC computes the DC model at internal junction voltages vbe, vbc.
func (q *BJT) evalDC(vbe, vbc float64) bjtDC {
	p := &q.P
	var d bjtDC
	d.cbe, d.gbe = junction(p.IS, vbe, p.NF*q.vt)
	cben, gben := junction(p.ISE, vbe, p.NE*q.vt)
	cben += q.gmin * vbe
	gben += q.gmin
	d.cbc, d.gbc = junction(p.IS, vbc, p.NR*q.vt)
	cbcn, gbcn := junction(p.ISC, vbc, p.NC*q.vt)
	cbcn += q.gmin * vbc
	gbcn += q.gmin

	var ivaf, ivar, oik, oikr float64
	if p.VAF > 0 {
		ivaf = 1 / p.VAF
	}
	if p.VAR > 0 {
		ivar = 1 / p.VAR
	}
	if p.IKF > 0 {
		oik = 1 / p.IKF
	}
	if p.IKR > 0 {
		oikr = 1 / p.IKR
	}
	q1 := 1 / (1 - ivaf*vbc - ivar*vbe)
	if oik == 0 && oikr == 0 {
		d.qb = q1
		d.dqbdve = q1 * d.qb * ivar
		d.dqbdvc = q1 * d.qb * ivaf
	} else {
		q2 := oik*d.cbe + oikr*d.cbc
		arg := math.Max(0, 1+4*q2)
		sq := 1.0
		if arg != 0 {
			sq = math.Sqrt(arg)
		}
		d.qb = q1 * (1 + sq) / 2
		d.dqbdve = q1 * (d.qb*ivar + oik*d.gbe/sq)
		d.dqbdvc = q1 * (d.qb*ivaf + oikr*d.gbc/sq)
	}
	d.it = (d.cbe - d.cbc) / d.qb
	d.itBE = (d.gbe - d.it*d.dqbdve) / d.qb
	d.itBC = (-d.gbc - d.it*d.dqbdvc) / d.qb
	d.ibe, d.gpi = d.cbe/p.BF+cben, d.gbe/p.BF+gben
	d.ibc, d.gmu = d.cbc/p.BR+cbcn, d.gbc/p.BR+gbcn

	gx := q.rbpr + q.rbpi/d.qb
	if p.IRB > 0 && q.rbpi != 0 {
		cb := d.ibe + d.ibc
		a1 := math.Max(cb/p.IRB, 1e-9)
		z := (-1 + math.Sqrt(1+14.59025*a1)) / 2.4317 / math.Sqrt(a1)
		tz := math.Tan(z)
		gx = q.rbpr + 3*q.rbpi*(tz-z)/z/tz/tz
	}
	if gx != 0 {
		gx = 1 / gx
	}
	d.gx = gx
	return d
}

// depletion returns a junction depletion charge and capacitance with
// SPICE's linear continuation above FC·VJ.
func depletion(v, cz, vj, m, fc, f1, f2, f3 float64) (qv, c float64) {
	if cz == 0 {
		return 0, 0
	}
	fcv := fc * vj
	if v < fcv {
		arg := 1 - v/vj
		sarg := math.Exp(-m * math.Log(arg))
		return vj * cz * (1 - arg*sarg) / (1 - m), cz * sarg
	}
	czf2 := cz / f2
	return cz*f1 + czf2*(f3*(v-fcv)+(m/(vj+vj))*(v*v-fcv*fcv)), czf2 * (f3 + m*v/vj)
}

// charges computes the four junction charges at the given voltages.
func (q *BJT) charges(vbe, vbc, vbx, vcs float64, d *bjtDC) bjtCharge {
	p := &q.P
	var ch bjtCharge
	cbe, gbe := d.cbe, d.gbe
	if p.TF != 0 && vbe > 0 {
		argtf, arg2, arg3 := 0.0, 0.0, 0.0
		if p.XTF != 0 {
			argtf = p.XTF
			if q.ovtf != 0 {
				argtf *= math.Exp(vbc * q.ovtf)
			}
			arg2 = argtf
			if p.ITF != 0 {
				t := cbe / (cbe + p.ITF)
				argtf *= t * t
				arg2 = argtf * (3 - t - t)
			}
			arg3 = cbe * argtf * q.ovtf
		}
		cbe = cbe * (1 + argtf) / d.qb
		gbe = (gbe*(1+arg2) - cbe*d.dqbdve) / d.qb
		ch.dqbeDvbc = p.TF * (arg3 - cbe*d.dqbdvc) / d.qb
	}
	qd, cd := depletion(vbe, p.CJE, p.VJE, p.MJE, p.FC, q.f1e, q.f2e, q.f3e)
	ch.q[qBE], ch.c[qBE] = p.TF*cbe+qd, p.TF*gbe+cd
	czbc := p.CJC * p.XCJC
	qd, cd = depletion(vbc, czbc, p.VJC, p.MJC, p.FC, q.f1c, q.f2c, q.f3c)
	ch.q[qBC], ch.c[qBC] = p.TR*d.cbc+qd, p.TR*d.gbc+cd
	ch.q[qBX], ch.c[qBX] = depletion(vbx, p.CJC-czbc, p.VJC, p.MJC, p.FC, q.f1c, q.f2c, q.f3c)
	if cz := p.CJS; cz != 0 {
		if vcs < 0 {
			arg := 1 - vcs/p.VJS
			sarg := math.Exp(-p.MJS * math.Log(arg))
			ch.q[qCS], ch.c[qCS] = p.VJS*cz*(1-arg*sarg)/(1-p.MJS), cz*sarg
		} else {
			ch.q[qCS], ch.c[qCS] = vcs*cz*(1+p.MJS*vcs/(2*p.VJS)), cz*(1+p.MJS*vcs/p.VJS)
		}
	}
	return ch
}

// voltages returns the polarity-corrected vbe, vbc, vbx, vcs at x.
func (q *BJT) voltages(x []float64) (vbe, vbc, vbx, vcs float64) {
	t := q.P.Type
	vb, vbp, vcp, vep := at(x, q.b), at(x, q.bp), at(x, q.cp), at(x, q.ep)
	return t * (vbp - vep), t * (vbp - vcp), t * (vb - vcp), t * (at(x, q.sa) - at(x, q.sb))
}

// Load linearizes the transistor at the current iterate and stamps it.
func (q *BJT) Load(ctx *LoadContext) error {
	q.gmin = ctx.Gmin
	vbe, vbc, vbx, vcs := q.voltages(ctx.X)
	q.limited = false
	if ctx.InitJct {
		vbe, vbc, q.limited = q.vcrit, 0, true
		vbx, vcs = 0, 0
	} else {
		var l1, l2 bool
		vbe, l1 = pnjlim(vbe, q.vbeAcc, q.vt, q.vcrit)
		vbc, l2 = pnjlim(vbc, q.vbcAcc, q.vt, q.vcrit)
		q.limited = l1 || l2
	}
	q.vbeAcc, q.vbcAcc = vbe, vbc // the next iteration limits relative to this one
	q.vbe, q.vbc, q.vbx, q.vcs = vbe, vbc, vbx, vcs
	q.dc = q.evalDC(vbe, vbc)
	q.isub, q.gsub = q.substrate(vcs)
	q.tran, q.ag = false, 0
	q.iq = [4]float64{}
	if ctx.Mode == ModeTran && q.hasCharge {
		q.ch = q.charges(vbe, vbc, vbx, vcs, &q.dc)
		q.tran = true
		if ctx.Method == Trapezoidal {
			q.ag = 2 / ctx.H
			for k := range q.iq {
				q.iq[k] = q.ag*(q.ch.q[k]-q.qPrev[k]) - q.iqPrev[k]
			}
		} else {
			q.ag = 1 / ctx.H
			for k := range q.iq {
				q.iq[k] = q.ag * (q.ch.q[k] - q.qPrev[k])
			}
		}
	}
	q.stamp(ctx)
	return nil
}

// branch stamps a current i0 + g1·Δv1 + g2·Δv2 flowing from node a through
// the device to node b, where v1 = T·(V(p1) − V(n1)) and likewise v2, and
// Δv is measured from the linearization point (v10, v20).
func (q *BJT) branch(ctx *LoadContext, a, b int, i0, g1 float64, p1, n1 int, v10, g2 float64, p2, n2 int, v20 float64) {
	m := ctx.M
	if g1 != 0 {
		m.Add(a, p1, g1)
		m.Add(a, n1, -g1)
		m.Add(b, p1, -g1)
		m.Add(b, n1, g1)
	}
	if g2 != 0 {
		m.Add(a, p2, g2)
		m.Add(a, n2, -g2)
		m.Add(b, p2, -g2)
		m.Add(b, n2, g2)
	}
	ieq := q.P.Type * (i0 - g1*v10 - g2*v20)
	ctx.AddRHS(a, -ieq)
	ctx.AddRHS(b, ieq)
}

func stampG(m Stamper, a, b int, g float64) {
	if g == 0 {
		return
	}
	m.Add(a, a, g)
	m.Add(b, b, g)
	m.Add(a, b, -g)
	m.Add(b, a, -g)
}

func (q *BJT) stamp(ctx *LoadContext) {
	d := &q.dc
	cp, bp, ep := q.cp, q.bp, q.ep
	q.branch(ctx, cp, ep, d.it, d.itBE, bp, ep, q.vbe, d.itBC, bp, cp, q.vbc)
	q.branch(ctx, bp, ep, d.ibe, d.gpi, bp, ep, q.vbe, 0, 0, 0, 0)
	q.branch(ctx, bp, cp, d.ibc, d.gmu, bp, cp, q.vbc, 0, 0, 0, 0)
	q.branch(ctx, q.sa, q.sb, q.isub, q.gsub, q.sa, q.sb, q.vcs, 0, 0, 0, 0)
	if q.tran {
		ag, ch := q.ag, &q.ch
		q.branch(ctx, bp, ep, q.iq[qBE], ag*ch.c[qBE], bp, ep, q.vbe, ag*ch.dqbeDvbc, bp, cp, q.vbc)
		q.branch(ctx, bp, cp, q.iq[qBC], ag*ch.c[qBC], bp, cp, q.vbc, 0, 0, 0, 0)
		q.branch(ctx, q.b, cp, q.iq[qBX], ag*ch.c[qBX], q.b, cp, q.vbx, 0, 0, 0, 0)
		q.branch(ctx, q.sa, q.sb, q.iq[qCS], ag*ch.c[qCS], q.sa, q.sb, q.vcs, 0, 0, 0, 0)
	}
	if q.bp != q.b {
		stampG(ctx.M, q.b, q.bp, d.gx)
	}
	if q.cp != q.u[0] {
		stampG(ctx.M, q.u[0], q.cp, q.gcpr)
	}
	if q.ep != q.u[1] {
		stampG(ctx.M, q.u[1], q.ep, q.gepr)
	}
}

// Limited reports whether the last Load limited a junction voltage.
func (q *BJT) Limited() bool { return q.limited }

// ConvTest compares the collector and base currents predicted by the
// linearization with the model at x (SPICE3 BJTconvTest).
func (q *BJT) ConvTest(x []float64, reltol, abstol float64) bool {
	vbe, vbc, _, _ := q.voltages(x)
	d := &q.dc
	dbe, dbc := vbe-q.vbe, vbc-q.vbc
	ichat := d.it - d.ibc + (d.itBE)*dbe + (d.itBC-d.gmu)*dbc
	ibhat := d.ibe + d.ibc + d.gpi*dbe + d.gmu*dbc
	n := q.evalDC(vbe, vbc)
	ic, ib := n.it-n.ibc, n.ibe+n.ibc
	return math.Abs(ichat-ic) <= reltol*math.Max(math.Abs(ichat), math.Abs(ic))+abstol &&
		math.Abs(ibhat-ib) <= reltol*math.Max(math.Abs(ibhat), math.Abs(ib))+abstol
}

// terminalCurrents returns the currents into the device at the collector,
// base, emitter and substrate terminals, from the linearization the solve
// used (so KCL holds exactly), and the power going into stored charge.
func (q *BJT) terminalCurrents(x []float64) (ic, ib, ie, is, pStore float64) {
	vbe, vbc, vbx, vcs := q.voltages(x)
	t := q.P.Type
	d := &q.dc
	dbe, dbc := vbe-q.vbe, vbc-q.vbc
	it := t * (d.it + d.itBE*dbe + d.itBC*dbc)
	ibe := t * (d.ibe + d.gpi*dbe)
	ibc := t * (d.ibc + d.gmu*dbc)
	var iqbe, iqbc, iqbx, iqcs float64
	if q.tran {
		ag, ch := q.ag, &q.ch
		iqbe = t * (q.iq[qBE] + ag*ch.c[qBE]*dbe + ag*ch.dqbeDvbc*dbc)
		iqbc = t * (q.iq[qBC] + ag*ch.c[qBC]*dbc)
		iqbx = t * (q.iq[qBX] + ag*ch.c[qBX]*(vbx-q.vbx))
		iqcs = t * (q.iq[qCS] + ag*ch.c[qCS]*(vcs-q.vcs))
		// Physical voltage × physical current: T² = 1.
		pStore = vbe*t*iqbe + vbc*t*iqbc + vbx*t*iqbx + vcs*t*iqcs
	}
	// Substrate junction: DC (gmin and ISS) plus charge, from sa to sb.
	isb := t*(q.isub+q.gsub*(vcs-q.vcs)) + iqcs
	ic = it - ibc - iqbc - iqbx
	ib = ibe + ibc + iqbe + iqbc + iqbx
	ie = -it - ibe - iqbe
	if q.P.SUBS > 0 {
		is, ic = isb, ic-isb
	} else {
		ib, is = ib+isb, -isb
	}
	q.vSol = [4]float64{vbe, vbc, vbx, vcs}
	return
}

// Operating reports V_CE, the collector current, and the total absorbed
// power Σ V_k·I_k over the terminals.
func (q *BJT) Operating(x []float64) OperatingPoint {
	ic, ib, ie, is, ps := q.terminalCurrents(x)
	vc, ve := at(x, q.u[0]), at(x, q.u[1])
	p := vc*ic + at(x, q.b)*ib + ve*ie + at(x, q.s)*is
	return OperatingPoint{Voltage: vc - ve, Current: ic, PowerAbsorbed: p, StoragePower: ps}
}

// Accept commits the converged junction state.
func (q *BJT) Accept(x []float64, t float64, mode Mode) {
	vbe, vbc, vbx, vcs := q.voltages(x)
	q.vbeAcc, q.vbcAcc, q.vbxAcc, q.vcsAcc = vbe, vbc, vbx, vcs
	d := q.evalDC(vbe, vbc)
	ch := q.charges(vbe, vbc, vbx, vcs, &d)
	for k := range q.qPrev {
		if mode == ModeTran && q.tran {
			dv := [4]float64{vbe - q.vbe, vbc - q.vbc, vbx - q.vbx, vcs - q.vcs}[k]
			i := q.iq[k] + q.ag*q.ch.c[k]*dv
			if k == qBE {
				i += q.ag * q.ch.dqbeDvbc * (vbc - q.vbc)
			}
			q.iqPrev[k] = i
		} else {
			q.iqPrev[k] = 0
		}
		if mode != ModeTran {
			q.hist[k].reset()
		}
		q.qPrev[k] = ch.q[k]
		q.hist[k].push(t, ch.q[k])
	}
}

// gl16 holds 16-point Gauss–Legendre nodes and weights on [−1, 1]
// (Abramowitz & Stegun, Table 25.4).
var gl16x = [8]float64{0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438,
	0.7554044083550030, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499}
var gl16w = [8]float64{0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
	0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541}

// integrate returns ∫ₐᵇ f by 16-point Gauss–Legendre on n equal panels.
func integrate(f func(float64) float64, a, b float64, n int) float64 {
	h := (b - a) / float64(n)
	s := 0.0
	for p := 0; p < n; p++ {
		mid, half := a+(float64(p)+0.5)*h, h/2
		for k := 0; k < 8; k++ {
			s += gl16w[k] * (f(mid-half*gl16x[k]) + f(mid+half*gl16x[k]))
		}
	}
	return s * h / 2
}

// StoredEnergy is the co-energy Σ (v·q(v) − ∫₀^v q(u) du) of the four
// charges at the junction voltages of the most recent Operating call, each
// integrated along its own voltage with the others held (quadrature on 8
// panels). For the charges that depend only on their own voltage this is
// the exact stored energy; the forward transit charge also depends on vbc
// (through qb and VTF), so its share is exact only where that dependence is
// weak — a mismatch shows in the transient ledger's integration mismatch.
func (q *BJT) StoredEnergy(OperatingPoint) float64 {
	if !q.hasCharge {
		return 0
	}
	v := q.vSol
	w := 0.0
	for k := 0; k < 4; k++ {
		qk := func(u float64) float64 {
			vv := v
			vv[k] = u
			d := q.evalDC(vv[qBE], vv[qBC])
			return q.charges(vv[qBE], vv[qBC], vv[qBX], vv[qCS], &d).q[k]
		}
		if v[k] == 0 {
			continue
		}
		w += v[k]*qk(v[k]) - integrate(qk, 0, v[k], 8)
	}
	return w
}

// Truncate proposes a step from the charges' truncation error.
func (q *BJT) Truncate(x []float64, tc *TruncContext) float64 {
	if !q.hasCharge || !q.tran {
		return inf
	}
	vbe, vbc, vbx, vcs := q.voltages(x)
	d := q.evalDC(vbe, vbc)
	ch := q.charges(vbe, vbc, vbx, vcs, &d)
	dv := [4]float64{vbe - q.vbe, vbc - q.vbc, vbx - q.vbx, vcs - q.vcs}
	h := inf
	for k := 0; k < 4; k++ {
		if ch.q[k] == 0 && q.qPrev[k] == 0 {
			continue
		}
		i := q.iq[k] + q.ag*q.ch.c[k]*dv[k]
		if k == qBE {
			i += q.ag * q.ch.dqbeDvbc * dv[qBC]
		}
		h = math.Min(h, lteStep(&q.hist[k], ch.q[k], i, q.iqPrev[k], tc.AbsTol, tc))
	}
	return h
}

// ExteriorPower: dissipated power (absorbed minus stored) leaves as heat.
func (q *BJT) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -(op.PowerAbsorbed - op.StoragePower)
	}
	return 0
}

// bjtAC holds the small-signal admittances at the operating point.
type bjtAC struct {
	gx, itBE, itBC, gpi, gmu  complex128
	ybe, ybeBC, ybc, ybx, ycs complex128
}

// LoadAC stamps the small-signal model at the accepted operating point.
func (q *BJT) LoadAC(ctx *ACContext) error {
	vbe, vbc := q.vbeAcc, q.vbcAcc
	vbx, vcs := q.vbxAcc, q.vcsAcc
	d := q.evalDC(vbe, vbc)
	ch := q.charges(vbe, vbc, vbx, vcs, &d)
	w := ctx.Omega
	q.omega = w
	a := &q.ac
	a.itBE, a.itBC = complex(d.itBE, 0), complex(d.itBC, 0)
	a.gpi, a.gmu, a.gx = complex(d.gpi, 0), complex(d.gmu, 0), complex(d.gx, 0)
	a.ybe, a.ybeBC = complex(0, w*ch.c[qBE]), complex(0, w*ch.dqbeDvbc)
	_, gsub := q.substrate(vcs)
	a.ybc, a.ybx, a.ycs = complex(0, w*ch.c[qBC]), complex(0, w*ch.c[qBX]), complex(gsub, w*ch.c[qCS])
	cp, bp, ep := q.cp, q.bp, q.ep
	m := ctx.M
	acBranch(m, cp, ep, a.itBE, bp, ep, a.itBC, bp, cp)
	acBranch(m, bp, ep, a.gpi+a.ybe, bp, ep, a.ybeBC, bp, cp)
	acBranch(m, bp, cp, a.gmu+a.ybc, bp, cp, 0, 0, 0)
	acBranch(m, q.b, cp, a.ybx, q.b, cp, 0, 0, 0)
	acBranch(m, q.sa, q.sb, a.ycs, q.sa, q.sb, 0, 0, 0)
	if q.bp != q.b {
		stampY(m, q.b, q.bp, a.gx)
	}
	if q.cp != q.u[0] {
		stampY(m, q.u[0], q.cp, complex(q.gcpr, 0))
	}
	if q.ep != q.u[1] {
		stampY(m, q.u[1], q.ep, complex(q.gepr, 0))
	}
	return nil
}

func acBranch(m CStamper, a, b int, y1 complex128, p1, n1 int, y2 complex128, p2, n2 int) {
	if y1 != 0 {
		m.Add(a, p1, y1)
		m.Add(a, n1, -y1)
		m.Add(b, p1, -y1)
		m.Add(b, n1, y1)
	}
	if y2 != 0 {
		m.Add(a, p2, y2)
		m.Add(a, n2, -y2)
		m.Add(b, p2, -y2)
		m.Add(b, n2, y2)
	}
}

// acCurrents returns the terminal current phasors into the device.
func (q *BJT) acCurrents(x []complex128) (ic, ib, ie, is complex128) {
	a := &q.ac
	vbe := catC(x, q.bp) - catC(x, q.ep)
	vbc := catC(x, q.bp) - catC(x, q.cp)
	vbx := catC(x, q.b) - catC(x, q.cp)
	vcs := catC(x, q.sa) - catC(x, q.sb)
	it := a.itBE*vbe + a.itBC*vbc
	ibe := (a.gpi+a.ybe)*vbe + a.ybeBC*vbc
	ibc := (a.gmu + a.ybc) * vbc
	ix, isub := a.ybx*vbx, a.ycs*vcs
	if q.P.SUBS > 0 {
		return it - ibc - ix - isub, ibe + ibc + ix, -it - ibe, isub
	}
	return it - ibc - ix, ibe + ibc + ix + isub, -it - ibe, -isub
}

// OperatingAC returns the V_CE and collector-current phasors.
func (q *BJT) OperatingAC(x []complex128) ACPoint {
	ic, _, _, _ := q.acCurrents(x)
	return ACPoint{Voltage: q.voltageAC(x), Current: ic}
}

// ACPower is ½·Σ V_k·I_k* over the four terminals.
func (q *BJT) ACPower(x []complex128) complex128 {
	ic, ib, ie, is := q.acCurrents(x)
	s := catC(x, q.u[0])*cmplx.Conj(ic) + catC(x, q.b)*cmplx.Conj(ib) +
		catC(x, q.u[1])*cmplx.Conj(ie) + catC(x, q.s)*cmplx.Conj(is)
	return s / 2
}
