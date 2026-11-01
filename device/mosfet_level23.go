// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import "math"

// Variables of the dual numbers in the drain-current equations.
const (
	dVGS = iota
	dVDS
	dVBS
)

// idsLevel2 is the SPICE2 level 2 (Grove–Frohman) drain current:
// A. Vladimirescu, S. Liu, "The simulation of MOS integrated circuits using
// SPICE2," UCB/ERL M80/7, 1980, §2.3 (threshold with short- and
// narrow-channel corrections, bulk-charge current, Grove–Frohman
// saturation voltage, Baum–Beneking velocity-saturated vdsat, mobility
// degradation by the critical field, channel-length modulation and
// punch-through, weak inversion with NFS); region structure and constants
// as in the SPICE3f5 / ngspice mos2load (BSD licence).
//
// The current is evaluated in dual numbers, so gm, gds and gmbs are its
// exact derivatives (SPICE's hand-written chain sets some terms, such as
// dudvds, to zero). Two SPICE special cases are kept as they are: for
// vds ≤ 1e−10 the current is 0 with the conductance of the vds → 0 limit,
// and without NFS the current is cut off below vbin rather than von.
// Where the velocity-saturation quartic has no real root SPICE would return
// NaN; here the Grove–Frohman vdsat is kept.
func (q *MOSFET) idsLevel2(vgsV, vdsV, vbsV float64) (cdrain, gm, gds, gmbs, von, vdsat float64) {
	p, lv := &q.p, &q.lv
	G, D, B := variable(vgsV, dVGS), variable(vdsV, dVDS), variable(vbsV, dVBS)
	phi := p.PHI
	sphi := math.Sqrt(phi)
	sphi3 := phi * sphi
	phiMinVbs := B.neg().addc(phi)

	// sarg = √(φ − vbs), barg = √(φ − vbs + vds), continued for forward
	// bias; their vbs derivatives enter the weak-inversion slope xn.
	var sarg, dsrgdb dual
	if B.v <= 0 {
		sarg = phiMinVbs.sqrt()
		dsrgdb = cst(-0.5).div(sarg)
	} else {
		sarg = cst(sphi).div(B.scale(0.5 / phi).addc(1))
		dsrgdb = sarg.mul(sarg).scale(-0.5 / sphi3)
	}
	var barg, dbrgdb dual
	if B.v-D.v <= 0 {
		barg = phiMinVbs.add(D).sqrt()
		dbrgdb = cst(-0.5).div(barg)
	} else {
		barg = cst(sphi).div(B.sub(D).scale(0.5 / phi).addc(1))
		dbrgdb = barg.mul(barg).scale(-0.5 / sphi3)
	}
	factor := lv.factor
	eta := 1 + factor
	vbin := phiMinVbs.scale(factor).addc(q.vbi)

	// Threshold with the short-channel (junction-depth) correction of γ.
	gamasd, dgddvb := cst(p.GAMMA), cst(0)
	if p.GAMMA > 0 || lv.hasNSUB {
		if p.XJ > 0 {
			xws, xwd := sarg.scale(lv.xd), barg.scale(lv.xd)
			args := xws.scale(2 / p.XJ).addc(1).sqrt()
			argd := xwd.scale(2 / p.XJ).addc(1).sqrt()
			t := 0.5 * p.XJ / q.leff
			argss, argsd := args.addc(-1).scale(t), argd.addc(-1).scale(t)
			gamasd = argss.add(argsd).neg().addc(1).scale(p.GAMMA)
			dbargs := dsrgdb.scale(lv.xd).div(args).scale(0.5 / q.leff)
			dbargd := dbrgdb.scale(lv.xd).div(argd).scale(0.5 / q.leff)
			dgddvb = dbargs.add(dbargd).scale(-p.GAMMA)
		}
	}
	vonD := vbin.add(gamasd.mul(sarg))
	nfs := p.NFS != 0
	var xn, argg dual
	if nfs {
		cdonco := gamasd.mul(dsrgdb).add(dgddvb.mul(sarg)).neg().addc(factor)
		xn = cdonco.addc(1 + lv.cfs)
		vtxn := xn.scale(q.vt)
		vonD = vonD.add(vtxn)
		argg = cst(1).div(vtxn)
	} else if G.v <= vbin.v {
		return 0, 0, 0, 0, vonD.v, 0 // cutoff
	}
	vgst := G.sub(vonD)

	// Mobility degradation: ufact = (UCRIT·ε_Si / (Cox·vgst))^UEXP above the
	// critical field.
	ufact := cst(1)
	if crit := p.UCRIT * 100 * epsSi / q.cox; vgst.v > crit {
		ufact = cst(crit).div(vgst).powc(p.UEXP)
	}
	ueff := ufact.scale(p.UO * 1e-4)

	// Saturation voltage (Grove–Frohman).
	vgsx := G
	if nfs {
		vgsx = dmax(G, vonD)
	}
	gammad := gamasd.scale(1 / eta)
	vgb := vgsx.sub(vbin).scale(1 / eta) // (vgsx − vbin)/η
	vdsatD := cst(0)
	if gammad.v > 0 {
		argv := vgb.add(phiMinVbs)
		if argv.v > 0 {
			gammd2 := gammad.mul(gammad)
			arg1 := argv.scale(4).div(gammd2).addc(1).sqrt()
			vdsatD = dmax(vgb.add(gammd2.mul(arg1.neg().addc(1)).scale(0.5)), cst(0))
		}
	} else {
		vdsatD = dmax(vgb, cst(0))
	}
	sarg3 := sarg.mul(sarg).mul(sarg)
	if p.VMAX > 0 {
		if x, ok := level2VelocitySat(gammad, vgb.add(phiMinVbs), phiMinVbs,
			ueff.scale(1/(p.VMAX*q.leff)), sarg3); ok {
			vdsatD = x.mul(x).sub(phiMinVbs)
		}
	}

	// Channel-length modulation: LAMBDA if given, else from the depletion
	// width of the pinched-off region.
	xlamda := cst(p.LAMBDA)
	var bsarg, bodys dual
	if D.v != 0 {
		if B.v-vdsatD.v <= 0 {
			bsarg = vdsatD.add(phiMinVbs).sqrt()
		} else {
			bsarg = cst(sphi).div(B.sub(vdsatD).scale(0.5 / phi).addc(1))
		}
		bodys = bsarg.mul(bsarg).mul(bsarg).sub(sarg3)
		if lv.hasNSUB && !(p.LAMBDA > 0) {
			if p.VMAX <= 0 {
				argv := D.sub(vdsatD).scale(0.25)
				sargv := argv.mul(argv).addc(1).sqrt()
				xlamda = argv.add(sargv).sqrt().mul(cst(lv.xd).div(D.scale(q.leff)))
			} else {
				xdv := lv.xd / math.Sqrt(p.NEFF)
				xlv := cst(p.VMAX * xdv * 0.5).div(ueff)
				argv := dmax(D.sub(vdsatD), cst(0))
				xls := xlv.mul(xlv).add(argv).sqrt()
				xlamda = xls.sub(xlv).mul(cst(xdv).div(D.scale(q.leff)))
			}
		}
	}
	// Punch-through limit of the channel shortening.
	xwb := lv.xd * lv.sqrtPB
	xld := q.leff - xwb
	clfact := xlamda.mul(D).neg().addc(1)
	if !lv.hasNSUB {
		xwb = 0.25e-6
	}
	if xleff := clfact.scale(q.leff); xleff.v < xwb {
		deltal := xlamda.mul(D).scale(q.leff)
		xleff = cst(xwb).div(deltal.addc(-xld).scale(1 / xwb).addc(1))
		clfact = xleff.scale(1 / q.leff)
	}
	beta1 := ufact.scale(q.beta).div(clfact)

	if D.v <= 1e-10 {
		// SPICE: zero current, conductance of the vds → 0 limit.
		var g dual
		switch {
		case G.v > vonD.v:
			g = beta1.mul(G.sub(vbin).sub(gamasd.mul(sarg)))
		case nfs:
			g = beta1.mul(vonD.sub(vbin).sub(gamasd.mul(sarg))).mul(argg.mul(G.sub(vonD)).exp())
		}
		return 0, 0, g.v, 0, vonD.v, vdsatD.v
	}
	var id dual
	if nfs && G.v <= vonD.v {
		// Weak inversion: the current at von, falling exponentially.
		if vdsatD.v <= 0 {
			return 0, 0, 0, 0, vonD.v, vdsatD.v
		}
		vdson := dmin(vdsatD, D)
		body := barg.mul(barg).mul(barg).sub(sarg3)
		if D.v > vdsatD.v {
			body = bodys
		}
		cdson := beta1.mul(vonD.sub(vbin).sub(vdson.scale(0.5 * eta)).mul(vdson).sub(gamasd.mul(body).scale(1 / 1.5)))
		id = cdson.mul(argg.mul(G.sub(vonD)).exp())
	} else if D.v <= vdsatD.v {
		body := barg.mul(barg).mul(barg).sub(sarg3)
		id = beta1.mul(G.sub(vbin).sub(D.scale(0.5 * eta)).mul(D).sub(gamasd.mul(body).scale(1 / 1.5)))
	} else {
		id = beta1.mul(G.sub(vbin).sub(vdsatD.scale(0.5 * eta)).mul(vdsatD).sub(gamasd.mul(bodys).scale(1 / 1.5)))
	}
	return fromDual(id, vonD, vdsatD)
}

// level2VelocitySat solves the Baum–Beneking quartic for x = √(vdsat + φ −
// vbs) (Vladimirescu & Liu 1980, eq. 2.3.22; mos2load's closed-form
// Ferrari solution): x⁴ + a1·x³ + b1·x² + c1·x + d1 = 0 with v1 = (vgs −
// vbin)/η + φ − vbs, v2 = φ − vbs and xv = VMAX·Leff/µeff (passed as its
// reciprocal invxv). The smallest positive root with |poly| ≤ 1e−6 is
// returned; ok is false when there is none.
func level2VelocitySat(gammad, v1, v2, invxv, sarg3 dual) (x dual, ok bool) {
	xv := cst(1).div(invxv)
	a1 := gammad.scale(1 / 0.75)
	b1 := v1.add(xv).scale(-2)
	c1 := gammad.mul(xv).scale(-2)
	d1 := v1.mul(v2.add(xv)).scale(2).sub(v2.mul(v2)).sub(gammad.mul(sarg3).scale(4.0 / 3.0))
	a := b1.neg()
	b := a1.mul(c1).sub(d1.scale(4))
	c := d1.neg().mul(a1.mul(a1).sub(b1.scale(4))).sub(c1.mul(c1))
	r := a.mul(a).scale(-1.0 / 3).add(b)
	s := a.mul(a).mul(a).scale(2.0 / 27).sub(a.mul(b).scale(1.0 / 3)).add(c)
	r3 := r.mul(r).mul(r)
	s2 := s.mul(s)
	pp := s2.scale(0.25).add(r3.scale(1.0 / 27))
	p0 := pp.abs()
	p2 := p0.sqrt()
	var y3 dual
	if pp.v < 0 {
		ro := s2.scale(0.25).add(p0).sqrt().log().scale(1.0 / 3).exp()
		fi := p2.scale(-2).div(s).atan()
		y3 = ro.mul(fi.scale(1.0 / 3).cos()).scale(2).sub(a.scale(1.0 / 3))
	} else {
		y3 = cbrtAbs(s.scale(-0.5).add(p2)).add(cbrtAbs(s.scale(-0.5).sub(p2))).sub(a.scale(1.0 / 3))
	}
	ra := a1.mul(a1).scale(0.25).sub(b1).add(y3)
	rb := y3.mul(y3).scale(0.25).sub(d1)
	if !(ra.v >= 0) || !(rb.v >= 0) {
		return dual{}, false
	}
	a3, b3 := ra.sqrt(), rb.sqrt()
	sig1 := [4]float64{1, -1, 1, -1}
	sig2 := [4]float64{1, 1, -1, -1}
	for i := 0; i < 4; i++ {
		a4 := a1.scale(0.5).add(a3.scale(sig1[i]))
		b4 := y3.scale(0.5).add(b3.scale(sig2[i]))
		delta4 := a4.mul(a4).scale(0.25).sub(b4)
		if !(delta4.v >= 0) {
			continue
		}
		t := delta4.sqrt()
		for _, r := range [2]dual{a4.scale(-0.5).add(t), a4.scale(-0.5).sub(t)} {
			if !(r.v > 0) {
				continue
			}
			x2 := r.v * r.v
			poly := x2*x2 + a1.v*x2*r.v + b1.v*x2 + c1.v*r.v + d1.v
			if !(math.Abs(poly) <= 1e-6) {
				continue
			}
			if !ok || r.v < x.v {
				x, ok = r, true
			}
		}
	}
	return x, ok
}

// cbrtAbs is |a|^(1/3) as mos2load computes it (exp(log|a|/3)); 0 at 0.
func cbrtAbs(a dual) dual {
	a = a.abs()
	if a.v == 0 {
		return dual{}
	}
	return a.powc(1.0 / 3)
}

// idsLevel3 is the SPICE2 level 3 semi-empirical drain current
// (Vladimirescu & Liu 1980, §2.4): threshold with the static-feedback
// (η), short-channel (XJ, LD) and narrow-channel (δ) corrections, mobility
// modulation by the gate field (θ), velocity saturation (VMAX),
// channel-length modulation with κ and punch-through, and weak inversion
// with NFS; region structure and constants as in the ngspice mos3load (BSD
// licence), including ngspice's default (non-"badmos3") channel-length
// modulation, which is continuous at vdsat. Evaluated in dual numbers, so
// gm, gds and gmbs are the exact derivatives of the current.
func (q *MOSFET) idsLevel3(vgsV, vdsV, vbsV float64) (cdrain, gm, gds, gmbs, von, vdsat float64) {
	const coeff0, coeff1, coeff2 = 0.0631353, 0.8013292, -0.01110777
	p, lv := &q.p, &q.lv
	G, D, B := variable(vgsV, dVGS), variable(vdsV, dVDS), variable(vbsV, dVBS)
	phi := p.PHI
	leff, weff := q.leff, q.weff
	var phibs, sqphbs dual
	if B.v <= 0 {
		phibs = B.neg().addc(phi)
		sqphbs = phibs.sqrt()
	} else {
		sqphbs = cst(math.Sqrt(phi)).div(B.scale(1 / (2 * phi)).addc(1))
		phibs = sqphbs.mul(sqphbs)
	}
	// Short-channel factor from the junction depth and depletion width.
	fshort := cst(1)
	if p.XJ != 0 && lv.xd != 0 {
		wponxj := sqphbs.scale(lv.xd / p.XJ)
		wconxj := wponxj.mul(wponxj).scale(coeff2).add(wponxj.scale(coeff1)).addc(coeff0)
		djonxj := p.LD / p.XJ
		arga := wconxj.addc(djonxj)
		argc := wponxj.div(wponxj.addc(1))
		argb := argc.mul(argc).neg().addc(1).sqrt()
		fshort = arga.mul(argb).addc(-djonxj).scale(-p.XJ / leff).addc(1)
	}
	gammas := fshort.scale(p.GAMMA)
	fbodys := gammas.div(sqphbs).scale(0.25)
	fbody := fbodys.addc(lv.narrow / weff)
	qbonco := gammas.mul(sqphbs).add(phibs.scale(lv.narrow / weff))
	vth := D.scale(-lv.eta3).addc(q.vbi).add(qbonco)
	vonD := vth
	nfs := p.NFS != 0
	var xn dual
	if nfs {
		xn = qbonco.div(phibs.scale(2)).addc(1 + lv.cfs)
		vonD = vth.add(xn.scale(q.vt))
	} else if G.v <= vonD.v {
		return 0, 0, 0, 0, vonD.v, 0
	}
	vgsx := dmax(G, vonD)
	onfg := vgsx.sub(vth).scale(p.THETA).addc(1)
	fgate := cst(1).div(onfg)
	us := fgate.scale(p.UO * 1e-4)
	onfbdy := cst(1).div(fbody.addc(1))
	arga := vgsx.sub(vth).mul(onfbdy)
	vdsatD := arga
	var vdsc dual
	if p.VMAX > 0 {
		vdsc = cst(leff * p.VMAX).div(us)
		argb := arga.mul(arga).add(vdsc.mul(vdsc)).sqrt()
		vdsatD = arga.add(vdsc).sub(argb)
	}
	vdsx := dmin(D, vdsatD)
	beta := fgate.scale(q.beta)
	if vdsx.v == 0 {
		// vds = 0: zero current, conductance β·(vgsx − vth).
		g := beta.mul(vgsx.sub(vth))
		if nfs && G.v < vonD.v {
			g = g.mul(G.sub(vonD).div(xn.scale(q.vt)).exp())
		}
		return 0, 0, g.v, 0, vonD.v, vdsatD.v
	}
	cdo := vgsx.sub(vth).sub(fbody.addc(1).mul(vdsx).scale(0.5))
	id := beta.mul(cdo).mul(vdsx)
	var fdrain dual
	if p.VMAX > 0 {
		fdrain = cst(1).div(vdsx.div(vdsc).addc(1))
		id = id.mul(fdrain)
	}
	// Channel-length modulation (with punch-through).
	clm := false
	var delxl dual
	switch {
	case D.v <= vdsatD.v:
		if p.VMAX <= 0 && lv.alpha != 0 {
			r := D.div(vdsatD)
			r2 := r.mul(r)
			delxl = vdsatD.scale(p.KAPPA * lv.alpha / 8).sqrt().mul(r2.mul(r2))
			clm = true
		}
	case p.VMAX <= 0:
		delxl = D.sub(vdsatD).add(vdsatD.scale(1.0 / 8)).scale(p.KAPPA * lv.alpha).sqrt()
		clm = true
	case lv.alpha != 0:
		gdsat := dmax(id.mul(fdrain.neg().addc(1)).div(vdsc), cst(1e-12))
		emax := id.div(gdsat).scale(p.KAPPA / leff)
		a := emax.scale(0.5 * lv.alpha)
		argb := a.mul(a).add(D.sub(vdsatD).scale(p.KAPPA * lv.alpha)).sqrt()
		delxl = argb.sub(a)
		clm = true
	}
	if clm {
		if delxl.v > 0.5*leff {
			delxl = cst(-leff * leff / 4).div(delxl).addc(leff)
		}
		id = id.div(delxl.scale(-1 / leff).addc(1))
	}
	if nfs && G.v < vonD.v {
		id = id.mul(G.sub(vonD).div(xn.scale(q.vt)).exp())
	}
	return fromDual(id, vonD, vdsatD)
}
