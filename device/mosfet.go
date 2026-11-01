// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// MOSFET is the SPICE3 MOS transistor, levels 1–3 (drain-current equations
// in mosfet_levels.go); everything else — terminals, junction diodes and
// charges, Meyer gate capacitances, series resistances, limiting,
// source/drain interchange — is shared by the levels.
//
// Maxwell derivation (CLAUDE.md §3c): the channel current is a field-
// controlled conduction current 𝔎 (gradual-channel approximation of the
// quasi-static field equations); the gate and junction charges add
// displacement-current terms ∂q/∂t, as for a capacitor.
//
// Sources (equations, not code): H. Shichman, D. A. Hodges, "Modeling and
// simulation of insulated-gate field-effect transistor switching circuits,"
// IEEE JSSC 3(3), 1968 (level 1); J. E. Meyer, "MOS models and circuit
// simulation," RCA Review 32, 1971 (gate capacitances); A. Vladimirescu,
// S. Liu, "The simulation of MOS integrated circuits using SPICE2," UCB/ERL
// M80/7, 1980 (levels 1–3, process parameters, junctions); the SPICE3f5 /
// ngspice mos1–3 equations (BSD licence) for region formulas, defaults and
// limiting (DEVfetlim, DEVlimvds, DEVpnjlim, DEVqmeyer).
//
// Terminals: Nodes = {drain, source} (the reported port: V_DS and the
// drain current), Extra = {gate, bulk}. Internal drain and source nodes sit
// behind RD and RS when they are non-zero.
//
// The Meyer gate capacitances are not charge-conserving: as in SPICE, the
// gate charge advances by the average capacitance times the voltage step.
// The ledger books the power into them exactly (v·i); the stored energy
// reported for them, ½·C·v² at the present voltage, is an estimate, so a
// transient's integration mismatch includes Meyer's inconsistency.
type MOSFET struct {
	Base
	level int
	typ   float64 // +1 NMOS, −1 PMOS
	p     mosParams

	// Derived (nominal temperature).
	vt               float64 // kT/q
	leff, weff       float64
	beta             float64 // KP·W/Leff·M
	cox              float64 // oxide capacitance per area [F/m²]
	oxideCap         float64 // Cox·Leff·W·M [F]
	vbi              float64 // VTO − type·γ·√φ (normalized)
	vto              float64
	gdpr, gspr       float64 // 1/RD, 1/RS
	isd, iss         float64 // drain and source junction saturation currents
	vcritD, vcritS   float64
	czbd, czbdsw     float64
	czbs, czbssw     float64
	fbt, fsw         [3]float64 // depletion continuation constants (bottom, sidewall)
	cgso, cgdo, cgbo float64    // overlap capacitances [F]
	hasJcap, hasGcap bool
	lv               levelConsts

	// Unknowns: d, s (Base.u), g, b; internal dp, sp.
	g, b, dp, sp int

	// Linearization of the most recent Load.
	vgs, vds, vbs, vbd, vgd, vgb float64
	mode                         float64 // +1 normal, −1 drain and source interchanged
	cdrain, gm, gds, gmbs        float64 // channel current (normal-mode frame) and derivatives
	von, vdsat                   float64
	cbs, gbs, cbd, gbd           float64
	limited, tran                bool
	ag                           float64
	gmin                         float64

	// Gate (Meyer, incl. overlap) and bulk charges: gs, gd, gb, bs, bd.
	capG [3]float64 // average gate capacitances of this step
	capB [2]float64 // bulk junction capacitances
	iq   [5]float64 // charge currents of this step
	dv   [5]float64 // voltage of each charge branch at the linearization point

	// Accepted state.
	acc           struct{ vgs, vds, vbs, vgd, vgb, vbd, von float64 }
	meyerAcc      [3]float64 // full Meyer capacitances at the accepted point
	qPrev, iqPrev [5]float64
	vPrev         [5]float64
	hist          [5]history
	ac            mosAC
	vSol          [5]float64
}

// mosParams are the model and instance parameters in SI units.
type mosParams struct {
	VTO, KP, GAMMA, PHI, LAMBDA, RD, RS, CBD, CBS, IS, PB float64
	CGSO, CGDO, CGBO, RSH, CJ, MJ, CJSW, MJSW, JS, TOX    float64
	NSUB, NSS, NFS, TPG, XJ, LD, UO, UCRIT, UEXP, UTRA    float64
	VMAX, NEFF, FC, DELTA, THETA, ETA, KAPPA, WD          float64
	L, W, AD, AS, PD, PS, NRD, NRS, M                     float64
	given                                                 map[string]bool
}

// Charge branch indices.
const (
	cGS = iota
	cGD
	cGB
	cBS
	cBD
)

// Physical constants for the process-parameter equations (as in SPICE3).
const (
	epsOx = 3.9 * 8.854214871e-12  // SiO₂ permittivity [F/m]
	epsSi = 11.7 * 8.854214871e-12 // silicon permittivity [F/m]
	niSi  = 1.45e16                // intrinsic carrier density at 300 K [m⁻³]
)

// DefaultMOSLength and DefaultMOSWidth are SPICE's DEFL and DEFW.
const (
	DefaultMOSLength = 100e-6
	DefaultMOSWidth  = 100e-6
)

// NewMOSFET builds a transistor from its element and model.
func NewMOSFET(e circuit.Element, m circuit.Model) (*MOSFET, error) {
	q := &MOSFET{Base: newBase(e), level: 1, typ: 1}
	if m.Type == "pmos" {
		q.typ = -1
	}
	if lv, ok := m.Params["level"]; ok {
		q.level = int(lv)
	}
	p := &q.p
	p.given = map[string]bool{}
	// Defaults: SPICE3f3 User's Manual §4.4.
	*p = mosParams{KP: 2e-5, PHI: 0.6, IS: 1e-14, PB: 0.8, MJ: 0.5, MJSW: 0.5, UO: 600, FC: 0.5,
		UCRIT: 1e4, NEFF: 1, DELTA: 0, given: p.given, TPG: 1,
		L: DefaultMOSLength, W: DefaultMOSWidth, M: 1}
	if q.level >= 2 {
		p.MJSW = 0.33 // mos2set, mos3set
	}
	if q.level >= 2 {
		p.TOX = 1e-7
		p.given["tox"] = true
	}
	set := func(dst *float64, keys ...string) {
		for _, k := range keys {
			if v, ok := m.Params[k]; ok {
				*dst = v
				p.given[keys[0]] = true
				return
			}
		}
	}
	set(&p.VTO, "vto", "vt0")
	set(&p.KP, "kp")
	set(&p.GAMMA, "gamma")
	set(&p.PHI, "phi")
	set(&p.LAMBDA, "lambda")
	set(&p.RD, "rd")
	set(&p.RS, "rs")
	set(&p.CBD, "cbd")
	set(&p.CBS, "cbs")
	set(&p.IS, "is")
	set(&p.PB, "pb")
	set(&p.CGSO, "cgso")
	set(&p.CGDO, "cgdo")
	set(&p.CGBO, "cgbo")
	set(&p.RSH, "rsh")
	set(&p.CJ, "cj")
	set(&p.MJ, "mj")
	set(&p.CJSW, "cjsw")
	set(&p.MJSW, "mjsw")
	set(&p.JS, "js")
	set(&p.TOX, "tox")
	set(&p.NSUB, "nsub")
	set(&p.NSS, "nss")
	set(&p.NFS, "nfs")
	set(&p.TPG, "tpg")
	set(&p.XJ, "xj")
	set(&p.LD, "ld")
	set(&p.UO, "uo", "u0")
	set(&p.UCRIT, "ucrit")
	set(&p.UEXP, "uexp")
	set(&p.UTRA, "utra")
	set(&p.VMAX, "vmax")
	set(&p.NEFF, "neff")
	set(&p.FC, "fc")
	set(&p.DELTA, "delta")
	set(&p.THETA, "theta")
	set(&p.ETA, "eta")
	set(&p.KAPPA, "kappa")
	set(&p.WD, "wd")
	if q.level == 3 && !p.given["kappa"] {
		p.KAPPA = 0.2
	}
	for k, dst := range map[string]*float64{"l": &p.L, "w": &p.W, "ad": &p.AD, "as": &p.AS, "pd": &p.PD,
		"ps": &p.PS, "nrd": &p.NRD, "nrs": &p.NRS, "m": &p.M} {
		if v, ok := e.Params[k]; ok {
			*dst = v
		}
	}
	if err := q.derive(); err != nil {
		return nil, fmt.Errorf("%w: %s: %v", ErrUnsupported, e.Name, err)
	}
	q.g, q.b = e.Extra[0], e.Extra[1]
	_ = q.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network; stores gate and junction charge")
	_ = q.SetIntrinsic(exchange.Thermal, exchange.Out, "channel, junction and series-resistance heating to surroundings")
	return q, nil
}

// derive computes the nominal-temperature quantities (SPICE3 mos1temp and
// Vladimirescu & Liu 1980 §2: process parameters → KP, φ, γ, VTO).
func (q *MOSFET) derive() error {
	p := &q.p
	q.vt = BoltzmanK * TNomK / ChargeQ
	q.leff = p.L - 2*p.LD
	q.weff = p.W
	if q.level == 3 {
		q.weff = p.W - 2*p.WD
	}
	if !(q.leff > 0) || !(q.weff > 0) {
		return fmt.Errorf("effective channel length %g or width %g is not positive", q.leff, q.weff)
	}
	if p.given["tox"] && p.TOX > 0 {
		q.cox = epsOx / p.TOX
		if !p.given["kp"] {
			p.KP = p.UO * q.cox * 1e-4
		}
		if p.given["nsub"] {
			nsub := p.NSUB * 1e6
			ni := niSi
			if q.level == 3 {
				// ngspice's level 3 (mos3temp.c), not levels 1–2, scales
				// ni from 300 K to TNOM: (T/300)^1.5·e^(Eg(T)/2k·(1/300 − 1/T)).
				// At 27 °C this is 1.0116, so φ is 0.6 mV lower. Checked
				// against ngspice-39 (testdata/netlists/mos3_subvt_dc.cir:
				// subthreshold current at VGS = 0 within 1e-6, against 1 %
				// with the unscaled 1.45e16).
				eg := 1.16 - 7.02e-4*TNomK*TNomK/(TNomK+1108)
				r := TNomK / 300
				ni *= r * math.Sqrt(r) * math.Exp(0.5*eg*(1.0/300-1/TNomK)/(BoltzmanK/ChargeQ))
			}
			if !(nsub > ni) {
				return fmt.Errorf("NSUB %g cm⁻³ is below the intrinsic density (SPICE: Nsub < Ni)", p.NSUB)
			}
			if !p.given["phi"] {
				p.PHI = math.Max(0.1, 2*q.vt*math.Log(nsub/ni))
			}
			if !p.given["gamma"] {
				p.GAMMA = math.Sqrt(2*epsSi*ChargeQ*nsub) / q.cox
			}
			if !p.given["vto"] {
				eg := 1.16 - 7.02e-4*TNomK*TNomK/(TNomK+1108)
				fermis := q.typ * 0.5 * p.PHI
				wkfng := 3.2
				if p.TPG != 0 {
					fermig := q.typ * p.TPG * 0.5 * eg
					wkfng = 3.25 + 0.5*eg - fermig
				}
				wkfngs := wkfng - (3.25 + 0.5*eg + fermis)
				vfb := wkfngs - p.NSS*1e4*ChargeQ/q.cox
				p.VTO = vfb + q.typ*(p.GAMMA*math.Sqrt(p.PHI)+p.PHI)
			}
			if !p.given["cj"] {
				p.CJ = math.Sqrt(epsSi * ChargeQ * nsub / (2 * p.PB))
			}
		}
	}
	q.oxideCap = q.cox * q.leff * q.weff * p.M
	q.beta = p.KP * q.weff / q.leff * p.M
	q.vto = q.typ * p.VTO // normalized threshold
	q.vbi = q.vto - p.GAMMA*math.Sqrt(p.PHI)
	rd, rs := p.RD, p.RS
	if !p.given["rd"] && p.RSH > 0 {
		rd = p.RSH * p.NRD
	}
	if !p.given["rs"] && p.RSH > 0 {
		rs = p.RSH * p.NRS
	}
	if rd > 0 {
		q.gdpr = p.M / rd
	}
	if rs > 0 {
		q.gspr = p.M / rs
	}
	// SPICE: JS·area for both junctions only when JS, AD and AS are all
	// non-zero, otherwise IS for both (mos1load, mos2load, mos3load).
	q.isd, q.iss = p.IS*p.M, p.IS*p.M
	if p.JS > 0 && p.AD > 0 && p.AS > 0 {
		q.isd, q.iss = p.JS*p.AD*p.M, p.JS*p.AS*p.M
	}
	q.vcritD = q.vt * math.Log(q.vt/(math.Sqrt2*math.Max(q.isd, 1e-300)))
	q.vcritS = q.vt * math.Log(q.vt/(math.Sqrt2*math.Max(q.iss, 1e-300)))
	if p.given["cbd"] {
		q.czbd = p.CBD * p.M
	} else {
		q.czbd = p.CJ * p.AD * p.M
	}
	if p.given["cbs"] {
		q.czbs = p.CBS * p.M
	} else {
		q.czbs = p.CJ * p.AS * p.M
	}
	q.czbdsw = p.CJSW * p.PD * p.M
	q.czbssw = p.CJSW * p.PS * p.M
	f1, f2, f3 := deplF(p.PB, p.MJ, p.FC)
	q.fbt = [3]float64{f1, f2, f3}
	f1, f2, f3 = deplF(p.PB, p.MJSW, p.FC)
	q.fsw = [3]float64{f1, f2, f3}
	q.cgso, q.cgdo, q.cgbo = p.CGSO*q.weff*p.M, p.CGDO*q.weff*p.M, p.CGBO*q.leff*p.M
	q.hasJcap = q.czbd+q.czbdsw+q.czbs+q.czbssw > 0
	q.hasGcap = q.oxideCap > 0 || q.cgso+q.cgdo+q.cgbo > 0
	if err := q.deriveLevel(); err != nil {
		return err
	}
	q.acc.von = q.vto
	return nil
}

// Setup maps terminals and allocates internal drain and source nodes.
func (q *MOSFET) Setup(ctx *SetupContext) error {
	q.setupNodes(ctx)
	q.g, q.b = ctx.Unknown(q.g), ctx.Unknown(q.b)
	q.dp, q.sp = q.u[0], q.u[1]
	if q.gdpr > 0 {
		q.dp = ctx.AllocUnknown()
	}
	if q.gspr > 0 {
		q.sp = ctx.AllocUnknown()
	}
	return nil
}

// voltages returns the type-normalized vgs, vds, vbs and vgb at x.
func (q *MOSFET) voltages(x []float64) (vgs, vds, vbs, vgb float64) {
	t := q.typ
	vg, vb, vdp, vsp := at(x, q.g), at(x, q.b), at(x, q.dp), at(x, q.sp)
	return t * (vg - vsp), t * (vdp - vsp), t * (vb - vsp), t * (vg - vb)
}

// fetlim limits a gate-voltage step (SPICE3 DEVfetlim).
func fetlim(vnew, vold, vto float64) float64 {
	vtsthi := math.Abs(2*(vold-vto)) + 2
	vtstlo := vtsthi/2 + 2
	vtox := vto + 3.5
	delv := vnew - vold
	if vold >= vto {
		if vold >= vtox {
			if delv <= 0 {
				if vnew >= vtox {
					if -delv > vtstlo {
						vnew = vold - vtstlo
					}
				} else {
					vnew = math.Max(vnew, vto+2)
				}
			} else if delv >= vtsthi {
				vnew = vold + vtsthi
			}
		} else if delv <= 0 {
			vnew = math.Max(vnew, vto-0.5)
		} else {
			vnew = math.Min(vnew, vto+4)
		}
	} else if delv <= 0 {
		if -delv > vtsthi {
			vnew = vold - vtsthi
		}
	} else {
		vtemp := vto + 0.5
		if vnew <= vtemp {
			if delv > vtstlo {
				vnew = vold + vtstlo
			}
		} else {
			vnew = vtemp
		}
	}
	return vnew
}

// limvds limits a drain–source step (SPICE3 DEVlimvds).
func limvds(vnew, vold float64) float64 {
	if vold >= 3.5 {
		if vnew > vold {
			return math.Min(vnew, 3*vold+2)
		}
		if vnew < 3.5 {
			return math.Max(vnew, 2)
		}
		return vnew
	}
	if vnew > vold {
		return math.Min(vnew, 4)
	}
	return math.Max(vnew, -0.5)
}

// bulkDiode returns a bulk–source or bulk–drain junction current and
// conductance, with gmin (SPICE3 mos1load): IS(e^(v/vt) − 1) + gmin·v, and
// below −3·vt the constant −IS + gmin·v (the exponential term is under 5 %
// of IS there).
func (q *MOSFET) bulkDiode(v, is float64) (i, g float64) {
	if v <= -3*q.vt {
		return q.gmin*v - is, q.gmin
	}
	e := math.Exp(math.Min(v/q.vt, 700))
	return is*(e-1) + q.gmin*v, is*e/q.vt + q.gmin
}

// meyer returns the full Meyer gate capacitances cgs, cgd, cgb for the
// (normal-mode) voltages, without overlap (Meyer 1971; SPICE3 DEVqmeyer
// returns half of these and sums two time points).
func meyer(vgs, vgd, von, vdsat, phi, cox float64) (cgs, cgd, cgb float64) {
	vgst := vgs - von
	vdsat = math.Max(vdsat, 0)
	switch {
	case vgst <= -phi:
		return 0, 0, cox
	case vgst <= -phi/2:
		return 0, 0, -vgst * cox / phi
	case vgst <= 0:
		return 2 * (vgst*cox/(1.5*phi) + cox/3), 0, -vgst * cox / phi
	}
	vds := vgs - vgd
	if vdsat <= vds {
		return 2 * cox / 3, 0, 0
	}
	vddif := 2*vdsat - vds
	vddif1 := vdsat - vds
	vddif2 := vddif * vddif
	return 2 * cox * (1 - vddif1*vddif1/vddif2) / 3, 2 * cox * (1 - vdsat*vdsat/vddif2) / 3, 0
}

// meyerCaps returns cgs, cgd, cgb in the terminal (not normal-mode) frame.
func (q *MOSFET) meyerCaps(vgs, vgd, von, vdsat, mode float64) [3]float64 {
	if q.oxideCap == 0 {
		return [3]float64{}
	}
	if mode > 0 {
		a, b, c := meyer(vgs, vgd, von, vdsat, q.p.PHI, q.oxideCap)
		return [3]float64{a, b, c}
	}
	a, b, c := meyer(vgd, vgs, von, vdsat, q.p.PHI, q.oxideCap)
	return [3]float64{b, a, c}
}

// bulkCharge returns a bulk junction's charge and capacitance (bottom and
// sidewall depletion charges, linear continuation above FC·PB).
func (q *MOSFET) bulkCharge(v, cz, czsw float64) (qv, c float64) {
	p := &q.p
	q1, c1 := depletion(v, cz, p.PB, p.MJ, p.FC, q.fbt[0], q.fbt[1], q.fbt[2])
	q2, c2 := depletion(v, czsw, p.PB, p.MJSW, p.FC, q.fsw[0], q.fsw[1], q.fsw[2])
	return q1 + q2, c1 + c2
}

// channel evaluates the drain current in the normal-mode frame for the
// terminal voltages, returning the mode as well.
func (q *MOSFET) channel(vgs, vds, vbs float64) (mode, cdrain, gm, gds, gmbs, von, vdsat float64) {
	mode = 1
	vg, vd, vb := vgs, vds, vbs
	if vds < 0 {
		mode = -1
		vg, vd, vb = vgs-vds, -vds, vbs-vds
	}
	cdrain, gm, gds, gmbs, von, vdsat = q.ids(vg, vd, vb)
	return
}

// Load linearizes the transistor at the current iterate and stamps it.
func (q *MOSFET) Load(ctx *LoadContext) error {
	q.gmin = ctx.Gmin
	vgs, vds, vbs, _ := q.voltages(ctx.X)
	q.limited = false
	if ctx.InitJct {
		vbs, vgs, vds = -1, q.vto, 0
		q.limited = true
	} else {
		a := &q.acc
		vgd := vgs - vds
		vgdOld := a.vgs - a.vds
		von := a.von
		if a.vds >= 0 {
			vgs = fetlim(vgs, a.vgs, von)
			vds = vgs - vgd
			vds = limvds(vds, a.vds)
			vgd = vgs - vds
		} else {
			vgd = fetlim(vgd, vgdOld, von)
			vds = vgs - vgd
			vds = -limvds(-vds, -a.vds)
			vgs = vgd + vds
		}
		var l bool
		if vds >= 0 {
			vbs, l = pnjlim(vbs, a.vbs, q.vt, q.vcritS)
		} else {
			vbd, l2 := pnjlim(vbs-vds, a.vbs-a.vds, q.vt, q.vcritD)
			vbs, l = vbd+vds, l2
		}
		q.limited = l
	}
	vbd, vgd := vbs-vds, vgs-vds
	vgb := vgs - vbs
	q.vgs, q.vds, q.vbs, q.vbd, q.vgd, q.vgb = vgs, vds, vbs, vbd, vgd, vgb
	// The next iteration limits relative to this one.
	q.acc.vgs, q.acc.vds, q.acc.vbs = vgs, vds, vbs
	q.cbs, q.gbs = q.bulkDiode(vbs, q.iss)
	q.cbd, q.gbd = q.bulkDiode(vbd, q.isd)
	q.mode, q.cdrain, q.gm, q.gds, q.gmbs, q.von, q.vdsat = q.channel(vgs, vds, vbs)
	q.acc.von = q.von

	q.tran, q.ag = false, 0
	q.iq = [5]float64{}
	q.dv = [5]float64{vgs, vgd, vgb, vbs, vbd}
	if ctx.Mode == ModeTran && (q.hasGcap || q.hasJcap) {
		q.tran = true
		if ctx.Method == Trapezoidal {
			q.ag = 2 / ctx.H
		} else {
			q.ag = 1 / ctx.H
		}
		mc := q.meyerCaps(vgs, vgd, q.von, q.vdsat, q.mode)
		ov := [3]float64{q.cgso, q.cgdo, q.cgbo}
		for k := 0; k < 3; k++ {
			q.capG[k] = (mc[k]+q.meyerAcc[k])/2 + ov[k]
		}
		qbs, cbs := q.bulkCharge(vbs, q.czbs, q.czbssw)
		qbd, cbd := q.bulkCharge(vbd, q.czbd, q.czbdsw)
		q.capB = [2]float64{cbs, cbd}
		qs := [5]float64{
			q.qPrev[cGS] + q.capG[0]*(vgs-q.vPrev[cGS]),
			q.qPrev[cGD] + q.capG[1]*(vgd-q.vPrev[cGD]),
			q.qPrev[cGB] + q.capG[2]*(vgb-q.vPrev[cGB]),
			qbs, qbd,
		}
		for k := range q.iq {
			q.iq[k] = q.ag * (qs[k] - q.qPrev[k])
			if ctx.Method == Trapezoidal {
				q.iq[k] -= q.iqPrev[k]
			}
		}
	}
	q.stamp(ctx)
	return nil
}

// cap returns the capacitance of charge branch k in this step.
func (q *MOSFET) cap(k int) float64 {
	if k < 3 {
		return q.capG[k]
	}
	return q.capB[k-3]
}

// branch stamps a current i0 + Σ g_k·Δv_k from a to b (type-normalized
// values, controls v_k = T·(V(p_k) − V(n_k))).
func (q *MOSFET) branch(ctx *LoadContext, a, b int, i0 float64, n int, g [3]float64, p, m [3]int, v0 [3]float64) {
	ieq := i0
	for k := 0; k < n; k++ {
		if g[k] == 0 {
			continue
		}
		ctx.M.Add(a, p[k], g[k])
		ctx.M.Add(a, m[k], -g[k])
		ctx.M.Add(b, p[k], -g[k])
		ctx.M.Add(b, m[k], g[k])
		ieq -= g[k] * v0[k]
	}
	ieq *= q.typ
	ctx.AddRHS(a, -ieq)
	ctx.AddRHS(b, ieq)
}

// channelBranch returns the channel current dp → sp and its derivatives on
// vgs, vds, vbs in the terminal frame.
func (q *MOSFET) channelBranch() (i0 float64, g [3]float64) {
	if q.mode > 0 {
		return q.cdrain, [3]float64{q.gm, q.gds, q.gmbs}
	}
	// Interchanged: the current flows source → drain, = cdrain(vgd, −vds, vbd).
	// d/dvgs = −gm, d/dvds = gm + gds + gmbs, d/dvbs = −gmbs.
	return -q.cdrain, [3]float64{-q.gm, q.gm + q.gds + q.gmbs, -q.gmbs}
}

func (q *MOSFET) stamp(ctx *LoadContext) {
	g, b, dp, sp := q.g, q.b, q.dp, q.sp
	i0, gc := q.channelBranch()
	q.branch(ctx, dp, sp, i0, 3, gc, [3]int{g, dp, b}, [3]int{sp, sp, sp}, [3]float64{q.vgs, q.vds, q.vbs})
	q.branch(ctx, b, sp, q.cbs, 1, [3]float64{q.gbs}, [3]int{b}, [3]int{sp}, [3]float64{q.vbs})
	q.branch(ctx, b, dp, q.cbd, 1, [3]float64{q.gbd}, [3]int{b}, [3]int{dp}, [3]float64{q.vbd})
	if q.tran {
		ends := [5][2]int{{g, sp}, {g, dp}, {g, b}, {b, sp}, {b, dp}}
		for k := 0; k < 5; k++ {
			c := q.ag * q.cap(k)
			e := ends[k]
			q.branch(ctx, e[0], e[1], q.iq[k], 1, [3]float64{c}, [3]int{e[0]}, [3]int{e[1]}, [3]float64{q.dv[k]})
		}
	}
	if q.dp != q.u[0] {
		stampG(ctx.M, q.u[0], q.dp, q.gdpr)
	}
	if q.sp != q.u[1] {
		stampG(ctx.M, q.u[1], q.sp, q.gspr)
	}
}

// Limited reports whether the last Load limited a voltage.
func (q *MOSFET) Limited() bool { return q.limited }

// ConvTest compares the drain and bulk currents predicted by the
// linearization with the model at x (SPICE3 MOS1convTest).
func (q *MOSFET) ConvTest(x []float64, reltol, abstol float64) bool {
	vgs, vds, vbs, _ := q.voltages(x)
	vbd := vbs - vds
	i0, gc := q.channelBranch()
	ich := i0 + gc[0]*(vgs-q.vgs) + gc[1]*(vds-q.vds) + gc[2]*(vbs-q.vbs)
	ibs := q.cbs + q.gbs*(vbs-q.vbs)
	ibd := q.cbd + q.gbd*(vbd-q.vbd)
	mode, cd, _, _, _, _, _ := q.channel(vgs, vds, vbs)
	ichN := mode * cd
	ibsN, _ := q.bulkDiode(vbs, q.iss)
	ibdN, _ := q.bulkDiode(vbd, q.isd)
	ok := func(a, b float64) bool { return math.Abs(a-b) <= reltol*math.Max(math.Abs(a), math.Abs(b))+abstol }
	return ok(ich-ibd, ichN-ibdN) && ok(ibs+ibd, ibsN+ibdN)
}

// terminalCurrents returns the currents into the device at drain, gate,
// source and bulk from the linearization, and the power into storage.
func (q *MOSFET) terminalCurrents(x []float64) (id, ig, is, ib, ps float64) {
	vgs, vds, vbs, vgb := q.voltages(x)
	vbd, vgd := vbs-vds, vgs-vds
	t := q.typ
	i0, gc := q.channelBranch()
	ich := t * (i0 + gc[0]*(vgs-q.vgs) + gc[1]*(vds-q.vds) + gc[2]*(vbs-q.vbs))
	ibs := t * (q.cbs + q.gbs*(vbs-q.vbs))
	ibd := t * (q.cbd + q.gbd*(vbd-q.vbd))
	var iq [5]float64
	if q.tran {
		v := [5]float64{vgs, vgd, vgb, vbs, vbd}
		for k := 0; k < 5; k++ {
			iq[k] = t * (q.iq[k] + q.ag*q.cap(k)*(v[k]-q.dv[k]))
			ps += v[k] * t * iq[k]
		}
	}
	id = ich - ibd - iq[cGD] - iq[cBD]
	is = -ich - ibs - iq[cGS] - iq[cBS]
	ig = iq[cGS] + iq[cGD] + iq[cGB]
	ib = ibs + ibd - iq[cGB] + iq[cBS] + iq[cBD]
	q.vSol = [5]float64{vgs, vgd, vgb, vbs, vbd}
	return
}

// Operating reports V_DS, the drain current and the total absorbed power.
func (q *MOSFET) Operating(x []float64) OperatingPoint {
	id, ig, is, ib, ps := q.terminalCurrents(x)
	vd, vs := at(x, q.u[0]), at(x, q.u[1])
	p := vd*id + at(x, q.g)*ig + vs*is + at(x, q.b)*ib
	return OperatingPoint{Voltage: vd - vs, Current: id, PowerAbsorbed: p, StoragePower: ps}
}

// Accept commits the converged state.
func (q *MOSFET) Accept(x []float64, t float64, mode Mode) {
	vgs, vds, vbs, vgb := q.voltages(x)
	vgd, vbd := vgs-vds, vbs-vds
	md, _, _, _, _, von, vdsat := q.channel(vgs, vds, vbs)
	q.acc.vgs, q.acc.vds, q.acc.vbs, q.acc.von = vgs, vds, vbs, von
	q.acc.vgd, q.acc.vgb, q.acc.vbd = vgd, vgb, vbd
	mc := q.meyerCaps(vgs, vgd, von, vdsat, md)
	v := [5]float64{vgs, vgd, vgb, vbs, vbd}
	qbs, _ := q.bulkCharge(vbs, q.czbs, q.czbssw)
	qbd, _ := q.bulkCharge(vbd, q.czbd, q.czbdsw)
	if mode == ModeTran && q.tran {
		var qs [5]float64
		for k := 0; k < 3; k++ {
			qs[k] = q.qPrev[k] + q.capG[k]*(v[k]-q.vPrev[k])
		}
		qs[cBS], qs[cBD] = qbs, qbd
		for k := 0; k < 5; k++ {
			q.iqPrev[k] = q.iq[k] + q.ag*q.cap(k)*(v[k]-q.dv[k])
			q.qPrev[k] = qs[k]
			q.hist[k].push(t, qs[k])
		}
	} else {
		// Operating point: gate charges start from C·v at the full
		// capacitance (SPICE MODEINITTRAN), junction charges exactly.
		ov := [3]float64{q.cgso, q.cgdo, q.cgbo}
		for k := 0; k < 3; k++ {
			q.qPrev[k] = (mc[k] + ov[k]) * v[k]
		}
		q.qPrev[cBS], q.qPrev[cBD] = qbs, qbd
		for k := 0; k < 5; k++ {
			q.iqPrev[k] = 0
			q.hist[k].reset()
			q.hist[k].push(t, q.qPrev[k])
		}
	}
	q.vPrev = v
	q.meyerAcc = mc
}

// StoredEnergy: exact co-energy of the junction charges, and ½·C·v² of the
// gate capacitances at the present voltages (an estimate; see the type
// comment on Meyer's model).
func (q *MOSFET) StoredEnergy(OperatingPoint) float64 {
	if !q.hasGcap && !q.hasJcap {
		return 0
	}
	v := q.vSol
	w := 0.0
	for k, cz := range [2][2]float64{{q.czbs, q.czbssw}, {q.czbd, q.czbdsw}} {
		vv := v[3+k]
		if vv == 0 || cz[0]+cz[1] == 0 {
			continue
		}
		qf := func(u float64) float64 { r, _ := q.bulkCharge(u, cz[0], cz[1]); return r }
		w += vv*qf(vv) - integrate(qf, 0, vv, 8)
	}
	md, _, _, _, _, von, vdsat := q.channel(v[0], v[0]-v[1], v[3])
	mc := q.meyerCaps(v[0], v[1], von, vdsat, md)
	ov := [3]float64{q.cgso, q.cgdo, q.cgbo}
	for k := 0; k < 3; k++ {
		w += 0.5 * (mc[k] + ov[k]) * v[k] * v[k]
	}
	return w
}

// Truncate proposes a step from the charges' truncation error.
func (q *MOSFET) Truncate(x []float64, tc *TruncContext) float64 {
	if !q.tran {
		return inf
	}
	vgs, vds, vbs, vgb := q.voltages(x)
	v := [5]float64{vgs, vgs - vds, vgb, vbs, vbs - vds}
	h := inf
	for k := 0; k < 5; k++ {
		var qk float64
		if k < 3 {
			qk = q.qPrev[k] + q.capG[k]*(v[k]-q.vPrev[k])
		} else if k == cBS {
			qk, _ = q.bulkCharge(v[k], q.czbs, q.czbssw)
		} else {
			qk, _ = q.bulkCharge(v[k], q.czbd, q.czbdsw)
		}
		if qk == 0 && q.qPrev[k] == 0 {
			continue
		}
		i := q.iq[k] + q.ag*q.cap(k)*(v[k]-q.dv[k])
		h = math.Min(h, lteStep(&q.hist[k], qk, i, q.iqPrev[k], tc.AbsTol, tc))
	}
	return h
}

// ExteriorPower: dissipated power (absorbed minus stored) leaves as heat.
func (q *MOSFET) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == exchange.Thermal {
		return -(op.PowerAbsorbed - op.StoragePower)
	}
	return 0
}

// mosAC holds the small-signal admittances.
type mosAC struct {
	gch      [3]float64 // channel derivatives on vgs, vds, vbs (terminal frame)
	gbs, gbd float64
	cap      [5]float64
	omega    float64
}

// LoadAC stamps the small-signal model at the accepted operating point:
// channel transconductances, junction conductances, and the Meyer
// capacitances (SPICE: twice the stored half-value, plus overlap) and
// junction capacitances.
func (q *MOSFET) LoadAC(ctx *ACContext) error {
	a := &q.acc
	md, cd, gm, gds, gmbs, von, vdsat := q.channel(a.vgs, a.vds, a.vbs)
	_ = cd
	if md > 0 {
		q.ac.gch = [3]float64{gm, gds, gmbs}
	} else {
		q.ac.gch = [3]float64{-gm, gm + gds + gmbs, -gmbs}
	}
	_, q.ac.gbs = q.bulkDiode(a.vbs, q.iss)
	_, q.ac.gbd = q.bulkDiode(a.vbd, q.isd)
	mc := q.meyerCaps(a.vgs, a.vgd, von, vdsat, md)
	_, cbs := q.bulkCharge(a.vbs, q.czbs, q.czbssw)
	_, cbd := q.bulkCharge(a.vbd, q.czbd, q.czbdsw)
	q.ac.cap = [5]float64{mc[0] + q.cgso, mc[1] + q.cgdo, mc[2] + q.cgbo, cbs, cbd}
	q.ac.omega = ctx.Omega
	g, b, dp, sp := q.g, q.b, q.dp, q.sp
	m := ctx.M
	ctrlP, ctrlN := [3]int{g, dp, b}, [3]int{sp, sp, sp}
	for k := 0; k < 3; k++ {
		acBranch(m, dp, sp, complex(q.ac.gch[k], 0), ctrlP[k], ctrlN[k], 0, 0, 0)
	}
	stampY(m, b, sp, complex(q.ac.gbs, 0))
	stampY(m, b, dp, complex(q.ac.gbd, 0))
	ends := [5][2]int{{g, sp}, {g, dp}, {g, b}, {b, sp}, {b, dp}}
	for k := 0; k < 5; k++ {
		stampY(m, ends[k][0], ends[k][1], complex(0, ctx.Omega*q.ac.cap[k]))
	}
	if q.dp != q.u[0] {
		stampY(m, q.u[0], q.dp, complex(q.gdpr, 0))
	}
	if q.sp != q.u[1] {
		stampY(m, q.u[1], q.sp, complex(q.gspr, 0))
	}
	return nil
}

// acCurrents returns the terminal current phasors into the device.
func (q *MOSFET) acCurrents(x []complex128) (id, ig, is, ib complex128) {
	g, b, dp, sp := catC(x, q.g), catC(x, q.b), catC(x, q.dp), catC(x, q.sp)
	a := &q.ac
	vgs, vds, vbs := g-sp, dp-sp, b-sp
	ich := complex(a.gch[0], 0)*vgs + complex(a.gch[1], 0)*vds + complex(a.gch[2], 0)*vbs
	ibs := complex(a.gbs, 0) * vbs
	ibd := complex(a.gbd, 0) * (b - dp)
	jw := complex(0, a.omega)
	iqs := [5]complex128{jw * complex(a.cap[0], 0) * vgs, jw * complex(a.cap[1], 0) * (g - dp),
		jw * complex(a.cap[2], 0) * (g - b), jw * complex(a.cap[3], 0) * vbs, jw * complex(a.cap[4], 0) * (b - dp)}
	id = ich - ibd - iqs[cGD] - iqs[cBD]
	is = -ich - ibs - iqs[cGS] - iqs[cBS]
	ig = iqs[cGS] + iqs[cGD] + iqs[cGB]
	ib = ibs + ibd - iqs[cGB] + iqs[cBS] + iqs[cBD]
	return
}

// OperatingAC returns the V_DS and drain-current phasors.
func (q *MOSFET) OperatingAC(x []complex128) ACPoint {
	id, _, _, _ := q.acCurrents(x)
	return ACPoint{Voltage: q.voltageAC(x), Current: id}
}

// ACPower is ½·Σ V_k·I_k* over the four terminals.
func (q *MOSFET) ACPower(x []complex128) complex128 {
	id, ig, is, ib := q.acCurrents(x)
	s := catC(x, q.u[0])*cmplx.Conj(id) + catC(x, q.g)*cmplx.Conj(ig) +
		catC(x, q.u[1])*cmplx.Conj(is) + catC(x, q.b)*cmplx.Conj(ib)
	return s / 2
}
