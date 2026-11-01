// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"
)

// levelConsts holds constants the level 2 and 3 currents derive once.
type levelConsts struct {
	xd      float64 // depletion-width coefficient √(2ε_Si / (q·NSUB)) [m/√V]; 0 without NSUB
	alpha   float64 // xd² (level 3 channel-length modulation)
	cfs     float64 // q·NFS / Cox: fast-surface-state term of the weak-inversion slope
	factor  float64 // level 2 narrow-channel factor δ·π·ε_Si / (4·Cox·W)
	narrow  float64 // level 3 narrow-channel factor δ·π·ε_Si / (2·Cox)
	eta3    float64 // level 3 static feedback η·8.15e−22 / (Cox·Leff³)
	sqrtPB  float64 // √PB (level 2 punch-through)
	hasNSUB bool
}

// deriveLevel checks the level and prepares its constants (SPICE3
// mos2temp, mos3temp).
func (q *MOSFET) deriveLevel() error {
	p := &q.p
	lv := &q.lv
	if p.given["nsub"] {
		lv.hasNSUB = true
		lv.alpha = 2 * epsSi / (ChargeQ * p.NSUB * 1e6)
		lv.xd = math.Sqrt(lv.alpha)
	}
	switch q.level {
	case 1:
		return nil
	case 2:
		lv.cfs = ChargeQ * p.NFS * 1e4 / q.cox
		lv.factor = p.DELTA * math.Pi * epsSi / (4 * q.cox * q.weff)
		lv.sqrtPB = math.Sqrt(p.PB)
		if p.VMAX > 0 && !(p.NEFF > 0) {
			return fmt.Errorf("NEFF must be positive")
		}
		return nil
	case 3:
		lv.cfs = ChargeQ * p.NFS * 1e4 / q.cox
		lv.narrow = p.DELTA * 0.5 * math.Pi * epsSi / q.cox
		lv.eta3 = p.ETA * 8.15e-22 / (q.cox * q.leff * q.leff * q.leff)
		return nil
	}
	return fmt.Errorf("MOSFET level %d is not implemented", q.level)
}

// ids is the drain current in the normal-mode frame (vds ≥ 0) and its
// derivatives on vgs, vds and vbs, with the threshold von and the
// saturation voltage vdsat.
func (q *MOSFET) ids(vgs, vds, vbs float64) (cdrain, gm, gds, gmbs, von, vdsat float64) {
	switch q.level {
	case 2:
		return q.idsLevel2(vgs, vds, vbs)
	case 3:
		return q.idsLevel3(vgs, vds, vbs)
	}
	return q.idsLevel1(vgs, vds, vbs)
}

// fromDual unpacks a dual drain current.
func fromDual(id, von, vdsat dual) (cdrain, gm, gds, gmbs, vn, vs float64) {
	return id.v, id.d[0], id.d[1], id.d[2], von.v, vdsat.v
}

// idsLevel1 is the Shichman–Hodges model (Shichman & Hodges 1968; SPICE3
// mos1load), with body effect and channel-length modulation:
//
//	von   = VBI + γ·√(φ − vbs)            (for vbs > 0 the square root is
//	                                       continued linearly, as in SPICE)
//	cutoff     vgs ≤ von:          id = 0
//	saturation vgs − von ≤ vds:    id = β/2·(vgs − von)²·(1 + λ·vds)
//	linear     otherwise:          id = β·vds·(vgs − von − vds/2)·(1 + λ·vds)
func (q *MOSFET) idsLevel1(vgs, vds, vbs float64) (cdrain, gm, gds, gmbs, von, vdsat float64) {
	p := &q.p
	var sarg float64
	if vbs <= 0 {
		sarg = math.Sqrt(p.PHI - vbs)
	} else {
		sarg = math.Sqrt(p.PHI)
		sarg = math.Max(0, sarg-vbs/(sarg+sarg))
	}
	von = q.vbi + p.GAMMA*sarg
	vgst := vgs - von
	vdsat = math.Max(vgst, 0)
	arg := 0.0
	if sarg > 0 {
		arg = p.GAMMA / (sarg + sarg)
	}
	if vgst <= 0 {
		return 0, 0, 0, 0, von, vdsat
	}
	betap := q.beta * (1 + p.LAMBDA*vds)
	if vgst <= vds {
		cdrain = betap * vgst * vgst * 0.5
		gm = betap * vgst
		gds = p.LAMBDA * q.beta * vgst * vgst * 0.5
	} else {
		cdrain = betap * vds * (vgst - 0.5*vds)
		gm = betap * vds
		gds = betap*(vgst-vds) + p.LAMBDA*q.beta*vds*(vgst-0.5*vds)
	}
	// For vbs > 0 this is SPICE's gmbs, γ/(2·sarg) with the continued
	// sarg, not the exact slope of the continued threshold (kept for
	// compatibility with SPICE's small-signal results).
	gmbs = gm * arg
	return cdrain, gm, gds, gmbs, von, vdsat
}
