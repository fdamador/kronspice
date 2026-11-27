// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strings"
)

// MagAmpSpec describes a half-wave self-saturating magnetic amplifier
// (P. Mali, Magnetic Amplifiers: Principles and Applications, Rider, 1960,
// pp. 34–35): a gate (load) winding of NG turns in series with a rectifier
// on one saturable Jiles–Atherton core, a control winding of NC turns and
// an optional feedback winding of NF turns carrying the load current.
//
// On each conducting half cycle the gate winding blocks the supply until
// the core saturates and then passes the rest of the half cycle; on the
// other half the rectifier blocks and the control ampere-turns set the
// flux the core is reset to. So the volt-seconds passed, and the mean
// output, follow the control current: the core blocks NG·Ae·ΔB, ΔB being
// the flux the gate has to drive back to saturation. With a filter
// capacitor CF the load current flows on through the feedback winding
// while the rectifier blocks; wound to aid the gate it opposes the reset,
// which is positive feedback: more gain, slower response (Mali
// pp. 54–55). The amplifier controls power; every watt delivered comes
// from the supply, and the ledger books control power, core loss and
// rectifier loss in their own elements.
type MagAmpSpec struct {
	// Core holds the Jiles–Atherton parameters by netlist name (ms, a, k,
	// c, alpha, area, len); turns are set from NG, NC, NF.
	Core       map[string]float64
	NG, NC, NF float64 // gate, control and feedback turns (NF = 0: none)
	CF         float64 // filter capacitance at the rectifier output [F] (0: none)
	// DiodeModel names a diode model defined with Builder.Model; empty
	// selects MagAmpDiode, defined here.
	DiodeModel string
}

// MagAmpDiode is the default rectifier model of a magnetic amplifier.
const MagAmpDiode = "magamp_d"

// MagAmp adds a magnetic amplifier fed from the AC node in (gate winding
// in → name_g, rectifier name_g → name_k, feedback winding name_k → out
// when NF > 0, otherwise name_k is out). The load goes from out to the
// supply's return. The control winding is ctlP → ctlN: control current
// entering ctlP aids saturation (more output), negative current resets
// the core (less output). Elements: U<name> (core), D<name> (rectifier),
// C<name> (filter, when CF > 0).
func (b *Builder) MagAmp(name, in, out, ctlP, ctlN, ret string, s MagAmpSpec) *Builder {
	n := func(x string) string { return strings.ToLower(name) + "_" + x }
	dio := s.DiodeModel
	if dio == "" {
		dio = MagAmpDiode
		b.Model(MagAmpDiode, "d", map[string]float64{"is": 1e-12, "rs": 0.05})
	}
	p := make(map[string]float64, len(s.Core)+3)
	for k, v := range s.Core {
		p[strings.ToLower(k)] = v
	}
	k := out
	if s.NF > 0 {
		k = n("k")
	}
	p["n1"], p["n2"] = s.NG, s.NC
	w := []string{in, n("g"), ctlP, ctlN}
	if s.NF > 0 {
		// Load current entering name_k aids the gate's ampere-turns.
		p["n3"] = s.NF
		w = append(w, k, out)
	}
	b.Core("U"+name, w, p)
	b.Diode("D"+name, n("g"), k, dio)
	if s.CF > 0 {
		b.Capacitor("C"+name, k, ret, s.CF)
	}
	return b
}

// ValidateMagAmp checks a MagAmpSpec before it is built.
func ValidateMagAmp(s MagAmpSpec) error {
	if !(s.NG > 0) || !(s.NC > 0) || math.IsInf(s.NG, 0) || math.IsInf(s.NC, 0) {
		return fmt.Errorf("%w: magnetic amplifier needs NG, NC > 0", ErrBadValue)
	}
	if !(s.NF >= 0) || !(s.CF >= 0) || math.IsInf(s.NF, 0) || math.IsInf(s.CF, 0) {
		return fmt.Errorf("%w: magnetic amplifier needs NF, CF ≥ 0", ErrBadValue)
	}
	return nil
}
