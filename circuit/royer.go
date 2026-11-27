// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strings"
)

// RoyerSpec describes a two-transistor Royer converter stage (G. H. Royer,
// "A switching transistor DC to AC converter having an output frequency
// proportional to the DC input voltage," AIEE Trans. 74, 322–326, 1955):
// a centre-tapped primary of 2·NP turns, a centre-tapped base (feedback)
// winding of 2·NB turns and an optional secondary of NS turns, all on one
// saturable Jiles–Atherton core (Kind Core). The transistors switch over
// each time the core saturates, so the stage runs at
//
//	f ≈ V_in / (4 · NP · B_sat · A_e)
//
// (RoyerFrequency). The stage is a composite of ordinary elements, so the
// ledger books every loss — core hysteresis, transistors, base resistors —
// in its own element.
type RoyerSpec struct {
	// Core holds the Jiles–Atherton parameters by netlist name (ms, a, k,
	// c, alpha, area, len); turns are set from NP, NB, NS.
	Core       map[string]float64
	NP, NB, NS float64 // turns: each primary half, each base half, secondary
	RB         float64 // base resistor of each transistor [Ω]
	RBias      float64 // start resistor from the supply to the base centre tap [Ω]
	// BJTModel and DiodeModel name NPN and diode models defined with
	// Builder.Model; empty selects RoyerNPN and RoyerDiode, defined here.
	BJTModel, DiodeModel string
	// Mismatch is the relative difference between the two base resistors
	// (Q1 gets RB·(1−Mismatch)). A perfectly symmetric stage sits on its
	// unstable balanced point; a real one starts from component mismatch.
	// 0 selects 0.01.
	Mismatch float64
}

// Default models for a Royer stage: a switching NPN with a collector–base
// capacitance (without one the switching edge is instantaneous and the
// transient cannot follow it) and a plain silicon diode for the base
// return.
const (
	RoyerNPN   = "royer_npn"
	RoyerDiode = "royer_d"
)

// RoyerFrequency is Royer's frequency for a primary half of np turns, core
// saturation flux density bsat [T] and core area area [m²].
func RoyerFrequency(vin, np, bsat, area float64) float64 {
	return vin / (4 * np * bsat * area)
}

// Royer adds a Royer stage fed from vin (emitters on gnd). The secondary,
// when NS > 0, is the winding outP → outN. Internal nodes are name_c1,
// name_c2 (collectors), name_b1, name_b2 (base-winding ends), name_x1,
// name_x2 (transistor bases) and name_bt (base centre tap); elements are
// U<name> (core), Q<name>1, Q<name>2, R<name>1, R<name>2 (base), R<name>s
// (start) and D<name> (base return from gnd to the tap).
func (b *Builder) Royer(name, vin, gnd, outP, outN string, s RoyerSpec) *Builder {
	n := func(x string) string { return strings.ToLower(name) + "_" + x }
	bjt, dio := s.BJTModel, s.DiodeModel
	if bjt == "" {
		bjt = RoyerNPN
		b.Model(RoyerNPN, "npn", map[string]float64{"bf": 100, "is": 1e-14, "cjc": 20e-12,
			"rb": 5, "rc": 0.5})
	}
	if dio == "" {
		dio = RoyerDiode
		b.Model(RoyerDiode, "d", map[string]float64{"is": 1e-14})
	}
	mm := s.Mismatch
	if mm == 0 {
		mm = 0.01
	}
	p := make(map[string]float64, len(s.Core)+5)
	for k, v := range s.Core {
		p[strings.ToLower(k)] = v
	}
	p["n1"], p["n2"], p["n3"], p["n4"] = s.NP, s.NP, s.NB, s.NB
	// Dot convention: Q1's half vin→c1 and Q2's half c2→vin are one
	// centre-tapped winding; the base winding b1→bt→b2 is wound the same
	// way, so Q1 pulling c1 down drives b1 up and b2 down.
	w := []string{vin, n("c1"), n("c2"), vin, n("b1"), n("bt"), n("bt"), n("b2")}
	if s.NS > 0 {
		p["n5"] = s.NS
		w = append(w, outP, outN)
	}
	b.Core("U"+name, w, p)
	b.Resistor("R"+name+"1", n("b1"), n("x1"), s.RB*(1-mm))
	b.Resistor("R"+name+"2", n("b2"), n("x2"), s.RB)
	b.BJT("Q"+name+"1", n("c1"), n("x1"), gnd, bjt)
	b.BJT("Q"+name+"2", n("c2"), n("x2"), gnd, bjt)
	b.Resistor("R"+name+"s", vin, n("bt"), s.RBias)
	return b.Diode("D"+name, gnd, n("bt"), dio)
}

// ValidateRoyer checks a RoyerSpec before it is built.
func ValidateRoyer(s RoyerSpec) error {
	for _, v := range []float64{s.NP, s.NB, s.RB, s.RBias} {
		if !(v > 0) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: Royer stage needs NP, NB, RB, RBias > 0", ErrBadValue)
		}
	}
	if !(s.NS >= 0) || !(s.Mismatch >= 0 && s.Mismatch < 0.5) {
		return fmt.Errorf("%w: Royer stage needs NS ≥ 0 and 0 ≤ Mismatch < 0.5", ErrBadValue)
	}
	return nil
}
