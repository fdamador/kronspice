// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strings"
)

// Arc is a dynamic arc or spark gap, kronspice A (CLAUDE.md §3g):
//
//	Aname n+ n- MAYR     TAU=τ P0=P0           [common]
//	Aname n+ n- CASSIE   TAU=τ U0=U0           [common]
//	Aname n+ n- HABEDANK TAUC=τc U0=U0 TAUM=τm P0=P0 [common]
//	common: [G0=g0] [VBD=Vbd] [GON=gon] [GOFF=goff] [GMAX=gmax]
//
// The element is a conductance g(t) with i = g·v whose state follows
// (in the current-driven form, equivalent to the published ones because
// v·i = i²/g):
//
//	Mayr (O. Mayr, Arch. Elektrotech. 37, 588, 1943; low current, near
//	current zero, constant cooling power P0):
//	    dg/dt = (g/τ)(v·i/P0 − 1) = (i²/P0 − g)/τ
//	Cassie (A. M. Cassie, CIGRE report 102, 1939; high current, constant
//	arc voltage U0):
//	    dg/dt = (g/τ)(v²/U0² − 1),  d(g²)/dt = 2(i²/U0² − g²)/τ
//	Habedank (U. Habedank, IEEE Trans. Power Delivery 8, 1993): Cassie
//	and Mayr conductances in series, 1/g = 1/g_C + 1/g_M, each driven by
//	the common current.
//
// For a constant current each law is linear in g (Mayr) or g² (Cassie),
// so the steady states are v = P0/|i| (Mayr, a falling characteristic: a
// negative differential resistance that sustains oscillation; Steinmetz,
// Lecture VIII §35) and v = U0 (Cassie). The arc is passive: v·i = g·v² ≥ 0,
// every joule taken from the circuit leaves as heat and light.
//
// G0 is the conductance at the operating point and t = 0 (default GOFF),
// GOFF the floor (the cold gap's insulation leakage, default 1e-9 S), GMAX
// the ceiling (default 1e6 S). VBD > 0 makes the element a spark gap: when
// an accepted transient point has |v| ≥ VBD with g < GON the gap breaks
// down and the state is set to GON (default 1 S); it then burns or
// extinguishes by its law. Breakdown here is a threshold; the statistical
// delay and the mist-dependent breakdown voltages of §3j enter later
// through fluid/mist.go as values of VBD and τ.
const Arc Kind = 'A'

// Arc model names.
const (
	ArcMayr     = "mayr"
	ArcCassie   = "cassie"
	ArcHabedank = "habedank"
)

// ArcParams are an arc element's parameters, defaults filled.
type ArcParams struct {
	Model        string
	Tau, P0, U0  float64 // Mayr / Cassie (TAU, P0, U0); Habedank uses TauC, TauM
	TauC, TauM   float64
	G0, Vbd, Gon float64
	Goff, Gmax   float64
}

// ArcParamsOf returns an arc element's parameters with defaults.
func ArcParamsOf(e Element) ArcParams {
	p := ArcParams{Model: e.Model, Tau: e.Params["tau"], P0: e.Params["p0"], U0: e.Params["u0"],
		TauC: e.Params["tauc"], TauM: e.Params["taum"], Vbd: e.Params["vbd"],
		Goff: 1e-9, Gmax: 1e6, Gon: 1}
	if v, ok := e.Params["goff"]; ok {
		p.Goff = v
	}
	if v, ok := e.Params["gmax"]; ok {
		p.Gmax = v
	}
	if v, ok := e.Params["gon"]; ok {
		p.Gon = v
	}
	p.G0 = p.Goff
	if v, ok := e.Params["g0"]; ok {
		p.G0 = v
	}
	return p
}

// Arc adds an arc between n+ and n−: model is "mayr", "cassie" or
// "habedank", params by lower-case name as in the netlist line.
func (b *Builder) Arc(name, np, nn, model string, params map[string]float64) *Builder {
	p := make(map[string]float64, len(params))
	for k, v := range params {
		p[strings.ToLower(k)] = v
	}
	return b.Add(Element{Name: name, Kind: Arc, Model: strings.ToLower(model), Params: p}, np, nn)
}

func validateArc(e Element) error {
	need := map[string][]string{
		ArcMayr:     {"tau", "p0"},
		ArcCassie:   {"tau", "u0"},
		ArcHabedank: {"tauc", "u0", "taum", "p0"},
	}
	req, ok := need[e.Model]
	if !ok {
		return fmt.Errorf("%w: %s arc model must be MAYR, CASSIE or HABEDANK, got %q", ErrBadValue, e.Name, e.Model)
	}
	allowed := map[string]bool{"g0": true, "vbd": true, "gon": true, "goff": true, "gmax": true}
	for _, k := range req {
		allowed[k] = true
		if v, ok := e.Params[k]; !ok || !(v > 0) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: %s %s arc needs %s > 0", ErrBadValue, e.Name, strings.ToUpper(e.Model), strings.ToUpper(k))
		}
	}
	for k, v := range e.Params {
		if !allowed[k] {
			return fmt.Errorf("%w: %s: parameter %q does not apply to a %s arc", ErrBadValue, e.Name, k, strings.ToUpper(e.Model))
		}
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: %s has a non-finite parameter", ErrBadValue, e.Name)
		}
	}
	p := ArcParamsOf(e)
	if !(p.Goff > 0) || !(p.Gmax > p.Goff) || !(p.Gon >= p.Goff && p.Gon <= p.Gmax) ||
		!(p.G0 >= p.Goff && p.G0 <= p.Gmax) {
		return fmt.Errorf("%w: %s needs 0 < GOFF ≤ G0, GON ≤ GMAX", ErrBadValue, e.Name)
	}
	if p.Vbd < 0 {
		return fmt.Errorf("%w: %s VBD must be ≥ 0", ErrBadValue, e.Name)
	}
	return nil
}
