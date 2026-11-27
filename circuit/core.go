// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strconv"
	"strings"
)

// Core is a saturable magnetic core with hysteresis carrying one or more
// windings, kronspice U (CLAUDE.md §3j):
//
//	Uname w1+ w1- [w2+ w2- …] MS=Ms A=a K=k [C=c] [ALPHA=α] AREA=Ae LEN=le N1=n1 [N2=n2 …]
//
// The magnetization follows the Jiles–Atherton model (D. C. Jiles,
// D. L. Atherton, "Theory of ferromagnetic hysteresis", J. Magn. Magn.
// Mater. 61, 48–60, 1986; with the zero-slope correction of Lederer et al.,
// IEEE Trans. Magn. 35, 1211, 1999):
//
//	He = H + α·M,  Man = Ms·(coth(He/a) − a/He)
//	dMirr/dHe = (Man − Mirr)/(δ·k)   (0 when δ·(Man − Mirr) < 0)
//	M = Mirr + c·(Man − Mirr),  δ = sign(dH/dt)
//
// with H = Σ N_j·i_j / le (Ampère's law round the magnetic path) and
// B = μ0·(H + M) (Maxwell's relation (D), Treatise Art. 619; §3c). Every
// winding sees v_k = N_k·Ae·dB/dt (Faraday). Winding resistance is a
// separate R element. Parameters are in SI: Ms, a, k [A/m], c and α
// dimensionless, Ae [m²], le [m], turns N_k > 0.
const Core Kind = 'U'

// CoreParams are a core element's parameters.
type CoreParams struct {
	Ms, A, K, C, Alpha float64
	Area, Len          float64
	Turns              []float64
}

// CoreParamsOf returns a core element's parameters.
func CoreParamsOf(e Element) CoreParams {
	p := CoreParams{Ms: e.Params["ms"], A: e.Params["a"], K: e.Params["k"], C: e.Params["c"],
		Alpha: e.Params["alpha"], Area: e.Params["area"], Len: e.Params["len"]}
	for w := 1; w <= 1+len(e.Extra)/2; w++ {
		p.Turns = append(p.Turns, e.Params["n"+strconv.Itoa(w)])
	}
	return p
}

// CoreWindings returns the node pairs of a core's windings.
func CoreWindings(e Element) [][2]int {
	w := [][2]int{e.Nodes}
	for i := 0; i+1 < len(e.Extra); i += 2 {
		w = append(w, [2]int{e.Extra[i], e.Extra[i+1]})
	}
	return w
}

// Core adds a Jiles–Atherton core with windings given as node pairs
// {w1+, w1−, w2+, w2−, …}; params by lower-case name as in the netlist line
// (n1, n2, … are the turns).
func (b *Builder) Core(name string, windings []string, params map[string]float64) *Builder {
	p := make(map[string]float64, len(params))
	for k, v := range params {
		p[strings.ToLower(k)] = v
	}
	if len(windings) < 2 {
		windings = append(windings, "0", "0")
	}
	return b.AddMulti(Element{Name: name, Kind: Core, Params: p}, windings[0], windings[1], windings[2:]...)
}

// MaxCoreWindings bounds the windings on one core.
const MaxCoreWindings = 64

func validateCore(e Element) error {
	if len(e.Extra)%2 != 0 || len(e.Extra)/2+1 > MaxCoreWindings {
		return fmt.Errorf("%w: %s windings must be node pairs (at most %d)", ErrBadValue, e.Name, MaxCoreWindings)
	}
	nw := 1 + len(e.Extra)/2
	for i, w := range CoreWindings(e) {
		if w[0] == w[1] {
			return fmt.Errorf("%w: %s winding %d is shorted on itself", ErrSelfConnected, e.Name, i+1)
		}
	}
	allowed := map[string]bool{"ms": true, "a": true, "k": true, "c": true, "alpha": true, "area": true, "len": true}
	for w := 1; w <= nw; w++ {
		allowed["n"+strconv.Itoa(w)] = true
	}
	for k, v := range e.Params {
		if !allowed[k] {
			return fmt.Errorf("%w: %s: unknown core parameter %q", ErrBadValue, e.Name, k)
		}
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: %s has a non-finite parameter", ErrBadValue, e.Name)
		}
	}
	for _, k := range []string{"ms", "a", "k", "area", "len"} {
		if !(e.Params[k] > 0) {
			return fmt.Errorf("%w: %s needs %s > 0 (no default material)", ErrBadValue, e.Name, strings.ToUpper(k))
		}
	}
	p := CoreParamsOf(e)
	if !(p.C >= 0 && p.C < 1) || p.Alpha < 0 {
		return fmt.Errorf("%w: %s needs 0 ≤ C < 1 and ALPHA ≥ 0", ErrBadValue, e.Name)
	}
	// The anhysteretic slope at the origin is Ms/(3a); the mean-field
	// feedback must not exceed it, or the demagnetized state is unstable.
	if p.Alpha*p.Ms/(3*p.A) >= 1 {
		return fmt.Errorf("%w: %s: ALPHA·MS/(3A) must be < 1", ErrBadValue, e.Name)
	}
	for i, n := range p.Turns {
		if !(n > 0) {
			return fmt.Errorf("%w: %s needs N%d > 0", ErrBadValue, e.Name, i+1)
		}
	}
	return nil
}
