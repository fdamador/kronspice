// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// Efficiency and Tellegen's theorem from the open-system ledger (CLAUDE.md
// §3c "Open systems and the energy ledger": efficiency = useful outflow /
// inflow as a GA fitness source).

// ErrNoInflow means the ledger shows no net energy entering, so an
// efficiency is undefined.
var ErrNoInflow = errors.New("em: no inflow (efficiency undefined)")

// OPEfficiency is Σ P(load) / net inflow at an operating point: the power
// the named elements absorb over the net power the exterior delivers
// through its non-thermal ports, Inflow − Outflow from the ledger (a source
// that absorbs power counts against it; the model-supplied power of
// idealized elements is never inflow).
func OPEfficiency(r *output.Result, loads ...string) (float64, error) {
	net := r.Ledger.Inflow - r.Ledger.Outflow
	if net <= 0 {
		return 0, ErrNoInflow
	}
	p := 0.0
	for _, n := range loads {
		el, ok := r.Element(n)
		if !ok {
			return 0, fmt.Errorf("em: OPEfficiency: no element %q", n)
		}
		p += el.PowerAbsorbed
	}
	return p / net, nil
}

// TranEfficiency is the energy the named elements absorb over [from, to]
// divided by the net energy the exterior delivered in the same interval,
// ∫(Inflow − Outflow) (trapezoidal over the solver's points, interpolated
// at the ends): energy a reactive load hands back to its source within a
// cycle is not counted as delivered. For a periodic steady state over
// whole periods the stored energy returns to its start, so this is the
// cycle efficiency. It also returns both energies [J].
func TranEfficiency(wf *output.Waveforms, from, to float64, loads ...string) (eta, eLoad, eIn float64, err error) {
	if !(to > from) {
		return 0, 0, 0, fmt.Errorf("em: TranEfficiency: empty window [%g, %g]", from, to)
	}
	idx := make([]int, 0, len(loads))
	for _, n := range loads {
		k := -1
		for i, e := range wf.ElementNames {
			if strings.EqualFold(e, n) {
				k = i
			}
		}
		if k < 0 {
			return 0, 0, 0, fmt.Errorf("em: TranEfficiency: no element %q", n)
		}
		idx = append(idx, k)
	}
	load := make([]float64, len(wf.Time))
	for _, k := range idx {
		for n := range load {
			load[n] += wf.P[k][n]
		}
	}
	eLoad = integrate(wf.Time, load, from, to)
	eIn = integrate(wf.Time, wf.Ledger.Inflow, from, to) - integrate(wf.Time, wf.Ledger.Outflow, from, to)
	if eIn <= 0 {
		return 0, eLoad, eIn, ErrNoInflow
	}
	return eLoad / eIn, eLoad, eIn, nil
}

func integrate(t, y []float64, from, to float64) float64 {
	s := 0.0
	for k := 1; k < len(t); k++ {
		t0, t1 := t[k-1], t[k]
		if t1 <= from || t0 >= to {
			continue
		}
		a, b := math.Max(t0, from), math.Min(t1, to)
		s += 0.5 * (output.Interp(t, y, a) + output.Interp(t, y, b)) * (b - a)
	}
	return s
}

// Tellegen returns Σ_k v_k·i_k over the elements of c, with v_k = V(n+) −
// V(n−) from the node voltages nodeV (indexed like c.NodeNames, ground
// first) and i_k from elemI (indexed like c.Elements), and the scale
// Σ|v_k·i_k|. Tellegen's theorem (B. D. H. Tellegen, "A general network
// theorem, with applications," Philips Res. Rep. 7, 1952) makes the sum
// zero for any voltages satisfying KVL and any currents satisfying KCL on
// the same graph — even from two different circuits, or two instants —
// because Σ v_k·i_k = Σ_n V_n·(net current leaving node n). It is defined
// here for circuits of two-terminal elements; a multi-terminal or coupled
// element gives an error.
func Tellegen(c *circuit.Circuit, nodeV, elemI []float64) (sum, scale float64, err error) {
	if len(nodeV) != c.NumNodes() || len(elemI) != len(c.Elements) {
		return 0, 0, fmt.Errorf("em: Tellegen: %d node voltages for %d nodes, %d currents for %d elements",
			len(nodeV), c.NumNodes(), len(elemI), len(c.Elements))
	}
	for k, e := range c.Elements {
		if len(e.Extra) > 0 || e.Kind == circuit.Mutual {
			return 0, 0, fmt.Errorf("em: Tellegen: %s is not a two-terminal element", e.Name)
		}
		p := (nodeV[e.Nodes[0]] - nodeV[e.Nodes[1]]) * elemI[k]
		sum += p
		scale += math.Abs(p)
	}
	return sum, scale, nil
}
