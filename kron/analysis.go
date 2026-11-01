// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// OP solves the operating point of a linear circuit in the chosen form and
// returns it as the same output.Result the MNA analysis returns, so the two
// can be compared directly (CLAUDE.md §3a rule 4). Node voltages, and every
// element's voltage, current and absorbed power, are filled in; the energy
// ledger is not (it is the subject of M9a.1).
//
// A nonlinear circuit is solved in the node-pair frame with state-
// dependent primitive tensors (NonlinearOP); the mesh frame cannot take
// a voltage-controlled nonlinearity and returns ErrNonlinear.
func OP(c *circuit.Circuit, form Form) (*output.Result, error) {
	n, err := Build(c)
	if errors.Is(err, ErrNonlinear) && form == NodeForm {
		return NonlinearOP(c, NonlinearOptions{})
	}
	if err != nil {
		return nil, err
	}
	sol, err := n.Solve(form, 0, true)
	if err != nil {
		return nil, err
	}
	return n.opResult(sol)
}

// TornOP is OP solved by tearing (diakoptics.go); cache may be nil.
func TornOP(c *circuit.Circuit, t Tearing, cache *TornCache) (*output.Result, error) {
	d, err := Tear(c, t, cache)
	if err != nil {
		return nil, err
	}
	sol, _, err := d.Solve(0, true)
	if err != nil {
		return nil, err
	}
	return d.n.opResult(sol)
}

func (n *Network) opResult(sol *Solution) (*output.Result, error) {
	c := n.c
	v, err := n.voltages(sol)
	if err != nil {
		return nil, err
	}
	res := &output.Result{Title: c.Title, Analysis: "op",
		NodeNames: append([]string(nil), c.NodeNames...), NodeVoltages: make([]float64, c.NumNodes())}
	for i := range res.NodeVoltages {
		res.NodeVoltages[i] = real(v[i])
	}
	for ei, e := range c.Elements {
		var E, I, P float64
		if k := n.elemBranch[ei]; k >= 0 {
			E, I = real(sol.E.At(k)), real(sol.I.At(k))
			P = real(n.absorbed(sol, k))
		}
		res.Elements = append(res.Elements, output.ElementResult{Name: e.Name, Kind: e.Kind.String(),
			Nodes:   [2]string{c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]},
			Voltage: E, Current: I, PowerAbsorbed: P})
	}
	return res, nil
}

// AC solves the small-signal response of a linear circuit at each
// frequency [Hz] in the chosen form, returning node and element phasors as
// output.ACResult (the ledger is left empty). A nonlinear circuit is
// linearized at its operating point in the node-pair frame (NonlinearAC).
func AC(c *circuit.Circuit, freqs []float64, form Form) (*output.ACResult, error) {
	n, err := Build(c)
	if errors.Is(err, ErrNonlinear) && form == NodeForm {
		return NonlinearAC(c, freqs, NonlinearOptions{})
	}
	if err != nil {
		return nil, err
	}
	return n.acResult(freqs, func(w float64) (*Solution, error) { return n.Solve(form, w, false) })
}

// TornAC is AC solved by tearing (diakoptics.go); cache may be nil and
// may be shared between goroutines.
func TornAC(c *circuit.Circuit, freqs []float64, t Tearing, cache *TornCache) (*output.ACResult, error) {
	d, err := Tear(c, t, cache)
	if err != nil {
		return nil, err
	}
	return d.n.acResult(freqs, func(w float64) (*Solution, error) {
		sol, _, err := d.Solve(w, false)
		return sol, err
	})
}

func (n *Network) acResult(freqs []float64, solve func(omega float64) (*Solution, error)) (*output.ACResult, error) {
	c := n.c
	res := &output.ACResult{Title: c.Title, Freq: append([]float64(nil), freqs...),
		NodeNames: append([]string(nil), c.NodeNames...), V: make([][]complex128, c.NumNodes())}
	for i := 1; i < c.NumNodes(); i++ {
		res.V[i] = make([]complex128, len(freqs))
	}
	for _, e := range c.Elements {
		res.ElementNames = append(res.ElementNames, e.Name)
		res.ElementKinds = append(res.ElementKinds, e.Kind.String())
		res.I = append(res.I, make([]complex128, len(freqs)))
	}
	for fi, f := range freqs {
		if !(f > 0) || math.IsInf(f, 0) {
			return nil, fmt.Errorf("%w: frequency %g", ErrUnsupported, f)
		}
		sol, err := solve(2 * math.Pi * f)
		if err != nil {
			return nil, fmt.Errorf("at %g Hz: %w", f, err)
		}
		v, err := n.voltages(sol)
		if err != nil {
			return nil, err
		}
		for i := 1; i < c.NumNodes(); i++ {
			res.V[i][fi] = v[i]
		}
		for ei, k := range n.elemBranch {
			if k >= 0 {
				res.I[ei][fi] = sol.I.At(k)
			}
		}
	}
	return res, nil
}

// absorbed is the instantaneous (DC) power absorbed by the element on
// branch k, E·i, including a transducer's second port.
func (n *Network) absorbed(sol *Solution, k int) complex128 {
	p := sol.E.At(k) * sol.I.At(k)
	if n.branches[k].kind == bGyr1 {
		k2 := n.branches[k].ctrl
		p += sol.E.At(k2) * sol.I.At(k2)
	}
	return p
}

// absorbedAC is the average power ½·Re Σ E·i* of the element on branch k.
func (n *Network) absorbedAC(sol *Solution, k int) float64 {
	p := 0.5 * real(sol.E.At(k)*cmplx.Conj(sol.I.At(k)))
	if n.branches[k].kind == bGyr1 {
		k2 := n.branches[k].ctrl
		p += 0.5 * real(sol.E.At(k2)*cmplx.Conj(sol.I.At(k2)))
	}
	return p
}

// voltages returns the node potentials of a solution.
func (n *Network) voltages(sol *Solution) ([]complex128, error) {
	return n.potentials(sol.E)
}
