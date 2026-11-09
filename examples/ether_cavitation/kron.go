// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Kron twin: the same discrete amplitude problem built as a kronspice
// circuit and solved by its operating-point analysis. Node i is a voltage
// a_i (1 V = the ether's amplitude), branch i a resistor H/(2κ2A_i) (Ω),
// node i a voltage-controlled nonlinear shunt i = W_i·V′(v) to ground
// (a piecewise-linear B element; its negative-slope parts are booked as
// model-supplied power in the ledger), fixed ends as voltage sources. The
// circuit's KCL at every node is the ladder's F_i = 0 (ladder.go), so the
// operating point equals the Go solution up to the PWL sampling of V′.
// NodeSet starts Newton at the Go profile: a saddle (the critical bubble)
// is an operating point like any other.

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// KronResult compares the circuit with the Go solution.
type KronResult struct {
	MaxDiff        float64 // max |V_node − a_i|
	Energy         float64 // ladder energy of the circuit's voltages
	LedgerResidual float64 // kronspice relative residual
	ModelSupplied  float64 // W
	Nodes          int
}

// KronTwin builds and solves the ladder circuit for profile a (the Go
// solution), with vref the energy reference density.
func KronTwin(l Ladder, e EOS, a []float64, vref float64) (KronResult, error) {
	var r KronResult
	n := len(a)
	b := circuit.NewBuilder("ether_cavitation: Kron ladder of the soft ether")
	nd := func(i int) string { return fmt.Sprintf("a%d", i) }
	// Each node's curve spans [−0.25, 1.25] at spacing dv with the node's
	// expected value a_i at the middle of a segment: Newton then starts
	// inside one linear piece. Without this, next to a saddle (the critical
	// bubble, one unstable mode) Newton cycles between neighbouring
	// segments and never meets its tolerance.
	const dv = 1.5 / 3000
	for i := 0; i < n; i++ {
		if i < n-1 {
			b.Resistor(fmt.Sprintf("R%d", i), nd(i), nd(i+1), l.H/(2*kappa2*l.A[i]))
		}
		fixed := (i == 0 && l.FixL) || (i == n-1 && l.FixR)
		if fixed {
			b.VoltageSource(fmt.Sprintf("V%d", i), nd(i), "0", a[i])
			continue
		}
		start := a[i] + dv/2 - dv*math.Ceil((a[i]+dv/2+0.25)/dv)
		cv := make([]float64, 0, 6006)
		for j := 0; j <= 3002; j++ {
			v := start + float64(j)*dv
			cv = append(cv, v, l.W[i]*e.D1(v))
		}
		b.NonlinearVC(fmt.Sprintf("B%d", i), nd(i), "0", cv)
	}
	c, err := b.Build()
	if err != nil {
		return r, err
	}
	ns := map[string]float64{}
	for i, v := range a {
		ns[nd(i)] = v
	}
	res, err := analysis.OP(c, analysis.Options{NodeSet: ns, RelTol: 1e-9, VnTol: 1e-12, ITL1: 100, NoGminStepping: true, NoSourceStepping: true, NoPseudoTransient: true})
	if err != nil {
		return r, err
	}
	got := make([]float64, n)
	for i := range a {
		v, _ := res.Voltage(nd(i))
		got[i] = v
		r.MaxDiff = math.Max(r.MaxDiff, math.Abs(v-a[i]))
	}
	r.Energy = l.Energy(e, got, vref)
	if res.Ledger != nil {
		r.LedgerResidual = res.Ledger.RelativeResidual()
		r.ModelSupplied = res.Ledger.ModelSupplied
	}
	r.Nodes = n
	return r, nil
}
