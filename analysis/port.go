// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"fmt"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// PortResult is the open-circuit voltage and short-circuit current of a
// terminal pair, each from a full (nonlinear) DC solve.
type PortResult struct {
	// Voc is V(plus) − V(minus) with nothing connected between them.
	Voc float64
	// Isc is the current the circuit drives out of plus, through an ideal
	// short, into minus.
	Isc float64
}

// OpenShort measures a terminal pair at its two extremes: open circuit and
// short circuit. For a linear circuit these fix its Thévenin/Norton
// equivalent (R = Voc/Isc); for a nonlinear one they are the two end points
// of its output characteristic, and are useful fitness objectives in their
// own right (for example a harvester's open-circuit voltage and short-
// circuit current). The short is a 0 V source added to a copy of c, so c is
// unchanged.
func OpenShort(c *circuit.Circuit, plus, minus string, opt Options) (PortResult, error) {
	var pr PortResult
	r, err := OP(c, opt)
	if err != nil {
		return pr, fmt.Errorf("open circuit: %w", err)
	}
	vp, ok1 := r.Voltage(plus)
	vm, ok2 := r.Voltage(minus)
	if !ok1 || !ok2 {
		return pr, fmt.Errorf("%w: %s or %s", ErrUnknownNode, plus, minus)
	}
	pr.Voc = vp - vm

	sc := c.Clone()
	name := "Vshort_probe"
	for sc.Index(name) >= 0 {
		name += "_"
	}
	ip, im := -1, -1
	for i, n := range sc.NodeNames {
		if strings.EqualFold(n, plus) || (i == 0 && (plus == "0" || strings.EqualFold(plus, "gnd"))) {
			ip = i
		}
		if strings.EqualFold(n, minus) || (i == 0 && (minus == "0" || strings.EqualFold(minus, "gnd"))) {
			im = i
		}
	}
	if ip == im {
		return pr, fmt.Errorf("%w: %s and %s are the same node", ErrUnknownNode, plus, minus)
	}
	sc.Elements = append(sc.Elements, circuit.Element{Name: name, Kind: circuit.VoltageSource, Nodes: [2]int{ip, im}, HasDC: true})
	rs, err := OP(sc, opt)
	if err != nil {
		return pr, fmt.Errorf("short circuit: %w", err)
	}
	e, _ := rs.Element(name)
	// The probe's current flows into its + terminal: out of the circuit at plus.
	pr.Isc = e.Current
	return pr, nil
}
