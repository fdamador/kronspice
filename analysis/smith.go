// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// ErrBadReference is returned for a Smith reference impedance that is not
// positive and finite.
var ErrBadReference = errors.New("analysis: Smith reference impedance must be positive and finite")

// Smith sweeps the driving-point impedance of the port (plus, minus) and
// maps it onto the Smith chart against the real reference z0:
// Γ = (Z − Z0)/(Z + Z0).
//
// The port is measured as a vector network analyser would: every
// independent AC source in a copy of c is set to zero (DC values are kept,
// so nonlinear devices are linearized about the same operating point) and
// a 1 A AC current probe is driven into plus and out of minus, so that
// Z = V(plus) − V(minus). c is unchanged.
//
// Negative resistance is not clipped. Where Re Z < 0 the point has
// |Γ| > 1 and Region SmithActive, and the point carries the AC ledger's
// ModelSupplied and Inflow columns at that frequency, which say where the
// returned power comes from (the ledger itself must close, as in AC).
// output.WriteSmithSVG draws such points on the extended chart.
func Smith(c *circuit.Circuit, plus, minus string, z0 float64, spec ACSpec, opt Options) (*output.SmithResult, error) {
	if !(z0 > 0) || math.IsInf(z0, 0) {
		return nil, fmt.Errorf("%w: %g", ErrBadReference, z0)
	}
	ip, im := c.NodeIndex(plus), c.NodeIndex(minus)
	if ip < 0 || im < 0 {
		return nil, fmt.Errorf("%w: %s or %s", ErrUnknownNode, plus, minus)
	}
	if ip == im {
		return nil, fmt.Errorf("%w: %s and %s are the same node", ErrUnknownNode, plus, minus)
	}
	pc := c.Clone()
	for i := range pc.Elements {
		pc.Elements[i].ACMag = 0
	}
	name := "Ismith_probe"
	for pc.Index(name) >= 0 {
		name += "_"
	}
	// SPICE convention: the current flows from n+ through the source to
	// n−, so Nodes = {minus, plus} delivers it into plus.
	pc.Elements = append(pc.Elements, circuit.Element{Name: name, Kind: circuit.CurrentSource, Nodes: [2]int{im, ip}, ACMag: 1})
	r, err := AC(pc, spec, opt)
	if err != nil {
		return nil, err
	}
	vp, _ := r.NodeVoltage(strings.ToLower(pc.NodeNames[ip]))
	vm, _ := r.NodeVoltage(strings.ToLower(pc.NodeNames[im]))
	res := &output.SmithResult{
		Title: c.Title, Plus: plus, Minus: minus, Z0: z0,
		Points:              make([]output.SmithPoint, len(r.Freq)),
		MaxRelativeResidual: r.Ledger.MaxRelativeResidual,
	}
	for k, f := range r.Freq {
		z := vp[k] - vm[k]
		res.Points[k] = output.SmithPoint{
			Freq: f, Z: z, Gamma: output.SmithGamma(z, z0),
			Region:        output.ClassifySmith(z, z0, 0),
			PortPower:     0.5 * real(z),
			ModelSupplied: r.Ledger.ModelSupplied[k],
			Inflow:        r.Ledger.Inflow[k],
		}
	}
	return res, nil
}
