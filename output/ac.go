// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"encoding/csv"
	"io"
	"math"
	"math/cmplx"
	"strconv"
	"strings"
)

// ACLedger is the average-power balance at each frequency. Phasors are peak
// amplitudes (as in SPICE), so average power is ½·Re(V·I*).
type ACLedger struct {
	// Per frequency [W] (reactive power in var), aligned with ACResult.Freq.
	Inflow, Outflow, Dissipation, Reactive, Residual []float64
	ModelSupplied                                    []float64 // average power from idealized active elements
	MaxRelativeResidual                              float64
}

// ACResult is the result of a small-signal AC analysis.
type ACResult struct {
	Title     string
	Freq      []float64 // [Hz]
	NodeNames []string
	// V[n][k] is node n's phasor at Freq[k]; V[0] (ground) is nil.
	V            [][]complex128
	ElementNames []string
	ElementKinds []string
	I            [][]complex128 // element currents n+ → n−
	Ledger       ACLedger
}

// NodeVoltage returns the phasor series of a named node.
func (r *ACResult) NodeVoltage(name string) ([]complex128, bool) {
	key := strings.ToLower(name)
	if key == "0" || key == "gnd" {
		return make([]complex128, len(r.Freq)), true
	}
	for i, n := range r.NodeNames {
		if i > 0 && n == key {
			return r.V[i], true
		}
	}
	return nil, false
}

// ElementCurrent returns the current phasor series of a named element.
func (r *ACResult) ElementCurrent(name string) ([]complex128, bool) {
	for i, n := range r.ElementNames {
		if strings.EqualFold(n, name) {
			return r.I[i], true
		}
	}
	return nil, false
}

// DB returns 20·log10|z|.
func DB(z complex128) float64 { return 20 * math.Log10(cmplx.Abs(z)) }

// PhaseDeg returns arg(z) in degrees.
func PhaseDeg(z complex128) float64 { return cmplx.Phase(z) * 180 / math.Pi }

// WriteACCSV writes freq, then vm/vp (magnitude, phase in degrees) for every
// node, then the ledger's average inflow and dissipation.
func WriteACCSV(w io.Writer, r *ACResult) error {
	cw := csv.NewWriter(w)
	head := []string{"freq"}
	for i, n := range r.NodeNames {
		if i > 0 {
			head = append(head, "vm("+n+")", "vp("+n+")")
		}
	}
	head = append(head, "ledger.inflow", "ledger.dissipation", "ledger.reactive")
	if err := cw.Write(head); err != nil {
		return err
	}
	f := func(x float64) string { return strconv.FormatFloat(x, 'g', 10, 64) }
	row := make([]string, len(head))
	for k, fr := range r.Freq {
		j := 0
		row[j] = f(fr)
		j++
		for i := 1; i < len(r.NodeNames); i++ {
			row[j], row[j+1] = f(cmplx.Abs(r.V[i][k])), f(PhaseDeg(r.V[i][k]))
			j += 2
		}
		row[j], row[j+1], row[j+2] = f(r.Ledger.Inflow[k]), f(r.Ledger.Dissipation[k]), f(r.Ledger.Reactive[k])
		if err := cw.Write(row); err != nil {
			return err
		}
	}
	cw.Flush()
	return cw.Error()
}

// DCSweepResult is the result of a .DC sweep: one operating point per value
// of the swept source(s).
type DCSweepResult struct {
	Title   string
	Sources []string    // swept sources, inner first
	Values  [][]float64 // Values[s][k] is source s's value at point k
	// V[n][k] is node n's voltage at point k; V[0] (ground) is nil.
	NodeNames           []string
	V                   [][]float64
	ElementNames        []string
	I                   [][]float64
	MaxRelativeResidual float64 // worst ledger residual over the sweep
}

// NodeVoltage returns the sweep series of a named node.
func (r *DCSweepResult) NodeVoltage(name string) ([]float64, bool) {
	key := strings.ToLower(name)
	for i, n := range r.NodeNames {
		if i > 0 && n == key {
			return r.V[i], true
		}
	}
	return nil, false
}

// WriteDCCSV writes the swept values, then v(node) and i(element) columns.
func WriteDCCSV(w io.Writer, r *DCSweepResult) error {
	cw := csv.NewWriter(w)
	var head []string
	for _, s := range r.Sources {
		head = append(head, strings.ToLower(s))
	}
	for i, n := range r.NodeNames {
		if i > 0 {
			head = append(head, "v("+n+")")
		}
	}
	for _, n := range r.ElementNames {
		head = append(head, "i("+strings.ToLower(n)+")")
	}
	if err := cw.Write(head); err != nil {
		return err
	}
	f := func(x float64) string { return strconv.FormatFloat(x, 'g', 10, 64) }
	if len(r.Values) == 0 {
		cw.Flush()
		return cw.Error()
	}
	for k := range r.Values[0] {
		row := make([]string, 0, len(head))
		for s := range r.Sources {
			row = append(row, f(r.Values[s][k]))
		}
		for i := 1; i < len(r.NodeNames); i++ {
			row = append(row, f(r.V[i][k]))
		}
		for i := range r.ElementNames {
			row = append(row, f(r.I[i][k]))
		}
		if err := cw.Write(row); err != nil {
			return err
		}
	}
	cw.Flush()
	return cw.Error()
}
