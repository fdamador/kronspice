// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"encoding/csv"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"
)

// TranLedger is the energy ledger of a transient analysis: a power balance at
// every time point plus energies integrated over the run (CLAUDE.md §3c).
type TranLedger struct {
	// Per time point [W]: aligned with Waveforms.Time.
	Inflow, Outflow, Dissipation, StorageRate, Residual []float64
	// ModelSupplied is power from idealized active elements (negative
	// resistances), kept apart from Inflow (see Ledger.ModelSupplied).
	ModelSupplied []float64
	// StoredEnergy is ½Cv² + ½Li² summed over the circuit at each point [J].
	StoredEnergy []float64

	// Energies over the recorded interval [J], integrated with the
	// trapezoidal rule over the time points.
	EnergyIn, EnergyOut, EnergyDissipated float64
	EnergyModelSupplied                   float64
	StoredStart, StoredEnd                float64
	// IntegrationMismatch = EnergyIn + EnergyModelSupplied − EnergyOut −
	// EnergyDissipated −
	// (StoredEnd − StoredStart). It measures time-integration error (for
	// example backward Euler's numerical damping), not missing physics.
	IntegrationMismatch float64
	// MaxRelativeResidual is the largest per-point |residual| divided by the
	// whole run's peak Σ|power|.
	MaxRelativeResidual float64
}

// TranStats describes how the time stepping went.
type TranStats struct {
	Accepted, Rejected int
	MinStep, MaxStep   float64
}

// Waveforms is the result of a transient analysis.
type Waveforms struct {
	Title     string
	Time      []float64
	NodeNames []string
	// V[n][k] is the voltage of node n at Time[k]; V[0] (ground) is nil.
	V            [][]float64
	ElementNames []string
	ElementKinds []string
	// I[e][k] and P[e][k] are element e's current (n+ → n−) and absorbed power.
	I, P   [][]float64
	Ledger TranLedger
	Stats  TranStats
}

// NodeVoltage returns the time series of a named node.
func (w *Waveforms) NodeVoltage(name string) ([]float64, bool) {
	key := strings.ToLower(name)
	if key == "gnd" || key == "0" {
		return make([]float64, len(w.Time)), true
	}
	for i, n := range w.NodeNames {
		if n == key && i > 0 {
			return w.V[i], true
		}
	}
	return nil, false
}

// ElementCurrent returns the current time series of a named element.
func (w *Waveforms) ElementCurrent(name string) ([]float64, bool) {
	for i, n := range w.ElementNames {
		if strings.EqualFold(n, name) {
			return w.I[i], true
		}
	}
	return nil, false
}

// Interp returns y at time tq by linear interpolation over the time points t.
func Interp(t, y []float64, tq float64) float64 {
	n := len(t)
	if n == 0 {
		return math.NaN()
	}
	if tq <= t[0] {
		return y[0]
	}
	if tq >= t[n-1] {
		return y[n-1]
	}
	k := sort.SearchFloat64s(t, tq)
	if t[k] == tq {
		return y[k]
	}
	t0, t1 := t[k-1], t[k]
	return y[k-1] + (y[k]-y[k-1])*(tq-t0)/(t1-t0)
}

// WriteTranCSV writes the waveforms resampled on a uniform grid of the given
// step (as SPICE prints them): time, v(node)…, i(element)…, then the ledger.
func WriteTranCSV(w io.Writer, wf *Waveforms, step float64) error {
	cw := csv.NewWriter(w)
	head := []string{"time"}
	for i, n := range wf.NodeNames {
		if i > 0 {
			head = append(head, "v("+n+")")
		}
	}
	for _, n := range wf.ElementNames {
		head = append(head, "i("+strings.ToLower(n)+")")
	}
	head = append(head, "ledger.inflow", "ledger.dissipation", "ledger.storage_rate", "ledger.stored_energy", "ledger.residual")
	if err := cw.Write(head); err != nil {
		return err
	}
	if len(wf.Time) == 0 {
		cw.Flush()
		return cw.Error()
	}
	f := func(x float64) string { return strconv.FormatFloat(x, 'g', 10, 64) }
	t0, t1 := wf.Time[0], wf.Time[len(wf.Time)-1]
	n := int(math.Floor((t1-t0)/step+1e-9)) + 1
	row := make([]string, len(head))
	for k := 0; k < n; k++ {
		tq := t0 + float64(k)*step
		if tq > t1 {
			tq = t1
		}
		j := 0
		row[j] = f(tq)
		j++
		for i := 1; i < len(wf.NodeNames); i++ {
			row[j] = f(Interp(wf.Time, wf.V[i], tq))
			j++
		}
		for i := range wf.ElementNames {
			row[j] = f(Interp(wf.Time, wf.I[i], tq))
			j++
		}
		l := &wf.Ledger
		for _, s := range [][]float64{l.Inflow, l.Dissipation, l.StorageRate, l.StoredEnergy, l.Residual} {
			row[j] = f(Interp(wf.Time, s, tq))
			j++
		}
		if err := cw.Write(row); err != nil {
			return err
		}
	}
	cw.Flush()
	return cw.Error()
}
