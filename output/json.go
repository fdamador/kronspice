// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"encoding/json"
	"io"
	"math/cmplx"
	"strings"
)

// Document is the machine-readable form of a simulation run, written by
// `kronspice -json`. It carries every analysis the netlist requested, in
// order, and on failure the error with its class so that a program (or an
// AI assistant) can act on it without parsing prose. Field names are stable
// and snake_case; new fields may be added, existing ones are not renamed.
type Document struct {
	Tool    string `json:"tool"`    // "kronspice"
	Version string `json:"version"` // build version, "dev" when unset
	OK      bool   `json:"ok"`
	// Error is set when OK is false. Analyses completed before the failure
	// are still listed.
	Error    *ErrorInfo    `json:"error,omitempty"`
	Title    string        `json:"title,omitempty"`
	Analyses []AnalysisDoc `json:"analyses"`
}

// ErrorInfo names the failure class (the same names evo uses for penalty
// classes: syntax, invalid, singular, no-convergence, ledger, timeout,
// error) and carries the error text.
type ErrorInfo struct {
	Class   string `json:"class"`
	Message string `json:"message"`
}

// AnalysisDoc is one analysis; exactly one of the pointers is set.
type AnalysisDoc struct {
	Kind string   `json:"kind"` // "op", "dc", "ac", "tran"
	OP   *OPDoc   `json:"op,omitempty"`
	DC   *DCDoc   `json:"dc,omitempty"`
	AC   *ACDoc   `json:"ac,omitempty"`
	Tran *TranDoc `json:"tran,omitempty"`
	// Smith is set on an "ac" analysis when a port was named (kronspice
	// -smith): the port's impedance and reflection over the same sweep.
	Smith *SmithDoc `json:"smith,omitempty"`
}

// OPDoc is an operating point.
type OPDoc struct {
	Nodes    map[string]float64 `json:"nodes"` // node name → voltage [V]; ground omitted
	Elements []ElementDoc       `json:"elements"`
	Ledger   *LedgerDoc         `json:"ledger,omitempty"`
}

// ElementDoc is one element at the operating point.
type ElementDoc struct {
	Name          string    `json:"name"`
	Kind          string    `json:"kind"`
	Nodes         [2]string `json:"nodes"`
	Voltage       float64   `json:"voltage"`        // V(n+) − V(n−) [V]
	Current       float64   `json:"current"`        // n+ → n− [A]
	PowerAbsorbed float64   `json:"power_absorbed"` // [W], negative when supplying
	StoredEnergy  float64   `json:"stored_energy"`  // [J]
	ModelSupplied float64   `json:"model_supplied"` // [W], idealized active elements only
}

// LedgerDoc is the operating-point energy ledger (CLAUDE.md §3c).
type LedgerDoc struct {
	Inflow           float64        `json:"inflow"`
	Outflow          float64        `json:"outflow"`
	Dissipation      float64        `json:"dissipation"`
	StorageRate      float64        `json:"storage_rate"`
	ModelSupplied    float64        `json:"model_supplied"`
	Hypothetical     float64        `json:"hypothetical"`
	Residual         float64        `json:"residual"`
	RelativeResidual float64        `json:"relative_residual"`
	Rows             []LedgerRowDoc `json:"rows"` // modeled and hypothetical ports only
	Supplied         []SuppliedRow  `json:"supplied,omitempty"`
}

// LedgerRowDoc is one element × channel entry.
type LedgerRowDoc struct {
	Element string  `json:"element"`
	Channel string  `json:"channel"`
	Status  string  `json:"status"`
	PowerW  float64 `json:"power_w"`
	Note    string  `json:"note,omitempty"`
}

// DCDoc is a .DC sweep.
type DCDoc struct {
	Sources             []string             `json:"sources"` // swept sources, inner first
	Values              [][]float64          `json:"values"`  // values[s][k]
	V                   map[string][]float64 `json:"v"`       // node → series
	I                   map[string][]float64 `json:"i"`       // element → series
	MaxRelativeResidual float64              `json:"max_relative_residual"`
}

// ACDoc is a small-signal sweep. Phasors are peak amplitudes, as in SPICE.
type ACDoc struct {
	Freq   []float64               `json:"freq"` // [Hz]
	V      map[string]PhasorSeries `json:"v"`
	I      map[string]PhasorSeries `json:"i"`
	Ledger ACLedgerDoc             `json:"ledger"`
}

// PhasorSeries is a complex series as magnitude and phase in degrees.
type PhasorSeries struct {
	Mag      []float64 `json:"mag"`
	PhaseDeg []float64 `json:"phase_deg"`
}

// ACLedgerDoc is the average-power balance per frequency.
type ACLedgerDoc struct {
	Inflow              []float64 `json:"inflow"`
	Outflow             []float64 `json:"outflow"`
	Dissipation         []float64 `json:"dissipation"`
	Reactive            []float64 `json:"reactive"`
	ModelSupplied       []float64 `json:"model_supplied"`
	Residual            []float64 `json:"residual"`
	MaxRelativeResidual float64   `json:"max_relative_residual"`
}

// TranDoc is a transient. When Step > 0 the series are resampled on that
// uniform grid (the .tran print step), as the CSV writer does; otherwise
// they are the solver's accepted time points.
type TranDoc struct {
	Step   float64              `json:"step"`
	Time   []float64            `json:"time"`
	V      map[string][]float64 `json:"v"`
	I      map[string][]float64 `json:"i"`
	Ledger TranLedgerDoc        `json:"ledger"`
	Stats  TranStatsDoc         `json:"stats"`
}

// TranLedgerDoc is the transient ledger: per-point powers on the same grid
// as TranDoc.Time, and energies over the run.
type TranLedgerDoc struct {
	Inflow        []float64 `json:"inflow"`
	Outflow       []float64 `json:"outflow"`
	Dissipation   []float64 `json:"dissipation"`
	StorageRate   []float64 `json:"storage_rate"`
	ModelSupplied []float64 `json:"model_supplied"`
	StoredEnergy  []float64 `json:"stored_energy"`
	Residual      []float64 `json:"residual"`

	EnergyIn            float64 `json:"energy_in"`
	EnergyOut           float64 `json:"energy_out"`
	EnergyDissipated    float64 `json:"energy_dissipated"`
	EnergyModelSupplied float64 `json:"energy_model_supplied"`
	StoredStart         float64 `json:"stored_start"`
	StoredEnd           float64 `json:"stored_end"`
	IntegrationMismatch float64 `json:"integration_mismatch"`
	MaxRelativeResidual float64 `json:"max_relative_residual"`
}

// TranStatsDoc describes the time stepping.
type TranStatsDoc struct {
	Accepted int     `json:"accepted"`
	Rejected int     `json:"rejected"`
	MinStep  float64 `json:"min_step"`
	MaxStep  float64 `json:"max_step"`
}

// OPDocument converts an operating-point result.
func OPDocument(r *Result) *OPDoc {
	d := &OPDoc{Nodes: map[string]float64{}}
	for i, n := range r.NodeNames {
		if i > 0 {
			d.Nodes[n] = r.NodeVoltages[i]
		}
	}
	for _, e := range r.Elements {
		d.Elements = append(d.Elements, ElementDoc{
			Name: e.Name, Kind: e.Kind, Nodes: e.Nodes,
			Voltage: e.Voltage, Current: e.Current, PowerAbsorbed: e.PowerAbsorbed,
			StoredEnergy: e.StoredEnergy, ModelSupplied: e.ModelSupplied,
		})
	}
	if l := r.Ledger; l != nil {
		ld := &LedgerDoc{
			Inflow: l.Inflow, Outflow: l.Outflow, Dissipation: l.Dissipation,
			StorageRate: l.StorageRate, ModelSupplied: l.ModelSupplied,
			Hypothetical: l.Hypothetical, Residual: l.Residual,
			RelativeResidual: l.RelativeResidual(), Supplied: l.Supplied,
			Rows: []LedgerRowDoc{},
		}
		for _, row := range l.Rows {
			if row.Status == 0 { // not modeled: exactly 0 W, listed by -all-channels only
				continue
			}
			ld.Rows = append(ld.Rows, LedgerRowDoc{
				Element: row.Element, Channel: row.Channel.String(), Status: row.Status.String(),
				PowerW: row.PowerW, Note: row.Note,
			})
		}
		d.Ledger = ld
	}
	return d
}

// DCDocument converts a sweep result.
func DCDocument(r *DCSweepResult) *DCDoc {
	d := &DCDoc{Sources: r.Sources, Values: r.Values, V: map[string][]float64{}, I: map[string][]float64{}, MaxRelativeResidual: r.MaxRelativeResidual}
	if d.Sources == nil {
		d.Sources = []string{}
	}
	if d.Values == nil {
		d.Values = [][]float64{}
	}
	for i, n := range r.NodeNames {
		if i > 0 {
			d.V[n] = r.V[i]
		}
	}
	for i, n := range r.ElementNames {
		d.I[strings.ToLower(n)] = r.I[i]
	}
	return d
}

// ACDocument converts an AC result.
func ACDocument(r *ACResult) *ACDoc {
	series := func(z []complex128) PhasorSeries {
		s := PhasorSeries{Mag: make([]float64, len(z)), PhaseDeg: make([]float64, len(z))}
		for k, v := range z {
			s.Mag[k], s.PhaseDeg[k] = cmplx.Abs(v), PhaseDeg(v)
		}
		return s
	}
	d := &ACDoc{Freq: r.Freq, V: map[string]PhasorSeries{}, I: map[string]PhasorSeries{}}
	for i, n := range r.NodeNames {
		if i > 0 {
			d.V[n] = series(r.V[i])
		}
	}
	for i, n := range r.ElementNames {
		d.I[strings.ToLower(n)] = series(r.I[i])
	}
	l := r.Ledger
	d.Ledger = ACLedgerDoc{
		Inflow: orEmpty(l.Inflow), Outflow: orEmpty(l.Outflow), Dissipation: orEmpty(l.Dissipation),
		Reactive: orEmpty(l.Reactive), ModelSupplied: orEmpty(l.ModelSupplied), Residual: orEmpty(l.Residual),
		MaxRelativeResidual: l.MaxRelativeResidual,
	}
	return d
}

// TranDocument converts waveforms, resampled on step when step > 0.
func TranDocument(wf *Waveforms, step float64) *TranDoc {
	d := &TranDoc{Step: step, V: map[string][]float64{}, I: map[string][]float64{},
		Stats: TranStatsDoc{Accepted: wf.Stats.Accepted, Rejected: wf.Stats.Rejected, MinStep: wf.Stats.MinStep, MaxStep: wf.Stats.MaxStep}}
	grid := wf.Time
	if step > 0 && len(wf.Time) > 1 {
		t0, t1 := wf.Time[0], wf.Time[len(wf.Time)-1]
		n := int((t1-t0)/step+1e-9) + 1
		grid = make([]float64, n)
		for k := range grid {
			tq := t0 + float64(k)*step
			if tq > t1 {
				tq = t1
			}
			grid[k] = tq
		}
	}
	sample := func(y []float64) []float64 {
		if y == nil {
			return []float64{}
		}
		if &grid[0] == &wf.Time[0] && len(grid) == len(wf.Time) {
			return y
		}
		out := make([]float64, len(grid))
		for k, tq := range grid {
			out[k] = Interp(wf.Time, y, tq)
		}
		return out
	}
	d.Time = grid
	for i, n := range wf.NodeNames {
		if i > 0 {
			d.V[n] = sample(wf.V[i])
		}
	}
	for i, n := range wf.ElementNames {
		d.I[strings.ToLower(n)] = sample(wf.I[i])
	}
	l := &wf.Ledger
	d.Ledger = TranLedgerDoc{
		Inflow: sample(l.Inflow), Outflow: sample(l.Outflow), Dissipation: sample(l.Dissipation),
		StorageRate: sample(l.StorageRate), ModelSupplied: sample(l.ModelSupplied),
		StoredEnergy: sample(l.StoredEnergy), Residual: sample(l.Residual),
		EnergyIn: l.EnergyIn, EnergyOut: l.EnergyOut, EnergyDissipated: l.EnergyDissipated,
		EnergyModelSupplied: l.EnergyModelSupplied, StoredStart: l.StoredStart, StoredEnd: l.StoredEnd,
		IntegrationMismatch: l.IntegrationMismatch, MaxRelativeResidual: l.MaxRelativeResidual,
	}
	return d
}

func orEmpty(s []float64) []float64 {
	if s == nil {
		return []float64{}
	}
	return s
}

// WriteJSON writes the document, indented, followed by a newline.
func WriteJSON(w io.Writer, d *Document) error {
	if d.Analyses == nil {
		d.Analyses = []AnalysisDoc{}
	}
	enc := json.NewEncoder(w)
	enc.SetIndent("", "  ")
	return enc.Encode(d)
}
