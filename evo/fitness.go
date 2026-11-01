// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

// Objective is one term of a weighted-sum fitness. Eval returns a
// non-negative error (0 = target met); the engine minimizes
// Σ weight·error. Objectives of one candidate share an Eval cache, so
// several measurements on one analysis run it once.
type Objective struct {
	Name   string
	Weight float64
	Eval   func(e *Eval) (float64, error)
	// Constraint makes the term a constraint under NSGA2: its error is a
	// violation (0 = satisfied), it is not one of the Pareto objectives,
	// and a candidate that violates any constraint is dominated by every
	// candidate that satisfies all of them (Deb et al. 2002, §VI). Under
	// WeightedSum a constraint is an ordinary weighted term.
	Constraint bool
}

// Eval is the per-candidate context handed to objectives: the decoded
// circuit, the evaluation deadline and a cache of analysis results.
type Eval struct {
	Ctx     context.Context
	Circuit *circuit.Circuit
	Values  []float64 // decoded gene values, in Space order (part values for topology genomes)
	// Parts is the number of evolved parts after pruning (topology
	// genomes), or −1 for value genomes.
	Parts int
	// Base options for every analysis (Context is filled in).
	Options analysis.Options

	op   *output.Result
	ac   map[analysis.ACSpec]*output.ACResult
	tran map[analysis.TranSpec]*output.Waveforms
	port map[[2]string]analysis.PortResult
	memo map[string]memoEntry
	tear *Tear

	// ModelSupplied is the largest power delivered by idealized active
	// elements (negative resistances, controlled sources) in any analysis
	// run for this candidate [W]. It is reported with every individual so a
	// design that relies on one is always visible (CLAUDE.md §3c).
	ModelSupplied float64
}

type memoEntry struct {
	v   any
	err error
}

// Memo runs f once per candidate and key, so several objectives can share
// an expensive measurement of their own (for example a resonance search).
func (e *Eval) Memo(key string, f func() (any, error)) (any, error) {
	if m, ok := e.memo[key]; ok {
		return m.v, m.err
	}
	v, err := f()
	if e.memo == nil {
		e.memo = map[string]memoEntry{}
	}
	e.memo[key] = memoEntry{v, err}
	return v, err
}

// AnalysisOptions returns the analysis options for this candidate, with
// its deadline set; custom objectives pass them to every analysis they run.
func (e *Eval) AnalysisOptions() analysis.Options { return e.opts() }

func (e *Eval) opts() analysis.Options {
	o := e.Options
	o.Context = e.Ctx
	return o
}

func (e *Eval) noteSupplied(p float64) {
	if p > e.ModelSupplied {
		e.ModelSupplied = p
	}
}

// OP returns the cached operating point.
func (e *Eval) OP() (*output.Result, error) {
	if e.op != nil {
		return e.op, nil
	}
	r, err := analysis.OP(e.Circuit, e.opts())
	if err != nil {
		return nil, err
	}
	e.noteSupplied(r.Ledger.ModelSupplied)
	e.op = r
	return r, nil
}

// AC returns the cached AC sweep for spec.
func (e *Eval) AC(spec analysis.ACSpec) (*output.ACResult, error) {
	if r, ok := e.ac[spec]; ok {
		return r, nil
	}
	r, err := e.tornAC(spec)
	if r == nil && err == nil {
		r, err = analysis.AC(e.Circuit, spec, e.opts())
	}
	if err != nil {
		return nil, err
	}
	for _, p := range r.Ledger.ModelSupplied {
		e.noteSupplied(p)
	}
	if e.ac == nil {
		e.ac = map[analysis.ACSpec]*output.ACResult{}
	}
	e.ac[spec] = r
	return r, nil
}

// tornAC solves an AC sweep by diakoptics when the problem tears and the
// candidate allows it; (nil, nil) means "use MNA" (see Tear).
func (e *Eval) tornAC(spec analysis.ACSpec) (*output.ACResult, error) {
	if e.tear == nil || !passiveLinear(e.Circuit) {
		return nil, nil
	}
	if err := e.Ctx.Err(); err != nil {
		return nil, err
	}
	freqs, err := spec.Frequencies()
	if err != nil {
		return nil, err
	}
	r, err := kron.TornAC(e.Circuit, freqs, kron.Tearing{Torn: e.tear.Torn}, e.tear.Cache)
	if err != nil {
		return nil, nil // MNA solves it or reports the reason
	}
	return r, nil
}

// passiveLinear reports whether every element is a non-negative R, a C,
// L or K, or an independent source: then no power is model-supplied and
// the torn result carries everything the MNA result would.
func passiveLinear(c *circuit.Circuit) bool {
	for _, el := range c.Elements {
		switch el.Kind {
		case circuit.Resistor:
			if el.Value < 0 {
				return false
			}
		case circuit.Capacitor, circuit.Inductor, circuit.Mutual, circuit.VoltageSource, circuit.CurrentSource:
		default:
			return false
		}
	}
	return true
}

// Tran returns the cached transient for spec.
func (e *Eval) Tran(spec analysis.TranSpec) (*output.Waveforms, error) {
	if r, ok := e.tran[spec]; ok {
		return r, nil
	}
	c := e.Circuit
	r, err := analysis.Tran(c, spec, analysis.TranOptions{Options: e.opts(), MaxSteps: 200000})
	if err != nil {
		return nil, err
	}
	for _, p := range r.Ledger.ModelSupplied {
		e.noteSupplied(p)
	}
	if e.tran == nil {
		e.tran = map[analysis.TranSpec]*output.Waveforms{}
	}
	e.tran[spec] = r
	return r, nil
}

// Port returns the cached open-circuit voltage and short-circuit current.
func (e *Eval) Port(plus, minus string) (analysis.PortResult, error) {
	k := [2]string{plus, minus}
	if r, ok := e.port[k]; ok {
		return r, nil
	}
	r, err := analysis.OpenShort(e.Circuit, plus, minus, e.opts())
	if err != nil {
		return r, err
	}
	if e.port == nil {
		e.port = map[[2]string]analysis.PortResult{}
	}
	e.port[k] = r
	return r, nil
}

// relErr is |got − want| / max(|want|, floor).
func relErr(got, want, floor float64) float64 {
	return math.Abs(got-want) / math.Max(math.Abs(want), floor)
}

// Cutoff scores the −3 dB frequency of V(node)/V(ref) against target [Hz]:
// error = |f − target| / target. high selects the upper edge (low-pass,
// band-pass upper), otherwise the lower edge.
func Cutoff(node, ref string, spec analysis.ACSpec, target float64, high bool) Objective {
	return Objective{Name: "cutoff " + node, Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.AC(spec)
		if err != nil {
			return 0, err
		}
		h, err := analysis.Transfer(r, node, ref)
		if err != nil {
			return 0, err
		}
		lo, hi, err := analysis.Cutoff3dB(r.Freq, h)
		if err != nil {
			return 0, err
		}
		f := lo
		if high {
			f = hi
		}
		if f == 0 || math.IsInf(f, 0) {
			return 0, analysis.ErrNoMeasurement
		}
		return relErr(f, target, 0), nil
	}}
}

// GainAt scores |H(f)| of V(node)/V(ref) against targetDB: error =
// |gain − target| in dB / 20 (one decade of gain counts 1).
func GainAt(node, ref string, spec analysis.ACSpec, f, targetDB float64) Objective {
	return Objective{Name: fmt.Sprintf("gain %s at %g Hz", node, f), Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.AC(spec)
		if err != nil {
			return 0, err
		}
		h, err := analysis.Transfer(r, node, ref)
		if err != nil {
			return 0, err
		}
		g, err := analysis.GainAt(r.Freq, h, f)
		if err != nil || g <= 0 {
			return 0, analysis.ErrNoMeasurement
		}
		return math.Abs(20*math.Log10(g)-targetDB) / 20, nil
	}}
}

// GainShortfall is a "more is better" gain objective for multi-objective
// runs: error = max(0, ceilingDB − gain at f) / 20, so it falls as the gain
// of V(node)/V(ref) rises and reaches 0 at the ceiling (one decade of gain
// counts 1).
func GainShortfall(node, ref string, spec analysis.ACSpec, f, ceilingDB float64) Objective {
	return Objective{Name: fmt.Sprintf("gain short %s at %g Hz", node, f), Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.AC(spec)
		if err != nil {
			return 0, err
		}
		h, err := analysis.Transfer(r, node, ref)
		if err != nil {
			return 0, err
		}
		g, err := analysis.GainAt(r.Freq, h, f)
		if err != nil || g <= 0 {
			return 0, analysis.ErrNoMeasurement
		}
		return math.Max(0, ceilingDB-20*math.Log10(g)) / 20, nil
	}}
}

// BandwidthShortfall is a "more is better" bandwidth objective: error =
// max(0, log10(ceilingHz / f₋₃dB)) in decades, with f₋₃dB the upper −3 dB
// edge of V(node)/V(ref) (the sweep's stop frequency when the response has
// not fallen by then).
func BandwidthShortfall(node, ref string, spec analysis.ACSpec, ceilingHz float64) Objective {
	return Objective{Name: "bandwidth short " + node, Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.AC(spec)
		if err != nil {
			return 0, err
		}
		h, err := analysis.Transfer(r, node, ref)
		if err != nil {
			return 0, err
		}
		_, hi, err := analysis.Cutoff3dB(r.Freq, h)
		if err != nil {
			return 0, err
		}
		hi = math.Min(hi, r.Freq[len(r.Freq)-1])
		if !(hi > 0) {
			return 0, analysis.ErrNoMeasurement
		}
		return math.Max(0, math.Log10(ceilingHz/hi)), nil
	}}
}

// SupplyPower is the power the named sources deliver at the operating
// point [W] (the sum of −V·I; a source that absorbs counts negative), as a
// "less is better" objective. A negative total scores 0.
func SupplyPower(sources ...string) Objective {
	return Objective{Name: "supply power", Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		p := 0.0
		for _, s := range sources {
			el, ok := r.Element(s)
			if !ok {
				return 0, fmt.Errorf("%w: no element %q", analysis.ErrNoMeasurement, s)
			}
			p -= el.PowerAbsorbed
		}
		return math.Max(0, p), nil
	}}
}

// OPEfficiency scores 1 − η at the operating point, η = Σ P(loads) /
// inflow from the ledger (em.OPEfficiency): "less is better", 0 for a
// lossless transfer. A negative η (loads delivering power) scores above 1.
func OPEfficiency(loads ...string) Objective {
	return Objective{Name: "1 − efficiency", Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		eta, err := em.OPEfficiency(r, loads...)
		if err != nil {
			return 0, fmt.Errorf("%w: %v", analysis.ErrNoMeasurement, err)
		}
		return 1 - eta, nil
	}}
}

// TranEfficiency scores 1 − η over [from, to] of a transient, η the energy
// absorbed by the loads over the energy that entered from the exterior
// (em.TranEfficiency); use whole periods of a steady state.
func TranEfficiency(spec analysis.TranSpec, from, to float64, loads ...string) Objective {
	return Objective{Name: "1 − efficiency", Weight: 1, Eval: func(e *Eval) (float64, error) {
		wf, err := e.Tran(spec)
		if err != nil {
			return 0, err
		}
		eta, _, _, err := em.TranEfficiency(wf, from, to, loads...)
		if err != nil {
			return 0, fmt.Errorf("%w: %v", analysis.ErrNoMeasurement, err)
		}
		return 1 - eta, nil
	}}
}

// Band is one band of a frequency-response mask: between From and To [Hz]
// the gain |V(node)/V(ref)| must lie within [MinDB, MaxDB]. Use
// math.Inf(-1) / math.Inf(1) for an open side.
type Band struct {
	From, To     float64
	MinDB, MaxDB float64
}

// MaskViolationFloor is the least error of a response that leaves its
// mask anywhere: a design that misses the spec by a hair always scores
// worse than one that meets it plus any parsimony penalty below the
// floor, so "fitness < MaskViolationFloor" means "spec met".
const MaskViolationFloor = 0.1

// FilterMask scores V(node)/V(ref) against a mask of bands, checked at
// every sweep point inside a band and at the band edges (interpolated).
// The error is 0 when the response stays inside the mask everywhere, and
// otherwise MaskViolationFloor + the fraction of points outside + the mean
// excess in dB/20 over all checked points (Koza et al. 1999 score filters by
// the sum of weighted deviations over frequency; this keeps that grading
// and adds the floor so a spec-meeting design is recognizable).
func FilterMask(node, ref string, spec analysis.ACSpec, bands []Band) Objective {
	return Objective{Name: "mask " + node, Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.AC(spec)
		if err != nil {
			return 0, err
		}
		h, err := analysis.Transfer(r, node, ref)
		if err != nil {
			return 0, err
		}
		n, bad, excess := 0, 0, 0.0
		check := func(f float64, b Band) error {
			// A band edge at the sweep's end may lie a rounding error
			// outside the generated frequencies.
			if lo, hi := r.Freq[0], r.Freq[len(r.Freq)-1]; f < lo && f > lo*(1-1e-9) {
				f = lo
			} else if f > hi && f < hi*(1+1e-9) {
				f = hi
			}
			g, err := analysis.GainAt(r.Freq, h, f)
			if err != nil {
				return err
			}
			db := math.Inf(-1)
			if g > 0 {
				db = 20 * math.Log10(g)
			}
			n++
			d := 0.0
			switch {
			case db > b.MaxDB:
				d = db - b.MaxDB
			case db < b.MinDB:
				d = math.Min(b.MinDB-db, 200) // a dead output counts 200 dB
			}
			if d > 0 {
				bad++
				excess += d
			}
			return nil
		}
		for _, b := range bands {
			if err := check(b.From, b); err != nil {
				return 0, err
			}
			for _, f := range r.Freq {
				if f > b.From && f < b.To {
					if err := check(f, b); err != nil {
						return 0, err
					}
				}
			}
			if err := check(b.To, b); err != nil {
				return 0, err
			}
		}
		if bad == 0 {
			return 0, nil
		}
		return MaskViolationFloor + float64(bad)/float64(n) + excess/20/float64(n), nil
	}}
}

// PartCount is the parsimony term of a topology genome: the number of
// evolved parts (weight it with the per-part penalty).
func PartCount() Objective {
	return Objective{Name: "parts", Weight: 1, Eval: func(e *Eval) (float64, error) {
		if e.Parts < 0 {
			return 0, fmt.Errorf("%w: parts objective needs a topology genome", analysis.ErrNoMeasurement)
		}
		return float64(e.Parts), nil
	}}
}

// OPVoltage scores a node's operating-point voltage against target [V].
func OPVoltage(node string, target float64) Objective {
	return Objective{Name: "v(" + node + ")", Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		v, ok := r.Voltage(node)
		if !ok {
			return 0, analysis.ErrNoMeasurement
		}
		return relErr(v, target, 1e-3), nil
	}}
}

// window returns the samples of y with from ≤ t ≤ to.
func window(t, y []float64, from, to float64) []float64 {
	var out []float64
	for k, tt := range t {
		if tt >= from && tt <= to {
			out = append(out, y[k])
		}
	}
	return out
}

// TranAverage scores the mean of V(node) over [from, to] against target.
// The mean is the time average (trapezoidal over the simulator's points).
func TranAverage(node string, spec analysis.TranSpec, from, to, target float64) Objective {
	return Objective{Name: "mean v(" + node + ")", Weight: 1, Eval: func(e *Eval) (float64, error) {
		wf, err := e.Tran(spec)
		if err != nil {
			return 0, err
		}
		y, ok := wf.NodeVoltage(node)
		if !ok {
			return 0, analysis.ErrNoMeasurement
		}
		area, span := 0.0, 0.0
		for k := 1; k < len(wf.Time); k++ {
			t0, t1 := wf.Time[k-1], wf.Time[k]
			if t1 <= from || t0 >= to {
				continue
			}
			a, b := math.Max(t0, from), math.Min(t1, to)
			ya := output.Interp(wf.Time, y, a)
			yb := output.Interp(wf.Time, y, b)
			area += 0.5 * (ya + yb) * (b - a)
			span += b - a
		}
		if span <= 0 {
			return 0, analysis.ErrNoMeasurement
		}
		return relErr(area/span, target, 1e-3), nil
	}}
}

// TranRipple is a constraint on the peak-to-peak value of V(node) over
// [from, to]: error = max(0, ripple − max) / max.
func TranRipple(node string, spec analysis.TranSpec, from, to, max float64) Objective {
	return Objective{Name: "ripple v(" + node + ")", Weight: 1, Eval: func(e *Eval) (float64, error) {
		wf, err := e.Tran(spec)
		if err != nil {
			return 0, err
		}
		y, ok := wf.NodeVoltage(node)
		if !ok {
			return 0, analysis.ErrNoMeasurement
		}
		w := window(wf.Time, y, from, to)
		if len(w) < 2 || max <= 0 {
			return 0, analysis.ErrNoMeasurement
		}
		lo, hi := w[0], w[0]
		for _, v := range w {
			lo, hi = math.Min(lo, v), math.Max(hi, v)
		}
		return math.Max(0, hi-lo-max) / max, nil
	}}
}

// PortVoc scores the open-circuit voltage between plus and minus.
func PortVoc(plus, minus string, target float64) Objective {
	return Objective{Name: "Voc " + plus + "," + minus, Weight: 1, Eval: func(e *Eval) (float64, error) {
		p, err := e.Port(plus, minus)
		if err != nil {
			return 0, err
		}
		return relErr(p.Voc, target, 1e-6), nil
	}}
}

// PortIsc scores the short-circuit current between plus and minus.
func PortIsc(plus, minus string, target float64) Objective {
	return Objective{Name: "Isc " + plus + "," + minus, Weight: 1, Eval: func(e *Eval) (float64, error) {
		p, err := e.Port(plus, minus)
		if err != nil {
			return 0, err
		}
		return relErr(p.Isc, target, 1e-12), nil
	}}
}

// ModelSuppliedLimit is a constraint on power from idealized active
// elements at the operating point: error = max(0, P − max) / max(max, 1e-12).
// With max = 0 any reliance on negative resistance or controlled-source
// power is penalized in proportion to the power.
func ModelSuppliedLimit(max float64) Objective {
	return Objective{Name: "model-supplied power", Weight: 1, Eval: func(e *Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		return math.Max(0, r.Ledger.ModelSupplied-max) / math.Max(max, 1e-12), nil
	}}
}

// Failure classes and their penalty fitness. A failed candidate always
// scores worse than any successful one with total error below 1e5, and the
// classes are ordered so evolution can still tell "almost works" (a
// measurement is undefined) from "broken" (singular matrix).
const (
	PenaltyNoMeasurement = 1e6
	PenaltyNoConvergence = 2e6
	PenaltySingular      = 3e6
	PenaltyInvalid       = 4e6
	PenaltyLedger        = 5e6
	PenaltyTimeout       = 6e6
	PenaltyOther         = 1e7
	PenaltyPanic         = 2e7
)

// Penalty maps an evaluation error to its penalty fitness and a short class
// name (CLAUDE.md §4). A ledger imbalance is penalized heavily: a candidate
// whose energy does not add up is a numerical artefact, never a discovery.
func Penalty(err error) (float64, string) {
	switch {
	case errors.Is(err, context.DeadlineExceeded), errors.Is(err, context.Canceled):
		return PenaltyTimeout, "timeout"
	case errors.Is(err, em.ErrLedgerResidual):
		return PenaltyLedger, "ledger"
	case errors.Is(err, analysis.ErrNoMeasurement):
		return PenaltyNoMeasurement, "no-measurement"
	case errors.Is(err, analysis.ErrNoConvergence), errors.Is(err, analysis.ErrTimestepTooSmall),
		errors.Is(err, analysis.ErrTooManySteps):
		return PenaltyNoConvergence, "no-convergence"
	case errors.Is(err, solver.ErrSingularMatrix), errors.Is(err, solver.ErrNonFinite), errors.Is(err, kron.ErrSingular):
		return PenaltySingular, "singular"
	case errors.Is(err, circuit.ErrFloatingNode), errors.Is(err, circuit.ErrBadValue), errors.Is(err, em.ErrOutsideModelValidity),
		errors.Is(err, em.ErrNotQuasiStatic),
		errors.Is(err, circuit.ErrBadModel), errors.Is(err, circuit.ErrSelfConnected):
		return PenaltyInvalid, "invalid"
	}
	return PenaltyOther, "error"
}
