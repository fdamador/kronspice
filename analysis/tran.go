// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/output"
)

// Errors from transient analysis (sentinels for GA penalties).
var (
	ErrTimestepTooSmall = errors.New("analysis: timestep too small")
	ErrTooManySteps     = errors.New("analysis: too many time steps")
	ErrBadTranSpec      = errors.New("analysis: invalid transient specification")
)

// TranSpec is a .tran request: print step, stop time, start of recorded
// output and maximum internal step (0 = SPICE default), and UIC.
type TranSpec struct {
	Step, Stop, Start, MaxStep float64
	UIC                        bool
}

// TranOptions adds time-stepping controls to Options. Zero values take the
// SPICE defaults noted on each field.
type TranOptions struct {
	Options
	Method   device.Method // Trapezoidal (default) or BackwardEuler
	ChgTol   float64       // default 1e-14 C
	TrTol    float64       // default 7
	MinStep  float64       // default 1e-12 × (stop time)
	MaxSteps int           // default 1,000,000 accepted steps
}

func (o *TranOptions) defaults(spec TranSpec) {
	if o.Method == 0 {
		o.Method = device.Trapezoidal
	}
	o.Options.defaults()
	if o.ChgTol <= 0 {
		o.ChgTol = 1e-14
	}
	if o.TrTol <= 0 {
		o.TrTol = 7
	}
	if o.MinStep <= 0 {
		o.MinStep = 1e-12 * spec.Stop
	}
	if o.MaxSteps <= 0 {
		o.MaxSteps = 1_000_000
	}
}

// Tran runs a transient analysis of c.
//
// Time stepping follows SPICE practice (Nagel 1975, ch. 6): the operating
// point at t = 0 starts the run (unless UIC), the first step and every step
// after a source breakpoint use backward Euler, later steps use the selected
// method (trapezoidal by default), each candidate step is checked against the
// local truncation error of every capacitor and inductor, and the engine
// lands exactly on every source breakpoint. Failures return
// ErrTimestepTooSmall or ErrTooManySteps instead of looping forever.
func Tran(c *circuit.Circuit, spec TranSpec, opt TranOptions) (*output.Waveforms, error) {
	if spec.Step <= 0 || spec.Stop <= 0 || spec.Start < 0 || spec.Start >= spec.Stop || spec.MaxStep < 0 ||
		math.IsNaN(spec.Step+spec.Stop+spec.Start+spec.MaxStep) || math.IsInf(spec.Stop, 0) {
		return nil, fmt.Errorf("%w: %+v", ErrBadTranSpec, spec)
	}
	var err error
	if spec.UIC {
		err = c.ValidateTransient()
	} else {
		err = c.Validate()
	}
	if err != nil {
		return nil, err
	}
	opt.defaults(spec)
	s, err := newSim(c, spec.Step, spec.Stop)
	if err != nil {
		return nil, err
	}
	if err := s.useSolver(&opt.Options); err != nil {
		return nil, err
	}
	return s.tran(spec, opt)
}

func (s *Sim) tran(spec TranSpec, opt TranOptions) (*output.Waveforms, error) {
	maxStep := spec.MaxStep
	if maxStep == 0 {
		// SPICE default: the smaller of the print step and (stop − start)/50.
		maxStep = math.Min(spec.Step, (spec.Stop-spec.Start)/50)
	}
	rec := newRecorder(s, spec.Start)

	// Starting point.
	if spec.UIC {
		for _, d := range s.devs {
			if ic, ok := d.(device.InitialConditioner); ok {
				ic.SetInitialCondition()
			}
		}
		s.lc.Mode, s.lc.Time = device.ModeUICInit, 0
		s.lc.InitJct = true
		if err := s.newton(opt.ITL1, &opt.Options); err != nil {
			return nil, fmt.Errorf("initial conditions: %w", err)
		}
		s.accept(0, device.ModeUICInit)
		rec.record(0, true) // imposed point: ledger not checked
	} else {
		s.lc.Mode, s.lc.Time = device.ModeTranOP, 0
		if err := s.dcSolve(&opt.Options); err != nil {
			return nil, fmt.Errorf("transient operating point: %w", err)
		}
		s.accept(0, device.ModeTranOP)
		if err := rec.record(0, false); err != nil && !opt.SkipEnergyCheck {
			return rec.wf, err
		}
	}

	tc := device.TruncContext{RelTol: opt.RelTol, AbsTol: opt.AbsTol, VnTol: opt.VnTol, ChgTol: opt.ChgTol, TrTol: opt.TrTol}
	t := 0.0
	h := math.Min(maxStep, math.Min(spec.Step, spec.Stop/100)) / 10
	order := device.BackwardEuler
	stats := &rec.wf.Stats
	stats.MinStep = math.Inf(1)
	eps := 1e-12 * spec.Stop
	xAcc := append([]float64(nil), s.x...) // last accepted solution
	reject := func(hNew float64) error {
		stats.Rejected++
		copy(s.x, xAcc)
		if hNew < opt.MinStep {
			return fmt.Errorf("%w: %.3g s at t = %g s", ErrTimestepTooSmall, hNew, t)
		}
		return nil
	}

	for t < spec.Stop-eps {
		if err := opt.interrupted(); err != nil {
			return rec.wf, err
		}
		if stats.Accepted >= opt.MaxSteps {
			return rec.wf, fmt.Errorf("%w: %d accepted steps at t = %g s", ErrTooManySteps, stats.Accepted, t)
		}
		nb := s.nextBreakpoint(t)
		tn := t + h
		hitBreak := false
		if tn >= nb-eps {
			tn, hitBreak = nb, true
		}
		if tn >= spec.Stop-eps {
			tn, hitBreak = spec.Stop, false
		}
		h = tn - t

		s.lc.Mode, s.lc.Time, s.lc.H, s.lc.Method = device.ModeTran, tn, h, order
		if err := s.newton(opt.ITL4, &opt.Options); err != nil {
			if !retryable(err) || len(s.nonlin) == 0 {
				return rec.wf, fmt.Errorf("t = %g s: %w", tn, err)
			}
			// Newton failed within itl4: retry with an eighth of the step
			// (SPICE3 DCtran).
			h /= 8
			if err := reject(h); err != nil {
				return rec.wf, fmt.Errorf("%w (Newton did not converge: %v)", err, ErrNoConvergence)
			}
			order = device.BackwardEuler
			continue
		}

		// Local truncation error control.
		tc.Time, tc.H, tc.Method = tn, h, order
		proposal := math.Inf(1)
		for _, i := range s.reactive {
			if p := s.devs[i].(device.Reactive).Truncate(s.x, &tc); p < proposal {
				proposal = p
			}
		}
		if proposal < 0.9*h {
			h = math.Max(proposal, h/8)
			if err := reject(h); err != nil {
				return rec.wf, err
			}
			continue
		}

		// Accept.
		s.accept(tn, device.ModeTran)
		if err := rec.record(tn, false); err != nil && !opt.SkipEnergyCheck {
			return rec.wf, err
		}
		copy(xAcc, s.x)
		stats.Accepted++
		stats.MinStep = math.Min(stats.MinStep, h)
		stats.MaxStep = math.Max(stats.MaxStep, h)
		t = tn

		// Next step.
		next := math.Min(proposal, 2*h)
		next = math.Min(next, maxStep)
		order = opt.Method
		if hitBreak {
			// Restart with backward Euler and a small step after a corner.
			order = device.BackwardEuler
			nb2 := s.nextBreakpoint(t)
			next = math.Min(next, 0.1*math.Min(maxStep, nb2-t))
		}
		h = math.Max(next, opt.MinStep)
	}
	if err := rec.finish(); err != nil && !opt.SkipEnergyCheck {
		return rec.wf, err
	}
	return rec.wf, nil
}

func (s *Sim) nextBreakpoint(t float64) float64 {
	nb := math.Inf(1)
	for _, d := range s.devs {
		if bp, ok := d.(device.Breakpointer); ok {
			if b := bp.NextBreakpoint(t); b < nb {
				nb = b
			}
		}
	}
	return nb
}

// recorder stores time points and the per-point energy ledger.
type recorder struct {
	s     *Sim
	start float64
	wf    *output.Waveforms
	tol   float64
	first bool
	// imposedFirst is true when the first recorded point is a UIC point:
	// its element powers belong to an imposed, not a solved, state (an
	// inductor's initial current can sit in a reference resistor until the
	// first step), so they are not integrated (finish).
	imposedFirst bool
	// peakScale is the largest Σ|power| seen so far. Round-off in the solved
	// powers scales with the largest quantities in the run, not with the
	// powers at one instant: when every power is momentarily near zero
	// (a ringing circuit at a zero crossing, or a sine drive just after
	// t = 0) the cancelling terms inside each companion model are still of
	// the run's size. Each point's residual is therefore judged against the
	// whole run's peak, after the run (finish); during the run only gross
	// imbalances (grossTol of the peak so far) stop it early.
	peakScale float64
	worstAbs  float64 // largest |residual| of a checked point [W]
	worstT    float64
}

// grossTol stops a transient early when one point is out of balance by
// this fraction of the peak so far: far beyond round-off, a model bug.
const grossTol = 1e-3

func newRecorder(s *Sim, start float64) *recorder {
	wf := &output.Waveforms{
		Title:     s.c.Title,
		NodeNames: append([]string(nil), s.c.NodeNames...),
		V:         make([][]float64, s.c.NumNodes()),
		I:         make([][]float64, len(s.devs)),
		P:         make([][]float64, len(s.devs)),
	}
	for _, d := range s.devs {
		wf.ElementNames = append(wf.ElementNames, d.Name())
		wf.ElementKinds = append(wf.ElementKinds, d.Kind().String())
	}
	return &recorder{s: s, start: start, wf: wf, tol: em.DefaultTolerance, first: true}
}

// record appends the current solution at time t (if t ≥ start) and checks the
// point's energy balance unless the point was imposed (UIC).
func (r *recorder) record(t float64, imposed bool) error {
	s := r.s
	var acc em.Accumulator
	stored := 0.0
	// Element values.
	ops := make([]device.OperatingPoint, len(s.devs))
	for i, d := range s.devs {
		ops[i] = d.Operating(s.x)
	}
	for _, p := range s.ports {
		acc.AddPort(p.ch, p.status, s.devs[p.dev].ExteriorPower(p.ch, ops[p.dev]))
	}
	for _, i := range s.supply {
		acc.AddModelSupplied(s.devs[i].(device.ModelSupplier).ModelSuppliedPower(ops[i]))
	}
	for _, i := range s.reactive {
		acc.AddStorage(ops[i].StoragePower)
		stored += s.devs[i].(device.Reactive).StoredEnergy(ops[i])
	}
	var err error
	if acc.Scale > r.peakScale {
		r.peakScale = acc.Scale
	}
	res := math.Abs(acc.Residual())
	if !imposed && res > math.Max(output.LedgerAbsFloor, s.roundoff()) {
		if res > r.worstAbs {
			r.worstAbs, r.worstT = res, t
		}
		if rr := res / r.peakScale; rr > grossTol {
			err = fmt.Errorf("%w at t = %g s: residual %.3g W (relative %.3g)", em.ErrLedgerResidual, t, acc.Residual(), rr)
		}
	}
	if t < r.start {
		return err
	}
	wf := r.wf
	if len(wf.Time) == 0 {
		r.imposedFirst = imposed
	}
	wf.Time = append(wf.Time, t)
	for n := 1; n < s.c.NumNodes(); n++ {
		wf.V[n] = append(wf.V[n], s.x[n-1])
	}
	for i := range s.devs {
		wf.I[i] = append(wf.I[i], ops[i].Current)
		wf.P[i] = append(wf.P[i], ops[i].PowerAbsorbed)
	}
	l := &wf.Ledger
	l.Inflow = append(l.Inflow, acc.Inflow)
	l.Outflow = append(l.Outflow, acc.Outflow)
	l.Dissipation = append(l.Dissipation, acc.Dissipation)
	l.StorageRate = append(l.StorageRate, acc.StorageRate)
	l.ModelSupplied = append(l.ModelSupplied, acc.ModelSupplied)
	l.Residual = append(l.Residual, acc.Residual())
	l.StoredEnergy = append(l.StoredEnergy, stored)
	return err
}

// finish integrates energies over the recorded points and checks every
// point's residual against the run's peak Σ|power|.
func (r *recorder) finish() error {
	l := &r.wf.Ledger
	var err error
	if r.peakScale > 0 {
		l.MaxRelativeResidual = r.worstAbs / r.peakScale
		if l.MaxRelativeResidual > r.tol {
			err = fmt.Errorf("%w at t = %g s: residual %.3g W (relative to the run's peak power %.3g)",
				em.ErrLedgerResidual, r.worstT, r.worstAbs, l.MaxRelativeResidual)
		}
	}
	t := r.wf.Time
	n := len(t)
	if n == 0 {
		return err
	}
	trap := func(y []float64) float64 {
		sum := 0.0
		for k := 1; k < n; k++ {
			if k == 1 && r.imposedFirst {
				// The first step is backward Euler, whose energy change is
				// booked at its right end; the imposed UIC point's powers
				// are not a solved state and are left out.
				sum += y[1] * (t[1] - t[0])
				continue
			}
			sum += 0.5 * (y[k] + y[k-1]) * (t[k] - t[k-1])
		}
		return sum
	}
	l.EnergyIn = trap(l.Inflow)
	l.EnergyOut = trap(l.Outflow)
	l.EnergyDissipated = trap(l.Dissipation)
	l.EnergyModelSupplied = trap(l.ModelSupplied)
	l.StoredStart, l.StoredEnd = l.StoredEnergy[0], l.StoredEnergy[n-1]
	l.IntegrationMismatch = l.EnergyIn + l.EnergyModelSupplied - l.EnergyOut - l.EnergyDissipated - (l.StoredEnd - l.StoredStart)
	return err
}
