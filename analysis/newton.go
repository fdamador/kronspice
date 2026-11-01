// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/device"

	"github.com/fdamador/kronspice/solver"
)

// ErrNoConvergence means Newton–Raphson did not converge within its
// iteration limit, including after every DC fallback (CLAUDE.md §3,
// "Nonlinear solve"). It is never an infinite loop.
var ErrNoConvergence = errors.New("analysis: no convergence")

// Bounds on the DC fallback loops, so a hopeless circuit always returns.
const (
	maxGminSteps   = 200
	maxSourceSteps = 1000
	maxPTSteps     = 2000
)

// newton iterates linearize → solve until the iterate converges or maxIter
// is reached. A circuit with no nonlinear device needs exactly one solve.
//
// Convergence follows SPICE3 (Nagel 1975, ch. 5; SPICE3f5 NIiter/NIconvTest):
// no device limited its junction voltages in the last load, every unknown
// changed by at most reltol·max(|x|, |x_old|) + vntol (voltages) or + abstol
// (currents), and every nonlinear device passes its own current check.
func (s *Sim) newton(maxIter int, o *Options) error {
	s.lc.Gmin = o.Gmin
	for it := 0; it < maxIter; it++ {
		if err := o.interrupted(); err != nil {
			return err
		}
		if err := s.solve(); err != nil {
			return err
		}
		if len(s.nonlin) == 0 {
			return nil
		}
		if s.lc.InitJct {
			s.lc.InitJct = false
			continue
		}
		if s.converged(o) {
			return nil
		}
	}
	return fmt.Errorf("%w after %d iterations", ErrNoConvergence, maxIter)
}

func (s *Sim) converged(o *Options) bool {
	for _, nl := range s.nonlin {
		if nl.Limited() {
			return false
		}
	}
	for i, x := range s.x {
		old := s.xOld[i]
		tol := o.AbsTol
		if s.isV[i] {
			tol = o.VnTol
		}
		tol += o.RelTol * math.Max(math.Abs(x), math.Abs(old))
		if math.Abs(x-old) > tol {
			return false
		}
	}
	for _, nl := range s.nonlin {
		if !nl.ConvTest(s.x, o.RelTol, o.AbsTol) {
			return false
		}
	}
	return true
}

// retryable reports whether a failed Newton solve may succeed with a
// fallback (as opposed to a structural error such as a bad matrix size).
func retryable(err error) bool {
	return errors.Is(err, ErrNoConvergence) || errors.Is(err, solver.ErrNonFinite) ||
		errors.Is(err, solver.ErrSingularMatrix)
}

// dcSolve finds a DC solution (operating point or transient starting point)
// with the fallback chain of CLAUDE.md §3: plain Newton, then gmin
// stepping, then source stepping, then pseudo-transient continuation
// (Kundert, The Designer's Guide to SPICE and Spectre, ch. 4).
func (s *Sim) dcSolve(o *Options) error {
	s.lc.InitJct = !s.warm
	if !s.warm && len(o.NodeSet) > 0 {
		if err := s.applyNodeSet(o.NodeSet); err != nil {
			return err
		}
		s.lc.InitJct = false
	}
	err := s.newton(o.ITL1, o)
	if err == nil || len(s.nonlin) == 0 || !retryable(err) {
		s.warm = s.warm || err == nil
		return err
	}
	first := err
	if !o.NoGminStepping {
		if err = s.gminStepping(o); err == nil {
			s.warm = true
			return nil
		}
	}
	if !o.NoSourceStepping {
		if err = s.sourceStepping(o); err == nil {
			s.warm = true
			return nil
		}
	}
	if !o.NoPseudoTransient {
		if err = s.pseudoTransient(o); err == nil {
			s.warm = true
			return nil
		}
	}
	return fmt.Errorf("%w (fallbacks failed: %v)", first, err)
}

// gminStepping solves with a large conductance from every node to ground and
// reduces it step by step to zero, each solve starting from the previous one
// (SPICE3f5 dynamic gmin stepping). A failed step is retried with a smaller
// reduction factor.
func (s *Sim) gminStepping(o *Options) error {
	defer func() { s.diagGmin = 0 }()
	save := make([]float64, len(s.x))
	for i := range s.x {
		s.x[i] = 0
	}
	s.lc.InitJct = true
	g, factor := 1e-2, 10.0
	ok := false // a step at some g has converged
	for n := 0; n < maxGminSteps; n++ {
		s.diagGmin = g
		if err := s.newton(o.ITL1, o); err == nil {
			copy(save, s.x)
			if g == 0 {
				return nil
			}
			ok = true
			factor = math.Min(factor*math.Sqrt(factor), 10)
			if g/factor < o.Gmin {
				g = 0
			} else {
				g /= factor
			}
			continue
		} else if !retryable(err) {
			return err
		}
		if !ok {
			return fmt.Errorf("%w: gmin stepping failed at its first step", ErrNoConvergence)
		}
		copy(s.x, save)
		prev := g * factor
		if g == 0 {
			prev = o.Gmin
		}
		factor = math.Sqrt(factor)
		if factor < 1.00005 {
			return fmt.Errorf("%w: gmin stepping stalled at %g S", ErrNoConvergence, prev)
		}
		g = prev / factor
	}
	return fmt.Errorf("%w: gmin stepping exceeded %d steps", ErrNoConvergence, maxGminSteps)
}

// sourceStepping ramps every independent source from 0 to its full value,
// each solve starting from the previous one; failed steps are halved.
func (s *Sim) sourceStepping(o *Options) error {
	defer func() { s.lc.SourceStepping, s.lc.SrcFact = false, 0 }()
	save := make([]float64, len(s.x))
	for i := range s.x {
		s.x[i] = 0
	}
	s.lc.SourceStepping = true
	done, step := 0.0, 0.1
	s.lc.SrcFact = 0
	if err := s.newton(o.ITL1, o); err != nil {
		return err
	}
	copy(save, s.x)
	for n := 0; n < maxSourceSteps; n++ {
		f := math.Min(done+step, 1)
		s.lc.SrcFact = f
		err := s.newton(o.ITL1, o)
		if err == nil {
			if f == 1 {
				return nil
			}
			copy(save, s.x)
			done, step = f, math.Min(step*1.5, 0.5)
			continue
		} else if !retryable(err) {
			return err
		}
		copy(s.x, save)
		step /= 2
		if step < 1e-3 {
			return fmt.Errorf("%w: source stepping stalled at %.4g of full value", ErrNoConvergence, done)
		}
	}
	return fmt.Errorf("%w: source stepping exceeded %d steps", ErrNoConvergence, maxSourceSteps)
}

// ErrUnknownNode is returned for a node name that is not in the circuit.
var ErrUnknownNode = errors.New("analysis: unknown node")

// applyNodeSet starts Newton from the given node voltages, and every device
// from its control value at that point.
func (s *Sim) applyNodeSet(ns map[string]float64) error {
	for i := range s.x {
		s.x[i] = 0
	}
	for name, v := range ns {
		k := -1
		for i, n := range s.c.NodeNames {
			if i > 0 && strings.EqualFold(n, name) {
				k = i
			}
		}
		if k < 0 || math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: nodeset %s = %g", ErrUnknownNode, name, v)
		}
		s.x[k-1] = v
	}
	for _, d := range s.devs {
		if _, ok := d.(device.Nonlinear); ok {
			d.Accept(s.x, 0, device.ModeOP) // seed limiting from the guess
		}
	}
	return nil
}

// pseudoTransient finds a DC solution by integrating a fictitious dynamic
// circuit to its steady state: a unit capacitance from every node to ground
// and a unit inductance in series with every branch-current unknown, with
// backward Euler and a growing step (pseudo-transient continuation; Kelley &
// Keyes, SIAM J. Numer. Anal. 35, 1998; Kundert ch. 4). At steady state the
// fictitious elements carry nothing, so the result solves the real circuit,
// and a final plain Newton solve confirms it.
//
// Unlike the stepping methods it can follow a nonlinear element past a
// turning point of its characteristic. An S-type element (glow lamp, arc)
// driven beyond its striking voltage has no solution on its low-current
// branch, and Newton cannot cross the negative-slope part of the curve
// towards the solution on the other side; the series pseudo-inductance lets
// the current run through that part, as the physical discharge does when it
// strikes.
func (s *Sim) pseudoTransient(o *Options) error {
	defer func() { s.ptG = 0 }()
	if s.ptPrev == nil {
		s.ptPrev = make([]float64, len(s.x))
	}
	for i := range s.x {
		s.x[i] = 0
	}
	copy(s.ptPrev, s.x)
	s.lc.InitJct = true
	h := 1e-6
	quiet := 0
	for n := 0; n < maxPTSteps; n++ {
		s.ptG = 1 / h
		err := s.newton(o.ITL1, o)
		if err != nil {
			if !retryable(err) {
				return err
			}
			copy(s.x, s.ptPrev)
			h /= 4
			if h < 1e-15 {
				return fmt.Errorf("%w: pseudo-transient step too small", ErrNoConvergence)
			}
			continue
		}
		// Steady state: the state stopped moving at a large step.
		moved := false
		for i, x := range s.x {
			tol := o.AbsTol
			if s.isV[i] {
				tol = o.VnTol
			}
			if math.Abs(x-s.ptPrev[i]) > o.RelTol*math.Abs(x)+tol {
				moved = true
				break
			}
		}
		copy(s.ptPrev, s.x)
		if !moved {
			quiet++
		} else {
			quiet = 0
		}
		if quiet >= 3 && h >= 1e6 {
			s.ptG = 0
			return s.newton(o.ITL1, o)
		}
		h = math.Min(h*2, 1e12)
	}
	return fmt.Errorf("%w: pseudo-transient exceeded %d steps", ErrNoConvergence, maxPTSteps)
}
