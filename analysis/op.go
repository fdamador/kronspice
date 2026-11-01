// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package analysis runs circuit analyses: the DC operating point (.op) and
// transient analysis (.tran).
package analysis

import (
	"context"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

// Options controls an analysis.
type Options struct {
	// SkipEnergyCheck disables the ledger residual check. The check is on by
	// default and always on in tests (CLAUDE.md §3c, rule 5).
	SkipEnergyCheck bool
	// LedgerTolerance overrides em.DefaultTolerance when > 0.
	LedgerTolerance float64

	// Newton–Raphson controls; zero values take the SPICE defaults.
	RelTol float64 // relative tolerance, default 1e-3
	AbsTol float64 // current tolerance, default 1e-12 A
	VnTol  float64 // voltage tolerance, default 1e-6 V
	Gmin   float64 // junction minimum conductance, default 1e-12 S
	ITL1   int     // DC iteration limit, default 100
	ITL4   int     // iteration limit per transient step, default 10
	// NoGminStepping and NoSourceStepping disable the DC fallbacks.
	NoGminStepping, NoSourceStepping, NoPseudoTransient bool

	// NodeSet is an initial guess for node voltages (SPICE .nodeset) used by
	// a cold DC solve; nodes not listed start at 0 V. With nonlinear
	// elements that have negative-slope regions a circuit can have several
	// operating points, and the guess decides which one Newton finds.
	NodeSet map[string]float64

	// Backend selects the linear solver: solver.Auto (default) uses the
	// dense LU below solver.DenseThreshold unknowns and the sparse LU
	// above; DenseLU or SparseLU force one. Results agree to round-off.
	Backend solver.Backend

	// Solver, when set, supplies the linear systems instead of Backend:
	// every real (DC, transient) and complex (AC) system of the analysis
	// is made by it. kron.JacobianTearing is one: it solves each Newton
	// Jacobian by Kron's diakoptics. Results must agree with the built-in
	// backends to round-off.
	Solver SolverFactory

	// Context, when set, bounds the analysis: it is checked at every Newton
	// iteration, time step, frequency and sweep point, and its error is
	// returned (errors.Is(err, context.DeadlineExceeded)) so one bad circuit
	// cannot stall a GA generation (CLAUDE.md §4).
	Context context.Context
}

// SolverFactory makes the linear systems of an analysis (Options.Solver).
// Each call returns a new, independent matrix of n unknowns; a factory may
// be shared by concurrent analyses.
type SolverFactory interface {
	Matrix(n int) (solver.Matrix, error)
	CMatrix(n int) (solver.CMatrix, error)
}

// cmatrix makes an AC system: the factory's, or the chosen backend's.
func (o *Options) cmatrix(n int) (solver.CMatrix, error) {
	if o.Solver != nil {
		return o.Solver.CMatrix(n)
	}
	return solver.NewCMatrix(n, o.Backend)
}

// interrupted returns the context's error, if the analysis must stop.
func (o *Options) interrupted() error {
	if o.Context == nil {
		return nil
	}
	if err := o.Context.Err(); err != nil {
		return fmt.Errorf("analysis interrupted: %w", err)
	}
	return nil
}

func (o *Options) defaults() {
	if o.RelTol <= 0 {
		o.RelTol = 1e-3
	}
	if o.AbsTol <= 0 {
		o.AbsTol = 1e-12
	}
	if o.VnTol <= 0 {
		o.VnTol = 1e-6
	}
	if o.Gmin <= 0 {
		o.Gmin = 1e-12
	}
	if o.ITL1 <= 0 {
		o.ITL1 = 100
	}
	if o.ITL4 <= 0 {
		o.ITL4 = 10
	}
}

func (o Options) ledgerTol() float64 {
	if o.LedgerTolerance > 0 {
		return o.LedgerTolerance
	}
	return em.DefaultTolerance
}

// portRef is one modeled (device, channel) pair, precomputed so the inner
// loop visits only active ports (CLAUDE.md §3.1, "zero cost when unused").
type portRef struct {
	dev    int
	ch     exchange.ChannelID
	status exchange.PortStatus
}

// Sim holds the devices and storage for repeated solves of one circuit. It is
// not safe for concurrent use; give each goroutine its own Sim.
type Sim struct {
	c        *circuit.Circuit
	devs     []device.Device
	reactive []int // indices of devices implementing device.Reactive
	nonlin   []device.Nonlinear
	supply   []int     // indices of devices implementing device.ModelSupplier
	isV      []bool    // unknown is a voltage (vntol) or a current (abstol)
	xOld     []float64 // previous Newton iterate
	diagGmin float64   // gmin-stepping conductance from every node to ground
	warm     bool      // x holds a converged DC solution (no InitJct needed)
	ptG      float64   // pseudo-transient conductance 1/h (0 when off)
	ptPrev   []float64 // pseudo-transient previous state
	ports    []portRef
	m        solver.Matrix
	backend  solver.Backend // resolved backend of m
	factory  SolverFactory  // maker of m when Options.Solver was used
	rhs      []float64
	x        []float64
	absBuf   []float64 // workspace for the round-off estimate
	lc       device.LoadContext
}

// New prepares a simulation of c for an operating point.
func New(c *circuit.Circuit) (*Sim, error) {
	if err := c.Validate(); err != nil {
		return nil, err
	}
	return newSim(c, 0, 0)
}

func newSim(c *circuit.Circuit, step, stop float64) (*Sim, error) {
	devs, err := device.FromCircuit(c)
	if err != nil {
		return nil, err
	}
	ctx := device.NewSetupContext(c.NumNodes())
	ctx.Step, ctx.Stop = step, stop
	s := &Sim{c: c, devs: devs}
	for i, d := range devs {
		if err := d.Setup(ctx); err != nil {
			return nil, fmt.Errorf("setup %s: %w", d.Name(), err)
		}
		if _, ok := d.(device.Reactive); ok {
			s.reactive = append(s.reactive, i)
		}
		if _, ok := d.(device.ModelSupplier); ok {
			s.supply = append(s.supply, i)
		}
		if nl, ok := d.(device.Nonlinear); ok {
			s.nonlin = append(s.nonlin, nl)
		}
		for _, p := range d.Ports() {
			if p.Status != exchange.NotModeled {
				s.ports = append(s.ports, portRef{dev: i, ch: p.Channel, status: p.Status})
			}
		}
	}
	for _, d := range devs {
		if l, ok := d.(device.Linker); ok {
			if err := l.Link(ctx); err != nil {
				return nil, fmt.Errorf("setup %s: %w", d.Name(), err)
			}
		}
	}
	s.backend = solver.Auto.Resolve(ctx.Size())
	m, err := solver.NewMatrix(ctx.Size(), s.backend)
	if err != nil {
		return nil, err
	}
	s.m, s.rhs, s.x = m, make([]float64, ctx.Size()), make([]float64, ctx.Size())
	s.xOld = make([]float64, ctx.Size())
	s.absBuf = make([]float64, ctx.Size())
	s.isV = ctx.VoltageUnknowns()
	s.lc = device.LoadContext{M: m, RHS: s.rhs, X: s.xOld, Gmin: 1e-12}
	return s, nil
}

// useBackend switches the system matrix to backend b (Auto resolves by
// size). Switching keeps no state but the matrix itself, which every solve
// restamps.
func (s *Sim) useBackend(b solver.Backend) error {
	b = b.Resolve(len(s.rhs))
	if b == s.backend && s.factory == nil {
		return nil
	}
	m, err := solver.NewMatrix(len(s.rhs), b)
	if err != nil {
		return err
	}
	s.m, s.lc.M, s.backend, s.factory = m, m, b, nil
	return nil
}

// useSolver selects the system matrix for an analysis: the factory's when
// opt.Solver is set, otherwise opt.Backend's.
func (s *Sim) useSolver(opt *Options) error {
	if opt.Solver == nil {
		return s.useBackend(opt.Backend)
	}
	if s.factory == opt.Solver {
		return nil
	}
	m, err := opt.Solver.Matrix(len(s.rhs))
	if err != nil {
		return err
	}
	s.m, s.lc.M, s.factory = m, m, opt.Solver
	return nil
}

// polish runs up to two more Newton iterations from a converged DC
// solution and accepts the result. Newton stops once its last step is below
// reltol·|x| + vntol, which can leave a junction voltage a few microvolts
// from the exact solution; the small-signal model inherits that as an error
// of e^(δv/vt) − 1 in every transconductance. Newton converges
// quadratically, so after one or two more iterations the operating point is
// exact to round-off and AC results no longer depend on the stopping
// tolerance. Linear circuits are already exact and are left alone.
func (s *Sim) polish() error {
	if len(s.nonlin) == 0 {
		return nil
	}
	for k := 0; k < 2; k++ {
		if err := s.solve(); err != nil {
			return err
		}
		done := true
		for i, v := range s.x {
			if math.Abs(v-s.xOld[i]) > 1e-12*math.Max(1, math.Abs(v)) {
				done = false
			}
		}
		if done {
			break
		}
	}
	s.accept(0, device.ModeOP)
	return nil
}

// roundoff estimates the ledger residual that floating-point round-off
// alone can produce at the current solution (see output.Ledger.Roundoff).
func (s *Sim) roundoff() float64 {
	s.m.AbsProduct(s.x, s.absBuf)
	sum := 0.0
	for i, v := range s.absBuf {
		sum += math.Abs(s.x[i]) * (v + math.Abs(s.rhs[i]))
	}
	return 64 * 0x1p-52 * sum
}

// Backend reports the linear solver the simulation uses.
func (s *Sim) Backend() solver.Backend { return s.backend }

// solve stamps every device, linearized around the current iterate, and
// solves for the next iterate.
func (s *Sim) solve() error {
	copy(s.xOld, s.x)
	s.m.Clear()
	for i := range s.rhs {
		s.rhs[i] = 0
	}
	for _, d := range s.devs {
		if err := d.Load(&s.lc); err != nil {
			return fmt.Errorf("load %s: %w", d.Name(), err)
		}
	}
	if g := s.ptG; g > 0 {
		// Pseudo-transient continuation: a unit capacitance from every
		// node to ground and a unit inductance in series with every
		// branch, integrated with backward Euler (see pseudoTransient).
		for i, v := range s.isV {
			if v {
				s.m.Add(i, i, g)
				s.rhs[i] += g * s.ptPrev[i]
			} else {
				s.m.Add(i, i, -g)
				s.rhs[i] -= g * s.ptPrev[i]
			}
		}
	}
	if g := s.diagGmin; g > 0 {
		for i, v := range s.isV {
			if v {
				s.m.Add(i, i, g)
			}
		}
	}
	return s.m.Solve(s.rhs, s.x)
}

func (s *Sim) accept(t float64, mode device.Mode) {
	for _, d := range s.devs {
		d.Accept(s.x, t, mode)
	}
}

// OP computes the DC operating point of c.
func OP(c *circuit.Circuit, opt Options) (*output.Result, error) {
	s, err := New(c)
	if err != nil {
		return nil, err
	}
	return s.OP(opt)
}

// OP computes the DC operating point.
func (s *Sim) OP(opt Options) (*output.Result, error) {
	opt.defaults()
	if err := s.useSolver(&opt); err != nil {
		return nil, err
	}
	s.lc.Mode, s.lc.Time = device.ModeOP, 0
	if err := s.dcSolve(&opt); err != nil {
		return nil, err
	}
	s.accept(0, device.ModeOP)
	r := s.result("op")
	if !opt.SkipEnergyCheck {
		if err := em.CheckLedger(r.Ledger, opt.ledgerTol()); err != nil {
			return r, err
		}
	}
	return r, nil
}

// result builds a full result, with the element × channel ledger, from the
// current solution.
func (s *Sim) result(kind string) *output.Result {
	r := &output.Result{
		Title:        s.c.Title,
		Analysis:     kind,
		NodeNames:    append([]string(nil), s.c.NodeNames...),
		NodeVoltages: make([]float64, s.c.NumNodes()),
	}
	for i := 1; i < s.c.NumNodes(); i++ {
		r.NodeVoltages[i] = s.x[i-1]
	}
	for _, d := range s.devs {
		op := d.Operating(s.x)
		nodes := d.Nodes()
		er := output.ElementResult{
			Name:          d.Name(),
			Kind:          d.Kind().String(),
			Nodes:         [2]string{s.c.NodeNames[nodes[0]], s.c.NodeNames[nodes[1]]},
			Voltage:       op.Voltage,
			Current:       op.Current,
			PowerAbsorbed: op.PowerAbsorbed,
			Ports:         d.Ports(),
		}
		if rd, ok := d.(device.Reactive); ok {
			er.Stores = true
			er.StoragePower = op.StoragePower
			er.StoredEnergy = rd.StoredEnergy(op)
		}
		if ms, ok := d.(device.ModelSupplier); ok {
			er.ModelSupplied = ms.ModelSuppliedPower(op)
		}
		for _, p := range er.Ports {
			if p.Status != exchange.NotModeled {
				er.ExteriorPower[p.Channel] = d.ExteriorPower(p.Channel, op)
			}
		}
		r.Elements = append(r.Elements, er)
	}
	r.Ledger = em.BuildLedger(r)
	r.Ledger.Roundoff = s.roundoff()
	return r
}
