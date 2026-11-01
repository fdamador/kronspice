// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"sort"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

// Nonlinear circuits in Kron's frame.
//
// In Kron's construction a network is a primitive network — every element
// on its own, in its own reference frame — joined by a connection tensor
// (Kron 1939 p. 133: separate networks, of any kind, interconnected by a
// C). For a nonlinear element the primitive tensor depends on the state:
// at a Newton iterate each device is a small primitive network whose
// admittance Y_d(x) and impressed currents J_d(x) are its linearization,
// written in the device's own local axes (its terminals and internal
// unknowns). The connection tensor C_d maps those axes to the network's
// node-pair (open-path) axes and does not change: the iron and copper do
// not move relative to each other. The transformed equations are
//
//	Y′(x) = Σ_d C_dᵀ·Y_d(x)·C_d,   I′(x) = Σ_d C_dᵀ·J_d(x),   Y′·x_new = I′
//
// Kron's holonomic frame with a state-dependent metric, solved by Newton's
// method. The node-pair frame is the open-path frame of CLAUDE.md §3c: a
// voltage-controlled (N-type) element enters directly, a current-
// controlled (S-type) one through its current as an axis, so it accepts
// every element the engine has. Machines (circuit.Machine) enter the same
// way, with their angle and winding currents as axes; Kron's own dynamical
// view of them is in machine.go.
//
// The mesh (closed-path) frame cannot take a voltage-controlled
// nonlinearity, which is why NonlinearOP and NonlinearAC work in the
// node-pair frame only (Kron 1939 p. 411: an axis is a junction pair when
// its E is unknown).

// ErrNoConvergence is returned when Newton's method in Kron's frame does
// not converge, after gmin and source stepping.
var ErrNoConvergence = errors.New("kron: no convergence")

// PrimitiveDevice is one device's primitive network at an iterate: its
// local axes (global unknown indices), admittance tensor and impressed
// currents in those axes.
type PrimitiveDevice struct {
	Name string
	Axes []int
	Y    [][]float64
	J    []float64
}

// recorder captures one device's stamps.
type recorder struct {
	i, j []int
	v    []float64
}

func (r *recorder) Add(i, j int, v float64) {
	if i < 0 || j < 0 {
		return
	}
	r.i, r.j, r.v = append(r.i, i), append(r.j, j), append(r.v, v)
}

// kronSim holds the devices of a nonlinear Kron solve.
type kronSim struct {
	c    *circuit.Circuit
	devs []device.Device
	n    int
	isV  []bool
	nl   []device.Nonlinear
	lc   device.LoadContext
	x    []float64
	xOld []float64
	gmin float64 // stepping conductance on every node axis
}

func newKronSim(c *circuit.Circuit) (*kronSim, error) {
	if err := c.Validate(); err != nil {
		return nil, err
	}
	devs, err := device.FromCircuit(c)
	if err != nil {
		return nil, err
	}
	ctx := device.NewSetupContext(c.NumNodes())
	s := &kronSim{c: c, devs: devs}
	for _, d := range devs {
		if err := d.Setup(ctx); err != nil {
			return nil, fmt.Errorf("setup %s: %w", d.Name(), err)
		}
		if nl, ok := d.(device.Nonlinear); ok {
			s.nl = append(s.nl, nl)
		}
	}
	for _, d := range devs {
		if l, ok := d.(device.Linker); ok {
			if err := l.Link(ctx); err != nil {
				return nil, fmt.Errorf("setup %s: %w", d.Name(), err)
			}
		}
	}
	s.n = ctx.Size()
	s.isV = ctx.VoltageUnknowns()
	s.x, s.xOld = make([]float64, s.n), make([]float64, s.n)
	s.lc = device.LoadContext{X: s.xOld, Gmin: 1e-12, Mode: device.ModeOP}
	return s, nil
}

// primitives linearizes every device at the present iterate and returns
// its primitive network (Y_d, J_d in the device's local axes).
func (s *kronSim) primitives() ([]PrimitiveDevice, error) {
	out := make([]PrimitiveDevice, 0, len(s.devs))
	rhs := make([]float64, s.n)
	for _, d := range s.devs {
		var r recorder
		for i := range rhs {
			rhs[i] = 0
		}
		s.lc.M, s.lc.RHS = &r, rhs
		if err := d.Load(&s.lc); err != nil {
			return nil, fmt.Errorf("load %s: %w", d.Name(), err)
		}
		axes := map[int]int{}
		var list []int
		note := func(g int) {
			if _, ok := axes[g]; !ok {
				axes[g] = -1
				list = append(list, g)
			}
		}
		for k := range r.v {
			note(r.i[k])
			note(r.j[k])
		}
		for g, v := range rhs {
			if v != 0 {
				note(g)
			}
		}
		sort.Ints(list)
		for li, g := range list {
			axes[g] = li
		}
		p := PrimitiveDevice{Name: d.Name(), Axes: list, Y: square(len(list)), J: make([]float64, len(list))}
		for k := range r.v {
			p.Y[axes[r.i[k]]][axes[r.j[k]]] += r.v[k]
		}
		for li, g := range list {
			p.J[li] = rhs[g]
		}
		out = append(out, p)
	}
	return out, nil
}

// connect forms Kron's transformed tensors Y′ = Σ Cᵀ·Y_d·C and
// I′ = Σ Cᵀ·J_d, C_d being the incidence of the device's local axes on
// the network's axes.
func connect(n int, prims []PrimitiveDevice) ([][]float64, []float64) {
	Y := square(n)
	I := make([]float64, n)
	for _, p := range prims {
		for a, ga := range p.Axes {
			I[ga] += p.J[a]
			for b, gb := range p.Axes {
				Y[ga][gb] += p.Y[a][b]
			}
		}
	}
	return Y, I
}

func (s *kronSim) step() error {
	copy(s.xOld, s.x)
	prims, err := s.primitives()
	if err != nil {
		return err
	}
	Y, I := connect(s.n, prims)
	if s.gmin > 0 {
		for i, v := range s.isV {
			if v {
				Y[i][i] += s.gmin
			}
		}
	}
	m, err := solver.NewMatrix(s.n, solver.Auto)
	if err != nil {
		return err
	}
	for i, row := range Y {
		for j, v := range row {
			if v != 0 {
				m.Add(i, j, v)
			}
		}
	}
	if err := m.Solve(I, s.x); err != nil {
		if errors.Is(err, solver.ErrSingularMatrix) {
			return fmt.Errorf("%w: %v", ErrSingular, err)
		}
		return err
	}
	return nil
}

// NonlinearOptions are the Newton controls (zero values: SPICE defaults).
type NonlinearOptions struct {
	RelTol, AbsTol, VnTol float64
	ITL1                  int
}

func (o *NonlinearOptions) defaults() {
	if o.RelTol <= 0 {
		o.RelTol = 1e-3
	}
	if o.AbsTol <= 0 {
		o.AbsTol = 1e-12
	}
	if o.VnTol <= 0 {
		o.VnTol = 1e-6
	}
	if o.ITL1 <= 0 {
		o.ITL1 = 100
	}
}

func (s *kronSim) newton(o NonlinearOptions) error {
	for it := 0; it < o.ITL1; it++ {
		if err := s.step(); err != nil {
			return err
		}
		if len(s.nl) == 0 {
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
	return fmt.Errorf("%w after %d iterations", ErrNoConvergence, o.ITL1)
}

func (s *kronSim) converged(o NonlinearOptions) bool {
	for _, nl := range s.nl {
		if nl.Limited() {
			return false
		}
	}
	for i, x := range s.x {
		tol := o.AbsTol
		if s.isV[i] {
			tol = o.VnTol
		}
		if math.Abs(x-s.xOld[i]) > tol+o.RelTol*math.Max(math.Abs(x), math.Abs(s.xOld[i])) {
			return false
		}
	}
	for _, nl := range s.nl {
		if !nl.ConvTest(s.x, o.RelTol, o.AbsTol) {
			return false
		}
	}
	return true
}

// solveDC finds the operating point: Newton, then gmin stepping, then
// source stepping (as the MNA engine, CLAUDE.md §3), then two polishing
// iterations so the result is exact to round-off.
func (s *kronSim) solveDC(o NonlinearOptions) error {
	retry := func(err error) bool { return errors.Is(err, ErrNoConvergence) || errors.Is(err, ErrSingular) }
	s.lc.InitJct = true
	err := s.newton(o)
	if err != nil && len(s.nl) > 0 && retry(err) {
		err = s.gminStep(o)
	}
	if err != nil && len(s.nl) > 0 && retry(err) {
		err = s.sourceStep(o)
	}
	if err != nil {
		return err
	}
	for k := 0; k < 2 && len(s.nl) > 0; k++ {
		if err := s.step(); err != nil {
			return err
		}
	}
	return nil
}

func (s *kronSim) gminStep(o NonlinearOptions) error {
	defer func() { s.gmin = 0 }()
	for i := range s.x {
		s.x[i] = 0
	}
	s.lc.InitJct = true
	for g := 1e-2; ; g /= 10 {
		if g < 1e-12 {
			g = 0
		}
		s.gmin = g
		if err := s.newton(o); err != nil {
			return err
		}
		if g == 0 {
			return nil
		}
	}
}

func (s *kronSim) sourceStep(o NonlinearOptions) error {
	defer func() { s.lc.SourceStepping, s.lc.SrcFact = false, 0 }()
	for i := range s.x {
		s.x[i] = 0
	}
	s.lc.SourceStepping = true
	save := make([]float64, s.n)
	done, step := 0.0, 0.1
	s.lc.SrcFact = 0
	if err := s.newton(o); err != nil {
		return err
	}
	copy(save, s.x)
	for n := 0; n < 1000; n++ {
		f := math.Min(done+step, 1)
		s.lc.SrcFact = f
		if err := s.newton(o); err == nil {
			if f == 1 {
				return nil
			}
			copy(save, s.x)
			done, step = f, math.Min(step*1.5, 0.5)
			continue
		}
		copy(s.x, save)
		step /= 2
		if step < 1e-3 {
			break
		}
	}
	return fmt.Errorf("%w: source stepping stalled", ErrNoConvergence)
}

// NonlinearOP solves the operating point of any circuit the engine
// accepts — diodes, transistors, B elements, machines — in Kron's
// node-pair frame (see the package comment above). The result has the same
// form as analysis.OP's: node voltages and element values.
func NonlinearOP(c *circuit.Circuit, opt NonlinearOptions) (*output.Result, error) {
	s, err := newKronSim(c)
	if err != nil {
		return nil, err
	}
	opt.defaults()
	if err := s.solveDC(opt); err != nil {
		return nil, err
	}
	for _, d := range s.devs {
		d.Accept(s.x, 0, device.ModeOP)
	}
	res := &output.Result{Title: c.Title, Analysis: "op",
		NodeNames: append([]string(nil), c.NodeNames...), NodeVoltages: make([]float64, c.NumNodes())}
	for i := 1; i < c.NumNodes(); i++ {
		res.NodeVoltages[i] = s.x[i-1]
	}
	for i, d := range s.devs {
		op := d.Operating(s.x)
		e := c.Elements[i]
		res.Elements = append(res.Elements, output.ElementResult{Name: d.Name(), Kind: d.Kind().String(),
			Nodes:   [2]string{c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]},
			Voltage: op.Voltage, Current: op.Current, PowerAbsorbed: op.PowerAbsorbed})
	}
	return res, nil
}

// NonlinearPrimitives returns the primitive networks of every device of c
// linearized at its operating point, and the operating point's unknowns:
// Kron's primitive frame of a nonlinear circuit, for inspection and for
// the twin test that their transformation equals the MNA Jacobian.
func NonlinearPrimitives(c *circuit.Circuit, opt NonlinearOptions) ([]PrimitiveDevice, []float64, error) {
	s, err := newKronSim(c)
	if err != nil {
		return nil, nil, err
	}
	opt.defaults()
	if err := s.solveDC(opt); err != nil {
		return nil, nil, err
	}
	copy(s.xOld, s.x)
	p, err := s.primitives()
	return p, append([]float64(nil), s.x...), err
}

// NonlinearAC is the small-signal response of any circuit in Kron's
// node-pair frame: the operating point by NonlinearOP's Newton method,
// then every device's small-signal primitive network Y_d(ω) at that point,
// transformed and solved at each frequency [Hz].
func NonlinearAC(c *circuit.Circuit, freqs []float64, opt NonlinearOptions) (*output.ACResult, error) {
	s, err := newKronSim(c)
	if err != nil {
		return nil, err
	}
	opt.defaults()
	if err := s.solveDC(opt); err != nil {
		return nil, err
	}
	for _, d := range s.devs {
		d.Accept(s.x, 0, device.ModeOP)
		d.Operating(s.x)
	}
	res := &output.ACResult{Title: c.Title, Freq: append([]float64(nil), freqs...),
		NodeNames: append([]string(nil), c.NodeNames...), V: make([][]complex128, c.NumNodes())}
	for i := 1; i < c.NumNodes(); i++ {
		res.V[i] = make([]complex128, len(freqs))
	}
	for _, d := range s.devs {
		res.ElementNames = append(res.ElementNames, d.Name())
		res.ElementKinds = append(res.ElementKinds, d.Kind().String())
		res.I = append(res.I, make([]complex128, len(freqs)))
	}
	x := make([]complex128, s.n)
	rhs := make([]complex128, s.n)
	for fi, f := range freqs {
		if !(f > 0) || math.IsInf(f, 0) {
			return nil, fmt.Errorf("%w: frequency %g", ErrUnsupported, f)
		}
		Y := zeros(s.n, s.n)
		for i := range rhs {
			rhs[i] = 0
		}
		for _, d := range s.devs {
			// Each device's primitive Y_d(ω), scattered by its incidence:
			// Σ C_dᵀ·Y_d·C_d.
			cr := &crecorder{}
			drhs := make([]complex128, s.n)
			if err := d.LoadAC(&device.ACContext{M: cr, RHS: drhs, Omega: 2 * math.Pi * f}); err != nil {
				return nil, fmt.Errorf("load AC %s: %w", d.Name(), err)
			}
			for k := range cr.v {
				Y[cr.i[k]][cr.j[k]] += cr.v[k]
			}
			for i, v := range drhs {
				rhs[i] += v
			}
		}
		sol, err := solveSystem(Y, rhs)
		if err != nil {
			return nil, fmt.Errorf("at %g Hz: %w", f, err)
		}
		copy(x, sol)
		for i := 1; i < c.NumNodes(); i++ {
			res.V[i][fi] = x[i-1]
		}
		for ei, d := range s.devs {
			res.I[ei][fi] = d.OperatingAC(x).Current
		}
	}
	return res, nil
}

type crecorder struct {
	i, j []int
	v    []complex128
}

func (r *crecorder) Add(i, j int, v complex128) {
	if i < 0 || j < 0 {
		return
	}
	r.i, r.j, r.v = append(r.i, i), append(r.j, j), append(r.v, v)
}
