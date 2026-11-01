// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package device implements circuit elements and how they stamp into the
// modified nodal analysis (MNA) system (CLAUDE.md §3).
//
// Every device embeds Base, which embeds exchange.Base, so every element
// carries one external energy port per registered channel (§3.1) from the
// moment it is created.
package device

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// ErrUnsupported is returned for circuit elements this package cannot build yet.
var ErrUnsupported = errors.New("device: unsupported element kind")

// Ground is the unknown index of the ground node; stamps to it are dropped.
const Ground = -1

// Stamper receives matrix contributions. solver.Matrix satisfies it, so this
// package does not depend on any particular solver backend.
type Stamper interface {
	Add(i, j int, v float64)
}

// SetupContext maps circuit nodes to MNA unknowns and allocates branch
// equations. Unknown layout: node k (k ≥ 1) → k−1; branch b → (numNodes−1)+b.
type SetupContext struct {
	numNodes int
	branches int
	internal []int          // extra unknowns that are internal node voltages
	named    map[string]int // branch unknown of each voltage source, by lower-cased name
	// Step and Stop are the transient print step and stop time (0 outside a
	// transient analysis). Sources use them for SPICE default arguments.
	Step, Stop float64
}

// NewSetupContext prepares a context for a circuit with numNodes nodes
// (including ground).
func NewSetupContext(numNodes int) *SetupContext {
	return &SetupContext{numNodes: numNodes}
}

// Unknown returns the MNA unknown index for a circuit node (Ground for node 0).
func (s *SetupContext) Unknown(node int) int { return node - 1 }

// AllocUnknown reserves an extra unknown for an internal node voltage, such
// as a diode's node behind its series resistance.
func (s *SetupContext) AllocUnknown() int {
	i := s.AllocBranch()
	s.internal = append(s.internal, i)
	return i
}

// RegisterBranch records the branch-current unknown of a named voltage
// source, so current-controlled elements can find it (Linker).
func (s *SetupContext) RegisterBranch(name string, i int) {
	if s.named == nil {
		s.named = map[string]int{}
	}
	s.named[strings.ToLower(name)] = i
}

// BranchOf returns the branch-current unknown of a voltage source.
func (s *SetupContext) BranchOf(name string) (int, bool) {
	i, ok := s.named[strings.ToLower(name)]
	return i, ok
}

// VoltageUnknowns reports, for every unknown, whether it is a node voltage
// (true) or a branch current (false). Newton convergence uses vntol for
// voltages and abstol for currents.
func (s *SetupContext) VoltageUnknowns() []bool {
	v := make([]bool, s.Size())
	for i := 0; i < s.numNodes-1; i++ {
		v[i] = true
	}
	for _, i := range s.internal {
		v[i] = true
	}
	return v
}

// AllocBranch reserves a branch-current unknown and returns its index.
func (s *SetupContext) AllocBranch() int {
	i := s.numNodes - 1 + s.branches
	s.branches++
	return i
}

// Size is the total number of unknowns after all devices are set up.
func (s *SetupContext) Size() int { return s.numNodes - 1 + s.branches }

// Mode says what kind of solve is being loaded.
type Mode uint8

const (
	// ModeOP is a stand-alone DC operating point: sources use their DC value.
	ModeOP Mode = iota
	// ModeTranOP is the operating point at t = 0 that starts a transient:
	// sources use their waveform value at t = 0.
	ModeTranOP
	// ModeUICInit imposes initial conditions (UIC) at t = 0 to report a
	// consistent starting point; reactive devices hold their IC values.
	ModeUICInit
	// ModeTran is one transient time step.
	ModeTran
)

// Method is the integration formula for a transient step.
type Method uint8

const (
	BackwardEuler Method = 1 // order 1; used for the first step and after breakpoints
	Trapezoidal   Method = 2 // order 2; default
)

// LoadContext receives stamps for one solve.
type LoadContext struct {
	M      Stamper
	RHS    []float64
	Mode   Mode
	Time   float64 // time of the solution being computed [s]
	H      float64 // step size from the last accepted point [s] (ModeTran)
	Method Method  // integration method (ModeTran)

	// X is the current Newton iterate; nonlinear devices linearize around it.
	X []float64
	// InitJct asks nonlinear devices to start junctions at their critical
	// voltage instead of X (SPICE MODEINITJCT), on the first iteration of a
	// cold operating point.
	InitJct bool
	// Gmin is the minimum conductance added across every pn junction
	// (SPICE gmin, default 1e-12 S).
	Gmin float64
	// SourceStepping scales every independent source by SrcFact during a DC
	// solve (SPICE source stepping).
	SourceStepping bool
	SrcFact        float64
}

// AddRHS adds v to right-hand-side entry i, ignoring ground.
func (l *LoadContext) AddRHS(i int, v float64) {
	if i >= 0 && i < len(l.RHS) {
		l.RHS[i] += v
	}
}

// CStamper receives complex matrix contributions; solver.CDense satisfies it.
type CStamper interface {
	Add(i, j int, v complex128)
}

// ACContext receives small-signal stamps for one frequency.
type ACContext struct {
	M     CStamper
	RHS   []complex128
	Omega float64 // angular frequency [rad/s]
}

// AddRHS adds v to right-hand-side entry i, ignoring ground.
func (a *ACContext) AddRHS(i int, v complex128) {
	if i >= 0 && i < len(a.RHS) {
		a.RHS[i] += v
	}
}

// ACPoint is a two-terminal element's small-signal phasors (peak amplitude,
// as in SPICE): voltage V(n+) − V(n−) and current from n+ to n−.
type ACPoint struct {
	Voltage, Current complex128
}

// OperatingPoint describes a two-terminal element at a solution.
type OperatingPoint struct {
	Voltage       float64 // V(n+) − V(n−) [V]
	Current       float64 // current through the element from n+ to n− [A]
	PowerAbsorbed float64 // Voltage·Current [W]; negative when the element supplies power
	// StoragePower is the part of PowerAbsorbed going into stored energy
	// (capacitor, inductor, junction charge); the rest is dissipated.
	StoragePower float64
}

// Device is one circuit element.
type Device interface {
	exchange.Exchanger
	Name() string
	Kind() circuit.Kind
	Nodes() [2]int
	// Setup maps nodes to unknowns and allocates branch equations.
	Setup(ctx *SetupContext) error
	// Load stamps the element for the current solve. It must not allocate.
	Load(ctx *LoadContext) error
	// Operating evaluates the element at solution vector x, using the
	// context of the most recent Load.
	Operating(x []float64) OperatingPoint
	// LoadAC stamps the small-signal model, linearized at the operating point
	// of the most recent Load/Accept, for one frequency.
	LoadAC(ctx *ACContext) error
	// OperatingAC evaluates the element's phasors at AC solution x.
	OperatingAC(x []complex128) ACPoint
	// Accept commits solution x at time t as the element's accepted state.
	Accept(x []float64, t float64, mode Mode)
	// ExteriorPower returns the power exchanged with the exterior through the
	// given channel's port at operating point op: + into the circuit, − out.
	// Unmodeled ports return exactly 0.
	ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64
}

// Reactive devices store energy and limit the transient step size.
type Reactive interface {
	Device
	// StoredEnergy returns the energy stored at operating point op [J].
	StoredEnergy(op OperatingPoint) float64
	// Truncate proposes a step size from the local truncation error of the
	// candidate solution x (CLAUDE.md §3, Transient). +Inf means no limit.
	Truncate(x []float64, tc *TruncContext) float64
}

// ACPowerer is implemented by elements with more than two terminals, whose
// complex power is not the product of one voltage and one current: it
// returns ½·Σ_k V_k·I_k* over the terminals (currents into the element) at
// AC solution x.
type ACPowerer interface {
	ACPower(x []complex128) complex128
}

// Linker is implemented by devices that refer to another device's unknowns
// (a current-controlled element needs the controlling source's branch).
// Link runs after every device's Setup.
type Linker interface {
	Link(ctx *SetupContext) error
}

// ModelSupplier is implemented by idealized active elements — a negative
// resistance now, controlled sources later — that can deliver power without
// a declared external port. That power is reported per element as
// "model-supplied" (CLAUDE.md §3c ledger table), so it is visible and never
// counted as inflow from the exterior. A real device of this kind draws the
// energy from a supply the idealized model leaves out.
type ModelSupplier interface {
	// ModelSuppliedPower is the power the element delivers to the circuit
	// [W], ≥ 0 when it supplies.
	ModelSuppliedPower(op OperatingPoint) float64
}

// Nonlinear devices are linearized around the Newton iterate.
type Nonlinear interface {
	Device
	// Limited reports whether the most recent Load limited a junction
	// voltage; a Newton iteration with limiting never counts as converged.
	Limited() bool
	// ConvTest is the device-level convergence check at the new iterate x:
	// the device current predicted by the linearization must agree with the
	// current at the linearization point (SPICE3 DIOconvTest).
	ConvTest(x []float64, reltol, abstol float64) bool
}

// Breakpointer devices have waveform corners the engine must land on.
type Breakpointer interface {
	NextBreakpoint(t float64) float64
}

// Base is embedded by every device: name, nodes and the full port set.
type Base struct {
	exchange.Base
	name  string
	kind  circuit.Kind
	nodes [2]int
	u     [2]int // unknown indices of the nodes
}

func newBase(e circuit.Element) Base {
	return Base{name: e.Name, kind: e.Kind, nodes: e.Nodes}
}

// Name returns the element name.
func (b *Base) Name() string { return b.name }

// Kind returns the element kind.
func (b *Base) Kind() circuit.Kind { return b.kind }

// Nodes returns the circuit node indices.
func (b *Base) Nodes() [2]int { return b.nodes }

// Accept is a no-op for memoryless devices.
func (b *Base) Accept([]float64, float64, Mode) {}

func (b *Base) setupNodes(ctx *SetupContext) {
	b.u[0] = ctx.Unknown(b.nodes[0])
	b.u[1] = ctx.Unknown(b.nodes[1])
}

// voltage returns V(n+) − V(n−) from solution x.
func (b *Base) voltage(x []float64) float64 {
	v := 0.0
	if b.u[0] >= 0 {
		v += x[b.u[0]]
	}
	if b.u[1] >= 0 {
		v -= x[b.u[1]]
	}
	return v
}

// voltageAC returns the phasor V(n+) − V(n−) from AC solution x.
func (b *Base) voltageAC(x []complex128) complex128 {
	var v complex128
	if b.u[0] >= 0 {
		v += x[b.u[0]]
	}
	if b.u[1] >= 0 {
		v -= x[b.u[1]]
	}
	return v
}

func stampY(m CStamper, a, b int, y complex128) {
	m.Add(a, a, y)
	m.Add(b, b, y)
	m.Add(a, b, -y)
	m.Add(b, a, -y)
}

// DCOverrider sources accept a temporary DC value, used by .DC sweeps.
type DCOverrider interface {
	SetDCOverride(v float64)
	ClearDCOverride()
}

// FromCircuit builds one device per circuit element.
func FromCircuit(c *circuit.Circuit) ([]Device, error) {
	devs := make([]Device, 0, len(c.Elements))
	for _, e := range c.Elements {
		var d Device
		// A source in a non-electrical domain draws its energy from that
		// domain's exterior (a force or velocity source from the mechanical
		// exterior, a pressure source from the acoustic one) unless its
		// channel was set.
		if (e.Kind == circuit.VoltageSource || e.Kind == circuit.CurrentSource) && e.Channel == exchange.Electrical {
			e.Channel = c.ElementDomain(e)
		}
		switch e.Kind {
		case circuit.Diode:
			m, ok := c.ModelFor(e)
			if !ok {
				return nil, fmt.Errorf("%w: %s has no model", ErrUnsupported, e.Name)
			}
			d = NewDiode(e, m)
		case circuit.Resistor:
			d = NewResistor(e)
		case circuit.Capacitor:
			d = NewCapacitor(e)
		case circuit.Inductor:
			d = NewInductor(e)
		case circuit.NonlinearCapacitor:
			d = NewNonlinearCapacitor(e)
		case circuit.NonlinearInductor:
			d = NewNonlinearInductor(e)
		case circuit.Gyrator:
			d = NewGyrator(e)
		case circuit.BearingMotor:
			d = NewBearing(e)
		case circuit.Machine:
			spec, err := c.MachineFor(e)
			if err != nil {
				return nil, err
			}
			d = NewMachine(e, spec)
		case circuit.Mutual:
			ls := [2]float64{}
			for i, n := range e.Couples {
				j := c.Index(n)
				if j < 0 || c.Elements[j].Kind != circuit.Inductor {
					return nil, fmt.Errorf("%w: %s: %s is not an inductor", ErrUnsupported, e.Name, n)
				}
				ls[i] = c.Elements[j].Value
			}
			d = NewMutual(e, ls[0], ls[1])
		case circuit.VoltageSource:
			d = NewVoltageSource(e)
		case circuit.CurrentSource:
			d = NewCurrentSource(e)
		case circuit.Behavioral:
			d = NewBehavioral(e, c.SelfControlled(e))
		case circuit.VCVS, circuit.CCCS, circuit.VCCS, circuit.CCVS:
			d = NewControlled(e)
		case circuit.MOSFET:
			m, ok := c.ModelFor(e)
			if !ok {
				return nil, fmt.Errorf("%w: %s has no model", ErrUnsupported, e.Name)
			}
			q, err := NewMOSFET(e, m)
			if err != nil {
				return nil, err
			}
			d = q
		case circuit.BJT:
			m, ok := c.ModelFor(e)
			if !ok {
				return nil, fmt.Errorf("%w: %s has no model", ErrUnsupported, e.Name)
			}
			d = NewBJT(e, m)
		default:
			return nil, fmt.Errorf("%w: %s (%v)", ErrUnsupported, e.Name, e.Kind)
		}
		devs = append(devs, d)
	}
	return devs, nil
}

// inf is +Inf, the "no limit" step proposal.
var inf = math.Inf(1)
