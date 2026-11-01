// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package circuit is the in-memory circuit model shared by the netlist
// parser, the genetic algorithm (through Builder) and the analyses.
//
// Node 0 is ground. The names "0" and "gnd" both refer to it.
package circuit

import (
	"errors"
	"fmt"
	"github.com/fdamador/kronspice/exchange"
	"math"
	"sort"
	"strings"
)

// Errors returned by validation. They are sentinel errors so the GA can map
// each class to a penalty with errors.Is.
var (
	ErrNoGround      = errors.New("circuit: no ground node (0 or gnd)")
	ErrFloatingNode  = errors.New("circuit: node has no path to ground")
	ErrDuplicateName = errors.New("circuit: duplicate element name")
	ErrBadValue      = errors.New("circuit: invalid element value")
	ErrBadWaveform   = errors.New("circuit: invalid source waveform")
	ErrUnknownKind   = errors.New("circuit: unknown element kind")
	ErrSelfConnected = errors.New("circuit: element connects a node to itself")
	ErrEmptyCircuit  = errors.New("circuit: no elements")
	ErrTooLarge      = errors.New("circuit: too many nodes or elements")
)

// Limits protect against runaway evolved circuits exhausting memory.
const (
	MaxNodes    = 100000
	MaxElements = 1000000
	MaxPWL      = 100000 // PWL points per source
)

// Kind is an element type.
type Kind byte

const (
	Resistor      Kind = 'R'
	Capacitor     Kind = 'C'
	Inductor      Kind = 'L'
	VoltageSource Kind = 'V'
	CurrentSource Kind = 'I'
	Diode         Kind = 'D'
	// Behavioral is a nonlinear element whose output (current or voltage)
	// is a piecewise-linear function of one controlling quantity (a node
	// voltage or a voltage source's current), as an ngspice B source.
	Behavioral Kind = 'B'
	// Linear controlled sources (SPICE E, F, G, H). Value is the gain.
	// E and G are controlled by the voltage between CtrlNodes; F and H by
	// the current of the voltage source CtrlSource.
	VCVS Kind = 'E' // voltage-controlled voltage source, gain [V/V]
	CCCS Kind = 'F' // current-controlled current source, gain [A/A]
	VCCS Kind = 'G' // voltage-controlled current source, transconductance [S]
	CCVS Kind = 'H' // current-controlled voltage source, transresistance [Ω]
	// BJT is a bipolar junction transistor (SPICE Q, Gummel–Poon). Nodes
	// holds {collector, emitter}; Extra holds {base} or {base, substrate}.
	BJT Kind = 'Q'
	// MOSFET is a MOS transistor (SPICE M, levels 1–3). Nodes holds {drain,
	// source}; Extra holds {gate, bulk}; Params the instance parameters
	// (l, w, ad, as, pd, ps, nrd, nrs, m).
	MOSFET Kind = 'M'
	// Mutual is a mutual inductance between two inductors (SPICE K):
	// Couples names them, Value is the coupling coefficient k (0 < k ≤ 1),
	// M = k·√(L1·L2). It has no terminals of its own (Nodes are 0, 0).
	Mutual Kind = 'K'
	// Gyrator is an ideal, lossless two-port transducer (kronspice Y):
	// Nodes = {p1+, p1−}, Extra = {p2+, p2−}, Value = g, with
	// i1 = g·V2 and i2 = −g·V1 (currents p+ → p− through each port), so
	// V1·i1 + V2·i2 = 0. Between an electrical and a mechanical port it is a
	// piezoelectric (or electrostatic) transducer of force factor g [N/V]
	// in the mobility analogy; between mechanical and acoustic ports, a
	// piston of area g [m²].
	Gyrator Kind = 'Y'
	// NonlinearCapacitor is a capacitor defined by its charge,
	// q(v, t) = s(t)·Q(v): Q is the piecewise-linear Curve v1, q1, v2, q2, …
	// (or Value·v when Curve is empty) and s the dimensionless Scale
	// waveform (1 when none). Written "C … Q=pwl(…) SCALE=…" in a netlist.
	NonlinearCapacitor Kind = 'c'
	// NonlinearInductor is an inductor defined by its flux linkage,
	// λ(i, t) = s(t)·Λ(i), with Curve i1, λ1, i2, λ2, … (or Value·i) and
	// Scale as for NonlinearCapacitor. Written "L … FLUX=pwl(…) SCALE=…".
	NonlinearInductor Kind = 'l'
)

// IsControlled reports whether k is a linear controlled source.
func (k Kind) IsControlled() bool { return k == VCVS || k == CCCS || k == VCCS || k == CCVS }

// VoltageControlled reports whether a controlled source senses a voltage.
func (k Kind) VoltageControlled() bool { return k == VCVS || k == VCCS }

// VoltageOutput reports whether a controlled source forces a voltage.
func (k Kind) VoltageOutput() bool { return k == VCVS || k == CCVS }

func (k Kind) String() string {
	switch k {
	case Resistor:
		return "resistor"
	case Capacitor:
		return "capacitor"
	case Inductor:
		return "inductor"
	case VoltageSource:
		return "voltage source"
	case CurrentSource:
		return "current source"
	case Diode:
		return "diode"
	case Behavioral:
		return "nonlinear element"
	case VCVS:
		return "VCVS"
	case CCCS:
		return "CCCS"
	case VCCS:
		return "VCCS"
	case CCVS:
		return "CCVS"
	case BJT:
		return "BJT"
	case MOSFET:
		return "MOSFET"
	case Mutual:
		return "K"
	case Gyrator:
		return "transducer"
	case Machine:
		return "machine"
	case BearingMotor:
		return "bearing motor"
	case NonlinearCapacitor:
		return "nonlinear capacitor"
	case NonlinearInductor:
		return "nonlinear inductor"
	default:
		return fmt.Sprintf("kind(%c)", byte(k))
	}
}

// WaveKind is a time-dependent source function.
type WaveKind byte

const (
	WaveNone  WaveKind = iota
	WavePulse          // PULSE(v1 v2 td tr tf pw per)
	WaveSin            // SIN(vo va freq td theta phase)
	WaveExp            // EXP(v1 v2 td1 tau1 td2 tau2)
	WavePWL            // PWL(t1 v1 t2 v2 ...)
)

func (w WaveKind) String() string {
	switch w {
	case WavePulse:
		return "PULSE"
	case WaveSin:
		return "SIN"
	case WaveExp:
		return "EXP"
	case WavePWL:
		return "PWL"
	default:
		return "none"
	}
}

// Waveform is a source's transient specification, as written in SPICE.
// Missing optional arguments take SPICE3 defaults at simulation time.
type Waveform struct {
	Kind WaveKind
	Args []float64
}

// Element is one circuit element. Nodes are indices into Circuit.NodeNames.
//
// For sources, Nodes[0] is the + terminal and Nodes[1] the − terminal. SPICE
// convention: a current source's current flows from n+ through the source to
// n−, i.e. it is drawn from the n+ node and delivered into the n− node; a
// voltage source's branch current is positive when it flows into n+.
type Element struct {
	Name  string
	Kind  Kind
	Nodes [2]int
	// Value is ohms (R), farads (C), henries (L), or the DC value of a source.
	Value float64
	// HasDC is true when a source's DC value was given explicitly.
	HasDC bool
	// Wave is a source's transient waveform (Kind WaveNone if constant).
	Wave Waveform
	// AC magnitude and phase (degrees) for small-signal analysis (used from M3).
	ACMag, ACPhase float64
	// IC is the initial condition: capacitor voltage or inductor current, used
	// when a transient analysis skips the operating point (UIC).
	IC    float64
	HasIC bool
	// Model names the .model card of a semiconductor element (lower-cased);
	// Area is its area factor (0 means 1).
	Model string
	Area  float64
	// Extra are terminals beyond Nodes for elements with more than two:
	// a BJT's base, and optionally its substrate (Nodes = {collector,
	// emitter}, so the element's reported voltage and current are V_CE
	// and the collector current).
	Extra []int
	// Params are instance parameters by lower-case name (a MOSFET's l, w,
	// ad, as, pd, ps, nrd, nrs, m), in SI units.
	Params map[string]float64
	// Couples are the names of the two inductors a Mutual element couples.
	Couples [2]string

	// Behavioral elements: output 'I' or 'V'; control 'V' (voltage between
	// CtrlNodes) or 'I' (current of the voltage source CtrlSource); Curve
	// holds x1, y1, x2, y2, … with strictly increasing x. Beyond the end
	// points the first and last segments are extended linearly. A
	// NonlinearCapacitor's Curve is its charge (v, q pairs) and a
	// NonlinearInductor's its flux linkage (i, λ pairs).
	BOut, BCtrl byte
	CtrlNodes   [2]int
	CtrlSource  string
	Curve       []float64
	// Scale is the time variation s(t) > 0 of a NonlinearCapacitor's charge
	// or a NonlinearInductor's flux (Kind WaveNone: constant 1). The power
	// that varying it exchanges, ṡ·∫₀ˣ Q(u) du, crosses the element's port
	// on Channel (Mechanical when Channel is Electrical): moving plates or
	// a moving core.
	Scale Waveform

	// Channel is the exterior channel an independent source draws its
	// energy from: Electrical (the default: a generator or supply outside
	// the circuit), or for example Thermal for a Seebeck EMF or Chemical
	// for a galvanic cell. The ledger books the source's power on that
	// channel, so energy is attributed to its physical origin (CLAUDE.md
	// §3.1). ChannelNote describes the source for reports.
	Channel     exchange.ChannelID
	ChannelNote string
}

// SelfControlled reports whether a behavioral element is a two-terminal
// nonlinear resistor: a current output controlled by its own terminal
// voltage (voltage-controlled, N-type), or a voltage output controlled by
// the current of a 0 V source in series with it that nothing else touches
// (current-controlled, S-type). Such an element dissipates or supplies its
// own power; any other behavioral element is a controlled source.
func (c *Circuit) SelfControlled(e Element) bool {
	if e.Kind != Behavioral {
		return false
	}
	if e.BOut == 'I' && e.BCtrl == 'V' {
		return (e.CtrlNodes == e.Nodes) || (e.CtrlNodes == [2]int{e.Nodes[1], e.Nodes[0]})
	}
	if e.BOut != 'V' || e.BCtrl != 'I' {
		return false
	}
	si := c.Index(e.CtrlSource)
	if si < 0 {
		return false
	}
	s := c.Elements[si]
	if s.Kind != VoltageSource || s.Value != 0 || s.Wave.Kind != WaveNone || s.ACMag != 0 {
		return false
	}
	// The source and the element share exactly one node, used by nothing else.
	shared := -1
	for _, a := range s.Nodes {
		for _, b := range e.Nodes {
			if a == b {
				if shared >= 0 && shared != a {
					return false
				}
				shared = a
			}
		}
	}
	if shared <= 0 {
		return false
	}
	n := 0
	for _, x := range c.Elements {
		for _, v := range x.Nodes {
			if v == shared {
				n++
			}
		}
		for _, v := range x.Extra {
			if v == shared {
				n++
			}
		}
		if (x.Kind == Behavioral && x.BCtrl == 'V') || x.Kind.VoltageControlled() {
			for _, v := range x.CtrlNodes {
				if v == shared && x.Name != e.Name {
					n += 10 // sensed elsewhere: not a private series node
				}
			}
		}
	}
	return n == 2
}

// Circuit is a validated circuit.
type Circuit struct {
	Title     string
	NodeNames []string // NodeNames[0] == "0"
	Elements  []Element
	Models    map[string]Model // by lower-cased name
	// NodeDomains is each node's energy domain, aligned with NodeNames: the
	// exchange channel whose across variable its potential is (voltage,
	// velocity, pressure, …). Nil means every node is electrical. Ground is
	// the common reference of every domain (zero velocity, ambient
	// pressure).
	NodeDomains []exchange.ChannelID
}

// NodeDomain returns node n's domain (electrical when none is set).
func (c *Circuit) NodeDomain(n int) exchange.ChannelID {
	if n <= 0 || n >= len(c.NodeDomains) {
		return exchange.Electrical
	}
	return c.NodeDomains[n]
}

// ElementDomain is the domain of an element's terminals (of its first port
// for a gyrator): the domain of its first non-ground node.
func (c *Circuit) ElementDomain(e Element) exchange.ChannelID {
	if e.Nodes[0] != 0 {
		return c.NodeDomain(e.Nodes[0])
	}
	return c.NodeDomain(e.Nodes[1])
}

// ModelFor returns the model of a semiconductor element.
func (c *Circuit) ModelFor(e Element) (Model, bool) {
	m, ok := c.Models[strings.ToLower(e.Model)]
	return m, ok
}

// NumNodes returns the number of nodes including ground.
func (c *Circuit) NumNodes() int { return len(c.NodeNames) }

// Validate checks structural soundness for a DC solution: every node needs a
// DC path to ground through resistors, inductors or voltage sources.
// Numerical problems visible only to the solver (for example a loop of
// voltage sources) are reported by the analysis as solver.ErrSingularMatrix.
func (c *Circuit) Validate() error {
	if err := c.validateElements(); err != nil {
		return err
	}
	return c.checkPaths(false)
}

// ValidateTransient is like Validate but also accepts capacitors as paths, as
// needed for transient analysis without an operating point (UIC).
func (c *Circuit) ValidateTransient() error {
	if err := c.validateElements(); err != nil {
		return err
	}
	return c.checkPaths(true)
}

func (c *Circuit) validateElements() error {
	if len(c.Elements) == 0 {
		return ErrEmptyCircuit
	}
	if len(c.NodeNames) > MaxNodes || len(c.Elements) > MaxElements {
		return ErrTooLarge
	}
	if len(c.Models) > MaxElements {
		return ErrTooLarge
	}
	mnames := make([]string, 0, len(c.Models))
	for name := range c.Models {
		mnames = append(mnames, name)
	}
	sort.Strings(mnames)
	for _, name := range mnames {
		if err := c.Models[name].validate(); err != nil {
			return err
		}
	}
	names := make(map[string]bool, len(c.Elements))
	coupled := map[[2]string]bool{}
	groundUsed := false
	for _, e := range c.Elements {
		key := strings.ToLower(e.Name)
		if names[key] {
			return fmt.Errorf("%w: %s", ErrDuplicateName, e.Name)
		}
		names[key] = true
		for _, n := range append(e.Nodes[:], e.Extra...) {
			if n < 0 || n >= len(c.NodeNames) {
				return fmt.Errorf("%w: %s has node index %d out of range", ErrBadValue, e.Name, n)
			}
			if n == 0 {
				groundUsed = true
			}
		}
		if (e.Kind == BJT || e.Kind == MOSFET || e.Kind == Gyrator || e.Kind == Machine || e.Kind == BearingMotor) != (len(e.Extra) > 0) {
			return fmt.Errorf("%w: %s has %d extra terminals", ErrBadValue, e.Name, len(e.Extra))
		}
		if e.Nodes[0] == e.Nodes[1] && e.Kind != BJT && e.Kind != MOSFET && e.Kind != Mutual && e.Kind != Machine {
			return fmt.Errorf("%w: %s", ErrSelfConnected, e.Name)
		}
		for _, v := range []float64{e.Value, e.IC, e.ACMag, e.ACPhase} {
			if math.IsNaN(v) || math.IsInf(v, 0) {
				return fmt.Errorf("%w: %s has a non-finite value", ErrBadValue, e.Name)
			}
		}
		switch e.Kind {
		case Resistor, Capacitor, Inductor:
			// A negative resistance is allowed (an idealized active element,
			// as in ngspice); its power is booked as model-supplied. Zero is
			// a short with no conductance to stamp, and C, L must be > 0.
			if e.Value == 0 || (e.Value < 0 && e.Kind != Resistor) {
				return fmt.Errorf("%w: %s %s must be > 0, got %g", ErrBadValue, e.Name, e.Kind, e.Value)
			}
			if e.Wave.Kind != WaveNone {
				return fmt.Errorf("%w: %s is not a source", ErrBadWaveform, e.Name)
			}
		case NonlinearCapacitor, NonlinearInductor:
			if err := validateStateLaw(e); err != nil {
				return err
			}
		case VoltageSource, CurrentSource:
			if err := validateWave(e.Name, e.Wave); err != nil {
				return err
			}
			if int(e.Channel) >= exchange.NumChannels || e.Channel == exchange.Neutrino ||
				e.Channel == exchange.HypotheticalMagneticCharge || e.Channel == exchange.Gravitational {
				// Neutrino and hypothetical channels need a cited coupling
				// model (§1 rule 3, §3f); gravity carries no circuit energy.
				return fmt.Errorf("%w: %s cannot draw on the %v channel", ErrBadValue, e.Name, e.Channel)
			}
		case Behavioral:
			if err := c.validateBehavioral(e); err != nil {
				return err
			}
		case Mutual:
			if err := c.validateMutual(e, coupled); err != nil {
				return err
			}
		case Gyrator:
			if len(e.Extra) != 2 || e.Extra[0] == e.Extra[1] {
				return fmt.Errorf("%w: %s needs two distinct nodes on each port", ErrBadValue, e.Name)
			}
			if e.Value == 0 {
				return fmt.Errorf("%w: %s gyration must not be zero", ErrBadValue, e.Name)
			}
		case VCVS, CCCS, VCCS, CCVS:
			if err := c.validateControl(e, e.Kind.VoltageControlled()); err != nil {
				return err
			}
		case Machine:
			if err := c.validateMachine(e); err != nil {
				return err
			}
		case BearingMotor:
			if err := validateBearing(e); err != nil {
				return err
			}
		case BJT:
			m, ok := c.ModelFor(e)
			if !ok {
				return fmt.Errorf("%w: %s uses undefined model %q", ErrBadModel, e.Name, e.Model)
			}
			if m.Type != "npn" && m.Type != "pnp" {
				return fmt.Errorf("%w: %s needs an NPN or PNP model, %q is type %s", ErrBadModel, e.Name,
					e.Model, strings.ToUpper(m.Type))
			}
			if len(e.Extra) > 2 {
				return fmt.Errorf("%w: %s has too many terminals", ErrBadValue, e.Name)
			}
			if e.Area < 0 || e.Area > 1e12 {
				return fmt.Errorf("%w: %s area %g", ErrBadValue, e.Name, e.Area)
			}
		case MOSFET:
			m, ok := c.ModelFor(e)
			if !ok {
				return fmt.Errorf("%w: %s uses undefined model %q", ErrBadModel, e.Name, e.Model)
			}
			if m.Type != "nmos" && m.Type != "pmos" {
				return fmt.Errorf("%w: %s needs an NMOS or PMOS model, %q is type %s", ErrBadModel, e.Name,
					e.Model, strings.ToUpper(m.Type))
			}
			if len(e.Extra) != 2 {
				return fmt.Errorf("%w: %s needs drain, gate, source and bulk", ErrBadValue, e.Name)
			}
			if lv, ok := m.Params["level"]; ok && lv != math.Trunc(lv) {
				return fmt.Errorf("%w: %s: LEVEL must be 1, 2 or 3", ErrBadModel, e.Model)
			}
			if err := validateMOSInstance(e); err != nil {
				return err
			}
		case Diode:
			m, ok := c.ModelFor(e)
			if !ok {
				return fmt.Errorf("%w: %s uses undefined model %q", ErrBadModel, e.Name, e.Model)
			}
			if m.Type != kindModelType[e.Kind] {
				return fmt.Errorf("%w: %s needs a %s model, %q is type %s", ErrBadModel, e.Name,
					strings.ToUpper(kindModelType[e.Kind]), e.Model, strings.ToUpper(m.Type))
			}
			if e.Area < 0 || e.Area > 1e12 {
				return fmt.Errorf("%w: %s area %g", ErrBadValue, e.Name, e.Area)
			}
		default:
			return fmt.Errorf("%w: %s", ErrUnknownKind, e.Name)
		}
	}
	if !groundUsed {
		return ErrNoGround
	}
	return c.validateDomains()
}

// validateDomains checks the energy domains: every node is in a domain with
// across and through variables; an element's terminals share one domain
// (each port of a gyrator its own; ground belongs to all); semiconductors
// and behavioral elements are electrical; a resistor's exterior channel
// can carry energy.
func (c *Circuit) validateDomains() error {
	if c.NodeDomains != nil && len(c.NodeDomains) != len(c.NodeNames) {
		return fmt.Errorf("%w: %d node domains for %d nodes", ErrBadValue, len(c.NodeDomains), len(c.NodeNames))
	}
	for i, d := range c.NodeDomains {
		if !d.Valid() || d.Kind() != exchange.AcrossThrough {
			return fmt.Errorf("%w: node %s: %v is not a domain with across and through variables", ErrBadValue, c.NodeNames[i], d)
		}
	}
	same := func(e Element, a, b int) error {
		if a != 0 && b != 0 && c.NodeDomain(a) != c.NodeDomain(b) {
			return fmt.Errorf("%w: %s joins %s (%v) and %s (%v): only a transducer crosses domains",
				ErrBadValue, e.Name, c.NodeNames[a], c.NodeDomain(a), c.NodeNames[b], c.NodeDomain(b))
		}
		return nil
	}
	for _, e := range c.Elements {
		switch e.Kind {
		case Mutual:
			continue
		case Gyrator:
			if err := same(e, e.Nodes[0], e.Nodes[1]); err != nil {
				return err
			}
			if err := same(e, e.Extra[0], e.Extra[1]); err != nil {
				return err
			}
			continue
		case BearingMotor:
			// The current path is electrical; the shaft port is mechanical.
			for _, n := range e.Nodes {
				if n != 0 && c.NodeDomain(n) != exchange.Electrical {
					return fmt.Errorf("%w: %s current node %s is %v, not electrical", ErrBadValue, e.Name, c.NodeNames[n], c.NodeDomain(n))
				}
			}
			for _, n := range e.Extra {
				if n != 0 && c.NodeDomain(n) != exchange.Mechanical {
					return fmt.Errorf("%w: %s shaft node %s is %v, not mechanical", ErrBadValue, e.Name, c.NodeNames[n], c.NodeDomain(n))
				}
			}
			continue
		case Machine:
			// Windings are electrical; the port is mechanical (angular
			// velocity across, torque through).
			w, port := MachineWindings(e)
			for _, p := range w {
				for _, n := range p {
					if c.NodeDomain(n) != exchange.Electrical {
						return fmt.Errorf("%w: %s winding node %s is %v, not electrical", ErrBadValue, e.Name, c.NodeNames[n], c.NodeDomain(n))
					}
				}
			}
			for _, n := range port {
				if n != 0 && c.NodeDomain(n) != exchange.Mechanical {
					return fmt.Errorf("%w: %s port node %s is %v, not mechanical", ErrBadValue, e.Name, c.NodeNames[n], c.NodeDomain(n))
				}
			}
			continue
		case Diode, BJT, MOSFET, Behavioral:
			for _, n := range append(e.Nodes[:], e.Extra...) {
				if c.NodeDomain(n) != exchange.Electrical {
					return fmt.Errorf("%w: %s is electrical, node %s is %v", ErrBadValue, e.Name, c.NodeNames[n], c.NodeDomain(n))
				}
			}
			continue
		case Resistor, NonlinearCapacitor, NonlinearInductor:
			switch ch := e.Channel; {
			case ch == exchange.Electrical:
			case !ch.Valid() || ch.Kind() == exchange.FrameOnly || ch == exchange.Neutrino || ch == exchange.HypotheticalMagneticCharge:
				return fmt.Errorf("%w: %s cannot exchange energy with the %v channel", ErrBadValue, e.Name, ch)
			}
		}
		if err := same(e, e.Nodes[0], e.Nodes[1]); err != nil {
			return err
		}
	}
	for _, e := range c.Elements {
		if e.Kind != Mutual {
			continue
		}
		a, b := c.Elements[c.Index(e.Couples[0])], c.Elements[c.Index(e.Couples[1])]
		if c.ElementDomain(a) != c.ElementDomain(b) {
			return fmt.Errorf("%w: %s couples inductors of different domains", ErrBadValue, e.Name)
		}
	}
	return nil
}

func validateWave(name string, w Waveform) error {
	for _, a := range w.Args {
		if math.IsNaN(a) || math.IsInf(a, 0) {
			return fmt.Errorf("%w: %s has a non-finite argument", ErrBadWaveform, name)
		}
	}
	n := len(w.Args)
	switch w.Kind {
	case WaveNone:
		if n != 0 {
			return fmt.Errorf("%w: %s", ErrBadWaveform, name)
		}
	case WavePulse:
		if n < 2 || n > 7 {
			return fmt.Errorf("%w: %s PULSE needs 2–7 arguments, got %d", ErrBadWaveform, name, n)
		}
		for i := 2; i < n; i++ { // td tr tf pw per must not be negative
			if w.Args[i] < 0 {
				return fmt.Errorf("%w: %s PULSE argument %d is negative", ErrBadWaveform, name, i+1)
			}
		}
	case WaveSin:
		if n < 2 || n > 6 {
			return fmt.Errorf("%w: %s SIN needs 2–6 arguments, got %d", ErrBadWaveform, name, n)
		}
		if n > 2 && w.Args[2] < 0 || n > 3 && w.Args[3] < 0 {
			return fmt.Errorf("%w: %s SIN frequency and delay must not be negative", ErrBadWaveform, name)
		}
	case WaveExp:
		if n < 2 || n > 6 {
			return fmt.Errorf("%w: %s EXP needs 2–6 arguments, got %d", ErrBadWaveform, name, n)
		}
		for i := 2; i < n; i++ {
			if w.Args[i] < 0 {
				return fmt.Errorf("%w: %s EXP argument %d is negative", ErrBadWaveform, name, i+1)
			}
		}
	case WavePWL:
		if n < 2 || n%2 != 0 || n/2 > MaxPWL {
			return fmt.Errorf("%w: %s PWL needs time/value pairs", ErrBadWaveform, name)
		}
		for i := 2; i < n; i += 2 {
			if w.Args[i] < w.Args[i-2] {
				return fmt.Errorf("%w: %s PWL times must not decrease", ErrBadWaveform, name)
			}
		}
		if w.Args[0] < 0 {
			return fmt.Errorf("%w: %s PWL times must not be negative", ErrBadWaveform, name)
		}
	default:
		return fmt.Errorf("%w: %s unknown waveform", ErrBadWaveform, name)
	}
	return nil
}

// checkPaths reports nodes with no path to ground. Current sources never
// provide a path; capacitors do only when capIsPath is set.
func (c *Circuit) checkPaths(capIsPath bool) error {
	n := len(c.NodeNames)
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	find := func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for _, e := range c.Elements {
		if e.Kind == CurrentSource || e.Kind == Mutual || ((e.Kind == Capacitor || e.Kind == NonlinearCapacitor) && !capIsPath) {
			continue
		}
		// A current-output behavioral element is a path only when it is a
		// nonlinear resistor (controlled by its own voltage).
		if e.Kind == Behavioral && e.BOut == 'I' && !c.SelfControlled(e) {
			continue
		}
		// Controlled current outputs (F, G) are no path; voltage outputs are.
		if e.Kind == VCCS || e.Kind == CCCS {
			continue
		}
		// A gyrator's ports are controlled current sources: no DC path
		// (a free mass under a constant force has no steady state).
		if e.Kind == Gyrator {
			continue
		}
		// A bearing motor conducts a → b (resistive path); its shaft port
		// drives torque, no path.
		if e.Kind == BearingMotor {
			if a, b := find(e.Nodes[0]), find(e.Nodes[1]); a != b {
				parent[a] = b
			}
			continue
		}
		// A machine's windings are inductors (shorts at DC); its mechanical
		// port drives torque, no path.
		if e.Kind == Machine {
			w, _ := MachineWindings(e)
			for _, p := range w {
				if a, b := find(p[0]), find(p[1]); a != b {
					parent[a] = b
				}
			}
			continue
		}
		a, b := find(e.Nodes[0]), find(e.Nodes[1])
		if a != b {
			parent[a] = b
		}
		// A MOSFET's bulk reaches source and drain through its junction
		// diodes (gmin keeps them conducting); the gate is insulated and is
		// no DC path.
		if e.Kind == MOSFET {
			if a, b := find(e.Extra[1]), find(e.Nodes[0]); a != b {
				parent[a] = b
			}
		}
		// A BJT's base reaches collector and emitter through its junctions
		// (gmin keeps them conducting); the substrate has only capacitance.
		if e.Kind == BJT {
			if a, b := find(e.Extra[0]), find(e.Nodes[0]); a != b {
				parent[a] = b
			}
		}
	}
	g := find(0)
	var floating []string
	for i := 1; i < n; i++ {
		if find(i) != g {
			floating = append(floating, c.NodeNames[i])
		}
	}
	if len(floating) > 0 {
		sort.Strings(floating)
		if len(floating) > 10 {
			floating = append(floating[:10], "…")
		}
		kind := "DC "
		if capIsPath {
			kind = ""
		}
		return fmt.Errorf("%w (%spath): %s", ErrFloatingNode, kind, strings.Join(floating, ", "))
	}
	return nil
}

// validateMutual checks a K element: two distinct existing inductors, a
// coupling coefficient 0 < k ≤ 1 (k = 1 is the ideal, ngspice's limit),
// and at most one K per pair.
func (c *Circuit) validateMutual(e Element, seen map[[2]string]bool) error {
	if !(e.Value > 0 && e.Value <= 1) {
		return fmt.Errorf("%w: %s: coupling coefficient %g must be in (0, 1]", ErrBadValue, e.Name, e.Value)
	}
	a, b := strings.ToLower(e.Couples[0]), strings.ToLower(e.Couples[1])
	if a == b {
		return fmt.Errorf("%w: %s couples %s with itself", ErrBadValue, e.Name, e.Couples[0])
	}
	for _, n := range []string{a, b} {
		i := c.Index(n)
		if i < 0 || c.Elements[i].Kind != Inductor {
			return fmt.Errorf("%w: %s: %q is not a linear inductor", ErrBadValue, e.Name, n)
		}
	}
	if a > b {
		a, b = b, a
	}
	if seen[[2]string{a, b}] {
		return fmt.Errorf("%w: %s: %s and %s are already coupled", ErrBadValue, e.Name, a, b)
	}
	seen[[2]string{a, b}] = true
	return nil
}

// mosInstanceParams are the MOSFET instance parameters (SPICE3f3 User's
// Manual §2.4.10): channel length and width, drain and source diffusion
// areas and perimeters, squares of drain and source diffusion, and the
// parallel multiplier.
var mosInstanceParams = map[string]bool{"l": true, "w": true, "ad": true, "as": true, "pd": true,
	"ps": true, "nrd": true, "nrs": true, "m": true}

// IsInstanceParam reports whether name (lower case) is an instance
// parameter of elements of kind k (MOSFET: l, w, ad, as, pd, ps, nrd, nrs,
// m).
func IsInstanceParam(k Kind, name string) bool {
	return k == MOSFET && mosInstanceParams[name]
}

func validateMOSInstance(e Element) error {
	keys := make([]string, 0, len(e.Params))
	for k := range e.Params {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		v := e.Params[k]
		if !mosInstanceParams[k] {
			return fmt.Errorf("%w: %s: unknown MOSFET instance parameter %q", ErrBadValue, e.Name, k)
		}
		if math.IsNaN(v) || math.IsInf(v, 0) || v < 0 || ((k == "l" || k == "w" || k == "m") && v == 0) || v > 1e6 {
			return fmt.Errorf("%w: %s: %s = %g", ErrBadValue, e.Name, k, v)
		}
	}
	return nil
}

// validateControl checks the control of a controlled source: two distinct
// control nodes, or an existing voltage source.
func (c *Circuit) validateControl(e Element, voltage bool) error {
	if voltage {
		for _, n := range e.CtrlNodes {
			if n < 0 || n >= len(c.NodeNames) {
				return fmt.Errorf("%w: %s control node index %d out of range", ErrBadValue, e.Name, n)
			}
		}
		if e.CtrlNodes[0] == e.CtrlNodes[1] {
			return fmt.Errorf("%w: %s is controlled by a voltage that is always 0", ErrBadValue, e.Name)
		}
		return nil
	}
	i := c.Index(e.CtrlSource)
	if i < 0 || c.Elements[i].Kind != VoltageSource {
		return fmt.Errorf("%w: %s is controlled by %q, which is not a voltage source", ErrBadValue, e.Name, e.CtrlSource)
	}
	return nil
}

// validateStateLaw checks a nonlinear capacitor or inductor. Its curve
// (charge against voltage, or flux against current) must rise strictly and
// pass through the origin: a positive incremental capacitance or
// inductance everywhere keeps the element passive, with stored energy
// ∫ v dq ≥ 0 measured from the uncharged state. Without a curve, Value is
// the linear C or L and must be > 0. The scale s(t) must stay > 0.
func validateStateLaw(e Element) error {
	if e.Wave.Kind != WaveNone {
		return fmt.Errorf("%w: %s is not a source", ErrBadWaveform, e.Name)
	}
	n := len(e.Curve)
	if n == 0 {
		if !(e.Value > 0) {
			return fmt.Errorf("%w: %s %s must be > 0, got %g", ErrBadValue, e.Name, e.Kind, e.Value)
		}
	} else {
		if e.Value != 0 {
			return fmt.Errorf("%w: %s has both a value and a curve", ErrBadValue, e.Name)
		}
		if n < 4 || n%2 != 0 || n > 2*MaxPWL {
			return fmt.Errorf("%w: %s needs 2 to %d (x, y) points", ErrBadValue, e.Name, MaxPWL)
		}
		ymax := 0.0
		for k, v := range e.Curve {
			if math.IsNaN(v) || math.IsInf(v, 0) {
				return fmt.Errorf("%w: %s curve has a non-finite value", ErrBadValue, e.Name)
			}
			if k >= 2 && !(v > e.Curve[k-2]) {
				return fmt.Errorf("%w: %s curve must strictly increase in x and y (positive incremental %s)",
					ErrBadValue, e.Name, strings.TrimPrefix(e.Kind.String(), "nonlinear "))
			}
			if k%2 == 1 {
				ymax = math.Max(ymax, math.Abs(v))
			}
		}
		if y0 := curveAt(e.Curve, 0); math.Abs(y0) > 1e-12*ymax {
			return fmt.Errorf("%w: %s curve must pass through the origin, has %g at 0", ErrBadValue, e.Name, y0)
		}
	}
	if e.Scale.Kind == WaveNone {
		return nil
	}
	if err := validateWave(e.Name, e.Scale); err != nil {
		return err
	}
	if e.Scale.Kind == WavePWL {
		// A jump in s would move charge or flux in zero time: an impulse
		// of exchanged power the ledger cannot book.
		for i := 2; i < len(e.Scale.Args); i += 2 {
			if !(e.Scale.Args[i] > e.Scale.Args[i-2]) {
				return fmt.Errorf("%w: %s scale PWL times must strictly increase", ErrBadWaveform, e.Name)
			}
		}
	}
	if lo, ok := scaleLowerBound(e.Scale); !ok || !(lo > 0) {
		return fmt.Errorf("%w: %s scale must stay > 0", ErrBadWaveform, e.Name)
	}
	return nil
}

// curveAt evaluates a piecewise-linear x, y curve at x, extending the end
// segments.
func curveAt(c []float64, x float64) float64 {
	n := len(c) / 2
	k := sort.Search(n, func(i int) bool { return c[2*i] >= x })
	switch {
	case k <= 0:
		k = 1
	case k >= n:
		k = n - 1
	}
	x0, y0, x1, y1 := c[2*(k-1)], c[2*(k-1)+1], c[2*k], c[2*k+1]
	return y0 + (y1-y0)*(x-x0)/(x1-x0)
}

// scaleLowerBound returns a lower bound of a waveform over all time that
// needs no analysis defaults; ok is false when there is none (a growing
// SIN). EXP gets min(v1, v2) here; the device checks the exact bound, which
// depends on the defaulted time constants.
func scaleLowerBound(w Waveform) (float64, bool) {
	a := func(i int) float64 {
		if i < len(w.Args) {
			return w.Args[i]
		}
		return 0
	}
	switch w.Kind {
	case WavePulse, WaveExp:
		return math.Min(a(0), a(1)), true
	case WaveSin:
		if a(4) < 0 {
			return 0, false
		}
		return a(0) - math.Abs(a(1)), true
	case WavePWL:
		lo := math.Inf(1)
		for i := 1; i < len(w.Args); i += 2 {
			lo = math.Min(lo, w.Args[i])
		}
		return lo, true
	}
	return 0, false
}

func (c *Circuit) validateBehavioral(e Element) error {
	if e.BOut != 'I' && e.BOut != 'V' {
		return fmt.Errorf("%w: %s output must be I or V", ErrBadValue, e.Name)
	}
	switch e.BCtrl {
	case 'V':
		for _, n := range e.CtrlNodes {
			if n < 0 || n >= len(c.NodeNames) {
				return fmt.Errorf("%w: %s control node index %d out of range", ErrBadValue, e.Name, n)
			}
		}
		if e.CtrlNodes[0] == e.CtrlNodes[1] {
			return fmt.Errorf("%w: %s is controlled by a voltage that is always 0", ErrBadValue, e.Name)
		}
	case 'I':
		i := c.Index(e.CtrlSource)
		if i < 0 || c.Elements[i].Kind != VoltageSource {
			return fmt.Errorf("%w: %s is controlled by %q, which is not a voltage source", ErrBadValue, e.Name, e.CtrlSource)
		}
	default:
		return fmt.Errorf("%w: %s control must be V or I", ErrBadValue, e.Name)
	}
	n := len(e.Curve)
	if n < 4 || n%2 != 0 || n > 2*MaxPWL {
		return fmt.Errorf("%w: %s needs 2 to %d (x, y) points", ErrBadValue, e.Name, MaxPWL)
	}
	for k, v := range e.Curve {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: %s curve has a non-finite value", ErrBadValue, e.Name)
		}
		if k >= 2 && k%2 == 0 && !(v > e.Curve[k-2]) {
			return fmt.Errorf("%w: %s curve x values must strictly increase", ErrBadValue, e.Name)
		}
	}
	return nil
}

// NodeIndex returns the index of a named node (case-insensitive; "0" and
// "gnd" are ground), or −1.
func (c *Circuit) NodeIndex(name string) int {
	key := strings.ToLower(strings.TrimSpace(name))
	if key == "gnd" {
		key = "0"
	}
	for i, n := range c.NodeNames {
		if n == key {
			return i
		}
	}
	return -1
}
