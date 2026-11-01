// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/netlist"
)

// Builder assembles a circuit programmatically. It is the API the genetic
// algorithm uses, so no netlist text is parsed per candidate.
type Builder struct {
	title  string
	index  map[string]int
	names  []string
	elems  []Element
	models map[string]Model
	domain map[int]exchange.ChannelID // non-electrical nodes
}

// NewBuilder returns a builder with only the ground node defined.
func NewBuilder(title string) *Builder {
	return &Builder{
		title: title,
		index: map[string]int{"0": 0, "gnd": 0},
		names: []string{"0"},
	}
}

// Node returns the index of a node, creating it if needed. Names are
// case-insensitive; "0" and "gnd" are ground.
func (b *Builder) Node(name string) int {
	key := strings.ToLower(strings.TrimSpace(name))
	if i, ok := b.index[key]; ok {
		return i
	}
	i := len(b.names)
	b.index[key] = i
	b.names = append(b.names, key)
	return i
}

// NonlinearVC adds a voltage-controlled nonlinear resistor (N-type, such as
// a tunnel diode): current from a to b is i = f(v), v = V(a) − V(b), with f
// the piecewise-linear curve v1, i1, v2, i2, … (ngspice
// "Bname a b I=pwl(V(a,b), ...)").
func (b *Builder) NonlinearVC(name, a, c string, curve []float64) *Builder {
	e := Element{Name: name, Kind: Behavioral, BOut: 'I', BCtrl: 'V', Curve: append([]float64(nil), curve...)}
	e.CtrlNodes = [2]int{b.Node(a), b.Node(c)}
	return b.Add(e, a, c)
}

// NonlinearCC adds a current-controlled nonlinear resistor (S-type, such as
// an arc or a neon lamp): v = g(i) for the current i from a to b, with g the
// curve i1, v1, i2, v2, … . The current is sensed by a 0 V source named
// "Vs_<name>" in series, so the netlist form is plain ngspice:
// "Vs_<name> a <name>_s 0" and "Bname <name>_s b V=pwl(I(Vs_<name>), ...)".
func (b *Builder) NonlinearCC(name, a, c string, curve []float64) *Builder {
	sense := "Vs_" + name
	mid := strings.ToLower(name) + "_s"
	b.Add(Element{Name: sense, Kind: VoltageSource, HasDC: true}, a, mid)
	e := Element{Name: name, Kind: Behavioral, BOut: 'V', BCtrl: 'I', CtrlSource: strings.ToLower(sense),
		Curve: append([]float64(nil), curve...)}
	return b.Add(e, mid, c)
}

// Add appends a fully specified element; n1 and n2 are node names.
func (b *Builder) Add(e Element, n1, n2 string) *Builder {
	e.Nodes = [2]int{b.Node(n1), b.Node(n2)}
	b.elems = append(b.elems, e)
	return b
}

// Mutual adds a mutual inductance (SPICE K) between inductors l1 and l2
// with coupling coefficient k.
func (b *Builder) Mutual(name, l1, l2 string, k float64) *Builder {
	return b.Add(Element{Name: name, Kind: Mutual, Value: k, Couples: [2]string{l1, l2}}, "0", "0")
}

// AddMulti appends an element with more than two terminals: Nodes gets the
// first two names, Extra the rest (for a BJT: collector, emitter, then
// base and optional substrate).
func (b *Builder) AddMulti(e Element, n1, n2 string, extra ...string) *Builder {
	e.Nodes = [2]int{b.Node(n1), b.Node(n2)}
	e.Extra = make([]int, len(extra))
	for i, n := range extra {
		e.Extra[i] = b.Node(n)
	}
	b.elems = append(b.elems, e)
	return b
}

// BJT adds a bipolar transistor with terminals in SPICE order (collector,
// base, emitter) and a model defined with Model (type "npn" or "pnp").
func (b *Builder) BJT(name, collector, base, emitter, model string) *Builder {
	return b.AddMulti(Element{Name: name, Kind: BJT, Model: strings.ToLower(model)}, collector, emitter, base)
}

// MOSFET adds a MOS transistor with terminals in SPICE order (drain,
// gate, source, bulk), a model of type "nmos" or "pmos", and instance
// parameters (l, w, ad, as, pd, ps, nrd, nrs, m; nil for the defaults).
func (b *Builder) MOSFET(name, drain, gate, source, bulk, model string, params map[string]float64) *Builder {
	var p map[string]float64
	if params != nil {
		p = make(map[string]float64, len(params))
		for k, v := range params {
			p[strings.ToLower(k)] = v
		}
	}
	return b.AddMulti(Element{Name: name, Kind: MOSFET, Model: strings.ToLower(model), Params: p}, drain, source, gate, bulk)
}

// VCVS adds E: V(p) − V(n) = gain·(V(cp) − V(cn)).
func (b *Builder) VCVS(name, p, n, cp, cn string, gain float64) *Builder {
	return b.Add(Element{Name: name, Kind: VCVS, Value: gain, CtrlNodes: [2]int{b.Node(cp), b.Node(cn)}}, p, n)
}

// VCCS adds G: a current gm·(V(cp) − V(cn)) from p through the source to n.
func (b *Builder) VCCS(name, p, n, cp, cn string, gm float64) *Builder {
	return b.Add(Element{Name: name, Kind: VCCS, Value: gm, CtrlNodes: [2]int{b.Node(cp), b.Node(cn)}}, p, n)
}

// CCCS adds F: a current gain·I(vsrc) from p through the source to n.
func (b *Builder) CCCS(name, p, n, vsrc string, gain float64) *Builder {
	return b.Add(Element{Name: name, Kind: CCCS, Value: gain, CtrlSource: strings.ToLower(vsrc)}, p, n)
}

// CCVS adds H: V(p) − V(n) = r·I(vsrc).
func (b *Builder) CCVS(name, p, n, vsrc string, r float64) *Builder {
	return b.Add(Element{Name: name, Kind: CCVS, Value: r, CtrlSource: strings.ToLower(vsrc)}, p, n)
}

// Model defines a .model card. Parameter names are case-insensitive.
func (b *Builder) Model(name, typ string, params map[string]float64) *Builder {
	if b.models == nil {
		b.models = map[string]Model{}
	}
	p := make(map[string]float64, len(params))
	for k, v := range params {
		p[strings.ToLower(k)] = v
	}
	key := strings.ToLower(name)
	b.models[key] = Model{Name: key, Type: strings.ToLower(typ), Params: p}
	return b
}

// Diode adds a diode from anode to cathode using the named model.
func (b *Builder) Diode(name, anode, cathode, model string) *Builder {
	return b.Add(Element{Name: name, Kind: Diode, Model: strings.ToLower(model)}, anode, cathode)
}

// Resistor adds a resistor of r ohms between nodes a and c.
func (b *Builder) Resistor(name, a, c string, ohms float64) *Builder {
	return b.Add(Element{Name: name, Kind: Resistor, Value: ohms}, a, c)
}

// Capacitor adds a capacitor of farads between nodes a and c.
func (b *Builder) Capacitor(name, a, c string, farads float64) *Builder {
	return b.Add(Element{Name: name, Kind: Capacitor, Value: farads}, a, c)
}

// CapacitorIC adds a capacitor with an initial voltage V(a) − V(c), used by UIC.
func (b *Builder) CapacitorIC(name, a, c string, farads, ic float64) *Builder {
	return b.Add(Element{Name: name, Kind: Capacitor, Value: farads, IC: ic, HasIC: true}, a, c)
}

// Inductor adds an inductor of henries between nodes a and c.
func (b *Builder) Inductor(name, a, c string, henries float64) *Builder {
	return b.Add(Element{Name: name, Kind: Inductor, Value: henries}, a, c)
}

// InductorIC adds an inductor with an initial current from a to c, used by UIC.
func (b *Builder) InductorIC(name, a, c string, henries, ic float64) *Builder {
	return b.Add(Element{Name: name, Kind: Inductor, Value: henries, IC: ic, HasIC: true}, a, c)
}

// NonlinearCapacitor adds a capacitor defined by its charge
// q(v, t) = s(t)·Q(v), with Q the piecewise-linear curve v1, q1, v2, q2, …
// (strictly rising, through the origin) and s the scale waveform (Kind
// WaveNone for none). The current is i = dq/dt = s·Q′(v)·dv/dt + ṡ·Q(v).
func (b *Builder) NonlinearCapacitor(name, a, c string, charge []float64, scale Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: NonlinearCapacitor, Curve: append([]float64(nil), charge...), Scale: scale}, a, c)
}

// VaryingCapacitor adds a linear capacitor whose capacitance varies in time,
// C(t) = farads·s(t): i = C·dv/dt + v·dC/dt.
func (b *Builder) VaryingCapacitor(name, a, c string, farads float64, scale Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: NonlinearCapacitor, Value: farads, Scale: scale}, a, c)
}

// NonlinearInductor adds an inductor defined by its flux linkage
// λ(i, t) = s(t)·Λ(i), with Λ the curve i1, λ1, i2, λ2, … (strictly rising,
// through the origin): v = dλ/dt = s·Λ′(i)·di/dt + ṡ·Λ(i).
func (b *Builder) NonlinearInductor(name, a, c string, flux []float64, scale Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: NonlinearInductor, Curve: append([]float64(nil), flux...), Scale: scale}, a, c)
}

// VaryingInductor adds a linear inductor whose inductance varies in time,
// L(t) = henries·s(t): v = L·di/dt + i·dL/dt.
func (b *Builder) VaryingInductor(name, a, c string, henries float64, scale Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: NonlinearInductor, Value: henries, Scale: scale}, a, c)
}

// VoltageSource adds a DC voltage source: V(plus) − V(minus) = volts.
func (b *Builder) VoltageSource(name, plus, minus string, volts float64) *Builder {
	return b.Add(Element{Name: name, Kind: VoltageSource, Value: volts, HasDC: true}, plus, minus)
}

// VoltageWave adds a voltage source with a transient waveform.
func (b *Builder) VoltageWave(name, plus, minus string, w Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: VoltageSource, Wave: w}, plus, minus)
}

// CurrentSource adds a DC current source. Following SPICE, amps flow from
// the plus node through the source into the minus node.
func (b *Builder) CurrentSource(name, plus, minus string, amps float64) *Builder {
	return b.Add(Element{Name: name, Kind: CurrentSource, Value: amps, HasDC: true}, plus, minus)
}

// CurrentWave adds a current source with a transient waveform.
func (b *Builder) CurrentWave(name, plus, minus string, w Waveform) *Builder {
	return b.Add(Element{Name: name, Kind: CurrentSource, Wave: w}, plus, minus)
}

func (b *Builder) circuit() *Circuit {
	models := make(map[string]Model, len(b.models))
	for k, v := range b.models {
		models[k] = v
	}
	c := &Circuit{
		Title:     b.title,
		NodeNames: append([]string(nil), b.names...),
		Elements:  append([]Element(nil), b.elems...),
		Models:    models,
	}
	if len(b.domain) > 0 {
		c.NodeDomains = make([]exchange.ChannelID, len(b.names))
		for n, d := range b.domain {
			c.NodeDomains[n] = d
		}
	}
	return c
}

// Domain puts nodes in an energy domain (exchange.Mechanical: across =
// velocity, through = force; exchange.Acoustic: pressure and volume
// velocity; …). Nodes not placed are electrical; ground belongs to every
// domain.
func (b *Builder) Domain(d exchange.ChannelID, nodes ...string) *Builder {
	if b.domain == nil {
		b.domain = map[int]exchange.ChannelID{}
	}
	for _, n := range nodes {
		if i := b.Node(n); i != 0 {
			b.domain[i] = d
		}
	}
	return b
}

// Gyrator adds an ideal lossless transducer between port 1 (p1, n1) and
// port 2 (p2, n2): i1 = g·V2, i2 = −g·V1 (see Kind Gyrator).
func (b *Builder) Gyrator(name, p1, n1, p2, n2 string, g float64) *Builder {
	return b.AddMulti(Element{Name: name, Kind: Gyrator, Value: g}, p1, n1, p2, n2)
}

// ResistorTo adds a resistor whose dissipated power leaves through channel
// ch instead of as heat: an acoustic radiation resistance (Acoustic), a
// mechanical load that does work outside the circuit (Mechanical), …
func (b *Builder) ResistorTo(name, a, c string, ohms float64, ch exchange.ChannelID) *Builder {
	return b.Add(Element{Name: name, Kind: Resistor, Value: ohms, Channel: ch}, a, c)
}

// Build validates (with DC paths required) and returns the circuit.
func (b *Builder) Build() (*Circuit, error) {
	c := b.circuit()
	if err := c.Validate(); err != nil {
		return nil, err
	}
	return c, nil
}

// BuildTransient validates for a transient analysis without an operating
// point (capacitors count as paths to ground) and returns the circuit.
func (b *Builder) BuildTransient() (*Circuit, error) {
	c := b.circuit()
	if err := c.ValidateTransient(); err != nil {
		return nil, err
	}
	return c, nil
}

var waveKinds = map[string]WaveKind{"pulse": WavePulse, "sin": WaveSin, "exp": WaveExp, "pwl": WavePWL}

// FromDeck converts a parsed netlist into a validated circuit. When the deck
// requests a transient analysis with UIC, capacitors count as paths.
func FromDeck(d *netlist.Deck) (*Circuit, error) {
	b := NewBuilder(d.Title)
	for _, m := range d.Models {
		if _, dup := b.models[m.Name]; dup {
			return nil, fmt.Errorf("%w: model %s defined twice (line %d)", ErrBadModel, m.Name, m.Line)
		}
		b.Model(m.Name, m.Type, m.Params)
	}
	if len(d.Domains) > 0 {
		names := make([]string, 0, len(d.Domains))
		for n := range d.Domains {
			names = append(names, n)
		}
		sort.Strings(names)
		for _, n := range names {
			ch, ok := exchange.ParseChannel(d.Domains[n])
			if !ok {
				return nil, fmt.Errorf("%w: .domain %s: unknown domain", ErrBadValue, d.Domains[n])
			}
			b.Domain(ch, n)
		}
	}
	for _, e := range d.Elements {
		if len(e.Nodes) != 2 && !(e.Kind == 'Q' && (len(e.Nodes) == 3 || len(e.Nodes) == 4)) && !((e.Kind == 'M' || e.Kind == 'Y' || e.Kind == 'P') && len(e.Nodes) == 4) && !(e.Kind == 'K' && len(e.Nodes) == 0) && !(e.Kind == 'N' && len(e.Nodes) >= 4 && len(e.Nodes)%2 == 0) {
			return nil, fmt.Errorf("%w: %s on line %d", ErrBadValue, e.Name, e.Line)
		}
		el := Element{
			Name: e.Name, Value: e.Value, HasDC: e.HasDC,
			ACMag: e.ACMag, ACPhase: e.ACPhase, IC: e.IC, HasIC: e.HasIC,
			Model: e.Model, Area: e.Area,
		}
		switch e.Kind {
		case 'R':
			el.Kind = Resistor
			if e.Ext != "" {
				ch, ok := exchange.ParseChannel(e.Ext)
				if !ok {
					return nil, fmt.Errorf("%w: %s on line %d: unknown channel %q", ErrBadValue, e.Name, e.Line, e.Ext)
				}
				el.Channel = ch
			}
		case 'C', 'L':
			el.Kind = map[byte]Kind{'C': Capacitor, 'L': Inductor}[e.Kind]
			if len(e.Curve) > 0 || e.ScaleKind != "" || e.Ext != "" {
				el.Kind = map[byte]Kind{'C': NonlinearCapacitor, 'L': NonlinearInductor}[e.Kind]
				el.Curve = append([]float64(nil), e.Curve...)
				if e.ScaleKind != "" {
					wk, ok := waveKinds[e.ScaleKind]
					if !ok {
						return nil, fmt.Errorf("%w: %s on line %d: %s", ErrBadWaveform, e.Name, e.Line, e.ScaleKind)
					}
					el.Scale = Waveform{Kind: wk, Args: append([]float64(nil), e.ScaleArgs...)}
				}
				if e.Ext != "" {
					ch, ok := exchange.ParseChannel(e.Ext)
					if !ok {
						return nil, fmt.Errorf("%w: %s on line %d: unknown channel %q", ErrBadValue, e.Name, e.Line, e.Ext)
					}
					el.Channel = ch
				}
			}
		case 'V':
			el.Kind = VoltageSource
		case 'I':
			el.Kind = CurrentSource
		case 'D':
			el.Kind = Diode
		case 'B':
			el.Kind = Behavioral
			el.BOut, el.BCtrl, el.CtrlSource = e.BOut, e.BCtrl, e.CtrlSource
			el.Curve = append([]float64(nil), e.Curve...)
			if len(e.CtrlNodes) == 2 {
				el.CtrlNodes = [2]int{b.Node(e.CtrlNodes[0]), b.Node(e.CtrlNodes[1])}
			}
		case 'E', 'F', 'G', 'H':
			el.Kind = Kind(e.Kind)
			el.CtrlSource = e.CtrlSource
			if len(e.CtrlNodes) == 2 {
				el.CtrlNodes = [2]int{b.Node(e.CtrlNodes[0]), b.Node(e.CtrlNodes[1])}
			}
		case 'M':
			el.Kind = MOSFET
			if len(e.Params) > 0 {
				el.Params = make(map[string]float64, len(e.Params))
				for k, v := range e.Params {
					el.Params[k] = v
				}
			}
			b.AddMulti(el, e.Nodes[0], e.Nodes[2], e.Nodes[1], e.Nodes[3])
			continue
		case 'Y':
			el.Kind = Gyrator
			b.AddMulti(el, e.Nodes[0], e.Nodes[1], e.Nodes[2], e.Nodes[3])
			continue
		case 'P':
			el.Kind = BearingMotor
			if len(e.Params) > 0 {
				el.Params = make(map[string]float64, len(e.Params))
				for k, v := range e.Params {
					el.Params[k] = v
				}
			}
			b.AddMulti(el, e.Nodes[0], e.Nodes[1], e.Nodes[2], e.Nodes[3])
			continue
		case 'N':
			el.Kind = Machine
			b.AddMulti(el, e.Nodes[0], e.Nodes[1], e.Nodes[2:]...)
			continue
		case 'K':
			el.Kind = Mutual
			if len(e.Couples) == 2 {
				el.Couples = [2]string{e.Couples[0], e.Couples[1]}
			}
			b.Add(el, "0", "0")
			continue
		case 'Q':
			el.Kind = BJT
			extra := append([]string{e.Nodes[1]}, e.Nodes[3:]...)
			b.AddMulti(el, e.Nodes[0], e.Nodes[2], extra...)
			continue
		default:
			return nil, fmt.Errorf("%w: %s on line %d", ErrUnknownKind, e.Name, e.Line)
		}
		if e.WaveKind != "" {
			wk, ok := waveKinds[e.WaveKind]
			if !ok {
				return nil, fmt.Errorf("%w: %s on line %d: %s", ErrBadWaveform, e.Name, e.Line, e.WaveKind)
			}
			el.Wave = Waveform{Kind: wk, Args: append([]float64(nil), e.WaveArgs...)}
		}
		b.Add(el, e.Nodes[0], e.Nodes[1])
	}
	for _, a := range d.Analyses {
		if a.Kind == "tran" && a.UIC {
			return b.BuildTransient()
		}
	}
	return b.Build()
}
