// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
)

// Errors.
var (
	// ErrNonlinear is returned for a circuit with nonlinear elements
	// (diodes, transistors, B elements): M9a covers linear DC and AC only
	// (CLAUDE.md §3a rule 3).
	ErrNonlinear = errors.New("kron: nonlinear element (linear circuits only)")
	// ErrDimension is a shape mismatch between tensors (a bug).
	ErrDimension = errors.New("kron: tensor dimensions do not match")
	// ErrSingular is returned when the transformed network equations have
	// no unique solution (a floating sub-network, a loop of voltage
	// sources, a cut set of current sources, …).
	ErrSingular = errors.New("kron: singular network equations")
	// ErrUnsupported is a construction this formulation cannot express.
	ErrUnsupported = errors.New("kron: unsupported")
)

// branchKind is what a primitive branch is.
type branchKind uint8

const (
	bR branchKind = iota
	bC
	bL
	bV
	bI
	bVCVS
	bCCCS
	bVCCS
	bCCVS
	bControl   // an open coil across a controlling node pair (no admittance)
	bRadiation // a radiation-resistance edge (§3a.1 item 5; chords.go)
	bGyr1      // a transducer's port 1: i1 = g·E2 (ctrl = port 2)
	bGyr2      // a transducer's port 2: i2 = −g·E1 (ctrl = port 1)
)

// branch is one primitive coil. Its current flows from → to through it
// and its voltage is E = V(from) − V(to), the SPICE sign of the element it
// comes from (n+ → n−).
type branch struct {
	kind     branchKind
	elem     int // circuit element index, −1 for a control coil
	from, to int
	value    float64 // R [Ω], C [F], L [H], gain, or a source's DC value
	ac       complex128
	ctrl     int // controlling branch (control coil, or the sensed voltage source)
}

// Network is a linear circuit as Kron's primitive network plus its
// interconnection: the branches on their own (the primitive frame), and
// the graph that connects them.
type Network struct {
	c        *circuit.Circuit
	branches []branch
	g        *graph
	nodes    int // junctions: the circuit's nodes plus internal ones (radiation edges)

	elemBranch []int           // circuit element → its branch, −1 for a K element
	devs       []device.Device // the elements' device models: their exterior ports
	chords     []chord
	radiation  []radiationEdge
}

// Counts are the frame dimensions of a network: b branches, N junctions
// (nodes including the datum), S sub-networks, m = b − N + S meshes and
// p = N − S junction pairs (Kron 1939 p. 75, eqs 3.1–3.2; Tensors for
// Circuits p. 45, eqs 6.12–6.13).
type Counts struct{ Branches, Junctions, SubNetworks, Meshes, JunctionPairs int }

// Counts returns the frame dimensions.
func (n *Network) Counts() Counts {
	b, N, S := len(n.branches), n.g.n, n.g.subNetworks()
	return Counts{Branches: b, Junctions: N, SubNetworks: S, Meshes: b - N + S, JunctionPairs: N - S}
}

// Build makes the primitive network of a linear circuit. Every element is
// one primitive branch; each voltage-controlled source also gets a control
// coil, a branch of zero admittance across its controlling node pair, so
// the controlling voltage is a branch voltage of the primitive frame (the
// coupling is then a mutual term between two primitive axes).
func Build(c *circuit.Circuit) (*Network, error) {
	if err := c.Validate(); err != nil {
		return nil, err
	}
	devs, err := device.FromCircuit(c)
	if err != nil {
		return nil, err
	}
	nw := &Network{c: c, devs: devs, nodes: c.NumNodes(), elemBranch: make([]int, len(c.Elements))}
	vsrc := map[string]int{} // voltage-source name → branch, for F and H
	ind := map[string]int{}  // inductor name → branch, for K
	var pending, mutual []int
	for ei, e := range c.Elements {
		nw.elemBranch[ei] = len(nw.branches)
		br := branch{elem: ei, from: e.Nodes[0], to: e.Nodes[1], value: e.Value, ctrl: -1,
			ac: cmplx.Rect(e.ACMag, e.ACPhase*math.Pi/180)}
		switch e.Kind {
		case circuit.Resistor:
			br.kind = bR
		case circuit.Capacitor:
			br.kind = bC
		case circuit.Inductor:
			br.kind = bL
			ind[strings.ToLower(e.Name)] = len(nw.branches)
		case circuit.Gyrator:
			// Two primitive axes coupled by an antisymmetric mutual
			// admittance, y12 = g, y21 = −g: a non-reciprocal primitive
			// network, as Kron treats active and electromechanical
			// elements.
			k1 := len(nw.branches)
			br.kind, br.ctrl = bGyr1, k1+1
			nw.branches = append(nw.branches, br,
				branch{kind: bGyr2, elem: ei, from: e.Extra[0], to: e.Extra[1], value: e.Value, ctrl: k1})
			continue
		case circuit.Mutual:
			// Not a branch: a chord between the two inductors' axes.
			nw.elemBranch[ei] = -1
			mutual = append(mutual, ei)
			continue
		case circuit.VoltageSource:
			br.kind, br.value = bV, dcValue(e)
			vsrc[strings.ToLower(e.Name)] = len(nw.branches)
		case circuit.CurrentSource:
			br.kind, br.value = bI, dcValue(e)
		case circuit.VCVS, circuit.VCCS, circuit.CCCS, circuit.CCVS:
			br.kind = map[circuit.Kind]branchKind{circuit.VCVS: bVCVS, circuit.VCCS: bVCCS,
				circuit.CCCS: bCCCS, circuit.CCVS: bCCVS}[e.Kind]
			if e.Kind == circuit.VCVS || e.Kind == circuit.VCCS {
				// Controlled by a node-pair voltage: add its control coil.
				br.ctrl = len(nw.branches) + 1
				nw.branches = append(nw.branches, br,
					branch{kind: bControl, elem: -1, from: e.CtrlNodes[0], to: e.CtrlNodes[1], ctrl: -1})
				continue
			}
			pending = append(pending, len(nw.branches))
		default:
			return nil, fmt.Errorf("%w: %s (%s)", ErrNonlinear, e.Name, e.Kind)
		}
		nw.branches = append(nw.branches, br)
	}
	for _, k := range pending {
		e := c.Elements[nw.branches[k].elem]
		vk, ok := vsrc[strings.ToLower(e.CtrlSource)]
		if !ok {
			return nil, fmt.Errorf("%w: %s: controlling source %s is not a voltage source", ErrUnsupported, e.Name, e.CtrlSource)
		}
		nw.branches[k].ctrl = vk
	}
	for _, ei := range mutual {
		e := c.Elements[ei]
		a, b := ind[strings.ToLower(e.Couples[0])], ind[strings.ToLower(e.Couples[1])]
		m := e.Value * math.Sqrt(nw.branches[a].value*nw.branches[b].value)
		// Mutual inductance is the classic chord: z_ab = z_ba = jωM, the
		// flux of one coil linking the other (CLAUDE.md §3c, inductor row).
		nw.chords = append(nw.chords, chord{a: a, b: b, z: func(w float64) complex128 { return complex(0, w*m) }, reciprocal: true})
	}
	nw.rebuildGraph()
	return nw, nil
}

// rebuildGraph rebuilds the interconnection after branches were added.
func (n *Network) rebuildGraph() {
	from := make([]int, len(n.branches))
	to := make([]int, len(n.branches))
	for k, b := range n.branches {
		from[k], to[k] = b.from, b.to
	}
	n.g = newGraph(n.nodes, from, to)
}

// dcValue is a source's operating-point value: its DC value, or its
// waveform at t = 0 when no DC value was given (SPICE3; as analysis.OP).
func dcValue(e circuit.Element) float64 {
	w := e.Wave
	if w.Kind == circuit.WaveNone || e.HasDC {
		return e.Value
	}
	arg := func(i int) float64 {
		if i < len(w.Args) {
			return w.Args[i]
		}
		return 0
	}
	switch w.Kind {
	case circuit.WavePulse, circuit.WaveExp:
		return arg(0)
	case circuit.WaveSin:
		return arg(0) + arg(1)*math.Sin(arg(5)*math.Pi/180)
	case circuit.WavePWL:
		p := w.Args
		if len(p) < 2 || 0 <= p[0] {
			return arg(1)
		}
		for i := 2; i+1 < len(p); i += 2 {
			if p[i] > 0 {
				t0, v0, t1, v1 := p[i-2], p[i-1], p[i], p[i+1]
				return v0 + (v1-v0)*(0-t0)/(t1-t0)
			}
		}
		return p[len(p)-1]
	}
	return e.Value
}

// term is one coefficient of a branch relation.
type term struct {
	k int // branch
	c complex128
}

// relation is a branch's constitutive equation, linear in the branch
// voltages E and currents i: Σ E-terms·E + Σ I-terms·i = s. Every element
// has one; together they are the primitive network's equations, in Kron's
// general form E + e = z·(i + I) (Kron 1939 p. 417; eq 19.3 p. 480)
// written so that zero impedances (voltage sources) and zero admittances
// (current sources, control coils) need no infinite entries.
type relation struct {
	e, i []term
	s    complex128
}

// own returns the relation's coefficients on its own branch's E and i.
func (r relation) own(k int) (ce, ci complex128) {
	for _, t := range r.e {
		if t.k == k {
			ce += t.c
		}
	}
	for _, t := range r.i {
		if t.k == k {
			ci += t.c
		}
	}
	return ce, ci
}

// relations returns every branch's relation at angular frequency omega;
// dc selects the operating point (capacitors open, inductors shorted,
// sources at their DC values) instead of the small-signal phasors.
func (n *Network) relations(omega float64, dc bool) []relation {
	out := make([]relation, len(n.branches))
	jw := complex(0, omega)
	for k, b := range n.branches {
		var r relation
		src := b.ac
		if dc {
			src = complex(b.value, 0)
		}
		g := complex(b.value, 0)
		switch b.kind {
		case bR: // E − R·i = 0 (Ohm's law, CLAUDE.md §3c)
			r = relation{e: []term{{k, 1}}, i: []term{{k, -g}}}
		case bC: // i − jωC·E = 0; open at DC
			if dc {
				r = relation{i: []term{{k, 1}}}
			} else {
				r = relation{e: []term{{k, -jw * g}}, i: []term{{k, 1}}}
			}
		case bL: // E − jωL·i = 0; shorted at DC
			if dc {
				r = relation{e: []term{{k, 1}}}
			} else {
				r = relation{e: []term{{k, 1}}, i: []term{{k, -jw * g}}}
			}
		case bV: // E = value
			r = relation{e: []term{{k, 1}}, s: src}
		case bI: // i = value
			r = relation{i: []term{{k, 1}}, s: src}
		case bVCVS: // E − µ·E_c = 0
			r = relation{e: []term{{k, 1}, {b.ctrl, -g}}}
		case bVCCS: // i − g·E_c = 0
			r = relation{e: []term{{b.ctrl, -g}}, i: []term{{k, 1}}}
		case bCCCS: // i − β·i_c = 0
			r = relation{i: []term{{k, 1}, {b.ctrl, -g}}}
		case bCCVS: // E − r·i_c = 0
			r = relation{e: []term{{k, 1}}, i: []term{{b.ctrl, -g}}}
		case bControl: // i = 0
			r = relation{i: []term{{k, 1}}}
		case bGyr1: // i1 − g·E2 = 0
			r = relation{e: []term{{b.ctrl, -g}}, i: []term{{k, 1}}}
		case bGyr2: // i2 + g·E1 = 0
			r = relation{e: []term{{b.ctrl, g}}, i: []term{{k, 1}}}
		case bRadiation: // E − R_rad(ω)·i = 0
			r = relation{e: []term{{k, 1}}, i: []term{{k, -complex(n.radiationR(k, omega), 0)}}}
		}
		out[k] = r
	}
	n.applyChords(out, omega)
	return out
}
