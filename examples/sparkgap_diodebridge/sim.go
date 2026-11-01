// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
)

// A spark gap is a switch with memory (it strikes at a threshold and stays
// on while its arc carries current), which kronspice's static elements
// cannot hold, so circuits with gaps are stepped by this small nodal solver:
// backward Euler for capacitors and inductors, Newton on the cells'
// piecewise-linear curves, gap states updated after each step. Circuits
// without gaps are checked against kronspice (TestSolverTwin).
//
// Energy: at every step Kirchhoff's current law holds exactly, so the sum
// of v·i over all elements is zero to round-off (Tellegen). The energy a
// capacitor or inductor takes in is split into its stored energy and
// backward Euler's numerical damping, which is reported separately.

// Net is a circuit for the solver.
type Net struct {
	nodes    map[string]int
	el       []*part
	last     []float64 // node voltages after the last step
	hasCells bool
	stored0  float64
}

type kind int

const (
	kR kind = iota
	kC
	kL
	kV // waveform source with series resistance
	kCell
	kGap
)

type part struct {
	name  string
	k     kind
	a, b  int
	val   float64                 // R [Ω], C [F], L [H], source resistance [Ω]
	src   func(t float64) float64 // kV
	cell  Cell
	gap   PolarGap
	group string // energy-report group

	// state
	v, i      float64 // last voltage a−b and current a→b
	on        bool
	armed     bool // a triggered gap fires on its next step if it holds voltage
	sgn       float64
	tOn, tOff float64
	strikesF  int
	strikesB  int
	heat      float64 // energy dissipated [J]
	heatBack  float64 // gaps: energy dissipated in backward arcs [J]
	numLoss   float64 // backward-Euler damping (C, L) [J]
	supplied  float64 // energy delivered by an ideal source [J]
}

// NewNet makes an empty circuit; node "0" is ground.
func NewNet() *Net { return &Net{nodes: map[string]int{"0": 0}} }

func (n *Net) node(s string) int {
	if k, ok := n.nodes[s]; ok {
		return k
	}
	k := len(n.nodes)
	n.nodes[s] = k
	return k
}

func (n *Net) add(p *part, a, b string) *part {
	p.a, p.b = n.node(a), n.node(b)
	p.tOff = -1
	if p.group == "" {
		p.group = p.name
	}
	n.el = append(n.el, p)
	return p
}

// R adds a resistor.
func (n *Net) R(name, a, b string, ohms float64) *part {
	return n.add(&part{name: name, k: kR, val: ohms}, a, b)
}

// C adds a capacitor with initial voltage v0.
func (n *Net) C(name, a, b string, f, v0 float64) *part {
	return n.add(&part{name: name, k: kC, val: f, v: v0}, a, b)
}

// L adds an inductor.
func (n *Net) L(name, a, b string, h float64) *part {
	return n.add(&part{name: name, k: kL, val: h}, a, b)
}

// V adds a voltage source v(t) (plus at a) with series resistance rs.
func (n *Net) V(name, a, b string, v func(float64) float64, rs float64) *part {
	return n.add(&part{name: name, k: kV, src: v, val: rs}, a, b)
}

// Cell adds a water cell (anode electrode at a) and, for valve metals, its
// oxide film capacitor.
func (n *Net) Cell(name, a, b string, c Cell, group string) *part {
	p := n.add(&part{name: name, k: kCell, cell: c, group: group}, a, b)
	if cf := c.FilmC(); cf > 0 {
		n.add(&part{name: "C" + name, k: kC, val: cf, group: group}, a, b)
	}
	return p
}

// Gap adds a polarity spark gap (point at a).
func (n *Net) Gap(name, a, b string, g PolarGap, group string) *part {
	return n.add(&part{name: name, k: kGap, gap: g, group: group}, a, b)
}

// slope returns the cell's local conductance and equivalent current at v.
func cellLin(c Cell, v float64) (g, ieq float64) {
	h := 1e-6 * math.Max(1, math.Abs(v))
	i0 := c.Current(v)
	// one-sided slope in the direction of v keeps piecewise segments exact
	g = (c.Current(v+h) - c.Current(v-h)) / (2 * h)
	if g < 1e-12 {
		g = 1e-12
	}
	return g, i0 - g*v
}

// ErrNoConvergence: Newton on the cells did not settle.
var ErrNoConvergence = errors.New("sparkgap_diodebridge: Newton did not converge")

// Step advances the circuit from t to t+dt.
func (n *Net) Step(t, dt float64) error {
	N := len(n.nodes) - 1
	x := make([]float64, N+1) // node voltages, x[0] = 0
	copy(x, n.last)           // starting guess: the last solution
	var err error
	for it := 0; it < 200; it++ {
		G := make([][]float64, N)
		for i := range G {
			G[i] = make([]float64, N)
		}
		I := make([]float64, N)
		stamp := func(a, b int, g, cur float64) { // i(a→b) = g·v + cur
			if a > 0 {
				G[a-1][a-1] += g
				I[a-1] -= cur
			}
			if b > 0 {
				G[b-1][b-1] += g
				I[b-1] += cur
			}
			if a > 0 && b > 0 {
				G[a-1][b-1] -= g
				G[b-1][a-1] -= g
			}
		}
		for _, p := range n.el {
			switch p.k {
			case kR:
				stamp(p.a, p.b, 1/p.val, 0)
			case kC:
				g := p.val / dt
				stamp(p.a, p.b, g, -g*p.v)
			case kL:
				g := dt / p.val
				stamp(p.a, p.b, g, p.i)
			case kV:
				stamp(p.a, p.b, 1/p.val, -p.src(t+dt)/p.val)
			case kCell:
				g, ieq := cellLin(p.cell, x[p.a]-x[p.b])
				stamp(p.a, p.b, g, ieq)
			case kGap:
				if p.on {
					stamp(p.a, p.b, 1/p.gap.Rarc, -p.sgn*p.gap.Varc/p.gap.Rarc)
				} else {
					stamp(p.a, p.b, 1/p.gap.Rleak(), 0)
				}
			}
		}
		for i := 0; i < N; i++ {
			G[i][i] += 1e-12 // keeps floating nodes defined
		}
		sol, e := gauss(G, I)
		if e != nil {
			return e
		}
		change := 0.0
		for i := 0; i < N; i++ {
			change = math.Max(change, math.Abs(sol[i]-x[i+1]))
			x[i+1] = sol[i]
		}
		if change < 1e-9*(1+maxAbs(x)) || !n.hasCells {
			err = nil
			break
		}
		err = ErrNoConvergence
	}
	if err != nil {
		return err
	}
	n.last = x
	// update states and energy
	for _, p := range n.el {
		v := x[p.a] - x[p.b]
		switch p.k {
		case kR:
			p.i = v / p.val
			p.heat += v * p.i * dt
		case kC:
			i := p.val * (v - p.v) / dt
			p.numLoss += 0.5 * p.val * (v - p.v) * (v - p.v) // v·i·dt − Δ(½Cv²)
			p.i = i
		case kL:
			i := p.i + dt/p.val*v
			p.numLoss += 0.5 * p.val * (i - p.i) * (i - p.i)
			p.i = i
		case kV:
			vs := p.src(t + dt)
			p.i = (v - vs) / p.val // into the + terminal
			p.supplied += -vs * p.i * dt
			p.heat += p.i * p.i * p.val * dt
		case kCell:
			p.i = p.cell.Current(v)
			p.heat += v * p.i * dt
		case kGap:
			g := p.gap
			if p.on {
				p.i = (v - p.sgn*g.Varc) / g.Rarc
			} else {
				p.i = v / g.Rleak()
			}
			p.heat += v * p.i * dt
			if p.on && p.sgn < 0 {
				p.heatBack += v * p.i * dt
			}
			if p.on {
				if t+dt-p.tOn >= g.Hold && p.i*p.sgn < g.Iext {
					p.on, p.tOff = false, t+dt
				}
			} else {
				k := 1.0
				if g.Recovery > 0 && p.tOff >= 0 {
					k = -math.Expm1(-(t + dt - p.tOff) / g.Recovery)
				}
				fwd, rev := g.restrike(k, g.Vpos()), g.restrike(k, g.Vneg())
				if p.armed && math.Abs(v) > 2*g.Varc { // trigger spark
					fwd, rev = 0, 0
				}
				p.armed = false
				if v >= fwd {
					p.on, p.sgn, p.tOn = true, 1, t+dt
					p.strikesF++
				} else if v <= -rev {
					p.on, p.sgn, p.tOn = true, -1, t+dt
					p.strikesB++
				}
			}
		}
		p.v = v
	}
	return nil
}

// Prepare must be called once after the circuit is built.
func (n *Net) Prepare() {
	for _, p := range n.el {
		if p.k == kCell {
			n.hasCells = true
		}
	}
	n.last = make([]float64, len(n.nodes))
	n.stored0 = n.stored()
	// capacitors with an initial voltage: start the guess there when they sit on ground
	for _, p := range n.el {
		if p.k == kC && p.v != 0 && p.b == 0 {
			n.last[p.a] = p.v
		}
	}
}

// Volt returns node s's voltage after the last step.
func (n *Net) Volt(s string) float64 { return n.last[n.nodes[s]] }

// Ledger is the circuit's energy account [J].
type Ledger struct {
	Supplied float64            // by the ideal sources
	Heat     map[string]float64 // dissipated, per group
	Stored   float64            // ½Cv² + ½Li² now
	Stored0  float64            // at the start
	NumLoss  float64            // backward-Euler damping
	Residual float64            // Supplied − heat − ΔStored − NumLoss
}

func (n *Net) stored() float64 {
	s := 0.0
	for _, p := range n.el {
		switch p.k {
		case kC:
			s += 0.5 * p.val * p.v * p.v
		case kL:
			s += 0.5 * p.val * p.i * p.i
		}
	}
	return s
}

// Ledger reports the energy account since Prepare.
func (n *Net) Ledger() Ledger {
	l := Ledger{Heat: map[string]float64{}, Stored: n.stored(), Stored0: n.stored0}
	heat := 0.0
	for _, p := range n.el {
		l.Supplied += p.supplied
		l.NumLoss += p.numLoss
		if p.k != kC && p.k != kL {
			l.Heat[p.group] += p.heat
			heat += p.heat
		}
	}
	l.Residual = l.Supplied - heat - (l.Stored - l.Stored0) - l.NumLoss
	return l
}

func maxAbs(x []float64) float64 {
	m := 0.0
	for _, v := range x {
		m = math.Max(m, math.Abs(v))
	}
	return m
}

// ErrSingular: the nodal matrix has no solution.
var ErrSingular = errors.New("sparkgap_diodebridge: singular nodal matrix")

func gauss(A [][]float64, b []float64) ([]float64, error) {
	n := len(b)
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if math.Abs(A[r][c]) > math.Abs(A[p][c]) {
				p = r
			}
		}
		if A[p][c] == 0 {
			return nil, ErrSingular
		}
		A[c], A[p] = A[p], A[c]
		b[c], b[p] = b[p], b[c]
		for r := c + 1; r < n; r++ {
			f := A[r][c] / A[c][c]
			if f == 0 {
				continue
			}
			for k := c; k < n; k++ {
				A[r][k] -= f * A[c][k]
			}
			b[r] -= f * b[c]
		}
	}
	x := make([]float64, n)
	for r := n - 1; r >= 0; r-- {
		s := b[r]
		for k := r + 1; k < n; k++ {
			s -= A[r][k] * x[k]
		}
		x[r] = s / A[r][r]
	}
	return x, nil
}
