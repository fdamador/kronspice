// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
)

// tcNet is a small nodal solver for the spark-gap Tesla coil: R, C (with an
// initial voltage), L, coupled pairs of inductors, voltage sources with a
// series resistance, and spark gaps (the PolarGap of gap.go). It steps by
// the trapezoidal rule, which keeps an LC ring-down's energy (the solver in
// sim.go uses backward Euler, whose numerical damping would eat a
// high-Q Tesla coil's ring-down), and updates each gap's state after the
// step, as sim.go does. Circuits without gaps are checked against kronspice
// (TestTeslaTwin). Energy: sources' work = heat in R and gaps + Δ(½Cv² +
// ½iᵀLi), with mid-step averages (closes to round-off).
type tcNet struct {
	nodes map[string]int
	el    []*tcPart
	v     []float64 // node voltages
	t     float64
	// ledger [J]
	supplied, heat, gapHeat float64
	stored0                 float64
}

type tcKind int

const (
	tcR tcKind = iota
	tcC
	tcL
	tcV
	tcGap
)

type tcPart struct {
	name string
	k    tcKind
	a, b int
	val  float64 // R, C, L, or source series R
	src  func(t float64) float64
	gap  PolarGap
	// coupling: partner inductor and mutual inductance
	pair *tcPart
	m    float64
	// state
	vab, i float64 // last branch voltage a−b and current a→b
	ieq    float64 // companion history current
	on     bool
	sgn    float64
	tOn    float64
	tOff   float64
	armed  bool
	strk   int
	forced bool // gap state set from outside (a timed switch)
}

func newTC() *tcNet { return &tcNet{nodes: map[string]int{"0": 0}} }

func (n *tcNet) node(s string) int {
	if k, ok := n.nodes[s]; ok {
		return k
	}
	k := len(n.nodes)
	n.nodes[s] = k
	return k
}

func (n *tcNet) add(p *tcPart, a, b string) *tcPart {
	p.a, p.b = n.node(a), n.node(b)
	p.tOff = -1
	n.el = append(n.el, p)
	return p
}

func (n *tcNet) R(name, a, b string, ohms float64) *tcPart {
	return n.add(&tcPart{name: name, k: tcR, val: ohms}, a, b)
}
func (n *tcNet) C(name, a, b string, f, v0 float64) *tcPart {
	return n.add(&tcPart{name: name, k: tcC, val: f, vab: v0}, a, b)
}
func (n *tcNet) L(name, a, b string, h float64) *tcPart {
	return n.add(&tcPart{name: name, k: tcL, val: h}, a, b)
}
func (n *tcNet) V(name, a, b string, v func(float64) float64, rs float64) *tcPart {
	return n.add(&tcPart{name: name, k: tcV, src: v, val: rs}, a, b)
}
func (n *tcNet) Gap(name, a, b string, g PolarGap) *tcPart {
	return n.add(&tcPart{name: name, k: tcGap, gap: g}, a, b)
}

// Couple couples two inductors with coefficient k.
func (n *tcNet) Couple(l1, l2 *tcPart, k float64) {
	m := k * math.Sqrt(l1.val*l2.val)
	l1.pair, l1.m = l2, m
	l2.pair, l2.m = l1, m
}

func (n *tcNet) stored() float64 {
	s := 0.0
	for _, p := range n.el {
		switch p.k {
		case tcC:
			s += 0.5 * p.val * p.vab * p.vab
		case tcL:
			s += 0.5 * p.val * p.i * p.i
			if p.pair != nil && p.name < p.pair.name {
				s += p.m * p.i * p.pair.i
			}
		}
	}
	return s
}

// Prepare finds a consistent starting state (capacitors at their initial
// voltage, inductors carrying no current, gaps in their state) and the
// ledger's start: every branch starts with the voltage and current of that
// solution, so Kirchhoff's laws hold from t = 0.
func (n *tcNet) Prepare() error {
	N := len(n.nodes) - 1
	A := make([][]float64, N)
	for i := range A {
		A[i] = make([]float64, N)
	}
	rhs := make([]float64, N)
	const gc = 1e3                          // a capacitor held at its voltage (stiff)
	stamp := func(a, b int, G, I float64) { // current a→b = G·v_ab + I
		if a > 0 {
			A[a-1][a-1] += G
			rhs[a-1] -= I
		}
		if b > 0 {
			A[b-1][b-1] += G
			rhs[b-1] += I
		}
		if a > 0 && b > 0 {
			A[a-1][b-1] -= G
			A[b-1][a-1] -= G
		}
	}
	gOf := func(p *tcPart) (G, I float64) {
		switch p.k {
		case tcR:
			return 1 / p.val, 0
		case tcC:
			return gc, -gc * p.vab
		case tcV:
			return 1 / p.val, -p.src(0) / p.val
		case tcGap:
			if p.on {
				return 1 / p.gap.Rarc, -p.sgn * p.gap.Varc / p.gap.Rarc
			}
			return 1 / p.gap.Rleak(), 0
		}
		return 0, 0 // inductor: no current at the start
	}
	for _, p := range n.el {
		if G, I := gOf(p); G != 0 {
			stamp(p.a, p.b, G, I)
		}
	}
	x, err := solveDense(A, rhs)
	if err != nil {
		return err
	}
	n.v = make([]float64, N+1)
	copy(n.v[1:], x)
	for _, p := range n.el {
		v := n.v[p.a] - n.v[p.b]
		G, I := gOf(p)
		if p.k == tcC {
			p.i = G*v + I // the current that holds it
			continue      // keep its initial voltage exactly
		}
		p.vab = v
		if p.k != tcL {
			p.i = G*v + I
		}
	}
	n.stored0 = n.stored()
	return nil
}

// Step advances by dt.
func (n *tcNet) Step(dt float64) error {
	N := len(n.nodes) - 1
	A := make([][]float64, N)
	for i := range A {
		A[i] = make([]float64, N)
	}
	rhs := make([]float64, N)
	g := func(a, b int, G float64) {
		if a > 0 {
			A[a-1][a-1] += G
		}
		if b > 0 {
			A[b-1][b-1] += G
		}
		if a > 0 && b > 0 {
			A[a-1][b-1] -= G
			A[b-1][a-1] -= G
		}
	}
	cur := func(a, b int, I float64) { // current I from a to b through the element (leaves a)
		if a > 0 {
			rhs[a-1] -= I
		}
		if b > 0 {
			rhs[b-1] += I
		}
	}
	gx := func(a, b, c, d int, G float64) { // current a→b = G·(v_c − v_d)
		for _, r := range [2][2]int{{a, 1}, {b, -1}} {
			if r[0] == 0 {
				continue
			}
			if c > 0 {
				A[r[0]-1][c-1] += float64(r[1]) * G
			}
			if d > 0 {
				A[r[0]-1][d-1] -= float64(r[1]) * G
			}
		}
	}
	t1 := n.t + dt
	for _, p := range n.el {
		switch p.k {
		case tcR:
			g(p.a, p.b, 1/p.val)
		case tcC:
			G := 2 * p.val / dt
			p.ieq = -G*p.vab - p.i // i_{n+1} = G·v_{n+1} + ieq
			g(p.a, p.b, G)
			cur(p.a, p.b, p.ieq)
		case tcV:
			g(p.a, p.b, 1/p.val)
			cur(p.a, p.b, -p.src(t1)/p.val)
		case tcGap:
			if p.on {
				g(p.a, p.b, 1/p.gap.Rarc)
				cur(p.a, p.b, -p.sgn*p.gap.Varc/p.gap.Rarc)
			} else {
				g(p.a, p.b, 1/p.gap.Rleak())
			}
		case tcL:
			// i_{n+1} = i_n + (dt/2)·Γ·(v_{n+1} + v_n), Γ = L⁻¹ (2×2 when coupled).
			if p.pair == nil {
				G := dt / (2 * p.val)
				g(p.a, p.b, G)
				cur(p.a, p.b, p.i+G*p.vab)
			} else {
				q := p.pair
				det := p.val*q.val - p.m*p.m
				gs, gm := dt/2*q.val/det, -dt/2*p.m/det
				g(p.a, p.b, gs)
				gx(p.a, p.b, q.a, q.b, gm)
				cur(p.a, p.b, p.i+gs*p.vab+gm*q.vab)
			}
		}
	}
	x, err := solveDense(A, rhs)
	if err != nil {
		return err
	}
	vn := make([]float64, N+1)
	copy(vn[1:], x)
	volt := func(k int) float64 { return vn[k] }
	// New branch quantities; ledger by the trapezoidal rule.
	newI := make([]float64, len(n.el))
	for idx, p := range n.el {
		v := volt(p.a) - volt(p.b)
		switch p.k {
		case tcL:
			if p.pair == nil {
				newI[idx] = p.i + dt/(2*p.val)*(v+p.vab)
			} else {
				q := p.pair
				det := p.val*q.val - p.m*p.m
				vq := volt(q.a) - volt(q.b)
				newI[idx] = p.i + dt/2*(q.val*(v+p.vab)-p.m*(vq+q.vab))/det
			}
		}
	}
	// Ledger with mid-step averages: v̄·ī·dt per element. The averaged
	// voltages and currents obey Kirchhoff's laws, so by Tellegen the
	// account closes to round-off, and for the trapezoidal L and C the
	// averaged power is exactly the change of ½Cv² and ½iᵀLi.
	for idx, p := range n.el {
		v := volt(p.a) - volt(p.b)
		var i float64
		switch p.k {
		case tcR:
			i = v / p.val
		case tcC:
			i = 2*p.val/dt*v + p.ieq
		case tcL:
			i = newI[idx]
		case tcV:
			i = (v - p.src(t1)) / p.val
		case tcGap:
			if p.on {
				i = (v - p.sgn*p.gap.Varc) / p.gap.Rarc
			} else {
				i = v / p.gap.Rleak()
			}
		}
		vb, ib := (v+p.vab)/2, (i+p.i)/2
		switch p.k {
		case tcR:
			n.heat += dt * vb * ib
		case tcV:
			vs := (p.src(t1) + p.src(n.t)) / 2
			n.supplied += -dt * ib * vs
			n.heat += dt * ib * ib * p.val
		case tcGap:
			n.heat += dt * vb * ib
			n.gapHeat += dt * vb * ib
		}
		p.vab, p.i = v, i
	}
	// Gap states after the step.
	for _, p := range n.el {
		if p.k != tcGap || p.forced {
			continue
		}
		gp := p.gap
		v := p.vab
		if p.on {
			if t1-p.tOn >= gp.Hold && p.i*p.sgn < gp.Iext {
				p.on, p.tOff = false, t1
			}
		} else {
			k := 1.0
			if gp.Recovery > 0 && p.tOff >= 0 {
				k = -math.Expm1(-(t1 - p.tOff) / gp.Recovery)
			}
			fwd, rev := gp.restrike(k, gp.Vpos()), gp.restrike(k, gp.Vneg())
			if p.armed && math.Abs(v) > 2*gp.Varc {
				fwd, rev = 0, 0
			}
			p.armed = false
			if v >= fwd {
				p.on, p.sgn, p.tOn = true, 1, t1
				p.strk++
			} else if v <= -rev {
				p.on, p.sgn, p.tOn = true, -1, t1
				p.strk++
			}
		}
	}
	n.v = vn
	n.t = t1
	return nil
}

// Volt is node s's voltage.
func (n *tcNet) Volt(s string) float64 { return n.v[n.nodes[s]] }

// Residual is supplied − heat − Δstored [J].
func (n *tcNet) Residual() float64 { return n.supplied - n.heat - (n.stored() - n.stored0) }

func solveDense(A [][]float64, b []float64) ([]float64, error) {
	n := len(b)
	M := make([][]float64, n)
	for i := range A {
		M[i] = append(append([]float64(nil), A[i]...), b[i])
	}
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if math.Abs(M[r][c]) > math.Abs(M[p][c]) {
				p = r
			}
		}
		if M[p][c] == 0 {
			return nil, errors.New("singular nodal matrix")
		}
		M[c], M[p] = M[p], M[c]
		for r := c + 1; r < n; r++ {
			f := M[r][c] / M[c][c]
			for k := c; k <= n; k++ {
				M[r][k] -= f * M[c][k]
			}
		}
	}
	x := make([]float64, n)
	for r := n - 1; r >= 0; r-- {
		s := M[r][n]
		for k := r + 1; k < n; k++ {
			s -= M[r][k] * x[k]
		}
		x[r] = s / M[r][r]
	}
	return x, nil
}
