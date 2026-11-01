// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// ErrBadTopology is returned for an invalid topology configuration.
var ErrBadTopology = errors.New("evo: invalid topology")

// Topology genome (CLAUDE.md §4, encoding 2): a variable-length list of
// two-terminal parts grown into a fixed embryo, as in Koza et al., Genetic
// Programming III (1999), ch. 25: the embryo is the test fixture — source,
// source resistance, load and the analyses — with nothing in the place the
// circuit will occupy, and evolution decides which parts exist, how they
// connect and their values.
//
// Unlike Koza's developmental (circuit-constructing program) encoding, the
// genome here is the netlist itself — a direct encoding: each Part names a
// kind, two nodes and a value. Nodes are ground, the embryo's terminals
// (for example "in" and "out") and a pool of internal nodes. Operators add,
// delete, rewire, re-kind and re-value parts, split a part into two in
// series through a new node, and cross two netlists over at one point.

// PartKind is a kind of part evolution may place, with its value range.
type PartKind struct {
	Kind     circuit.Kind // Resistor, Capacitor or Inductor
	Min, Max float64      // magnitude range (0 < Min ≤ Max)
	// Series snaps values to a preferred-number series ("" for none).
	Series string
}

// Topology configures a topology genome.
type Topology struct {
	// Terminals are embryo nodes parts may connect to; ground is always
	// available.
	Terminals []string
	// Internal is the number of internal nodes parts may use (default 8).
	Internal int
	// MaxParts bounds the netlist size (default 20).
	MaxParts int
	// InitParts: each initial individual gets 1 to InitParts random parts
	// (default 3). The embryo itself contributes none.
	InitParts int
	Kinds     []PartKind
}

// Part is one gene of a topology genome.
type Part struct {
	Kind byte    `json:"k"` // 'R', 'C', 'L'
	A    int     `json:"a"` // node: 0 ground, 1..T terminals, T+1.. internal
	B    int     `json:"b"`
	X    float64 `json:"x"` // log10 of the magnitude
}

// internalName is the name of internal node k (1-based).
func internalName(k int) string { return fmt.Sprintf("g_n%d", k) }

func (t *Topology) defaults() {
	if t.Internal == 0 {
		t.Internal = 8
	}
	if t.MaxParts == 0 {
		t.MaxParts = 20
	}
	if t.InitParts == 0 {
		t.InitParts = 3
	}
}

// nodes is the number of node indices parts may use (ground included).
func (t *Topology) nodes() int { return 1 + len(t.Terminals) + t.Internal }

func (t *Topology) kind(k byte) (PartKind, bool) {
	for _, pk := range t.Kinds {
		if byte(pk.Kind) == k {
			return pk, true
		}
	}
	return PartKind{}, false
}

// Validate checks the configuration against the embryo.
func (t *Topology) Validate(embryo *circuit.Circuit) error {
	t.defaults()
	if len(t.Terminals) == 0 {
		return fmt.Errorf("%w: no terminals", ErrBadTopology)
	}
	if t.Internal < 0 || t.Internal > 1000 || t.MaxParts < 1 || t.MaxParts > 10000 || t.InitParts < 1 || t.InitParts > t.MaxParts {
		return fmt.Errorf("%w: internal %d, max_parts %d, init_parts %d", ErrBadTopology, t.Internal, t.MaxParts, t.InitParts)
	}
	have := map[string]bool{}
	for _, n := range embryo.NodeNames {
		have[strings.ToLower(n)] = true
	}
	for _, n := range t.Terminals {
		if !have[strings.ToLower(n)] || n == "0" || strings.EqualFold(n, "gnd") {
			return fmt.Errorf("%w: terminal %q is not a node of the embryo", ErrBadTopology, n)
		}
	}
	for k := 1; k <= t.Internal; k++ {
		if have[internalName(k)] {
			return fmt.Errorf("%w: the embryo already uses node name %s", ErrBadTopology, internalName(k))
		}
	}
	if len(t.Kinds) == 0 {
		return fmt.Errorf("%w: no part kinds", ErrBadTopology)
	}
	seen := map[circuit.Kind]bool{}
	for _, pk := range t.Kinds {
		if pk.Kind != circuit.Resistor && pk.Kind != circuit.Capacitor && pk.Kind != circuit.Inductor {
			return fmt.Errorf("%w: part kind %v (R, C or L)", ErrBadTopology, pk.Kind)
		}
		if seen[pk.Kind] {
			return fmt.Errorf("%w: part kind %v listed twice", ErrBadTopology, pk.Kind)
		}
		seen[pk.Kind] = true
		if !(pk.Min > 0) || !(pk.Max >= pk.Min) || math.IsInf(pk.Max, 0) {
			return fmt.Errorf("%w: %v range [%g, %g]", ErrBadTopology, pk.Kind, pk.Min, pk.Max)
		}
		if _, err := series(pk.Series); err != nil {
			return fmt.Errorf("%w: %v", ErrBadTopology, err)
		}
	}
	return nil
}

// value returns a part's decoded magnitude.
func (t *Topology) value(p Part) float64 {
	pk, _ := t.kind(p.Kind)
	ser, _ := series(pk.Series)
	return snap(math.Pow(10, p.X), ser, pk.Min, pk.Max)
}

// prune returns the parts that matter: a part connecting a node to itself
// is dropped, and so, repeatedly, is any part with an end on an internal
// node that nothing else touches (it carries no current, so removing it
// leaves every node voltage unchanged). Evolution may leave such debris;
// pruning keeps it from causing floating-node errors.
func (t *Topology) prune(parts []Part) []Part {
	term := 1 + len(t.Terminals)
	keep := make([]bool, len(parts))
	for i, p := range parts {
		keep[i] = p.A != p.B
	}
	deg := make([]int, t.nodes())
	for changed := true; changed; {
		changed = false
		for i := range deg {
			deg[i] = 0
		}
		for i, p := range parts {
			if keep[i] {
				deg[p.A]++
				deg[p.B]++
			}
		}
		for i, p := range parts {
			if keep[i] && ((p.A >= term && deg[p.A] == 1) || (p.B >= term && deg[p.B] == 1)) {
				keep[i], changed = false, true
			}
		}
	}
	var out []Part
	for i, p := range parts {
		if keep[i] {
			out = append(out, p)
		}
	}
	return out
}

// Decode builds the circuit: the embryo plus the pruned parts, named
// R_g1, C_g2, … in genome order. It returns the part values in that order.
func (t *Topology) Decode(embryo *circuit.Circuit, g Genome) (*circuit.Circuit, []float64) {
	c := embryo.Clone()
	idx := map[string]int{}
	for i, n := range c.NodeNames {
		idx[strings.ToLower(n)] = i
	}
	node := func(k int) int {
		var name string
		switch {
		case k == 0:
			return 0
		case k <= len(t.Terminals):
			name = strings.ToLower(t.Terminals[k-1])
		default:
			name = internalName(k - len(t.Terminals))
		}
		if i, ok := idx[name]; ok {
			return i
		}
		c.NodeNames = append(c.NodeNames, name)
		idx[name] = len(c.NodeNames) - 1
		return idx[name]
	}
	parts := t.prune(g.Parts)
	vals := make([]float64, len(parts))
	for i, p := range parts {
		v := t.value(p)
		vals[i] = v
		c.Elements = append(c.Elements, circuit.Element{
			Name: fmt.Sprintf("%c_g%d", p.Kind, i+1), Kind: circuit.Kind(p.Kind),
			Nodes: [2]int{node(p.A), node(p.B)}, Value: v,
		})
	}
	return c, vals
}

// randomPart draws a part whose nodes come from those already in use plus
// one fresh internal node, so new parts mostly attach to the existing
// circuit instead of floating off on their own.
func (t *Topology) randomPart(r *rand.Rand, parts []Part) Part {
	pk := t.Kinds[r.IntN(len(t.Kinds))]
	pool := t.nodePool(parts)
	a := pool[r.IntN(len(pool))]
	b := pool[r.IntN(len(pool))]
	for tries := 0; b == a && tries < 8; tries++ {
		b = pool[r.IntN(len(pool))]
	}
	lo, hi := math.Log10(pk.Min), math.Log10(pk.Max)
	return Part{Kind: byte(pk.Kind), A: a, B: b, X: lo + r.Float64()*(hi-lo)}
}

// nodePool lists ground, the terminals, the internal nodes in use and the
// lowest unused internal node.
func (t *Topology) nodePool(parts []Part) []int {
	term := 1 + len(t.Terminals)
	used := make([]bool, t.nodes())
	for _, p := range parts {
		used[p.A], used[p.B] = true, true
	}
	pool := make([]int, 0, t.nodes())
	for k := 0; k < term; k++ {
		pool = append(pool, k)
	}
	fresh := false
	for k := term; k < t.nodes(); k++ {
		if used[k] {
			pool = append(pool, k)
		} else if !fresh {
			pool = append(pool, k)
			fresh = true
		}
	}
	return pool
}

// randomGenome draws 1..InitParts random parts.
func (t *Topology) randomGenome(r *rand.Rand) Genome {
	n := 1 + r.IntN(t.InitParts)
	var parts []Part
	for i := 0; i < n; i++ {
		parts = append(parts, t.randomPart(r, parts))
	}
	return Genome{Parts: parts}
}

// crossover is one-point crossover on the part lists (Koza-style subtree
// exchange in netlist form): the children are a[:i]+b[j:] and b[:j]+a[i:].
func (t *Topology) crossover(r *rand.Rand, a, b Genome) (Genome, Genome) {
	i, j := r.IntN(len(a.Parts)+1), r.IntN(len(b.Parts)+1)
	c1 := append(append([]Part(nil), a.Parts[:i]...), b.Parts[j:]...)
	c2 := append(append([]Part(nil), b.Parts[:j]...), a.Parts[i:]...)
	return Genome{Parts: c1}, Genome{Parts: c2}
}

// Mutation operators and their relative weights.
const (
	opValue = iota
	opAdd
	opDelete
	opRewire
	opKind
	opSeries
	opParallel
	numOps
)

var opWeight = [numOps]float64{opValue: 0.35, opAdd: 0.15, opDelete: 0.12, opRewire: 0.13, opKind: 0.05,
	opSeries: 0.1, opParallel: 0.1}

// mutate applies one operator, then each further one with probability
// 0.3 (a geometric number of operators, mean 1.43). sigma is the value
// step in decades.
func (t *Topology) mutate(r *rand.Rand, g *Genome, sigma float64) {
	for first := true; first || r.Float64() < 0.3; first = false {
		u := r.Float64()
		op := 0
		for ; op < numOps-1; op++ {
			if u < opWeight[op] {
				break
			}
			u -= opWeight[op]
		}
		t.apply(r, g, op, sigma)
	}
	for len(g.Parts) > t.MaxParts {
		k := r.IntN(len(g.Parts))
		g.Parts = append(g.Parts[:k], g.Parts[k+1:]...)
	}
}

func (t *Topology) apply(r *rand.Rand, g *Genome, op int, sigma float64) {
	n := len(g.Parts)
	if n == 0 {
		op = opAdd
	}
	switch op {
	case opValue:
		p := &g.Parts[r.IntN(n)]
		p.X += sigma * r.NormFloat64()
	case opAdd:
		g.Parts = append(g.Parts, t.randomPart(r, g.Parts))
	case opDelete:
		k := r.IntN(n)
		g.Parts = append(g.Parts[:k], g.Parts[k+1:]...)
	case opRewire:
		p := &g.Parts[r.IntN(n)]
		pool := t.nodePool(g.Parts)
		if r.IntN(2) == 0 {
			p.A = pool[r.IntN(len(pool))]
		} else {
			p.B = pool[r.IntN(len(pool))]
		}
	case opKind:
		p := &g.Parts[r.IntN(n)]
		old, _ := t.kind(p.Kind)
		nk := t.Kinds[r.IntN(len(t.Kinds))]
		// Keep the relative position in the range.
		lo, hi := math.Log10(old.Min), math.Log10(old.Max)
		f := 0.5
		if hi > lo {
			f = (p.X - lo) / (hi - lo)
		}
		p.Kind = byte(nk.Kind)
		nlo, nhi := math.Log10(nk.Min), math.Log10(nk.Max)
		p.X = nlo + f*(nhi-nlo)
	case opSeries:
		// a–b becomes a–m–b with a new part from m to b.
		pool := t.nodePool(g.Parts)
		m := pool[len(pool)-1]
		k := r.IntN(n)
		if m < 1+len(t.Terminals) || g.Parts[k].A == m || g.Parts[k].B == m {
			return // no free internal node
		}
		np := t.randomPart(r, g.Parts)
		np.A, np.B = m, g.Parts[k].B
		g.Parts[k].B = m
		g.Parts = append(g.Parts, np)
	case opParallel:
		np := t.randomPart(r, g.Parts)
		k := r.IntN(n)
		np.A, np.B = g.Parts[k].A, g.Parts[k].B
		g.Parts = append(g.Parts, np)
	}
}

// repair clamps values into range, maps unknown kinds and out-of-range
// nodes back into the configuration, so every genome decodes.
func (t *Topology) repair(g *Genome) {
	nn := t.nodes()
	for i := range g.Parts {
		p := &g.Parts[i]
		pk, ok := t.kind(p.Kind)
		if !ok {
			pk = t.Kinds[0]
			p.Kind = byte(pk.Kind)
		}
		lo, hi := math.Log10(pk.Min), math.Log10(pk.Max)
		if math.IsNaN(p.X) {
			p.X = (lo + hi) / 2
		}
		p.X = math.Min(math.Max(p.X, lo), hi)
		p.A = ((p.A % nn) + nn) % nn
		p.B = ((p.B % nn) + nn) % nn
	}
	if len(g.Parts) > t.MaxParts {
		g.Parts = g.Parts[:t.MaxParts]
	}
}
