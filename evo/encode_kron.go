// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"fmt"
	"math"
	"math/rand/v2"

	"github.com/fdamador/kronspice/circuit"
)

// Connection-tensor encoding (CLAUDE.md §3a, "Why this matters for the
// genetic algorithm": "the connection tensor gives a clean, algebraic genome
// for topology evolution: mutating C rewires the circuit while the primitive
// branch library stays fixed").
//
// Kron starts from the primitive network — every branch on its own,
// unconnected — and obtains the circuit by a change of reference frame
// that identifies the branches' junctions with the network's (Kron, Tensor
// Analysis of Networks, 1939, Ch. III pp. 75–76; Ch. V p. 133: "separate
// networks are interconnected by a C"). Here the primitive network is a
// fixed library of slots (Topology.Library), one branch each, and the
// genome is that interconnection: for every slot the network junction each
// of its terminals is identified with (the rows of the 0/1 junction
// tensor P, one 1 per primitive junction, so the branch–node incidence is
// A = P_a − P_b), plus the branch's value. A slot can also be open: its
// junction pair inactive, present in the frame and carrying no current
// (Kron 1939 p. 485; p. 429: inactive junction pairs may be ignored), so it
// leaves the circuit without leaving the genome.
//
// Because every genome has the same slots in the same order, crossover is
// aligned: a child takes each slot's row of P from one parent or the
// other, the branch keeping its kind. The direct encoding's one-point
// crossover splices unrelated lists. The mutations are Kron's own
// transformations of the frame:
//
//	value     change one branch's value
//	connect   open or close a slot's junction pair
//	rewire    move one terminal to another junction (one row of P)
//	identify  join two junctions into one (interconnection, Ch. V p. 133)
//	tear      split one junction into two, its terminals divided between
//	          them (Kron, Diakoptics, 1963; tearing as a transformation,
//	          Quart. Appl. Math. 13, 147, 1955)
//	flip      reverse a polarized branch (sources, diodes, transistors,
//	          ECE parts with κ_a ≠ κ_b)
//
// Internal junction labels are arbitrary; repair renumbers them in order
// of first use over the slots, so parents of common descent label the same
// junction alike and crossover keeps it whole.

// slotKinds returns the library: Topology.Library, or MaxParts slots
// dealt round-robin over Kinds.
func (t *Topology) slotKinds() []circuit.Kind {
	if len(t.Library) > 0 {
		return t.Library
	}
	lib := make([]circuit.Kind, t.MaxParts)
	for i := range lib {
		lib[i] = t.Kinds[i%len(t.Kinds)].Kind
	}
	return lib
}

func (t *Topology) validateKron() error {
	if t.Encoding != Connection {
		if len(t.Library) > 0 {
			return fmt.Errorf("%w: library is for the kron (connection) encoding", ErrBadTopology)
		}
		return nil
	}
	if len(t.Library) > 10000 {
		return fmt.Errorf("%w: library of %d slots", ErrBadTopology, len(t.Library))
	}
	for _, k := range t.Library {
		if _, ok := t.kind(byte(k)); !ok {
			return fmt.Errorf("%w: library slot %v is not a listed part kind", ErrBadTopology, k)
		}
	}
	return nil
}

// kronActive is the probability that a slot starts connected.
func (t *Topology) kronActive() float64 {
	return math.Min(1, float64(t.InitParts)/float64(len(t.slotKinds())))
}

// kronSlot draws terminals and a value for slot kind k.
func (t *Topology) kronSlot(r *rand.Rand, k circuit.Kind, parts []Part) Part {
	pk, _ := t.kind(byte(k))
	pool := t.nodePool(parts)
	lo, hi := math.Log10(pk.Min), math.Log10(pk.Max)
	p := Part{Kind: byte(k), X: lo + r.Float64()*(hi-lo)}
	t.kronTerminals(r, &p, pool)
	return p
}

// kronTerminals draws distinct-where-possible junctions for p from pool.
func (t *Topology) kronTerminals(r *rand.Rand, p *Part, pool []int) {
	p.A = pool[r.IntN(len(pool))]
	p.B = pool[r.IntN(len(pool))]
	for tries := 0; p.B == p.A && tries < 8; tries++ {
		p.B = pool[r.IntN(len(pool))]
	}
	p.C = 0
	if threeTerminal(p.Kind) {
		p.C = pool[r.IntN(len(pool))]
		for tries := 0; (p.C == p.A || p.C == p.B) && tries < 8; tries++ {
			p.C = pool[r.IntN(len(pool))]
		}
	}
}

// active returns the connected slots (for node pools).
func activeParts(parts []Part) []Part {
	var out []Part
	for _, p := range parts {
		if !p.Off {
			out = append(out, p)
		}
	}
	return out
}

func (t *Topology) kronRandomGenome(r *rand.Rand) Genome {
	lib := t.slotKinds()
	pOn := t.kronActive()
	parts := make([]Part, len(lib))
	on := 0
	for i, k := range lib {
		parts[i] = t.kronSlot(r, k, activeParts(parts[:i]))
		parts[i].Off = r.Float64() >= pOn
		if !parts[i].Off {
			on++
		}
	}
	if on == 0 {
		parts[r.IntN(len(parts))].Off = false
	}
	g := Genome{Parts: parts}
	t.kronCanonical(&g)
	return g
}

// kronCrossover is uniform crossover on aligned slots.
func (t *Topology) kronCrossover(r *rand.Rand, a, b Genome) (Genome, Genome) {
	c1, c2 := a.Clone(), b.Clone()
	n := min(len(c1.Parts), len(c2.Parts))
	for i := 0; i < n; i++ {
		if r.IntN(2) == 0 {
			c1.Parts[i], c2.Parts[i] = c2.Parts[i], c1.Parts[i]
		}
	}
	t.kronCanonical(&c1)
	t.kronCanonical(&c2)
	return c1, c2
}

// Kron mutation operators and their weights.
const (
	kmValue = iota
	kmConnect
	kmRewire
	kmIdentify
	kmTear
	kmNum
	kmFlip = kmNum // drawn only when a part kind is polarized
)

var kmWeight = [kmNum]float64{kmValue: 0.35, kmConnect: 0.2, kmRewire: 0.2, kmIdentify: 0.12, kmTear: 0.13}

// KronOpNames names the connection-tensor mutations (for logs).
var KronOpNames = [kmNum + 1]string{"value", "connect", "rewire", "identify", "tear", "flip"}

func (t *Topology) kronMutate(r *rand.Rand, g *Genome, sigma float64) {
	pol := t.polarized()
	for first := true; first || r.Float64() < 0.3; first = false {
		if pol && r.Float64() < 0.1 {
			t.kronApply(r, g, kmFlip, sigma)
			continue
		}
		u := r.Float64()
		op := 0
		for ; op < kmNum-1; op++ {
			if u < kmWeight[op] {
				break
			}
			u -= kmWeight[op]
		}
		t.kronApply(r, g, op, sigma)
	}
	t.kronCanonical(g)
}

// onSlots lists the indices of connected slots.
func onSlots(parts []Part) []int {
	var on []int
	for i, p := range parts {
		if !p.Off {
			on = append(on, i)
		}
	}
	return on
}

// kronApply applies one operator. Operators that need a connected slot
// fall back to connect when none is.
func (t *Topology) kronApply(r *rand.Rand, g *Genome, op int, sigma float64) {
	if len(g.Parts) == 0 {
		return
	}
	on := onSlots(g.Parts)
	if len(on) == 0 {
		op = kmConnect
	}
	switch op {
	case kmValue:
		g.Parts[on[r.IntN(len(on))]].X += sigma * r.NormFloat64()
	case kmConnect:
		k := r.IntN(len(g.Parts))
		p := &g.Parts[k]
		if !p.Off {
			if len(on) > 1 {
				p.Off = true
			}
			return
		}
		p.Off = false
		// Close the junction pair where it was, or attach it afresh to the
		// junctions in use.
		if r.IntN(2) == 0 {
			t.kronTerminals(r, p, t.nodePool(activeParts(g.Parts)))
		}
	case kmRewire:
		p := &g.Parts[on[r.IntN(len(on))]]
		pool := t.nodePool(activeParts(g.Parts))
		ends := 2
		if threeTerminal(p.Kind) {
			ends = 3
		}
		n := pool[r.IntN(len(pool))]
		switch r.IntN(ends) {
		case 0:
			p.A = n
		case 1:
			p.B = n
		default:
			p.C = n
		}
	case kmIdentify:
		// Join junction j into i; j is internal, i any junction in use.
		used := usedNodes(activeParts(g.Parts), t.nodes())
		term := 1 + len(t.Terminals)
		var internal, all []int
		for n, u := range used {
			if !u {
				continue
			}
			all = append(all, n)
			if n >= term {
				internal = append(internal, n)
			}
		}
		if len(internal) == 0 || len(all) < 2 {
			return
		}
		j := internal[r.IntN(len(internal))]
		i := all[r.IntN(len(all))]
		for tries := 0; i == j && tries < 8; tries++ {
			i = all[r.IntN(len(all))]
		}
		if i == j {
			return
		}
		for k := range g.Parts {
			for _, e := range []*int{&g.Parts[k].A, &g.Parts[k].B, &g.Parts[k].C} {
				if *e == j {
					*e = i
				}
			}
		}
	case kmTear:
		// Split junction n: move a random part (≥ 1, not all) of the
		// terminals on it to a fresh internal junction.
		act := activeParts(g.Parts)
		pool := t.nodePool(act)
		fresh := pool[len(pool)-1]
		if fresh < 1+len(t.Terminals) || usedNodes(act, t.nodes())[fresh] {
			return // no free internal junction
		}
		type end struct{ slot, which int }
		var ends map[int][]end = map[int][]end{}
		for _, k := range on {
			p := g.Parts[k]
			ts := p.terminals()
			for w, n := range ts {
				ends[n] = append(ends[n], end{k, w})
			}
		}
		var cand []int
		for n := 0; n < t.nodes(); n++ {
			if len(ends[n]) >= 2 {
				cand = append(cand, n)
			}
		}
		if len(cand) == 0 {
			return
		}
		n := cand[r.IntN(len(cand))]
		es := ends[n]
		move := make([]bool, len(es))
		moved := 0
		for i := range es {
			if r.IntN(2) == 0 {
				move[i] = true
				moved++
			}
		}
		if moved == 0 {
			move[r.IntN(len(es))], moved = true, 1
		}
		if moved == len(es) {
			move[r.IntN(len(es))] = false
		}
		for i, e := range es {
			if !move[i] {
				continue
			}
			p := &g.Parts[e.slot]
			switch e.which {
			case 0:
				p.A = fresh
			case 1:
				p.B = fresh
			default:
				p.C = fresh
			}
		}
	case kmFlip:
		p := &g.Parts[on[r.IntN(len(on))]]
		p.A, p.B = p.B, p.A
	}
}

// usedNodes marks the junctions the parts touch.
func usedNodes(parts []Part, nodes int) []bool {
	used := make([]bool, nodes)
	for _, p := range parts {
		for _, n := range p.terminals() {
			if n >= 0 && n < nodes {
				used[n] = true
			}
		}
	}
	return used
}

// kronCanonical renumbers internal junctions in order of first use over
// the slots, connected slots first, so the same junction gets the same
// label in related genomes.
func (t *Topology) kronCanonical(g *Genome) {
	term := 1 + len(t.Terminals)
	label := map[int]int{}
	next := term
	visit := func(off bool) {
		for i := range g.Parts {
			p := &g.Parts[i]
			if p.Off != off {
				continue
			}
			for _, e := range []*int{&p.A, &p.B, &p.C} {
				if *e < term {
					continue
				}
				if _, ok := label[*e]; !ok {
					label[*e] = next
					next++
				}
			}
		}
	}
	visit(false)
	visit(true)
	for i := range g.Parts {
		p := &g.Parts[i]
		for _, e := range []*int{&p.A, &p.B, &p.C} {
			if l, ok := label[*e]; ok {
				*e = l
			}
		}
	}
	// Relabelling keeps the number of distinct internal junctions, so the
	// labels stay within t.nodes().
}

// kronRepair keeps the slots aligned with the library: missing slots are
// added open, extra ones dropped, kinds reset to the library's.
func (t *Topology) kronRepair(g *Genome) {
	lib := t.slotKinds()
	if len(g.Parts) > len(lib) {
		g.Parts = g.Parts[:len(lib)]
	}
	for len(g.Parts) < len(lib) {
		k := lib[len(g.Parts)]
		pk, _ := t.kind(byte(k))
		g.Parts = append(g.Parts, Part{Kind: byte(k), X: (math.Log10(pk.Min) + math.Log10(pk.Max)) / 2, Off: true})
	}
	for i := range g.Parts {
		if g.Parts[i].Kind != byte(lib[i]) {
			pk, _ := t.kind(byte(lib[i]))
			g.Parts[i].Kind = byte(lib[i])
			g.Parts[i].X = (math.Log10(pk.Min) + math.Log10(pk.Max)) / 2
		}
	}
}
