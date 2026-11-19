// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"fmt"
	"math"
	"math/rand/v2"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// Developmental encoding (CLAUDE.md §4, encoding 2, second form), after
// J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, Genetic
// Programming III: Darwinian Invention and Problem Solving (Morgan
// Kaufmann, 1999), ch. 25: the genome is not the netlist but a program
// that builds it. The embryo holds the fixture and one or more modifiable
// wires; each wire carries a writing head and is developed by its own
// program tree (Koza's result-producing branches). Running the programs
// rewrites the wires into components.
//
// Function set (construction-continuing positions):
//
//	R, L, C (value, cc)        turn the component at the head into a part
//	                            whose value the arithmetic subtree sets
//	SERIES (cc, cc, cc)         Z(a→b) → Z(a→n1), wire(n1→n2), Z′(n2→b)
//	PARALLEL (cc, cc, cc, cc)   keeps Z(a→b), adds wire(a→n1), Z′(n1→n2),
//	                            wire(n2→b)
//	VIA_k (cc, cc, cc)          Z(a→b) → Z(a→m), Z′(m→b), wire(m→k), k =
//	                            ground or an embryo terminal (Koza's
//	                            TWO_GROUND and TWO_VIA)
//	V, I, D (value, cc)         the same for a voltage source, current
//	                            source or diode (value: DC magnitude or
//	                            area); orientation gives polarity
//	Q, M (value, cc, cc, cc)    replace Z(a→b) by a transistor and three
//	                            lead wires: a→collector, emitter→b and
//	                            base→a (drain, source, gate for M); the
//	                            base lead can grow into a bias part or be
//	                            cut; value: area (Q) or width W (M)
//	FLIP (cc)                   reverses the component's terminals
//	NOP (cc)                    passes the head on
//	END                         the head stops; the component stays
//	SAFE_CUT                    the component is removed (open circuit)
//
// Arithmetic (value) positions: +, − and constants. Koza interprets the
// result on a logarithmic scale; here it is decades away from the middle
// of the part kind's range, then clamped and snapped like every other
// value (an adaptation: GP III fixes the unit per kind, nF, µH, kΩ).
//
// FLIP reverses the component under the head, as Koza's does for diodes
// and transistors. A plain R, L or C has no polarity, so FLIP is drawn
// into random programs only when a part kind is polarized (PartKind.KappaA
// ≠ KappaB: the hypothetical ECE scalar terms, register entry H48); then
// the orientation decides which terminal's potential the body sees scaled
// by 1 + s·κ. Programs may contain FLIP either way (it is then a no-op).
// Koza's automatically defined functions and architecture-altering
// operations are left out (the number of branches is fixed by the
// embryo's wires).
//
// The heads are executed breadth-first, so no branch of the tree is
// developed far ahead of the others (GP III §25.4). Wires still in place
// at the end are wires: their nodes are merged; a wire that joins two
// embryo nodes becomes a 0 V source between them. Then pruning drops
// self-loops and dangling parts as for the direct encoding.

// Encoding selects the topology genome's form.
type Encoding uint8

const (
	// Direct: the genome is the part list (encode_topology.go).
	Direct Encoding = iota
	// Developmental: the genome is a circuit-constructing program.
	Developmental
)

func (e Encoding) String() string {
	if e == Developmental {
		return "developmental"
	}
	return "direct"
}

// ParseEncoding reads "direct" (or "") and "developmental".
func ParseEncoding(s string) (Encoding, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "direct", "netlist":
		return Direct, nil
	case "developmental", "develop", "koza":
		return Developmental, nil
	}
	return 0, fmt.Errorf("%w: encoding %q (direct or developmental)", ErrBadTopology, s)
}

// Op is one function or terminal of a developmental program, stored in
// prefix order.
type Op struct {
	F uint8   `json:"f"`           // function code (dv… constants)
	K int     `json:"k,omitempty"` // VIA target: 0 ground, 1..T terminals
	V float64 `json:"v,omitempty"` // constant (decades) for dvConst
}

// Function codes. The numbering is part of the checkpoint format: append
// new codes, never renumber.
const (
	dvEnd uint8 = iota
	dvCut
	dvNop
	dvR
	dvL
	dvC
	dvSeries
	dvParallel
	dvVia
	dvAdd
	dvSub
	dvConst
	dvFlip
	dvV
	dvI
	dvD
	dvQ
	dvM
	dvNumOps
)

var dvName = [dvNumOps]string{"END", "SAFE_CUT", "NOP", "R", "L", "C", "SERIES", "PARALLEL", "VIA", "+", "-", "CONST", "FLIP",
	"V", "I", "D", "Q", "M"}

// dvArity is the number of subtrees of each function.
var dvArity = [dvNumOps]int{dvEnd: 0, dvCut: 0, dvNop: 1, dvR: 2, dvL: 2, dvC: 2,
	dvSeries: 3, dvParallel: 4, dvVia: 3, dvAdd: 2, dvSub: 2, dvConst: 0, dvFlip: 1,
	dvV: 2, dvI: 2, dvD: 2, dvQ: 4, dvM: 4}

// isTransistor reports whether f inserts a three-terminal part.
func isTransistor(f uint8) bool { return f == dvQ || f == dvM }

// isArith reports whether a function belongs in a value (arithmetic)
// position rather than a construction-continuing one.
func isArith(f uint8) bool { return f == dvAdd || f == dvSub || f == dvConst }

func dvKind(f uint8) byte {
	switch f {
	case dvR:
		return 'R'
	case dvL:
		return 'L'
	case dvC:
		return 'C'
	case dvV:
		return 'V'
	case dvI:
		return 'I'
	case dvD:
		return 'D'
	case dvQ:
		return 'Q'
	case dvM:
		return 'M'
	}
	return 0
}

// subtreeEnd returns the index just past the subtree that starts at i, or
// len(prog)+1 if the program ends first.
func subtreeEnd(prog []Op, i int) int {
	need := 1
	j := i
	for need > 0 {
		if j >= len(prog) || prog[j].F >= dvNumOps {
			return len(prog) + 1
		}
		need += dvArity[prog[j].F] - 1
		j++
	}
	return j
}

// FormatProgram writes a program in Koza's LISP style, for logs and tests.
func FormatProgram(prog []Op) string {
	var b strings.Builder
	var walk func(i int) int
	walk = func(i int) int {
		if i >= len(prog) || prog[i].F >= dvNumOps {
			b.WriteString("?")
			return len(prog)
		}
		op := prog[i]
		switch {
		case op.F == dvConst:
			fmt.Fprintf(&b, "%.4g", op.V)
			return i + 1
		case dvArity[op.F] == 0:
			b.WriteString(dvName[op.F])
			return i + 1
		}
		b.WriteString("(" + dvName[op.F])
		if op.F == dvVia {
			fmt.Fprintf(&b, "%d", op.K)
		}
		j := i + 1
		for a := 0; a < dvArity[op.F]; a++ {
			b.WriteString(" ")
			j = walk(j)
		}
		b.WriteString(")")
		return j
	}
	walk(0)
	return b.String()
}

// validateDev checks the developmental fields against the embryo.
func (t *Topology) validateDev(embryo *circuit.Circuit) error {
	if t.Encoding > Developmental {
		return fmt.Errorf("%w: encoding %d", ErrBadTopology, t.Encoding)
	}
	if t.Encoding != Developmental {
		return nil
	}
	if len(t.Wires) == 0 {
		return fmt.Errorf("%w: a developmental encoding needs at least one modifiable wire", ErrBadTopology)
	}
	if t.MaxSize < 1 || t.MaxSize > 100000 || t.InitDepth < 1 || t.InitDepth > 30 {
		return fmt.Errorf("%w: max_size %d, init_depth %d", ErrBadTopology, t.MaxSize, t.InitDepth)
	}
	for _, w := range t.Wires {
		a, okA := t.devNode(w[0])
		b, okB := t.devNode(w[1])
		if !okA || !okB {
			return fmt.Errorf("%w: wire %s–%s: ends must be terminals or ground", ErrBadTopology, w[0], w[1])
		}
		if a == b {
			return fmt.Errorf("%w: wire %s–%s joins a node to itself", ErrBadTopology, w[0], w[1])
		}
	}
	// Development creates internal nodes without a fixed bound.
	for _, n := range embryo.NodeNames {
		if strings.HasPrefix(strings.ToLower(n), "g_n") {
			return fmt.Errorf("%w: the embryo uses node name %s, reserved for grown nodes", ErrBadTopology, n)
		}
	}
	return nil
}

// devNode maps a wire end to its node index (0 ground, 1..T terminals).
func (t *Topology) devNode(name string) (int, bool) {
	if name == "0" || strings.EqualFold(name, "gnd") {
		return 0, true
	}
	for i, n := range t.Terminals {
		if strings.EqualFold(n, name) {
			return i + 1, true
		}
	}
	return 0, false
}

// devKinds lists the component-creating functions of the configured kinds.
func (t *Topology) devKinds() []uint8 {
	var out []uint8
	// R, L, C first: a passive-only configuration draws exactly as before.
	for _, f := range []uint8{dvR, dvL, dvC, dvV, dvI, dvD, dvQ, dvM} {
		if _, ok := t.kind(dvKind(f)); ok {
			out = append(out, f)
		}
	}
	return out
}

// randomTree appends a random subtree (construction position unless arith)
// by Koza's grow method: below depth 1 only terminals are drawn.
func (t *Topology) randomTree(r *rand.Rand, prog []Op, arith bool, depth int) []Op {
	if arith {
		if depth <= 1 || r.Float64() < 0.5 {
			return append(prog, Op{F: dvConst, V: 2*r.Float64() - 1})
		}
		f := dvAdd
		if r.IntN(2) == 1 {
			f = dvSub
		}
		prog = append(prog, Op{F: f})
		prog = t.randomTree(r, prog, true, depth-1)
		return t.randomTree(r, prog, true, depth-1)
	}
	if depth <= 1 || r.Float64() < 0.2 {
		if r.Float64() < 0.85 {
			return append(prog, Op{F: dvEnd})
		}
		return append(prog, Op{F: dvCut})
	}
	kinds := t.devKinds()
	// Component-creating functions are drawn as often as all the
	// topology-modifying ones together.
	var op Op
	if r.IntN(2) == 0 {
		op.F = kinds[r.IntN(len(kinds))]
	} else {
		mods := []uint8{dvSeries, dvParallel, dvVia, dvVia, dvNop}
		if t.polarized() {
			mods = append(mods, dvFlip)
		}
		op.F = mods[r.IntN(len(mods))]
		if op.F == dvVia {
			op.K = r.IntN(1 + len(t.Terminals))
		}
	}
	prog = append(prog, op)
	if k := dvKind(op.F); k != 0 {
		prog = t.randomTree(r, prog, true, 3)
		for a := 1; a < dvArity[op.F]; a++ {
			prog = t.randomTree(r, prog, false, depth-1)
		}
		return prog
	}
	for a := 0; a < dvArity[op.F]; a++ {
		prog = t.randomTree(r, prog, false, depth-1)
	}
	return prog
}

// randomProgram draws one program of at most MaxSize ops, with depth
// ramped over 2..InitDepth (Koza's ramped initial population).
func (t *Topology) randomProgram(r *rand.Rand) []Op {
	for tries := 0; tries < 20; tries++ {
		d := 2
		if t.InitDepth > 2 {
			d += r.IntN(t.InitDepth - 1)
		}
		p := t.randomTree(r, nil, false, d)
		if len(p) <= t.MaxSize {
			return p
		}
	}
	return []Op{{F: dvEnd}}
}

func (t *Topology) devRandomGenome(r *rand.Rand) Genome {
	prog := make([][]Op, len(t.Wires))
	for i := range prog {
		prog[i] = t.randomProgram(r)
	}
	return Genome{Program: prog}
}

// pickPoint chooses a point of prog of the given type (arith or not),
// favouring functions over terminals 90 % to 10 % as Koza does; -1 if
// the program has no point of that type.
func pickPoint(r *rand.Rand, prog []Op, arith bool) int {
	var inner, leaf []int
	for i, op := range prog {
		if isArith(op.F) != arith {
			continue
		}
		if dvArity[op.F] > 0 {
			inner = append(inner, i)
		} else {
			leaf = append(leaf, i)
		}
	}
	switch {
	case len(inner) == 0 && len(leaf) == 0:
		return -1
	case len(leaf) == 0 || (len(inner) > 0 && r.Float64() < 0.9):
		return inner[r.IntN(len(inner))]
	}
	return leaf[r.IntN(len(leaf))]
}

// splice returns prog with its subtree at i replaced by sub.
func splice(prog []Op, i int, sub []Op) []Op {
	end := subtreeEnd(prog, i)
	out := make([]Op, 0, len(prog)-(end-i)+len(sub))
	out = append(out, prog[:i]...)
	out = append(out, sub...)
	return append(out, prog[end:]...)
}

// devCrossover swaps subtrees of the same type between the two parents'
// programs for one embryo wire (Koza's subtree crossover with point
// typing). A child longer than MaxSize is replaced by a copy of its
// parent.
func (t *Topology) devCrossover(r *rand.Rand, a, b Genome) (Genome, Genome) {
	c1, c2 := a.Clone(), b.Clone()
	if len(a.Program) == 0 || len(a.Program) != len(b.Program) {
		return c1, c2
	}
	w := r.IntN(len(a.Program))
	pa, pb := a.Program[w], b.Program[w]
	arith := r.Float64() < 0.2
	i, j := pickPoint(r, pa, arith), pickPoint(r, pb, arith)
	if i < 0 || j < 0 {
		arith = !arith
		i, j = pickPoint(r, pa, arith), pickPoint(r, pb, arith)
		if i < 0 || j < 0 {
			return c1, c2
		}
	}
	subA := pa[i:subtreeEnd(pa, i)]
	subB := pb[j:subtreeEnd(pb, j)]
	if n := splice(pa, i, subB); len(n) <= t.MaxSize {
		c1.Program[w] = n
	}
	if n := splice(pb, j, subA); len(n) <= t.MaxSize {
		c2.Program[w] = n
	}
	return c1, c2
}

// devMutate applies one mutation, then each further one with probability
// 0.3: half of them move a constant by sigma decades (when there is one),
// the rest replace a random subtree by a new random one of the same type.
func (t *Topology) devMutate(r *rand.Rand, g *Genome, sigma float64) {
	if len(g.Program) == 0 {
		return
	}
	for first := true; first || r.Float64() < 0.3; first = false {
		w := r.IntN(len(g.Program))
		prog := g.Program[w]
		if r.IntN(2) == 0 {
			if i := pickPoint(r, prog, true); i >= 0 && prog[i].F == dvConst {
				p := append([]Op(nil), prog...)
				p[i].V += sigma * r.NormFloat64()
				g.Program[w] = p
				continue
			}
		}
		arith := r.Float64() < 0.15
		i := pickPoint(r, prog, arith)
		if i < 0 {
			arith = !arith
			if i = pickPoint(r, prog, arith); i < 0 {
				continue
			}
		}
		for tries := 0; tries < 5; tries++ {
			n := splice(prog, i, t.randomTree(r, nil, arith, 4))
			if len(n) <= t.MaxSize {
				g.Program[w] = n
				break
			}
		}
	}
}

// validProgram reports whether prog is one well-formed, well-typed tree
// of construction type using only configured kinds.
func (t *Topology) validProgram(prog []Op) bool {
	if len(prog) == 0 || len(prog) > t.MaxSize || subtreeEnd(prog, 0) != len(prog) {
		return false
	}
	ok := true
	var walk func(i int, arith bool) int
	walk = func(i int, arith bool) int {
		op := prog[i]
		if isArith(op.F) != arith {
			ok = false
			return subtreeEnd(prog, i)
		}
		if k := dvKind(op.F); k != 0 {
			if _, has := t.kind(k); !has {
				ok = false
			}
			j := walk(i+1, true)
			for a := 1; a < dvArity[op.F]; a++ {
				j = walk(j, false)
			}
			return j
		}
		if op.F == dvConst && (math.IsNaN(op.V) || math.IsInf(op.V, 0)) {
			ok = false
		}
		if op.F == dvVia && (op.K < 0 || op.K > len(t.Terminals)) {
			ok = false
		}
		j := i + 1
		for a := 0; a < dvArity[op.F]; a++ {
			j = walk(j, arith)
		}
		return j
	}
	walk(0, false)
	return ok
}

// devRepair makes every genome decodable: one program per wire, each a
// valid tree; an invalid program becomes END (the wire stays a wire).
func (t *Topology) devRepair(g *Genome) {
	g.Parts = nil
	if len(g.Program) != len(t.Wires) {
		p := make([][]Op, len(t.Wires))
		copy(p, g.Program)
		g.Program = p
	}
	for i, prog := range g.Program {
		if !t.validProgram(prog) {
			g.Program[i] = []Op{{F: dvEnd}}
		}
	}
}

// dvComp is a component while the circuit develops: a wire (kind 0) or
// a part.
type dvComp struct {
	kind    byte
	a, b, c int // c: a transistor's base or gate
	x       float64
	cut     bool
}

// develop runs the programs from the embryo's wires. It returns the
// components and the number of nodes used (0 ground, 1..T terminals,
// then the grown internal nodes).
func (t *Topology) develop(g Genome) ([]dvComp, int) {
	type head struct{ comp, prog, pos int }
	var comps []dvComp
	var queue []head
	for w, ends := range t.Wires {
		a, _ := t.devNode(ends[0])
		b, _ := t.devNode(ends[1])
		comps = append(comps, dvComp{a: a, b: b})
		if w < len(g.Program) && len(g.Program[w]) > 0 {
			queue = append(queue, head{comp: w, prog: w, pos: 0})
		}
	}
	next := 1 + len(t.Terminals)
	parts := 0
	// Bound on all components (wires included), so a program of MaxSize
	// ops cannot build an unbounded circuit.
	maxComps := len(t.Wires) + 8*t.MaxParts
	newNode := func() int { next++; return next - 1 }
	add := func(c dvComp) int {
		comps = append(comps, c)
		if c.kind != 0 {
			parts++
		}
		return len(comps) - 1
	}
	for len(queue) > 0 {
		h := queue[0]
		queue = queue[1:]
		prog := g.Program[h.prog]
		if h.pos >= len(prog) {
			continue
		}
		op := prog[h.pos]
		c := h.comp
		push := func(comp, pos int) { queue = append(queue, head{comp, h.prog, pos}) }
		// children returns the start of each of op's n subtrees.
		children := func(n int) []int {
			s := make([]int, n)
			j := h.pos + 1
			for a := 0; a < n; a++ {
				s[a] = j
				j = subtreeEnd(prog, j)
			}
			return s
		}
		switch op.F {
		case dvEnd:
		case dvCut:
			comps[c].cut = true
		case dvNop:
			push(c, h.pos+1)
		case dvFlip:
			comps[c].a, comps[c].b = comps[c].b, comps[c].a
			push(c, h.pos+1)
		case dvR, dvL, dvC, dvV, dvI, dvD:
			v, j := evalArith(prog, h.pos+1)
			k := dvKind(op.F)
			if comps[c].kind == 0 && parts >= t.MaxParts {
				push(c, j) // part budget spent: the wire stays a wire
				break
			}
			if comps[c].kind == 0 {
				parts++
			}
			comps[c].kind = k
			comps[c].x = t.devValue(k, v)
			push(c, j)
		case dvQ, dvM:
			v, j := evalArith(prog, h.pos+1)
			s1 := j
			s2 := subtreeEnd(prog, s1)
			s3 := subtreeEnd(prog, s2)
			wasPart := comps[c].kind != 0
			if len(comps)+3 > maxComps || (!wasPart && parts >= t.MaxParts) {
				push(c, s1) // budget spent: the component goes on unchanged
				break
			}
			// Z(a→b) becomes the collector lead a→n1; the transistor sits
			// between n1 (collector), n2 (emitter) and n3 (base).
			z := comps[c]
			if wasPart {
				parts-- // Z is replaced by a lead wire
			}
			n1, n2, n3 := newNode(), newNode(), newNode()
			comps[c] = dvComp{a: z.a, b: n1}
			k := dvKind(op.F)
			add(dvComp{kind: k, a: n1, b: n2, c: n3, x: t.devValue(k, v)})
			we := add(dvComp{a: n2, b: z.b})
			wb := add(dvComp{a: n3, b: z.a})
			push(c, s1)
			push(we, s2)
			push(wb, s3)
		case dvSeries, dvParallel, dvVia:
			ch := children(dvArity[op.F])
			copies := comps[c].kind != 0
			if len(comps)+3 > maxComps || (copies && parts >= t.MaxParts) {
				push(c, ch[0]) // budget spent: only the first branch goes on
				break
			}
			z := comps[c]
			z.cut = false
			switch op.F {
			case dvSeries:
				n1, n2 := newNode(), newNode()
				comps[c].b = n1
				w := add(dvComp{a: n1, b: n2})
				z.a = n2 // z.b is still the old far end
				zc := add(z)
				push(c, ch[0])
				push(w, ch[1])
				push(zc, ch[2])
			case dvParallel:
				n1, n2 := newNode(), newNode()
				w1 := add(dvComp{a: z.a, b: n1})
				zb := z.b
				z.a, z.b = n1, n2
				zc := add(z)
				w2 := add(dvComp{a: n2, b: zb})
				push(c, ch[0])
				push(w1, ch[1])
				push(zc, ch[2])
				push(w2, ch[3])
			case dvVia:
				m := newNode()
				comps[c].b = m
				z.a = m
				zc := add(z)
				w := add(dvComp{a: m, b: op.K})
				push(c, ch[0])
				push(zc, ch[1])
				push(w, ch[2])
			}
		}
	}
	return comps, next
}

// devValue maps an arithmetic result v [decades from the middle of the
// kind's range] to the clamped log10 value.
func (t *Topology) devValue(k byte, v float64) float64 {
	pk, _ := t.kind(k)
	lo, hi := math.Log10(pk.Min), math.Log10(pk.Max)
	return math.Min(math.Max((lo+hi)/2+v, lo), hi)
}

// evalArith evaluates the value subtree at i and returns its result and
// the index just past it.
func evalArith(prog []Op, i int) (float64, int) {
	if i >= len(prog) {
		return 0, i
	}
	switch prog[i].F {
	case dvAdd, dvSub:
		a, j := evalArith(prog, i+1)
		b, k := evalArith(prog, j)
		if prog[i].F == dvSub {
			b = -b
		}
		return a + b, k
	case dvConst:
		return prog[i].V, i + 1
	}
	return 0, subtreeEnd(prog, i)
}

// devParts grows the circuit and reduces it to parts on t's node indices
// (internal nodes renumbered compactly) and shorts between embryo nodes.
func (t *Topology) devParts(g Genome) (parts []Part, shorts [][2]int, internal int) {
	comps, n := t.develop(g)
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	// Merge wire ends; the lower index represents the class, so ground and
	// the terminals keep their own numbers.
	for _, c := range comps {
		if c.kind == 0 && !c.cut {
			ra, rb := find(c.a), find(c.b)
			if ra < rb {
				parent[rb] = ra
			} else if rb < ra {
				parent[ra] = rb
			}
		}
	}
	term := 1 + len(t.Terminals)
	for k := 1; k < term; k++ {
		if r := find(k); r != k {
			shorts = append(shorts, [2]int{r, k})
		}
	}
	renum := map[int]int{}
	node := func(k int) int {
		r := find(k)
		if r < term {
			return r
		}
		if m, ok := renum[r]; ok {
			return m
		}
		renum[r] = term + len(renum)
		return renum[r]
	}
	for _, c := range comps {
		if c.kind != 0 && !c.cut {
			p := Part{Kind: c.kind, A: node(c.a), B: node(c.b), X: c.x}
			if threeTerminal(c.kind) {
				p.C = node(c.c)
			}
			parts = append(parts, p)
		}
	}
	return parts, shorts, len(renum)
}

// devDecode builds the developed circuit.
func (t *Topology) devDecode(embryo *circuit.Circuit, g Genome) (*circuit.Circuit, []float64) {
	parts, shorts, internal := t.devParts(g)
	tt := *t
	tt.Internal = internal
	return tt.emit(embryo, tt.prune(parts), shorts)
}
