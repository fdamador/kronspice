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
// Two encodings share this type. Direct (here): the genome is the netlist
// itself, each Part naming a kind, two nodes and a value. Developmental
// (encode_develop.go): Koza's circuit-constructing program, which grows
// the parts from the embryo's modifiable wires. In the direct encoding,
// nodes are ground, the embryo's terminals (for example "in" and "out")
// and a pool of internal nodes; operators add, delete, rewire, re-kind and
// re-value parts, split a part into two in series through a new node, and
// cross two netlists over at one point.

// PartKind is a kind of part evolution may place, with its value range.
//
// Passive kinds R, C, L take their value in Ω, F, H. Sources V and I take
// a DC magnitude [V], [A]; the part's orientation gives its polarity. D, Q
// and M name a .model of the embryo (Model): D and Q take an area factor
// as value (Min = Max = 1 for the model as given), M takes the channel
// width W [m] with the length Length [m] (0: the model default). A Q part
// is collector A, emitter B, base C; an M part drain A, source B, gate C,
// with its bulk tied to the source.
type PartKind struct {
	Kind     circuit.Kind // R, C, L, V, I, D, Q or M
	Min, Max float64      // value range (0 < Min ≤ Max)
	// Series snaps values to a preferred-number series ("" for none).
	Series string
	// KappaA and KappaB give every part of this kind the hypothetical ECE
	// scalar terms κ_a, κ_b at its two terminals (element parameters
	// kappa_a, kappa_b; circuit/ece.go, register entry H48), so the body
	// sees (1 + s·κ_a)·φ_a − (1 + s·κ_b)·φ_b. They need a .ece declaration
	// in the embryo. With κ_a ≠ κ_b the part is polarized: its orientation
	// changes the circuit, and FLIP (encode_develop.go) and the direct
	// encoding's flip mutation can reverse it. With ECE2 (s = −2) a κ above
	// ½ inverts that terminal's potential, a negative conductance whose
	// power is booked as model-supplied, never as inflow.
	KappaA, KappaB float64
	// Model names the embryo's .model card (D, Q and M only).
	Model string
	// Length is a MOSFET part's channel length L [m] (0: the model default).
	Length float64
}

// threeTerminal reports whether parts of kind k use the third node C.
func threeTerminal(k byte) bool { return k == byte(circuit.BJT) || k == byte(circuit.MOSFET) }

// polarKind reports whether a part of kind k has a polarity of its own
// (sources, diodes, transistors).
func polarKind(k circuit.Kind) bool {
	switch k {
	case circuit.VoltageSource, circuit.CurrentSource, circuit.Diode, circuit.BJT, circuit.MOSFET:
		return true
	}
	return false
}

// polarized reports whether any part kind's orientation matters: a source,
// diode or transistor, or a part with κ_a ≠ κ_b. For plain R, L and C
// (Maxwell) it does not: reversing one changes nothing.
func (t *Topology) polarized() bool {
	for _, pk := range t.Kinds {
		if pk.KappaA != pk.KappaB || polarKind(pk.Kind) {
			return true
		}
	}
	return false
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

	// Encoding is Direct (the default: the genome is the part list) or
	// Developmental (the genome is a Koza circuit-constructing program;
	// encode_develop.go). Internal and InitParts apply to Direct only.
	Encoding Encoding
	// Wires are the embryo's modifiable wires (Developmental): each joins
	// two nodes (terminals or "0") and carries one writing head, developed
	// by its own program (Koza's result-producing branch).
	Wires [][2]string
	// MaxSize bounds each program's length in functions and terminals
	// (Developmental; default 200).
	MaxSize int
	// InitDepth is the largest depth of the initial random programs
	// (Developmental; default 6).
	InitDepth int
	// Library is the primitive network of the kron (connection) encoding:
	// one slot per branch, each of a listed kind (encode_kron.go). Empty:
	// MaxParts slots dealt round-robin over Kinds.
	Library []circuit.Kind
}

// Part is one gene of a topology genome.
type Part struct {
	Kind byte    `json:"k"` // 'R', 'C', 'L', 'V', 'I', 'D', 'Q', 'M'
	A    int     `json:"a"` // node: 0 ground, 1..T terminals, T+1.. internal
	B    int     `json:"b"`
	X    float64 `json:"x"` // log10 of the magnitude
	// C is the third node of a transistor (base or gate); 0 otherwise.
	C int `json:"c,omitempty"`
	// Off marks an open slot of the kron (connection) encoding: its
	// junction pair is inactive and the part is left out of the circuit.
	Off bool `json:"off,omitempty"`
}

// terminals lists the nodes a part touches.
func (p Part) terminals() []int {
	if threeTerminal(p.Kind) {
		return []int{p.A, p.B, p.C}
	}
	return []int{p.A, p.B}
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
	if t.MaxSize == 0 {
		t.MaxSize = 200
	}
	if t.InitDepth == 0 {
		t.InitDepth = 6
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
	if err := t.validateDev(embryo); err != nil {
		return err
	}
	if len(t.Kinds) == 0 {
		return fmt.Errorf("%w: no part kinds", ErrBadTopology)
	}
	seen := map[circuit.Kind]bool{}
	for _, pk := range t.Kinds {
		switch pk.Kind {
		case circuit.Resistor, circuit.Capacitor, circuit.Inductor, circuit.VoltageSource, circuit.CurrentSource:
			if pk.Model != "" {
				return fmt.Errorf("%w: part kind %v takes no model", ErrBadTopology, pk.Kind)
			}
		case circuit.Diode, circuit.BJT, circuit.MOSFET:
			m, ok := embryo.Models[strings.ToLower(pk.Model)]
			want := map[circuit.Kind][]string{circuit.Diode: {"d"}, circuit.BJT: {"npn", "pnp"}, circuit.MOSFET: {"nmos", "pmos"}}[pk.Kind]
			if !ok || (m.Type != want[0] && (len(want) < 2 || m.Type != want[1])) {
				return fmt.Errorf("%w: part kind %v needs a .model of type %s in the embryo, got %q", ErrBadTopology,
					pk.Kind, strings.Join(want, " or "), pk.Model)
			}
		default:
			return fmt.Errorf("%w: part kind %v (R, C, L, V, I, D, Q or M)", ErrBadTopology, pk.Kind)
		}
		if pk.Length < 0 || math.IsNaN(pk.Length) || math.IsInf(pk.Length, 0) || (pk.Length != 0 && pk.Kind != circuit.MOSFET) {
			return fmt.Errorf("%w: %v length %g (MOSFET parts only)", ErrBadTopology, pk.Kind, pk.Length)
		}
		if (pk.KappaA != 0 || pk.KappaB != 0) && pk.Kind != circuit.Resistor && pk.Kind != circuit.Capacitor && pk.Kind != circuit.Inductor {
			return fmt.Errorf("%w: %v: kappa_a, kappa_b are for R, C and L", ErrBadTopology, pk.Kind)
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
		if math.IsNaN(pk.KappaA) || math.IsInf(pk.KappaA, 0) || math.IsNaN(pk.KappaB) || math.IsInf(pk.KappaB, 0) {
			return fmt.Errorf("%w: %v kappa_a, kappa_b must be finite", ErrBadTopology, pk.Kind)
		}
		if (pk.KappaA != 0 || pk.KappaB != 0) && (embryo.ECE == nil || embryo.ECE.Form == circuit.ECENone) {
			return fmt.Errorf("%w: %v kappa_a, kappa_b need a .ece declaration (form ece or ece2) in the embryo", ErrBadTopology, pk.Kind)
		}
	}
	return t.validateKron()
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
		// A two-terminal part on one node is a self-loop; a transistor is
		// one only with all three terminals on the same node.
		keep[i] = !p.Off && (p.A != p.B || (threeTerminal(p.Kind) && p.C != p.A))
	}
	deg := make([]int, t.nodes())
	for changed := true; changed; {
		changed = false
		for i := range deg {
			deg[i] = 0
		}
		for i, p := range parts {
			if keep[i] {
				for _, n := range p.terminals() {
					deg[n]++
				}
			}
		}
		for i, p := range parts {
			if !keep[i] {
				continue
			}
			for _, n := range p.terminals() {
				if n >= term && deg[n] == 1 {
					keep[i], changed = false, true
					break
				}
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
// A developmental genome (Encoding Developmental) is first grown into its
// parts by running its program (encode_develop.go).
func (t *Topology) Decode(embryo *circuit.Circuit, g Genome) (*circuit.Circuit, []float64) {
	if t.Encoding == Developmental {
		return t.devDecode(embryo, g)
	}
	return t.emit(embryo, t.prune(g.Parts), nil)
}

// emit adds parts (already pruned) and shorts to a copy of the embryo.
// Node indices are those of t: 0 ground, 1..T terminals, then internal
// nodes g_n1, …. Each short is a 0 V source V_w1, V_w2, … joining two
// embryo nodes (a wire the developmental encoding left in place).
func (t *Topology) emit(embryo *circuit.Circuit, parts []Part, shorts [][2]int) (*circuit.Circuit, []float64) {
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
	vals := make([]float64, len(parts))
	for i, p := range parts {
		v := t.value(p)
		vals[i] = v
		pk, _ := t.kind(p.Kind)
		e := circuit.Element{
			Name: fmt.Sprintf("%c_g%d", p.Kind, i+1), Kind: circuit.Kind(p.Kind),
			Nodes: [2]int{node(p.A), node(p.B)}, Value: v,
		}
		switch circuit.Kind(p.Kind) {
		case circuit.VoltageSource, circuit.CurrentSource:
			e.HasDC = true // DC value v from A (+) to B
		case circuit.Diode:
			e.Value, e.Model, e.Area = 0, strings.ToLower(pk.Model), v
		case circuit.BJT: // collector A, emitter B, base C
			e.Value, e.Model, e.Area = 0, strings.ToLower(pk.Model), v
			e.Extra = []int{node(p.C)}
		case circuit.MOSFET: // drain A, source B, gate C, bulk on the source
			e.Value, e.Model = 0, strings.ToLower(pk.Model)
			e.Extra = []int{node(p.C), node(p.B)}
			e.Params = map[string]float64{"w": v}
			if pk.Length > 0 {
				e.Params["l"] = pk.Length
			}
		}
		// Terminal A of the part is the element's first node, so a flipped
		// part carries κ_a at what was its B end.
		if pk.KappaA != 0 || pk.KappaB != 0 {
			e.Params = map[string]float64{"kappa_a": pk.KappaA, "kappa_b": pk.KappaB}
		}
		c.Elements = append(c.Elements, e)
	}
	for i, w := range shorts {
		c.Elements = append(c.Elements, circuit.Element{
			Name: fmt.Sprintf("V_w%d", i+1), Kind: circuit.VoltageSource, HasDC: true,
			Nodes: [2]int{node(w[0]), node(w[1])},
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
	p := Part{Kind: byte(pk.Kind), A: a, B: b, X: lo + r.Float64()*(hi-lo)}
	if threeTerminal(p.Kind) {
		p.C = pool[r.IntN(len(pool))]
		for tries := 0; (p.C == a || p.C == b) && tries < 8; tries++ {
			p.C = pool[r.IntN(len(pool))]
		}
	}
	return p
}

// nodePool lists ground, the terminals, the internal nodes in use and the
// lowest unused internal node.
func (t *Topology) nodePool(parts []Part) []int {
	term := 1 + len(t.Terminals)
	used := make([]bool, t.nodes())
	for _, p := range parts {
		for _, n := range p.terminals() {
			used[n] = true
		}
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
	if t.Encoding == Developmental {
		return t.devRandomGenome(r)
	}
	if t.Encoding == Connection {
		return t.kronRandomGenome(r)
	}
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
	if t.Encoding == Developmental {
		return t.devCrossover(r, a, b)
	}
	if t.Encoding == Connection {
		return t.kronCrossover(r, a, b)
	}
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
	// opFlip reverses a part's terminals. It is outside the weight table and
	// drawn only when a part kind is polarized, so for plain R, L, C the
	// operator sequence (and every earlier run) is unchanged.
	opFlip = numOps
)

var opWeight = [numOps]float64{opValue: 0.35, opAdd: 0.15, opDelete: 0.12, opRewire: 0.13, opKind: 0.05,
	opSeries: 0.1, opParallel: 0.1}

// mutate applies one operator, then each further one with probability
// 0.3 (a geometric number of operators, mean 1.43). sigma is the value
// step in decades.
func (t *Topology) mutate(r *rand.Rand, g *Genome, sigma float64) {
	if t.Encoding == Developmental {
		t.devMutate(r, g, sigma)
		return
	}
	if t.Encoding == Connection {
		t.kronMutate(r, g, sigma)
		return
	}
	pol := t.polarized()
	for first := true; first || r.Float64() < 0.3; first = false {
		if pol && r.Float64() < 0.1 {
			t.apply(r, g, opFlip, sigma)
			continue
		}
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
		ends := 2
		if threeTerminal(p.Kind) {
			ends = 3
		}
		switch r.IntN(ends) {
		case 0:
			p.A = pool[r.IntN(len(pool))]
		case 1:
			p.B = pool[r.IntN(len(pool))]
		default:
			p.C = pool[r.IntN(len(pool))]
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
		if threeTerminal(p.Kind) && !threeTerminal(byte(old.Kind)) {
			pool := t.nodePool(g.Parts)
			p.C = pool[r.IntN(len(pool))]
		}
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
	case opFlip:
		p := &g.Parts[r.IntN(n)]
		p.A, p.B = p.B, p.A
	}
}

// repair clamps values into range, maps unknown kinds and out-of-range
// nodes back into the configuration, so every genome decodes.
func (t *Topology) repair(g *Genome) {
	if t.Encoding == Developmental {
		t.devRepair(g)
		return
	}
	if t.Encoding == Connection {
		t.kronRepair(g)
	}
	nn := t.nodes()
	for i := range g.Parts {
		p := &g.Parts[i]
		if t.Encoding != Connection {
			p.Off = false
		}
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
		p.C = ((p.C % nn) + nn) % nn
		if !threeTerminal(p.Kind) {
			p.C = 0
		}
	}
	if t.Encoding == Connection {
		t.kronCanonical(g)
	} else if len(g.Parts) > t.MaxParts {
		g.Parts = g.Parts[:t.MaxParts]
	}
}
