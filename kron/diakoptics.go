// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"sync"
	"sync/atomic"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/solver"
)

// Diakoptics: Kron's piecewise solution by tearing (G. Kron, Diakoptics,
// Macdonald 1963; H. H. Happ, Diakoptics and Networks, 1971, Ch. XII —
// Happ's own factorized algorithms are in happ.go).
//
// A set of two-terminal branches is cut out of the network. What remains
// falls apart into pieces, and in the node-pair frame its matrix B is
// block-diagonal, one block per piece. The torn branches put back a
// low-rank term, A = B + U·D·Uᵀ, where D holds their admittances and each
// column of U is the incidence (+1, −1) of one torn branch. Then
//
//	z = B⁻¹·c                      every piece on its own
//	(I + D·Uᵀ·B⁻¹·U)·w = D·Uᵀ·z    the interconnection system, order t
//	x = z − B⁻¹·U·w
//
// and w is exactly the vector of torn-branch currents: the unknowns of the
// interconnection are Kron's "intersection" currents, which flow through
// the open paths that tearing created (CLAUDE.md §3a.1; Kron 1939 p. 133,
// parts analysed separately and recombined by a C). This is the matrix
// form of Aitchison's algorithm (P. W. Aitchison, "Diakoptics as a general
// approach in engineering," J. Eng. Math. 1983, Sec. 2) with the coupling
// entries carried by U and D; the interconnection matrix is non-singular
// exactly when A is (his Lemma 1), but every block of B must be too (Lemma
// 2): a piece with no path to the datum floats, and is reported.
//
// B⁻¹ is applied block by block, and each block's factors and its columns
// B⁻¹·U are cached by the block's exact content, so when a genetic
// algorithm changes one part of a circuit only the piece holding it is
// refactored (CLAUDE.md §3a, caching for the GA).

// Tearing names the elements to cut: resistors, capacitors or inductors,
// each torn as a whole branch. An inductor cannot be torn for the
// operating point, where it is a short circuit.
type Tearing struct {
	Torn []string
}

// Diakoptics is a network prepared for torn solution.
type Diakoptics struct {
	n     *Network
	torn  []int  // torn branches, in Tearing order
	mask  []bool // branch → torn
	cache *TornCache
}

// TearStats describes one torn solve.
type TearStats struct {
	// Pieces is the number of diagonal blocks of B; Sizes their orders.
	Pieces int
	Sizes  []int
	// Interconnection is the order of the interconnection system: the
	// torn branches with non-zero admittance at this frequency.
	Interconnection int
	// Hits and Misses count pieces whose factors came from the cache or
	// were computed.
	Hits, Misses int
	// Refinements is the number of refinement steps taken.
	Refinements int
	// BackwardError is the componentwise backward error of the solution
	// in the untorn system, max |A·x − c|ᵢ / (|A|·|x| + |c|)ᵢ.
	BackwardError float64
}

// Refinement limits: at most maxRefine passes, stopping below refineDone;
// a solve whose backward error stays above refineAccept is rejected.
const (
	maxRefine    = 12
	refineDone   = 1e-14
	refineAccept = 1e-12
)

// Tear builds the network of a linear circuit and checks the cut. The
// cache may be nil (nothing is reused between solves) and may be shared by
// any number of goroutines.
func Tear(c *circuit.Circuit, t Tearing, cache *TornCache) (*Diakoptics, error) {
	n, err := Build(c)
	if err != nil {
		return nil, err
	}
	d := &Diakoptics{n: n, mask: make([]bool, len(n.branches)), cache: cache}
	for _, name := range t.Torn {
		k, err := n.branchOf(name)
		if err != nil {
			return nil, err
		}
		switch n.branches[k].kind {
		case bR, bC, bL:
		default:
			return nil, fmt.Errorf("%w: tearing %s: only R, C and L branches can be torn", ErrUnsupported, name)
		}
		if n.branches[k].kind == bR && n.branches[k].value == 0 {
			return nil, fmt.Errorf("%w: tearing %s: a zero resistance has no admittance", ErrUnsupported, name)
		}
		if d.mask[k] {
			return nil, fmt.Errorf("%w: %s is torn twice", ErrUnsupported, name)
		}
		for _, ch := range n.chords {
			if ch.a == k || ch.b == k {
				return nil, fmt.Errorf("%w: tearing %s: it carries a chord (mutual coupling); tear both coils' piece elsewhere", ErrUnsupported, name)
			}
		}
		for _, r := range n.radiation {
			if r.elem == n.branches[k].elem {
				return nil, fmt.Errorf("%w: tearing %s: it has a radiation edge", ErrUnsupported, name)
			}
		}
		d.mask[k] = true
		d.torn = append(d.torn, k)
	}
	return d, nil
}

// Network returns the untorn network (for comparison and the ledger).
func (d *Diakoptics) Network() *Network { return d.n }

// Counts returns the frame dimensions of the torn network: the torn
// branches removed, S counting the sub-networks that remain (Kron 1939
// p. 75, eqs 3.1–3.2).
func (d *Diakoptics) Counts() Counts {
	var from, to []int
	for k, b := range d.n.branches {
		if !d.mask[k] {
			from, to = append(from, b.from), append(to, b.to)
		}
	}
	g := newGraph(d.n.nodes, from, to)
	b, N, S := len(from), g.n, g.subNetworks()
	return Counts{Branches: b, Junctions: N, SubNetworks: S, Meshes: b - N + S, JunctionPairs: N - S}
}

// tornAdmittance is a torn branch's admittance at omega.
func (d *Diakoptics) tornAdmittance(k int, omega float64, dc bool) (complex128, error) {
	b := d.n.branches[k]
	switch b.kind {
	case bR:
		return complex(1/b.value, 0), nil
	case bC:
		if dc {
			return 0, nil
		}
		return complex(0, omega*b.value), nil
	case bL:
		if dc || omega == 0 {
			return 0, fmt.Errorf("%w: torn inductor %s is a short circuit at the operating point", ErrUnsupported, d.n.c.Elements[b.elem].Name)
		}
		return 1 / complex(0, omega*b.value), nil
	}
	return 0, fmt.Errorf("%w: torn branch %d", ErrUnsupported, k)
}

// piece is one diagonal block of B in a solve.
type piece struct {
	idx   []int    // global unknowns, ascending
	ports []int    // local positions of the torn-branch terminals in it
	block []lentry // its entries in local indices, row-major
	e     *pieceFactor
}

// lentry is a block entry in local indices.
type lentry struct {
	i, j int
	v    complex128
}

// Solve solves the torn network at angular frequency omega (dc for the
// operating point), in the node-pair frame.
func (d *Diakoptics) Solve(omega float64, dc bool) (*Solution, TearStats, error) {
	var st TearStats
	n := d.n
	rel := n.relations(omega, dc)
	y := make([]complex128, len(d.torn))
	for p, k := range d.torn {
		v, err := d.tornAdmittance(k, omega, dc)
		if err != nil {
			return nil, st, err
		}
		y[p] = v
		rel[k] = relation{i: []term{{k, 1}}} // cut: carries no current in B
	}
	ns, err := n.assembleNode(rel, d.mask)
	if err != nil {
		return nil, st, err
	}
	size := ns.nv + ns.na
	B := ns.rows
	// The active torn branches (non-zero admittance) and their terminals.
	var act []int // index into d.torn
	for p := range d.torn {
		if y[p] != 0 {
			act = append(act, p)
		}
	}
	t := len(act)
	st.Interconnection = t
	isPort := make([]bool, size)
	for _, p := range act {
		k := d.torn[p]
		for _, c := range [2]int{ns.from[k], ns.to[k]} {
			if c >= 0 {
				isPort[c] = true
			}
		}
	}
	// Pieces: connected components of B's pattern (union by least index,
	// so the numbering is deterministic).
	parent := make([]int, size)
	for i := range parent {
		parent[i] = i
	}
	find := func(i int) int {
		for parent[i] != i {
			parent[i] = parent[parent[i]]
			i = parent[i]
		}
		return i
	}
	for i, row := range B {
		for _, e := range row {
			if e.v != 0 && e.j != i {
				if a, b := find(i), find(e.j); a != b {
					if a < b {
						parent[b] = a
					} else {
						parent[a] = b
					}
				}
			}
		}
	}
	pieceOf := make([]int, size)
	local := make([]int, size)
	root := make([]int, size)
	for i := range root {
		root[i] = -1
	}
	var pieces []*piece
	for i := 0; i < size; i++ {
		r := find(i)
		if root[r] < 0 {
			root[r] = len(pieces)
			pieces = append(pieces, &piece{})
		}
		pc := pieces[root[r]]
		pieceOf[i], local[i] = root[r], len(pc.idx)
		if isPort[i] {
			pc.ports = append(pc.ports, len(pc.idx))
		}
		pc.idx = append(pc.idx, i)
	}
	for i, row := range B {
		pc := pieces[pieceOf[i]]
		for _, e := range row {
			if e.v != 0 {
				pc.block = append(pc.block, lentry{local[i], local[e.j], e.v})
			}
		}
	}
	st.Pieces = len(pieces)
	for _, pc := range pieces {
		st.Sizes = append(st.Sizes, len(pc.idx))
		e, hit, err := d.factorPiece(pc)
		if err != nil {
			return nil, st, fmt.Errorf("%w (Aitchison 1983, Lemma 2: every torn piece needs its own reference; move the cut)", d.singularPiece(ns, pc, err))
		}
		if hit {
			st.Hits++
		} else {
			st.Misses++
		}
		pc.e = e
	}
	// B⁻¹ applied piece by piece; a piece with a zero right-hand side has a
	// zero solution and is skipped.
	bl := make([]complex128, size)
	binv := func(r []complex128) ([]complex128, error) {
		z := make([]complex128, size)
		for _, pc := range pieces {
			nz := false
			b := bl[:len(pc.idx)]
			for li, g := range pc.idx {
				b[li] = r[g]
				nz = nz || r[g] != 0
			}
			if !nz {
				continue
			}
			xl, err := pc.e.solve(b)
			if err != nil {
				return nil, err
			}
			for li, g := range pc.idx {
				z[g] = xl[li]
			}
		}
		return z, nil
	}
	// colAt(c, j): entry j of B⁻¹·e_c for a port c (zero across pieces).
	portSlot := make([]int, size)
	for _, pc := range pieces {
		for sl, li := range pc.ports {
			portSlot[pc.idx[li]] = sl
		}
	}
	colAt := func(c, j int) complex128 {
		if c < 0 || j < 0 || pieceOf[c] != pieceOf[j] {
			return 0
		}
		return pieces[pieceOf[c]].e.w[portSlot[c]][local[j]]
	}
	// The interconnection matrix I + D·Uᵀ·B⁻¹·U, factored once.
	var ic *pieceFactor
	if t > 0 {
		M := make([]lentry, 0, t*t)
		for a, p := range act {
			ka := d.torn[p]
			for b, q := range act {
				kb := d.torn[q]
				g := colAt(ns.from[kb], ns.from[ka]) - colAt(ns.to[kb], ns.from[ka]) -
					colAt(ns.from[kb], ns.to[ka]) + colAt(ns.to[kb], ns.to[ka])
				v := y[p] * g
				if a == b {
					v++
				}
				if v != 0 {
					M = append(M, lentry{a, b, v})
				}
			}
		}
		ic, err = factorBlock(t, M, nil)
		if err != nil {
			return nil, st, fmt.Errorf("%w: interconnection system of the torn network: %v (singular exactly when the untorn network is, Aitchison 1983 Lemma 1; otherwise a piece is nearly floating)", ErrSingular, err)
		}
	}
	pot := func(x []complex128, c int) complex128 {
		if c < 0 {
			return 0
		}
		return x[c]
	}
	// apply is the torn operator: z = B⁻¹r, w = (I + D·G)⁻¹·D·Uᵀz (the
	// torn-branch currents), x = z − B⁻¹U·w.
	apply := func(r []complex128) ([]complex128, error) {
		z, err := binv(r)
		if err != nil {
			return nil, err
		}
		if t == 0 {
			return z, nil
		}
		h := make([]complex128, t)
		for a, p := range act {
			k := d.torn[p]
			h[a] = y[p] * (pot(z, ns.from[k]) - pot(z, ns.to[k]))
		}
		w, err := ic.solve(h)
		if err != nil {
			return nil, err
		}
		for b, q := range act {
			if w[b] == 0 {
				continue
			}
			k := d.torn[q]
			for _, c := range [2]int{ns.from[k], ns.to[k]} {
				if c < 0 {
					continue
				}
				sgn := w[b]
				if c == ns.to[k] {
					sgn = -sgn
				}
				pc := pieces[pieceOf[c]]
				col := pc.e.w[portSlot[c]]
				for li, g := range pc.idx {
					z[g] -= sgn * col[li]
				}
			}
		}
		return z, nil
	}
	x, err := apply(ns.rhs)
	if err != nil {
		return nil, st, err
	}
	// Iterative refinement against the whole A = B + U·D·Uᵀ. The torn
	// operator is exact in exact arithmetic but loses accuracy when a piece
	// is nearly floating (a small capacitance its only reference at low
	// frequency): B is then ill-conditioned even though A is not.
	// Refinement recovers the accuracy while the error contracts; the
	// componentwise backward error decides, and a solve that does not
	// reach it is reported instead of returned.
	r := make([]complex128, size)
	den := make([]float64, size)
	berr := math.Inf(1)
	for pass := 0; pass <= maxRefine; pass++ {
		for i, row := range B {
			s, m := ns.rhs[i], cmplx.Abs(ns.rhs[i])
			for _, e := range row {
				p := e.v * x[e.j]
				s -= p
				m += cmplx.Abs(p)
			}
			r[i], den[i] = s, m
		}
		for _, p := range act {
			k := d.torn[p]
			f, o := ns.from[k], ns.to[k]
			i := y[p] * (pot(x, f) - pot(x, o))
			a := cmplx.Abs(y[p]) * (cmplx.Abs(pot(x, f)) + cmplx.Abs(pot(x, o)))
			if f >= 0 {
				r[f] -= i
				den[f] += a
			}
			if o >= 0 {
				r[o] += i
				den[o] += a
			}
		}
		berr = 0
		for i := range r {
			if den[i] > 0 {
				berr = math.Max(berr, cmplx.Abs(r[i])/den[i])
			} else if r[i] != 0 {
				berr = math.Inf(1)
			}
		}
		if berr <= refineDone || pass == maxRefine {
			break
		}
		dx, err := apply(r)
		if err != nil {
			return nil, st, err
		}
		for j := range x {
			x[j] += dx[j]
		}
		st.Refinements++
	}
	if !(berr <= refineAccept) {
		return nil, st, fmt.Errorf("%w: torn solve did not converge (backward error %.2g): a torn piece is nearly floating, so its block is ill-conditioned (Aitchison 1983, Lemma 2); move the cut", ErrSingular, berr)
	}
	st.BackwardError = berr
	for p, k := range d.torn {
		ns.yrow[k] = []term{{k, y[p]}}
	}
	return ns.finish(x), st, nil
}

// singularPiece names the nodes of a piece whose block cannot be factored.
func (d *Diakoptics) singularPiece(ns *nodeSystem, pc *piece, err error) error {
	var names []string
	for _, g := range pc.idx {
		if g < ns.nv && g+1 < len(d.n.c.NodeNames) {
			names = append(names, d.n.c.NodeNames[g+1])
		}
		if len(names) == 8 {
			names = append(names, "…")
			break
		}
	}
	return fmt.Errorf("%w: torn piece {%s} (%d unknowns) is singular: %v", ErrSingular, strings.Join(names, " "), len(pc.idx), err)
}

// pieceFactor is a factored block of B with its columns B⁻¹·e_port. It is
// immutable once built, so cached factors are shared between goroutines.
type pieceFactor struct {
	m      *solver.CDense // equilibrated block
	rs, cs []float64
	w      [][]complex128 // per port: B⁻¹·e_port, local
	// The block and ports it was built from: a cache hit is confirmed
	// against them, so a hash collision can never return wrong factors.
	n     int
	block []lentry
	ports []int
}

func (e *pieceFactor) solve(b []complex128) ([]complex128, error) {
	n := len(b)
	bb := make([]complex128, n)
	for i := range b {
		bb[i] = b[i] * complex(e.rs[i], 0)
	}
	x := make([]complex128, n)
	if err := e.m.SolveFactored(bb, x); err != nil {
		return nil, err
	}
	for j := range x {
		x[j] *= complex(e.cs[j], 0)
	}
	return x, nil
}

// cost is the storage of a factor in complex values.
func (e *pieceFactor) cost() int { return e.n*e.n + len(e.w)*e.n + len(e.block) }

// cacheBudget bounds a TornCache's storage (complex values).
const cacheBudget = 1 << 23

func (e *pieceFactor) same(n int, block []lentry, ports []int) bool {
	if e.n != n || len(e.block) != len(block) || len(e.ports) != len(ports) {
		return false
	}
	for i, b := range block {
		if e.block[i] != b {
			return false
		}
	}
	for i, p := range ports {
		if e.ports[i] != p {
			return false
		}
	}
	return true
}

// factorPiece returns the factors of a piece, from the cache when its
// block and port pattern were seen before (compared exactly).
func (d *Diakoptics) factorPiece(pc *piece) (*pieceFactor, bool, error) {
	var key uint64
	if d.cache != nil {
		key = pieceKey(len(pc.idx), pc.block, pc.ports)
		if e := d.cache.get(key, len(pc.idx), pc.block, pc.ports); e != nil {
			return e, true, nil
		}
	}
	e, err := factorBlock(len(pc.idx), pc.block, pc.ports)
	if err != nil {
		return nil, false, err
	}
	if d.cache != nil {
		e = d.cache.put(key, e)
	}
	return e, false, nil
}

// factorBlock equilibrates (rows and columns to unit maximum, as
// solveSystem) and factors an nl×nl block, and solves for B⁻¹·e at each
// local port position.
func factorBlock(nl int, block []lentry, ports []int) (*pieceFactor, error) {
	rs, cs := make([]float64, nl), make([]float64, nl)
	for _, e := range block {
		rs[e.i] = math.Max(rs[e.i], cmplx.Abs(e.v))
	}
	for i, m := range rs {
		if m == 0 {
			return nil, fmt.Errorf("equation %d is empty", i)
		}
		rs[i] = 1 / m
	}
	for _, e := range block {
		cs[e.j] = math.Max(cs[e.j], cmplx.Abs(e.v)*rs[e.i])
	}
	for j, m := range cs {
		if m == 0 {
			return nil, fmt.Errorf("unknown %d appears in no equation", j)
		}
		cs[j] = 1 / m
	}
	m, err := solver.NewCDense(nl)
	if err != nil {
		return nil, err
	}
	for _, e := range block {
		m.Add(e.i, e.j, e.v*complex(rs[e.i]*cs[e.j], 0))
	}
	if err := m.Factor(); err != nil {
		return nil, err
	}
	f := &pieceFactor{m: m, rs: rs, cs: cs, n: nl, block: block, ports: ports}
	for _, lp := range ports {
		u := make([]complex128, nl)
		u[lp] = 1
		x, err := f.solve(u)
		if err != nil {
			return nil, err
		}
		f.w = append(f.w, x)
	}
	return f, nil
}

// pieceKey hashes a block's exact entries and port positions (a
// multiply–xorshift mix of the bit patterns). Results never depend on it:
// a hit is confirmed by comparing the stored block entry by entry.
func pieceKey(n int, block []lentry, ports []int) uint64 {
	h := uint64(0x9e3779b97f4a7c15) ^ uint64(n)
	mix := func(u uint64) {
		h ^= u
		h *= 0xbf58476d1ce4e5b9
		h ^= h >> 31
	}
	for _, e := range block {
		mix(uint64(e.i)<<32 | uint64(uint32(e.j)))
		mix(math.Float64bits(real(e.v)))
		mix(math.Float64bits(imag(e.v)))
	}
	mix(math.MaxUint64)
	for _, p := range ports {
		mix(uint64(p))
	}
	return h
}

// TornCache holds factored torn pieces by content, safe for concurrent
// use. It is bounded by a number of pieces and by the storage of their
// factors (cacheBudget complex values, 128 MiB); when either is reached,
// the entries hit since the last sweep are kept (up to half of each bound)
// and the rest dropped (the fixed pieces a GA keeps reusing survive;
// the one-off pieces of individual candidates do not). Results never
// depend on what is cached, only the time taken.
type TornCache struct {
	mu   sync.Mutex
	m    map[uint64][]*cacheEntry
	n    int
	max  int
	size int // stored complex values
	hits atomic.Int64
	miss atomic.Int64
}

type cacheEntry struct {
	f    *pieceFactor
	used atomic.Bool
}

// NewTornCache returns a cache of at most max pieces (max ≤ 0: 4096).
func NewTornCache(max int) *TornCache {
	if max <= 0 {
		max = 4096
	}
	return &TornCache{m: map[uint64][]*cacheEntry{}, max: max}
}

func (c *TornCache) get(k uint64, n int, block []lentry, ports []int) *pieceFactor {
	c.mu.Lock()
	list := c.m[k]
	c.mu.Unlock()
	for _, e := range list {
		if e.f.same(n, block, ports) {
			e.used.Store(true)
			c.hits.Add(1)
			return e.f
		}
	}
	c.miss.Add(1)
	return nil
}

// put stores f unless an equal piece is already stored, and returns the
// stored one.
func (c *TornCache) put(k uint64, f *pieceFactor) *pieceFactor {
	c.mu.Lock()
	defer c.mu.Unlock()
	for _, e := range c.m[k] {
		if e.f.same(f.n, f.block, f.ports) {
			return e.f
		}
	}
	cost := f.cost()
	if cost > cacheBudget/2 {
		return f // too large to keep
	}
	if c.n >= c.max || c.size+cost > cacheBudget {
		kept := map[uint64][]*cacheEntry{}
		c.n, c.size = 0, 0
		for key, list := range c.m {
			for _, e := range list {
				if e.used.Load() && c.n < c.max/2 && c.size+e.f.cost() <= cacheBudget/2 {
					e.used.Store(false)
					kept[key] = append(kept[key], e)
					c.n++
					c.size += e.f.cost()
				}
			}
		}
		c.m = kept
	}
	c.m[k] = append(c.m[k], &cacheEntry{f: f})
	c.n++
	c.size += cost
	return f
}

// Stats returns the cache's hits, misses and current size.
func (c *TornCache) Stats() (hits, misses int64, entries int) {
	c.mu.Lock()
	defer c.mu.Unlock()
	return c.hits.Load(), c.miss.Load(), c.n
}
