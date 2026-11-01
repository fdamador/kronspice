// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"sync/atomic"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/solver"
)

// Tearing of the Newton Jacobian (milestone M9c; CLAUDE.md §3a rule 3).
//
// A nonlinear circuit is solved by Newton–Raphson: at every iteration (and
// every transient step) its devices are linearized at the present iterate
// and the linear system J·x = b is solved. J is an MNA matrix stamped by
// the devices, and tearing applies to it unchanged. JacobianTearing is a
// linear-solver backend (analysis.Options.Solver) that solves every such
// system by diakoptics, so the whole Newton machinery of the analysis —
// junction limiting, gmin and source stepping, pseudo-transient
// continuation, transient integration — runs on torn systems.
//
// The cut is given as elements whose node-to-node coupling is torn:
// resistors and capacitors, and also nonlinear devices — a diode, a
// transistor — whose linearized conductances then become coupling entries
// between pieces (tearTerminals lists what is accepted). The unknowns fall into pieces: the
// connected components of J's pattern with the torn elements' node pairs
// removed; an unknown of a device (a branch current, an internal node)
// joins the piece of the nodes it couples to. Every entry of J that joins
// two pieces is a coupling entry, and with B the block-diagonal part and
// E = P·M·Q the coupling entries (P selecting their rows, Q their
// columns),
//
//	z = B⁻¹·b                         each piece on its own
//	(I + M·Q·B⁻¹·P)·w = M·Q·z         interconnection, order = coupling rows
//	x = z − B⁻¹·P·w
//
// This is Aitchison's algorithm for a matrix block-diagonal apart from a
// few coupling entries (P. W. Aitchison, J. Eng. Math. 1983, Sec. 2, as
// set out in CLAUDE.md §3a). A torn resistor keeps its diagonal stamps in
// the blocks of its two pieces: each piece sees it as a conductance to the
// datum, so a piece never floats because of the cut (Aitchison's Lemma 2
// holds wherever the untorn circuit has a DC path to ground).
//
// Reuse between iterations. A piece whose block is bit-identical to the
// previous solve keeps its factors: the linear pieces of a nonlinear
// circuit are factored once per analysis, and only the pieces holding
// devices that moved are refactored — Kron's reason for tearing, applied
// to Newton's method. With ReuseTol > 0 a piece whose block changed by
// less than that fraction also keeps its (now approximate) factors; the
// refinement below, against the exact J, removes the difference.
//
// Accuracy. Every solve is refined against the whole J and accepted only
// at a componentwise backward error of 1e−12 or better (it stops at 1e−14,
// the level of a backward-stable LU). A solve that does not reach it, or a
// singular piece, is solved by an ordinary LU of the whole J instead and
// counted as a fallback: tearing never changes which circuits converge.

// JacobianTearing is an analysis.SolverFactory that tears every linear
// system of an analysis. It may be shared by concurrent analyses.
type JacobianTearing struct {
	cuts  [][2]int // node unknowns of each torn element (−1: ground)
	names []string
	// ReuseTol lets a piece keep its factors while its block changes by
	// at most this fraction of its largest entry (0: only when unchanged).
	ReuseTol float64

	solves, factors, reuses, stale, fallbacks, refines atomic.Int64
	maxPieces                                          atomic.Int64
}

// JacobianStats counts what the torn solves did.
type JacobianStats struct {
	// Solves is the number of linear systems solved.
	Solves int64
	// Factors, Reuses and StaleReuses count pieces factored, pieces whose
	// unchanged factors were reused, and pieces reused within ReuseTol.
	Factors, Reuses, StaleReuses int64
	// Fallbacks counts systems solved by an untorn LU instead.
	Fallbacks int64
	// Refinements is the total number of refinement steps.
	Refinements int64
	// MaxPieces is the largest number of pieces seen in one system.
	MaxPieces int64
}

// NewJacobianTearing prepares a cut of circuit c's Jacobian along the
// named elements, which must be resistors or capacitors between two
// non-ground nodes.
func NewJacobianTearing(c *circuit.Circuit, t Tearing) (*JacobianTearing, error) {
	if err := c.Validate(); err != nil {
		return nil, err
	}
	j := &JacobianTearing{}
	seen := map[string]bool{}
	for _, name := range t.Torn {
		key := strings.ToLower(name)
		if seen[key] {
			return nil, fmt.Errorf("%w: %s is torn twice", ErrUnsupported, name)
		}
		seen[key] = true
		var e *circuit.Element
		for i := range c.Elements {
			if strings.EqualFold(c.Elements[i].Name, name) {
				e = &c.Elements[i]
				break
			}
		}
		if e == nil {
			return nil, fmt.Errorf("%w: no element %q", ErrUnsupported, name)
		}
		terms, err := tearTerminals(e)
		if err != nil {
			return nil, err
		}
		if len(terms) < 2 {
			return nil, fmt.Errorf("%w: tearing %s: it connects to ground, so tearing it separates nothing", ErrUnsupported, name)
		}
		for x := 0; x < len(terms); x++ {
			for y := x + 1; y < len(terms); y++ {
				if terms[x] != terms[y] {
					j.cuts = append(j.cuts, [2]int{terms[x] - 1, terms[y] - 1})
				}
			}
		}
		j.names = append(j.names, e.Name)
	}
	return j, nil
}

// tearTerminals returns the non-ground terminals of an element whose
// Jacobian coupling may be torn: any element whose stamps couple node
// voltages only — resistors, capacitors, diodes, transistors, VCCS (G) and
// voltage-controlled B elements. Elements with a branch-current unknown
// (voltage sources, inductors, E, F, H, current-controlled B) are refused:
// their coupling runs through that unknown, which the cut does not name.
func tearTerminals(e *circuit.Element) ([]int, error) {
	switch e.Kind {
	case circuit.VoltageSource, circuit.Inductor, circuit.NonlinearInductor, circuit.VCVS, circuit.CCCS, circuit.CCVS, circuit.Mutual:
		return nil, fmt.Errorf("%w: tearing %s: its coupling runs through a branch-current unknown; tear a neighbouring element", ErrUnsupported, e.Name)
	}
	var out []int
	add := func(n int) {
		if n == 0 {
			return
		}
		for _, m := range out {
			if m == n {
				return
			}
		}
		out = append(out, n)
	}
	add(e.Nodes[0])
	add(e.Nodes[1])
	for _, n := range e.Extra {
		add(n)
	}
	if e.Kind == circuit.VCCS {
		add(e.CtrlNodes[0])
		add(e.CtrlNodes[1])
	}
	return out, nil
}

// Stats returns the counters accumulated over every system made so far.
func (j *JacobianTearing) Stats() JacobianStats {
	return JacobianStats{Solves: j.solves.Load(), Factors: j.factors.Load(), Reuses: j.reuses.Load(),
		StaleReuses: j.stale.Load(), Fallbacks: j.fallbacks.Load(), Refinements: j.refines.Load(),
		MaxPieces: j.maxPieces.Load()}
}

// Matrix makes a torn real system (DC and transient Jacobians).
func (j *JacobianTearing) Matrix(n int) (solver.Matrix, error) {
	return newTorn[float64](n, j, math.Abs, func(f float64) float64 { return f }, func(n int) (fallbackSolver[float64], error) {
		return solver.NewMatrix(n, solver.Auto)
	})
}

// CMatrix makes a torn complex system (AC).
func (j *JacobianTearing) CMatrix(n int) (solver.CMatrix, error) {
	return newTorn[complex128](n, j, cmplx.Abs, func(f float64) complex128 { return complex(f, 0) }, func(n int) (fallbackSolver[complex128], error) {
		return solver.NewCMatrix(n, solver.Auto)
	})
}

type num interface{ float64 | complex128 }

// fallbackSolver is the untorn LU used when a torn solve is refused.
type fallbackSolver[T num] interface {
	Clear()
	Add(i, j int, v T)
	Solve(b, x []T) error
}

// tornMat is one torn linear system. Its pattern is learned from the
// stamps; the tearing plan is rebuilt only when a new entry appears.
type tornMat[T num] struct {
	n    int
	jt   *JacobianTearing
	abs  func(T) float64
	real func(float64) T

	slot     map[int64]int32
	ri, ci   []int32
	val      []T
	plan     *tornPlan[T]
	newFB    func(int) (fallbackSolver[T], error)
	fallback fallbackSolver[T]
}

func newTorn[T num](n int, j *JacobianTearing, abs func(T) float64, real func(float64) T, fb func(int) (fallbackSolver[T], error)) (*tornMat[T], error) {
	if n <= 0 {
		return nil, fmt.Errorf("%w: %d unknowns", solver.ErrBadSize, n)
	}
	return &tornMat[T]{n: n, jt: j, abs: abs, real: real, slot: map[int64]int32{}, newFB: fb}, nil
}

func (m *tornMat[T]) Size() int { return m.n }

func (m *tornMat[T]) Clear() {
	var zero T
	for i := range m.val {
		m.val[i] = zero
	}
}

func (m *tornMat[T]) Add(i, j int, v T) {
	if i < 0 || j < 0 || i >= m.n || j >= m.n {
		return
	}
	k := int64(i)*int64(m.n) + int64(j)
	s, ok := m.slot[k]
	if !ok {
		s = int32(len(m.val))
		m.slot[k] = s
		m.ri, m.ci = append(m.ri, int32(i)), append(m.ci, int32(j))
		var zero T
		m.val = append(m.val, zero)
		m.plan = nil
	}
	m.val[s] += v
}

func (m *tornMat[T]) AbsProduct(x []T, y []float64) {
	for i := range y {
		y[i] = 0
	}
	for s, v := range m.val {
		y[m.ri[s]] += m.abs(v) * m.abs(x[m.ci[s]])
	}
}

// tornPlan is the tearing of one pattern.
type tornPlan[T num] struct {
	pieces  []*tpiece[T]
	pieceOf []int32
	local   []int32
	off     []int32 // coupling slots
	rows    []int32 // distinct rows of the coupling slots
	rowIdx  []int32 // global row → index in rows, −1
	icN     int
	ic      *denseLU[T]
	icEnt   []lent[T]
	// work space
	z, x, r, h, w, dx []T
	den               []float64
}

// tpiece is one diagonal block.
type tpiece[T num] struct {
	idx   []int32 // global unknowns
	slots []int32 // block slots
	li    []int32 // local row and column of each block slot
	lj    []int32
	ports []int32 // indices into plan.rows of the coupling rows in this piece
	lu    *denseLU[T]
	last  []T
	have  bool
	W     [][]T // per port: B⁻¹·e_row, local
	ent   []lent[T]
	b, xl []T
}

type lent[T num] struct {
	i, j int32
	v    T
}

func (m *tornMat[T]) buildPlan() error {
	n := m.n
	cut := map[int64]bool{}
	for _, c := range m.jt.cuts {
		if c[0] < n && c[1] < n {
			cut[int64(c[0])*int64(n)+int64(c[1])] = true
			cut[int64(c[1])*int64(n)+int64(c[0])] = true
		}
	}
	parent := make([]int32, n)
	for i := range parent {
		parent[i] = int32(i)
	}
	find := func(i int32) int32 {
		for parent[i] != i {
			parent[i] = parent[parent[i]]
			i = parent[i]
		}
		return i
	}
	for s := range m.val {
		i, j := m.ri[s], m.ci[s]
		if i == j || cut[int64(i)*int64(n)+int64(j)] {
			continue
		}
		if a, b := find(i), find(j); a != b {
			if a < b {
				parent[b] = a
			} else {
				parent[a] = b
			}
		}
	}
	p := &tornPlan[T]{pieceOf: make([]int32, n), local: make([]int32, n), rowIdx: make([]int32, n)}
	root := make([]int32, n)
	for i := range root {
		root[i] = -1
		p.rowIdx[i] = -1
	}
	for i := int32(0); i < int32(n); i++ {
		r := find(i)
		if root[r] < 0 {
			root[r] = int32(len(p.pieces))
			p.pieces = append(p.pieces, &tpiece[T]{})
		}
		pc := p.pieces[root[r]]
		p.pieceOf[i], p.local[i] = root[r], int32(len(pc.idx))
		pc.idx = append(pc.idx, i)
	}
	for s := range m.val {
		i, j := m.ri[s], m.ci[s]
		if p.pieceOf[i] == p.pieceOf[j] {
			pc := p.pieces[p.pieceOf[i]]
			pc.slots = append(pc.slots, int32(s))
			pc.li, pc.lj = append(pc.li, p.local[i]), append(pc.lj, p.local[j])
			continue
		}
		p.off = append(p.off, int32(s))
		if p.rowIdx[i] < 0 {
			p.rowIdx[i] = int32(len(p.rows))
			p.rows = append(p.rows, i)
			pc := p.pieces[p.pieceOf[i]]
			pc.ports = append(pc.ports, p.rowIdx[i])
		}
	}
	for _, pc := range p.pieces {
		nl := len(pc.idx)
		if nl > solver.MaxDense {
			return fmt.Errorf("%w: torn piece of %d unknowns", solver.ErrBadSize, nl)
		}
		pc.lu = newDenseLU[T](nl, m.abs, m.real)
		pc.last = make([]T, len(pc.slots))
		pc.ent = make([]lent[T], len(pc.slots))
		pc.W = make([][]T, len(pc.ports))
		for k := range pc.W {
			pc.W[k] = make([]T, nl)
		}
		pc.b, pc.xl = make([]T, nl), make([]T, nl)
	}
	p.icN = len(p.rows)
	if p.icN > 0 {
		p.ic = newDenseLU[T](p.icN, m.abs, m.real)
	}
	p.z, p.x, p.r, p.dx = make([]T, n), make([]T, n), make([]T, n), make([]T, n)
	p.h, p.w = make([]T, p.icN), make([]T, p.icN)
	p.den = make([]float64, n)
	m.plan = p
	if np := int64(len(p.pieces)); np > m.jt.maxPieces.Load() {
		m.jt.maxPieces.Store(np) // a statistic: a lost race only under-reports
	}
	return nil
}

// Solve solves J·x = b by tearing, or by an untorn LU when the torn solve
// is refused.
func (m *tornMat[T]) Solve(b, x []T) error {
	if len(b) != m.n || len(x) != m.n {
		return fmt.Errorf("%w: rhs %d, x %d, matrix %d", solver.ErrBadSize, len(b), len(x), m.n)
	}
	m.jt.solves.Add(1)
	for _, v := range m.val {
		if f := m.abs(v); math.IsNaN(f) || math.IsInf(f, 0) {
			return solver.ErrNonFinite
		}
	}
	if m.plan == nil {
		if err := m.buildPlan(); err != nil {
			return m.untorn(b, x)
		}
	}
	ok, stale := m.torn(b, x, true)
	if !ok && stale {
		// Stale factors did not refine to the exact solution: refactor
		// the pieces that changed and try once more.
		ok, _ = m.torn(b, x, false)
	}
	if !ok {
		return m.untorn(b, x)
	}
	return nil
}

func (m *tornMat[T]) untorn(b, x []T) error {
	m.jt.fallbacks.Add(1)
	if m.fallback == nil {
		f, err := m.newFB(m.n)
		if err != nil {
			return err
		}
		m.fallback = f
	}
	m.fallback.Clear()
	for s, v := range m.val {
		m.fallback.Add(int(m.ri[s]), int(m.ci[s]), v)
	}
	return m.fallback.Solve(b, x)
}

// torn reports whether the torn solve succeeded (x is then the solution),
// and whether it used stale factors (allowed by allowStale).
func (m *tornMat[T]) torn(b, x []T, allowStale bool) (ok, usedStale bool) {
	p := m.plan
	var zero T
	// Factor (or reuse) each piece.
	for _, pc := range p.pieces {
		same, small := pc.have, pc.have && allowStale && m.jt.ReuseTol > 0
		mx, dmax := 0.0, 0.0
		for k, s := range pc.slots {
			v := m.val[s]
			if v != pc.last[k] {
				same = false
				if small {
					dmax = math.Max(dmax, m.abs(v-pc.last[k]))
				}
			}
			if small {
				mx = math.Max(mx, m.abs(v))
			}
		}
		switch {
		case same:
			m.jt.reuses.Add(1)
			continue
		case small && dmax <= m.jt.ReuseTol*mx:
			m.jt.stale.Add(1)
			usedStale = true
			continue
		}
		for k, s := range pc.slots {
			pc.ent[k] = lent[T]{pc.li[k], pc.lj[k], m.val[s]}
			pc.last[k] = m.val[s]
		}
		pc.have = false
		if err := pc.lu.factor(pc.ent); err != nil {
			return false, usedStale
		}
		m.jt.factors.Add(1)
		for k, rk := range pc.ports {
			u := pc.b
			for i := range u {
				u[i] = zero
			}
			u[p.local[p.rows[rk]]] = 1
			pc.lu.solve(u, pc.W[k])
		}
		pc.have = true
	}
	// Interconnection: I + M·Q·B⁻¹·P.
	if p.icN > 0 {
		ent := p.icEnt[:0]
		for a := 0; a < p.icN; a++ {
			ent = append(ent, lent[T]{int32(a), int32(a), 1})
		}
		for _, s := range p.off {
			r, c, v := m.ri[s], m.ci[s], m.val[s]
			if v == zero {
				continue
			}
			pc := p.pieces[p.pieceOf[c]]
			lc := p.local[c]
			for k, rk := range pc.ports {
				if g := pc.W[k][lc]; g != zero {
					ent = append(ent, lent[T]{p.rowIdx[r], rk, v * g})
				}
			}
		}
		p.icEnt = ent
		if err := p.ic.factor(ent); err != nil {
			return false, usedStale
		}
	}
	apply := func(rhs, out []T) {
		z := out
		for _, pc := range p.pieces {
			nz := false
			for li, g := range pc.idx {
				pc.b[li] = rhs[g]
				nz = nz || rhs[g] != zero
			}
			if !nz {
				for _, g := range pc.idx {
					z[g] = zero
				}
				continue
			}
			pc.lu.solve(pc.b, pc.xl)
			for li, g := range pc.idx {
				z[g] = pc.xl[li]
			}
		}
		if p.icN == 0 {
			return
		}
		for a := range p.h {
			p.h[a] = zero
		}
		for _, s := range p.off {
			p.h[p.rowIdx[m.ri[s]]] += m.val[s] * z[m.ci[s]]
		}
		p.ic.solve(p.h, p.w)
		for _, pc := range p.pieces {
			for k, rk := range pc.ports {
				wk := p.w[rk]
				if wk == zero {
					continue
				}
				col := pc.W[k]
				for li, g := range pc.idx {
					z[g] -= wk * col[li]
				}
			}
		}
	}
	apply(b, x)
	berr := math.Inf(1)
	for pass := 0; pass <= maxRefine; pass++ {
		for i := range p.r {
			p.r[i], p.den[i] = b[i], m.abs(b[i])
		}
		for s, v := range m.val {
			i := m.ri[s]
			t := v * x[m.ci[s]]
			p.r[i] -= t
			p.den[i] += m.abs(t)
		}
		berr = 0
		for i, r := range p.r {
			if p.den[i] > 0 {
				berr = math.Max(berr, m.abs(r)/p.den[i])
			} else if r != zero {
				berr = math.Inf(1)
			}
		}
		if berr <= refineDone || pass == maxRefine {
			break
		}
		apply(p.r, p.dx)
		for i := range x {
			x[i] += p.dx[i]
		}
		m.jt.refines.Add(1)
	}
	if !(berr <= refineAccept) {
		return false, usedStale
	}
	for _, v := range x {
		if f := m.abs(v); math.IsNaN(f) || math.IsInf(f, 0) {
			return false, usedStale
		}
	}
	return true, usedStale
}

// denseLU is an equilibrated dense LU with partial pivoting over float64
// or complex128, allocated once (Golub & Van Loan, Alg. 3.4.1).
type denseLU[T num] struct {
	n      int
	lu     []T
	piv    []int
	rs, cs []float64
	rt, ct []T // rs, cs as T
	col    []float64
	abs    func(T) float64
	real   func(float64) T
}

func newDenseLU[T num](n int, abs func(T) float64, real func(float64) T) *denseLU[T] {
	return &denseLU[T]{n: n, lu: make([]T, n*n), piv: make([]int, n), rs: make([]float64, n),
		cs: make([]float64, n), rt: make([]T, n), ct: make([]T, n), col: make([]float64, n), abs: abs, real: real}
}

// factor equilibrates (rows and columns to unit maximum) and factors the
// matrix given by its entries.
func (d *denseLU[T]) factor(ent []lent[T]) error {
	n := d.n
	var zero T
	for i := range d.lu {
		d.lu[i] = zero
	}
	for i := range d.rs {
		d.rs[i], d.cs[i], d.col[i] = 0, 0, 0
	}
	for _, e := range ent {
		d.rs[e.i] = math.Max(d.rs[e.i], d.abs(e.v))
	}
	for i, v := range d.rs {
		if v == 0 {
			return fmt.Errorf("%w: equation %d is empty", solver.ErrSingularMatrix, i)
		}
		d.rs[i] = 1 / v
	}
	for _, e := range ent {
		d.cs[e.j] = math.Max(d.cs[e.j], d.abs(e.v)*d.rs[e.i])
	}
	for j, v := range d.cs {
		if v == 0 {
			return fmt.Errorf("%w: unknown %d appears in no equation", solver.ErrSingularMatrix, j)
		}
		d.cs[j] = 1 / v
	}
	for i := range d.rs {
		d.rt[i], d.ct[i] = d.real(d.rs[i]), d.real(d.cs[i])
	}
	for _, e := range ent {
		d.lu[int(e.i)*n+int(e.j)] += e.v * d.rt[e.i] * d.ct[e.j]
	}
	for i := 0; i < n; i++ {
		for j, v := range d.lu[i*n : (i+1)*n] {
			d.col[j] = math.Max(d.col[j], d.abs(v))
		}
	}
	lu := d.lu
	for k := 0; k < n; k++ {
		p := k
		maxv := d.abs(lu[k*n+k])
		for i := k + 1; i < n; i++ {
			if v := d.abs(lu[i*n+k]); v > maxv {
				maxv, p = v, i
			}
		}
		if maxv <= d.col[k]*float64(n)*1e-14 {
			return fmt.Errorf("%w: pivot %d is %.3g", solver.ErrSingularMatrix, k, maxv)
		}
		d.piv[k] = p
		if p != k {
			rk, rp := lu[k*n:(k+1)*n], lu[p*n:(p+1)*n]
			for j := range rk {
				rk[j], rp[j] = rp[j], rk[j]
			}
		}
		pivot := lu[k*n+k]
		for i := k + 1; i < n; i++ {
			f := lu[i*n+k] / pivot
			lu[i*n+k] = f
			if f == zero {
				continue
			}
			ri, rk := lu[i*n:(i+1)*n], lu[k*n:(k+1)*n]
			for j := k + 1; j < n; j++ {
				ri[j] -= f * rk[j]
			}
		}
	}
	return nil
}

// solve solves with the factors; x is its working storage (b is not
// changed, and b and x must differ).
func (d *denseLU[T]) solve(b, x []T) {
	n := d.n
	lu := d.lu
	for i := range x {
		x[i] = b[i] * d.rt[i]
	}
	for k := 0; k < n; k++ {
		if p := d.piv[k]; p != k {
			x[k], x[p] = x[p], x[k]
		}
	}
	for i := 0; i < n; i++ {
		s := x[i]
		for j := 0; j < i; j++ {
			s -= lu[i*n+j] * x[j]
		}
		x[i] = s
	}
	for i := n - 1; i >= 0; i-- {
		s := x[i]
		for j := i + 1; j < n; j++ {
			s -= lu[i*n+j] * x[j]
		}
		x[i] = s / lu[i*n+i]
	}
	for j := range x {
		x[j] *= d.ct[j]
	}
}
