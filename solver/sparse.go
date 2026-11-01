// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package solver

import (
	"fmt"
	"math"
	"math/cmplx"
	"sort"
)

// Sparse LU backend (CLAUDE.md §3, §5): Gilbert–Peierls left-looking LU with
// threshold partial pivoting, after a minimum-degree symmetric ordering.
//
// References:
//   - J. R. Gilbert, T. Peierls, "Sparse partial pivoting in time
//     proportional to arithmetic operations," SIAM J. Sci. Stat. Comput. 9(5),
//     1988 — column-by-column factorization; the nonzero pattern of each
//     column of L\A(:,k) is the set reachable in the graph of L (Sec. 3).
//   - T. A. Davis, Direct Methods for Sparse Linear Systems, SIAM 2006:
//     cs_reach / cs_dfs (non-recursive depth-first search, Sec. 3.2),
//     cs_lu (Sec. 6.3), minimum degree (Ch. 7). Reimplemented from the text.
//   - T. A. Davis, E. Palamadai Natarajan, "Algorithm 907: KLU," ACM TOMS
//     37(3), 2010 — diagonal preference with a relative pivot tolerance,
//     which keeps the fill-reducing order where it is numerically safe.
//
// The matrix pattern is discovered while stamping: the first time an entry
// (i, j) is added it gets a slot; later stampings reuse it. When the pattern
// grows, the ordering is recomputed at the next Solve. Once the pattern is
// stable (after the first Newton iteration), Clear, Add and Solve do not
// allocate (§8 rule 7): L and U keep their capacity between factorizations.

// DiagPivotTol is the relative tolerance for keeping the diagonal as the
// pivot: it is kept when |a_kk| ≥ DiagPivotTol·max|a_ik| (KLU's default).
const DiagPivotTol = 1e-3

// MaxSparse bounds sparse systems.
const MaxSparse = 1 << 22

type scalar interface{ float64 | complex128 }

// sparseLU is the factorization engine shared by the real and complex
// backends.
type sparseLU[T scalar] struct {
	n    int
	abs  func(T) float64
	fin  func(T) bool
	slot map[int64]int32 // (i, j) → index into val
	ri   []int32         // slot → row
	ci   []int32         // slot → column
	val  []T
	// dirty means the pattern grew since the ordering was computed.
	dirty bool
	// Stamp-order cache: devices stamp in the same order every iteration,
	// so the k-th Add after Clear usually hits the same entry as last
	// time. Replaying that sequence avoids a map lookup per stamp.
	seqKey  []int64
	seqSlot []int32
	cursor  int

	// Ordering: permuted index k is original unknown q[k]; qinv inverts it.
	q, qinv []int32
	// Permuted matrix in compressed columns: column k holds rows ap[k]:ap[k+1]
	// as permuted row indices arow with values from val[aslot].
	ap, arow, aslot []int32

	// L (unit lower, diagonal not stored) and U (upper incl. diagonal) by
	// columns. L's row indices are permuted row ids; U's are pivot steps.
	lp, li []int32
	lx     []T
	up, ui []int32
	ux     []T
	udiag  []T
	pinv   []int32 // permuted row id → pivot step, −1 if not yet pivotal
	prow   []int32 // pivot step → permuted row id

	// Workspace.
	x      []T     // dense column, indexed by permuted row id
	xi     []int32 // reach stack (2n)
	pstack []int32
	mark   []int32
	gen    int32
	w      []T
	valid  bool // a factorization is held
}

func newSparseLU[T scalar](n int, abs func(T) float64, fin func(T) bool) (*sparseLU[T], error) {
	if n <= 0 || n > MaxSparse {
		return nil, fmt.Errorf("%w: %d", ErrBadSize, n)
	}
	return &sparseLU[T]{
		n: n, abs: abs, fin: fin,
		slot:  make(map[int64]int32, 8*n),
		dirty: true,
		q:     make([]int32, n), qinv: make([]int32, n),
		ap: make([]int32, n+1),
		lp: make([]int32, n+1), up: make([]int32, n+1),
		udiag: make([]T, n),
		pinv:  make([]int32, n), prow: make([]int32, n),
		x: make([]T, n), xi: make([]int32, 2*n), pstack: make([]int32, n),
		mark: make([]int32, n), w: make([]T, n),
	}, nil
}

func (s *sparseLU[T]) clear() {
	for i := range s.val {
		s.val[i] = 0
	}
	s.valid = false
	s.cursor = 0
}

func (s *sparseLU[T]) add(i, j int, v T) {
	if i < 0 || j < 0 || i >= s.n || j >= s.n {
		return
	}
	key := int64(i)*int64(s.n) + int64(j)
	c := s.cursor
	s.cursor++
	if c < len(s.seqKey) && s.seqKey[c] == key {
		s.val[s.seqSlot[c]] += v
		s.valid = false
		return
	}
	k, ok := s.slot[key]
	if !ok {
		k = int32(len(s.val))
		s.slot[key] = k
		s.ri = append(s.ri, int32(i))
		s.ci = append(s.ci, int32(j))
		s.val = append(s.val, 0)
		s.dirty = true
	}
	if c < len(s.seqKey) {
		s.seqKey[c], s.seqSlot[c] = key, k
	} else {
		s.seqKey = append(s.seqKey, key)
		s.seqSlot = append(s.seqSlot, k)
	}
	s.val[k] += v
	s.valid = false
}

func (s *sparseLU[T]) at(i, j int) T {
	if k, ok := s.slot[int64(i)*int64(s.n)+int64(j)]; ok {
		return s.val[k]
	}
	return 0
}

// analyze computes the ordering and the permuted column structure. It runs
// only when the pattern changed, so its allocations are setup cost.
func (s *sparseLU[T]) analyze() {
	n := s.n
	order := minimumDegree(n, s.ri, s.ci)
	for k, v := range order {
		s.q[k] = int32(v)
		s.qinv[v] = int32(k)
	}
	// Count entries per permuted column, then fill (rows sorted for a
	// deterministic traversal order).
	cnt := make([]int32, n+1)
	for k := range s.val {
		cnt[s.qinv[s.ci[k]]+1]++
	}
	for k := 0; k < n; k++ {
		cnt[k+1] += cnt[k]
	}
	copy(s.ap, cnt)
	nz := len(s.val)
	s.arow = make([]int32, nz)
	s.aslot = make([]int32, nz)
	next := append([]int32(nil), cnt[:n]...)
	for k := range s.val {
		c := s.qinv[s.ci[k]]
		p := next[c]
		next[c]++
		s.arow[p] = s.qinv[s.ri[k]]
		s.aslot[p] = int32(k)
	}
	for c := 0; c < n; c++ {
		lo, hi := s.ap[c], s.ap[c+1]
		sort.Sort(colSorter{s.arow[lo:hi], s.aslot[lo:hi]})
	}
	s.dirty = false
}

type colSorter struct{ r, s []int32 }

func (c colSorter) Len() int           { return len(c.r) }
func (c colSorter) Less(i, j int) bool { return c.r[i] < c.r[j] }
func (c colSorter) Swap(i, j int) {
	c.r[i], c.r[j] = c.r[j], c.r[i]
	c.s[i], c.s[j] = c.s[j], c.s[i]
}

// factor runs Gilbert–Peierls on the permuted matrix.
func (s *sparseLU[T]) factor() error {
	n := s.n
	if s.dirty {
		s.analyze()
	}
	scale := 0.0
	for _, v := range s.val {
		if !s.fin(v) {
			return ErrNonFinite
		}
		if a := s.abs(v); a > scale {
			scale = a
		}
	}
	if scale == 0 {
		return ErrSingularMatrix
	}
	for i := range s.pinv {
		s.pinv[i] = -1
	}
	s.li, s.lx = s.li[:0], s.lx[:0]
	s.ui, s.ux = s.ui[:0], s.ux[:0]
	x := s.x
	for k := 0; k < n; k++ {
		s.lp[k] = int32(len(s.li))
		s.up[k] = int32(len(s.ui))
		// Reach of column k's pattern in the graph of L: top..n-1 of xi
		// in topological order.
		top := s.reach(k)
		for p := top; p < n; p++ {
			x[s.xi[p]] = 0
		}
		// Same column-relative singularity criterion as the dense backends.
		colMax := 0.0
		for p := s.ap[k]; p < s.ap[k+1]; p++ {
			v := s.val[s.aslot[p]]
			x[s.arow[p]] = v
			if a := s.abs(v); a > colMax {
				colMax = a
			}
		}
		tiny := colMax * float64(n) * 1e-14
		// Sparse triangular solve x = L \ A(:,k).
		for p := top; p < n; p++ {
			r := s.xi[p]
			j := s.pinv[r]
			if j < 0 {
				continue
			}
			xj := x[r]
			if xj == 0 {
				continue
			}
			for q := s.lp[j]; q < s.lp[j+1]; q++ {
				x[s.li[q]] -= s.lx[q] * xj
			}
		}
		// Pivot among the rows not yet pivotal.
		piv, best := int32(-1), -1.0
		for p := top; p < n; p++ {
			r := s.xi[p]
			if s.pinv[r] >= 0 {
				// Entry of U.
				s.ui = append(s.ui, s.pinv[r])
				s.ux = append(s.ux, x[r])
				continue
			}
			if a := s.abs(x[r]); a > best || (a == best && r < piv) {
				best, piv = a, r
			}
		}
		if piv < 0 || best <= tiny {
			s.lp[k+1], s.up[k+1] = int32(len(s.li)), int32(len(s.ui))
			return fmt.Errorf("%w: pivot %d is %.3g", ErrSingularMatrix, k, math.Max(best, 0))
		}
		// Diagonal preference (KLU): keep row k if it is large enough.
		if s.pinv[k] < 0 && s.marked(int32(k)) && s.abs(x[k]) >= DiagPivotTol*best && s.abs(x[k]) > tiny {
			piv = int32(k)
		}
		pv := x[piv]
		s.udiag[k] = pv
		s.pinv[piv] = int32(k)
		s.prow[k] = piv
		for p := top; p < n; p++ {
			r := s.xi[p]
			if s.pinv[r] < 0 {
				s.li = append(s.li, r)
				s.lx = append(s.lx, x[r]/pv)
			}
		}
	}
	s.lp[n], s.up[n] = int32(len(s.li)), int32(len(s.ui))
	s.valid = true
	return nil
}

// marked reports whether permuted row r is in the current reach.
func (s *sparseLU[T]) marked(r int32) bool { return s.mark[r] == s.gen }

// reach performs the depth-first search of Gilbert & Peierls (Davis 2006,
// cs_reach) from every row of A(:,k). Pivotal rows lead on to the rows of
// their L column. It returns top; the reach is xi[top:n].
func (s *sparseLU[T]) reach(k int) int {
	n := s.n
	s.gen++
	if s.gen == math.MaxInt32 {
		for i := range s.mark {
			s.mark[i] = 0
		}
		s.gen = 1
	}
	top := n
	stack := s.xi // uses xi[0:n) as the DFS stack, xi[top:n) as output
	for p := s.ap[k]; p < s.ap[k+1]; p++ {
		r := s.arow[p]
		if s.mark[r] == s.gen {
			continue
		}
		// Non-recursive DFS from r.
		head := 0
		stack[0] = r
		for head >= 0 {
			v := stack[head]
			j := s.pinv[v]
			if s.mark[v] != s.gen {
				s.mark[v] = s.gen
				if j >= 0 {
					s.pstack[head] = s.lp[j]
				}
			}
			done := true
			if j >= 0 {
				end := s.lp[j+1]
				for q := s.pstack[head]; q < end; q++ {
					w := s.li[q]
					if s.mark[w] == s.gen {
						continue
					}
					s.pstack[head] = q + 1
					head++
					stack[head] = w
					done = false
					break
				}
			}
			if done {
				head--
				top--
				s.xi[top] = v
			}
		}
	}
	return top
}

// solve solves A·x = b with the held factorization.
func (s *sparseLU[T]) solve(b, x []T) error {
	n := s.n
	if len(b) != n || len(x) != n {
		return fmt.Errorf("%w: rhs %d, x %d, matrix %d", ErrBadSize, len(b), len(x), n)
	}
	if !s.valid {
		if err := s.factor(); err != nil {
			return err
		}
	}
	// w indexed by permuted row id: w = b permuted.
	w := s.w
	for i := 0; i < n; i++ {
		w[s.qinv[i]] = b[i]
	}
	// Forward: L y = P w, y stored in x (pivot-step order) via s.x.
	y := s.x
	for k := 0; k < n; k++ {
		yk := w[s.prow[k]]
		y[k] = yk
		if yk == 0 {
			continue
		}
		for p := s.lp[k]; p < s.lp[k+1]; p++ {
			w[s.li[p]] -= s.lx[p] * yk
		}
	}
	// Backward: U z = y, column-oriented.
	for k := n - 1; k >= 0; k-- {
		zk := y[k] / s.udiag[k]
		y[k] = zk
		if zk == 0 {
			continue
		}
		for p := s.up[k]; p < s.up[k+1]; p++ {
			y[s.ui[p]] -= s.ux[p] * zk
		}
	}
	for k := 0; k < n; k++ {
		v := y[k]
		if !s.fin(v) {
			return ErrNonFinite
		}
		x[s.q[k]] = v
	}
	return nil
}

// nnz returns the number of stored entries in A, L and U (diagonal of U
// included).
func (s *sparseLU[T]) nnz() (a, l, u int) { return len(s.val), len(s.li), len(s.ui) + s.n }

// minimumDegree returns a fill-reducing elimination order for the pattern
// of A + Aᵀ: repeatedly eliminate the node of least degree (ties to the
// lowest index, so the order is deterministic) and join its neighbours into
// a clique (Tinney scheme 2; Davis 2006, Ch. 7). Simple and exact, not the
// approximate-degree variant: it runs once per pattern.
func minimumDegree(n int, ri, ci []int32) []int {
	adj := make([]map[int32]struct{}, n)
	for i := range adj {
		adj[i] = map[int32]struct{}{}
	}
	for k := range ri {
		i, j := ri[k], ci[k]
		if i != j {
			adj[i][j] = struct{}{}
			adj[j][i] = struct{}{}
		}
	}
	gone := make([]bool, n)
	order := make([]int, 0, n)
	nbr := make([]int32, 0, 16)
	for len(order) < n {
		v, best := -1, math.MaxInt
		for i := 0; i < n; i++ {
			if !gone[i] && len(adj[i]) < best {
				v, best = i, len(adj[i])
			}
		}
		gone[v] = true
		order = append(order, v)
		nbr = nbr[:0]
		for u := range adj[v] {
			nbr = append(nbr, u)
		}
		for _, u := range nbr {
			delete(adj[u], int32(v))
			for _, w := range nbr {
				if w != u {
					adj[u][w] = struct{}{}
				}
			}
		}
		adj[v] = nil
	}
	return order
}

// Sparse is the real sparse backend for MNA (Matrix interface).
type Sparse struct{ lu *sparseLU[float64] }

// NewSparse allocates an n×n sparse system with an empty pattern.
func NewSparse(n int) (*Sparse, error) {
	lu, err := newSparseLU(n, math.Abs, func(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) })
	if err != nil {
		return nil, err
	}
	return &Sparse{lu}, nil
}

// Size returns the number of unknowns.
func (s *Sparse) Size() int { return s.lu.n }

// Clear zeroes every stored entry, keeping the pattern.
func (s *Sparse) Clear() { s.lu.clear() }

// Add adds v to entry (i, j); out-of-range indices (ground) are ignored.
func (s *Sparse) Add(i, j int, v float64) { s.lu.add(i, j, v) }

// At returns entry (i, j) of the stamped matrix.
func (s *Sparse) At(i, j int) float64 { return s.lu.at(i, j) }

// Solve factors the matrix (if it changed) and solves A·x = b.
func (s *Sparse) Solve(b, x []float64) error { return s.lu.solve(b, x) }

// AbsProduct sets y_i = Σ_j |A_ij|·|x_j|.
func (s *Sparse) AbsProduct(x, y []float64) {
	for i := range y {
		y[i] = 0
	}
	l := s.lu
	for k, v := range l.val {
		y[l.ri[k]] += math.Abs(v) * math.Abs(x[l.ci[k]])
	}
}

// NNZ reports stored entries of A, L and U.
func (s *Sparse) NNZ() (a, l, u int) { return s.lu.nnz() }

// CSparse is the complex sparse backend for AC analysis.
type CSparse struct{ lu *sparseLU[complex128] }

// NewCSparse allocates an n×n complex sparse system.
func NewCSparse(n int) (*CSparse, error) {
	lu, err := newSparseLU(n, cmplx.Abs, cfinite)
	if err != nil {
		return nil, err
	}
	return &CSparse{lu}, nil
}

// Size returns the number of unknowns.
func (s *CSparse) Size() int { return s.lu.n }

// Clear zeroes every stored entry, keeping the pattern.
func (s *CSparse) Clear() { s.lu.clear() }

// Add adds v to entry (i, j); out-of-range indices (ground) are ignored.
func (s *CSparse) Add(i, j int, v complex128) { s.lu.add(i, j, v) }

// At returns entry (i, j) of the stamped matrix.
func (s *CSparse) At(i, j int) complex128 { return s.lu.at(i, j) }

// AbsProduct sets y_i = Σ_j |A_ij|·|x_j|.
func (s *CSparse) AbsProduct(x []complex128, y []float64) {
	for i := range y {
		y[i] = 0
	}
	l := s.lu
	for k, v := range l.val {
		y[l.ri[k]] += cmplx.Abs(v) * cmplx.Abs(x[l.ci[k]])
	}
}

// Solve factors the matrix (if it changed) and solves A·x = b.
func (s *CSparse) Solve(b, x []complex128) error { return s.lu.solve(b, x) }
