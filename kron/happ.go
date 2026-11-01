// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math/cmplx"
)

// Happ's factorized piecewise algorithms (H. H. Happ, Diakoptics and
// Networks, Academic Press, 1971, Ch. XII; reference/
// happ_1971_diakoptics_and_networks.pdf, PDF page = printed page + 17).
//
// Both functions work on explicit matrices rather than on a circuit: they
// are Happ's algorithms as he states them, so the book's numerical
// examples can be reproduced step by step, and they are cross-checked
// against MNA and against the Aitchison-form tearing of diakoptics.go in
// happ_test.go. Nothing is ever inverted (CLAUDE.md §3a rule 2): Z_cc and
// each subdivision's nodal admittance are factored and solved.

// PiecewiseSolution is the result of Happ's Table 12.1.
type PiecewiseSolution struct {
	// Ic are the closed-path currents i^c (Step 6).
	Ic []complex128
	// Vb are the branch voltages V_b = Z_bb (I^b + i^b) (Step 9), and
	// Jb = I^b + i^b the total branch currents.
	Vb, Jb []complex128
	// Eo are the open-path (junction-pair) voltages E_o = V_o − e_o
	// (Step 12).
	Eo []complex128
	// Power is Kron's invariant e_b·J_b, equal to e_c·i^c + e_o·I° in the
	// orthogonal frame (Happ Ch. IV §5).
	Power complex128
}

// PiecewiseZ is the factorized solution of the orthogonal network in Happ
// Table 12.1 (p. 237, the unconstrained case; derivation pp. 234–236,
// from the orthogonal equations (4.14)–(4.15)). The primitive network is
// zbb (b×b, mutual impedances allowed) with branch voltage sources eb; the
// connection tensor is split into its open-path columns co (b×p, carrying
// the impressed currents io) and its closed-path columns cc (b×m). The
// twelve steps are:
//
//	1  I^b   = C_o I°              7  i^b   = C_c i^c
//	2  V_b⁰  = Z_bb I^b            8  V_b¹  = Z_bb i^b
//	3  e_c   = C_cᵀ e_b            9  V_b   = V_b⁰ + V_b¹
//	4  V_c⁰  = C_cᵀ V_b⁰           10 V_o   = C_oᵀ V_b
//	5  e_c'  = e_c − V_c⁰          11 e_o   = C_oᵀ e_b
//	6  i^c   = Z_cc⁻¹ e_c'          12 E_o   = V_o − e_o
//
// with Z_cc = C_cᵀ Z_bb C_c factored, not inverted. Steps 1–2 are the
// solution with every closed path open (Happ's "torn" state), 3–6 the
// interconnection, 7–9 its contribution.
func PiecewiseZ(zbb, co, cc [][]complex128, eb, io []complex128) (*PiecewiseSolution, error) {
	b := len(zbb)
	if len(eb) != b || len(co) != b || len(cc) != b {
		return nil, fmt.Errorf("%w: Z_bb is %d×%d; e_b, C_o, C_c have %d, %d, %d rows", ErrDimension, b, b, len(eb), len(co), len(cc))
	}
	for i, row := range zbb {
		if len(row) != b {
			return nil, fmt.Errorf("%w: Z_bb row %d has %d entries, want %d", ErrDimension, i, len(row), b)
		}
	}
	p, m := len(io), 0
	if b > 0 {
		m = len(cc[0])
	}
	for i := 0; i < b; i++ {
		if len(co[i]) != p || len(cc[i]) != m {
			return nil, fmt.Errorf("%w: C_o row %d has %d columns (want %d), C_c row %d has %d (want %d)", ErrDimension, i, len(co[i]), p, i, len(cc[i]), m)
		}
	}
	if p+m != b {
		return nil, fmt.Errorf("%w: %d open + %d closed paths for %d branches (Kron's frame must be square, p + m = b)", ErrDimension, p, m, b)
	}
	mulZ := func(x []complex128) []complex128 {
		y := make([]complex128, b)
		for i, row := range zbb {
			for j, v := range row {
				y[i] += v * x[j]
			}
		}
		return y
	}
	// C x and Cᵀ y for a b×k column block.
	mulC := func(c [][]complex128, x []complex128) []complex128 {
		y := make([]complex128, b)
		for i, row := range c {
			for j, v := range row {
				y[i] += v * x[j]
			}
		}
		return y
	}
	tmulC := func(c [][]complex128, k int, y []complex128) []complex128 {
		x := make([]complex128, k)
		for i, row := range c {
			for j, v := range row {
				x[j] += v * y[i]
			}
		}
		return x
	}
	ib0 := mulC(co, io)          // 1
	vb0 := mulZ(ib0)             // 2
	ec := tmulC(cc, m, eb)       // 3
	vc0 := tmulC(cc, m, vb0)     // 4
	ecp := make([]complex128, m) // 5
	for k := range ecp {
		ecp[k] = ec[k] - vc0[k]
	}
	zcc := make([][]complex128, m) // C_cᵀ Z_bb C_c
	for j := 0; j < m; j++ {
		col := make([]complex128, b)
		for i := 0; i < b; i++ {
			col[i] = cc[i][j]
		}
		zc := tmulC(cc, m, mulZ(col))
		for i := 0; i < m; i++ {
			if zcc[i] == nil {
				zcc[i] = make([]complex128, m)
			}
			zcc[i][j] = zc[i]
		}
	}
	ic, err := solveSystem(zcc, ecp) // 6
	if err != nil {
		return nil, fmt.Errorf("%w: closed-path impedance Z_cc of Happ Table 12.1: %v", ErrSingular, err)
	}
	ib := mulC(cc, ic) // 7
	vb1 := mulZ(ib)    // 8
	out := &PiecewiseSolution{Ic: ic, Vb: make([]complex128, b), Jb: make([]complex128, b)}
	for i := range out.Vb { // 9
		out.Vb[i] = vb0[i] + vb1[i]
		out.Jb[i] = ib0[i] + ib[i]
		out.Power += eb[i] * out.Jb[i]
	}
	vo := tmulC(co, p, out.Vb) // 10
	eo := tmulC(co, p, eb)     // 11
	out.Eo = make([]complex128, p)
	for k := range out.Eo { // 12
		out.Eo[k] = vo[k] - eo[k]
	}
	return out, nil
}

// Happ's two-level diakoptics with a single common bus (Ch. XII §4,
// pp. 250–264; Figs. 12.9–12.14; algorithm Fig. 12.12, the rightmost
// column of Table 12.5 with e_T = 0 and I^L = 0; H. H. Happ, "Diakoptics
// and piecewise methods," IEEE Trans. PAS 89, 1373–1382, 1970).
//
// The network is cut along its intersubdivision branches (links). What
// remains is a set of subdivisions, each hinged radially at the common bus
// (ground): its tree branches to that bus are its radial equivalent, with
// impedance Z_k = Y_k⁻¹ where Y_k is the subdivision's nodal admittance
// (Fig. 12.10a). With no mutual coupling between subdivisions Z_TT is
// block-diagonal (eq. 12.104), so each block is solved on its own. Each
// link closes one path through two subdivisions and the common bus
// (Fig. 12.10b), and Z_cc = Uᵀ Z_TT U + Z_LL (eq. 12.78 with C_Lc = 1).
// The interconnection system has the order of the number of links.
//
// Steps (Fig. 12.12, with U the link incidence, +1 at a link's From bus
// and −1 at its To bus):
//
//	a  V_T⁰  = Z_TT I^T                torn subdivisions
//	b  e_c'  = e_L + Uᵀ V_T⁰            voltage across each cut, plus its source
//	c  i^c   = Z_cc⁻¹ e_c'               closed-path (link) currents
//	d  I^T'  = −U i^c                    injected into the subdivisions
//	e  V_T¹  = Z_TT I^T'
//	f  V_T   = V_T⁰ + V_T¹
//
// A link's current flows from its From bus to its To bus, and its source
// e_L acts in that direction: z_ℓ i_ℓ = V_From − V_To + e_ℓ.

// BusRef names bus Bus of subdivision Sub.
type BusRef struct{ Sub, Bus int }

// Link is an intersubdivision branch.
type Link struct {
	From, To BusRef
	Z        complex128 // link impedance; zero is allowed (Happ's Z_b2 = 0)
	E        complex128 // series voltage source, acting From → To
}

// TwoLevel is a network prepared for Happ's two-level solution.
type TwoLevel struct {
	subs  []*denseLU[complex128]
	sizes []int
	links []Link
	// zcol[k][ℓ]: column of Z_TT·U for link ℓ restricted to subdivision k
	// (nil when ℓ does not touch k).
	zcol [][][]complex128
	zcc  [][]complex128
}

// NewTwoLevel factors every subdivision's nodal admittance y[k] (its buses
// against the common bus) and builds Z_cc. A subdivision whose buses do
// not all reach the common bus has a singular Y_k; it is reported, not
// solved (Aitchison 1983 Lemma 2).
func NewTwoLevel(y [][][]complex128, links []Link) (*TwoLevel, error) {
	t := &TwoLevel{subs: make([]*denseLU[complex128], len(y)), sizes: make([]int, len(y)),
		links: append([]Link(nil), links...), zcol: make([][][]complex128, len(y))}
	for k := range y {
		if err := t.setSub(k, y[k]); err != nil {
			return nil, err
		}
	}
	for i, l := range t.links {
		for _, r := range []BusRef{l.From, l.To} {
			if r.Sub < 0 || r.Sub >= len(y) || r.Bus < 0 || r.Bus >= t.sizes[r.Sub] {
				return nil, fmt.Errorf("%w: link %d ends at bus %d of subdivision %d, which does not exist", ErrDimension, i, r.Bus, r.Sub)
			}
		}
		if l.From.Sub == l.To.Sub {
			return nil, fmt.Errorf("%w: link %d joins two buses of subdivision %d; it belongs inside that subdivision", ErrDimension, i, l.From.Sub)
		}
	}
	for k := range y {
		t.columns(k)
	}
	t.buildZcc()
	return t, nil
}

func (t *TwoLevel) setSub(k int, y [][]complex128) error {
	n := len(y)
	ent := make([]lent[complex128], 0, n*n)
	for i, row := range y {
		if len(row) != n {
			return fmt.Errorf("%w: Y of subdivision %d is not square (row %d has %d entries, want %d)", ErrDimension, k, i, len(row), n)
		}
		for j, v := range row {
			if v != 0 {
				ent = append(ent, lent[complex128]{int32(i), int32(j), v})
			}
		}
	}
	lu := newDenseLU[complex128](n, cmplx.Abs, func(x float64) complex128 { return complex(x, 0) })
	if err := lu.factor(ent); err != nil {
		return fmt.Errorf("%w: subdivision %d (%d buses) is singular — a bus with no path to the common bus: %v", ErrSingular, k, n, err)
	}
	t.subs[k], t.sizes[k] = lu, n
	return nil
}

// columns computes Z_k·U for the links touching subdivision k.
func (t *TwoLevel) columns(k int) {
	n := t.sizes[k]
	t.zcol[k] = make([][]complex128, len(t.links))
	u := make([]complex128, n)
	for l, ln := range t.links {
		touch := false
		for i := range u {
			u[i] = 0
		}
		if ln.From.Sub == k {
			u[ln.From.Bus] += 1
			touch = true
		}
		if ln.To.Sub == k {
			u[ln.To.Bus] -= 1
			touch = true
		}
		if !touch {
			continue
		}
		x := make([]complex128, n)
		t.subs[k].solve(u, x)
		t.zcol[k][l] = x
	}
}

// uT returns (Uᵀ x)_ℓ for a vector given per subdivision.
func (t *TwoLevel) uT(l int, x [][]complex128) complex128 {
	ln := t.links[l]
	return x[ln.From.Sub][ln.From.Bus] - x[ln.To.Sub][ln.To.Bus]
}

func (t *TwoLevel) buildZcc() {
	nl := len(t.links)
	t.zcc = make([][]complex128, nl)
	for i := range t.zcc {
		t.zcc[i] = make([]complex128, nl)
	}
	for j := 0; j < nl; j++ {
		for k := range t.zcol {
			col := t.zcol[k][j]
			if col == nil {
				continue
			}
			for i, ln := range t.links {
				if ln.From.Sub == k {
					t.zcc[i][j] += col[ln.From.Bus]
				}
				if ln.To.Sub == k {
					t.zcc[i][j] -= col[ln.To.Bus]
				}
			}
		}
		t.zcc[j][j] += t.links[j].Z
	}
}

// SetSubdivision replaces subdivision k's admittance (a genetic-algorithm
// mutation inside it): only that block is refactored and only its
// columns of Z_TT·U recomputed; Z_cc is rebuilt from the cached columns.
func (t *TwoLevel) SetSubdivision(k int, y [][]complex128) error {
	if k < 0 || k >= len(t.subs) {
		return fmt.Errorf("%w: no subdivision %d", ErrDimension, k)
	}
	if len(y) != t.sizes[k] {
		return fmt.Errorf("%w: subdivision %d has %d buses, new Y has %d", ErrDimension, k, t.sizes[k], len(y))
	}
	if err := t.setSub(k, y); err != nil {
		return err
	}
	t.columns(k)
	t.buildZcc()
	return nil
}

// Zcc returns a copy of the interconnection matrix Uᵀ Z_TT U + Z_LL.
func (t *TwoLevel) Zcc() [][]complex128 {
	out := make([][]complex128, len(t.zcc))
	for i, r := range t.zcc {
		out[i] = append([]complex128(nil), r...)
	}
	return out
}

// RadialZ returns subdivision k's radial-equivalent impedance Z_k = Y_k⁻¹,
// column by column from its factors (for reports and checks against
// Happ's Fig. 12.13; the solution itself never needs it).
func (t *TwoLevel) RadialZ(k int) [][]complex128 {
	n := t.sizes[k]
	z := make([][]complex128, n)
	e := make([]complex128, n)
	x := make([]complex128, n)
	for i := range z {
		z[i] = make([]complex128, n)
	}
	for j := 0; j < n; j++ {
		for i := range e {
			e[i] = 0
		}
		e[j] = 1
		t.subs[k].solve(e, x)
		for i := range x {
			z[i][j] = x[i]
		}
	}
	return z
}

// TwoLevelSolution holds every intermediate of Fig. 12.12, per subdivision
// where the quantity lives on buses and per link where it lives on cuts.
type TwoLevelSolution struct {
	V0, V1, V [][]complex128 // steps a, e, f
	Ecp       []complex128   // step b, e_c'
	Ic        []complex128   // step c, link currents i^c
	Inj       [][]complex128 // step d, I^T' injected into each subdivision
	// Power: Σ V·I^T delivered by the bus injections plus Σ e_ℓ i_ℓ by the
	// link sources, equal to the power absorbed by subdivisions and links.
	Supplied, Absorbed complex128
}

// Solve runs Fig. 12.12 for bus injections inj[k][bus] (current into each
// bus from the common bus).
func (t *TwoLevel) Solve(inj [][]complex128) (*TwoLevelSolution, error) {
	if len(inj) != len(t.subs) {
		return nil, fmt.Errorf("%w: injections for %d subdivisions, network has %d", ErrDimension, len(inj), len(t.subs))
	}
	s := &TwoLevelSolution{}
	alloc := func() [][]complex128 {
		v := make([][]complex128, len(t.subs))
		for k, n := range t.sizes {
			v[k] = make([]complex128, n)
		}
		return v
	}
	s.V0, s.V1, s.V, s.Inj = alloc(), alloc(), alloc(), alloc()
	for k, lu := range t.subs { // a
		if len(inj[k]) != t.sizes[k] {
			return nil, fmt.Errorf("%w: subdivision %d has %d buses, %d injections given", ErrDimension, k, t.sizes[k], len(inj[k]))
		}
		lu.solve(inj[k], s.V0[k])
	}
	nl := len(t.links)
	s.Ecp = make([]complex128, nl) // b
	for l, ln := range t.links {
		s.Ecp[l] = ln.E + t.uT(l, s.V0)
	}
	ic, err := solveSystem(t.zcc, s.Ecp) // c
	if err != nil {
		return nil, fmt.Errorf("%w: interconnection matrix Z_cc of the links (singular exactly when the whole network is, Aitchison 1983 Lemma 1): %v", ErrSingular, err)
	}
	s.Ic = ic
	for l, ln := range t.links { // d
		s.Inj[ln.From.Sub][ln.From.Bus] -= ic[l]
		s.Inj[ln.To.Sub][ln.To.Bus] += ic[l]
	}
	for k := range t.subs { // e: V_T¹ = Z_TT I^T' = −Σ_ℓ (Z_k U)_ℓ i_ℓ, from the cached columns
		for l := range t.links {
			col := t.zcol[k][l]
			if col == nil {
				continue
			}
			for i, v := range col {
				s.V1[k][i] -= v * ic[l]
			}
		}
		for i := range s.V[k] { // f
			s.V[k][i] = s.V0[k][i] + s.V1[k][i]
			s.Supplied += s.V[k][i] * cmplx.Conj(inj[k][i])
		}
	}
	for l, ln := range t.links {
		s.Supplied += ln.E * cmplx.Conj(ic[l])
		s.Absorbed += ln.Z * ic[l] * cmplx.Conj(ic[l])
	}
	// Power absorbed inside the subdivisions: Σ V·(I^T + I^T')* = Σ V·(Y V)*.
	for k := range t.subs {
		for i := range s.V[k] {
			s.Absorbed += s.V[k][i] * cmplx.Conj(inj[k][i]+s.Inj[k][i])
		}
	}
	return s, nil
}

// ClosedPathsZ is Happ's mesh solution with the open paths ignored (Ch.
// XI §2, steps 1–6 on p. 213, eqs. 11.4–11.6): current sources I^b enter
// as equivalent voltage sources −Z_bb I^b (eq. 4.33), so only the
// closed-path columns cc (b×m) of the connection tensor are needed:
//
//	e_c' = C_cᵀ (e_b − Z_bb I^b),  Z_cc = C_cᵀ Z_bb C_c,  i^c = Z_cc⁻¹ e_c',
//	i^b = C_c i^c,  J^b = I^b + i^b,  E_b = Z_bb J^b − e_b.
//
// It returns J^b and E_b.
func ClosedPathsZ(zbb, cc [][]complex128, eb, ib []complex128) (jb, ebOut []complex128, err error) {
	b := len(zbb)
	if len(cc) != b || len(eb) != b || len(ib) != b {
		return nil, nil, fmt.Errorf("%w: Z_bb has %d rows; C_c, e_b, I^b have %d, %d, %d", ErrDimension, b, len(cc), len(eb), len(ib))
	}
	m := 0
	if b > 0 {
		m = len(cc[0])
	}
	zmul := func(x []complex128) []complex128 {
		y := make([]complex128, b)
		for i, row := range zbb {
			if len(row) != b {
				continue
			}
			for j, v := range row {
				y[i] += v * x[j]
			}
		}
		return y
	}
	zi := zmul(ib)
	ecp := make([]complex128, m)
	zcc := make([][]complex128, m)
	for i := range zcc {
		zcc[i] = make([]complex128, m)
	}
	col := make([]complex128, b)
	for j := 0; j < m; j++ {
		for i := 0; i < b; i++ {
			if len(cc[i]) != m {
				return nil, nil, fmt.Errorf("%w: C_c row %d has %d columns, want %d", ErrDimension, i, len(cc[i]), m)
			}
			col[i] = cc[i][j]
			ecp[j] += cc[i][j] * (eb[i] - zi[i])
		}
		zc := zmul(col)
		for r := 0; r < m; r++ {
			for i := 0; i < b; i++ {
				zcc[r][j] += cc[i][r] * zc[i]
			}
		}
	}
	ic, err := solveSystem(zcc, ecp)
	if err != nil {
		return nil, nil, fmt.Errorf("%w: closed-path impedance Z_cc (Happ eq. 11.5): %v", ErrSingular, err)
	}
	jb = make([]complex128, b)
	for i := 0; i < b; i++ {
		jb[i] = ib[i]
		for k := 0; k < m; k++ {
			jb[i] += cc[i][k] * ic[k]
		}
	}
	ebOut = zmul(jb)
	for i := range ebOut {
		ebOut[i] -= eb[i]
	}
	return jb, ebOut, nil
}

// OpenPathsY is the dual, Happ's nodal solution with the closed paths
// ignored (Ch. XI §3, steps 1–6 on p. 218, eqs. 11.22–11.24): voltage
// sources e_b enter as equivalent current sources −Y_bb e_b (eq. 4.38), so
// only the open-path rows ao (p×b) of the inverse connection tensor A^s_b
// are needed — the junction pairs, not the whole open paths:
//
//	I^o' = A_o (I^b − Y_bb e_b),  Y^oo = A_o Y_bb A_oᵀ,  E_o = (Y^oo)⁻¹ I^o',
//	E_b = A_oᵀ E_o,  J^b = Y_bb (E_b + e_b).
//
// It returns E_o, E_b and J^b.
func OpenPathsY(ybb, ao [][]complex128, eb, ib []complex128) (eo, ebOut, jb []complex128, err error) {
	b := len(ybb)
	p := len(ao)
	if len(eb) != b || len(ib) != b {
		return nil, nil, nil, fmt.Errorf("%w: Y_bb has %d rows; e_b, I^b have %d, %d", ErrDimension, b, len(eb), len(ib))
	}
	for k, row := range ao {
		if len(row) != b {
			return nil, nil, nil, fmt.Errorf("%w: A_o row %d has %d columns, want %d", ErrDimension, k, len(row), b)
		}
	}
	ymul := func(x []complex128) []complex128 {
		y := make([]complex128, b)
		for i, row := range ybb {
			if len(row) != b {
				continue
			}
			for j, v := range row {
				y[i] += v * x[j]
			}
		}
		return y
	}
	ye := ymul(eb)
	io := make([]complex128, p)
	for k, row := range ao {
		for i, v := range row {
			io[k] += v * (ib[i] - ye[i])
		}
	}
	yoo := make([][]complex128, p)
	col := make([]complex128, b)
	for j := 0; j < p; j++ {
		copy(col, ao[j])
		yc := ymul(col)
		for r := 0; r < p; r++ {
			if yoo[r] == nil {
				yoo[r] = make([]complex128, p)
			}
			for i := 0; i < b; i++ {
				yoo[r][j] += ao[r][i] * yc[i]
			}
		}
	}
	eo, err = solveSystem(yoo, io)
	if err != nil {
		return nil, nil, nil, fmt.Errorf("%w: open-path admittance Y^oo (Happ eq. 11.23): %v", ErrSingular, err)
	}
	ebOut = make([]complex128, b)
	for k, row := range ao {
		for i, v := range row {
			ebOut[i] += v * eo[k]
		}
	}
	x := make([]complex128, b)
	for i := range x {
		x[i] = ebOut[i] + eb[i]
	}
	jb = ymul(x)
	return eo, ebOut, jb, nil
}
