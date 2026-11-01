// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
)

// H. H. Happ, Diakoptics and Networks, Academic Press, 1971, Appendix I
// (pp. 275–283; reference/happ_1971_diakoptics_and_networks.pdf, PDF
// pages 292–300): the orthogonal network of Chapter III worked
// numerically. A six-branch primitive network (eq. 1.3, p. 276) with
// mutual impedances and voltage sources, three impressed junction
// currents I° = (1, 1, 1) into junctions 1, 2, 3 against the datum 4, and
// the open- and closed-path contours of Fig. I.2 (connection tensor C^b_s
// and its inverse A^s_b, Table I.1, p. 278). The published results:
// Z_ss = Cᵀ Z C (p. 279), the closed-path currents i^c from eq. (3.76)
// (p. 280), the branch currents J^b = C J^s (eq. 1.6, p. 280), and
// i^s = A i^b with zero open-path parts (eq. 3.24, p. 283).
func TestHappAppendixI(t *testing.T) {
	zbb := [6][6]float64{
		{1, 0, 0, 0, 0, -2},
		{0, 2, 1, 0, 0, 0},
		{0, 1, 2, 0, 0, 0},
		{0, 0, 0, 4, -2, 0},
		{0, 0, 0, -2, 2, 0},
		{-2, 0, 0, 0, 0, 5},
	}
	eb := [6]float64{1, 2, 3, 1, 4, 2}
	// C^b_s: rows b1…b6, columns 1, 2, 3 (open paths), c1, c2, c3 (closed).
	cbs := [6][6]float64{
		{0, 0, 0, 1, 1, 0},
		{0, 0, 0, 1, 0, -1},
		{0, -1, 0, 0, 1, 1},
		{0, 0, 1, 0, 1, 1},
		{-1, 0, 0, 0, 1, 0},
		{0, -1, 0, 0, 0, 1},
	}
	// A^s_b: rows 1, 2, 3, c1, c2, c3; columns b1…b6.
	asb := [6][6]float64{
		{0, 0, 1, 0, -1, -1},
		{1, -1, -1, 0, 0, 0},
		{-1, 1, 0, 1, 0, 0},
		{1, 0, -1, 0, 0, 1},
		{0, 0, 1, 0, 0, -1},
		{1, -1, -1, 0, 0, 1},
	}
	z := NewMat[Primitive, Primitive](6, 6)
	c := NewMat[Primitive, Orthogonal](6, 6)
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			z.Set(i, j, complex(zbb[i][j], 0))
			c.Set(i, j, complex(cbs[i][j], 0))
		}
	}
	// A is the inverse of C (Happ eq. 2.9 with Table I.1).
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			var s float64
			for k := 0; k < 6; k++ {
				s += asb[i][k] * cbs[k][j]
			}
			want := 0.0
			if i == j {
				want = 1
			}
			if s != want {
				t.Fatalf("(A·C)[%d][%d] = %g, want %g", i, j, s, want)
			}
		}
	}
	zss, err := Transform(c, z)
	if err != nil {
		t.Fatal(err)
	}
	wantZ := [6][6]float64{ // p. 279
		{2, 0, 2, 0, 0, 2},
		{0, 7, 0, 1, 0, -6},
		{2, 0, 4, 0, 2, 4},
		{0, 1, 0, 3, 2, -3},
		{0, 0, 2, 2, 5, 1},
		{2, -6, 4, -3, 1, 11},
	}
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			if real(zss.At(i, j)) != wantZ[i][j] || imag(zss.At(i, j)) != 0 {
				t.Errorf("Z_ss[%d][%d] = %v, Happ %g", i, j, zss.At(i, j), wantZ[i][j])
			}
		}
	}
	// e_s = Cᵀ e_b: open parts (−4, −5, 1), closed parts (3, 9, 4) (p. 277).
	eVec := NewVec[Primitive](6)
	for i, v := range eb {
		eVec.Set(i, complex(v, 0))
	}
	es, err := TMulVec(c, eVec)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range []float64{-4, -5, 1, 3, 9, 4} {
		if real(es.At(i)) != w {
			t.Errorf("e_s[%d] = %v, Happ %g", i, es.At(i), w)
		}
	}
	// Eq. (3.76): i^c = Z_cc⁻¹ (e_c − Z_co I°).
	io := []float64{1, 1, 1}
	a := make([][]complex128, 3)
	rhs := make([]complex128, 3)
	for r := 0; r < 3; r++ {
		a[r] = make([]complex128, 3)
		rhs[r] = es.At(3 + r)
		for k := 0; k < 3; k++ {
			a[r][k] = zss.At(3+r, 3+k)
			rhs[r] -= zss.At(3+r, k) * complex(io[k], 0)
		}
	}
	ic, err := solveSystem(a, rhs)
	if err != nil {
		t.Fatal(err)
	}
	// Happ multiplies his six-digit (Z_cc)⁻¹ by the right-hand side, so his
	// printed values carry that rounding (up to 6e−6). The exact solution
	// is (1, 82, 15)/61 (det Z_cc = 61): check it to round-off, and Happ's
	// figures to 1e−5.
	for i, w := range []float64{0.016390, 1.344262, 0.245904} { // p. 280
		if math.Abs(real(ic[i])-w) > 1e-5 {
			t.Errorf("i^c[%d] = %.6f, Happ %.6f", i, real(ic[i]), w)
		}
		if exact := []float64{1, 82, 15}[i] / 61; math.Abs(real(ic[i])-exact) > 1e-13 {
			t.Errorf("i^c[%d] = %.15f, exact %.15f", i, real(ic[i]), exact)
		}
	}
	// Eq. (1.6): J^b = C J^s, J^s = (I°, i^c).
	js := NewVec[Orthogonal](6)
	for k := 0; k < 3; k++ {
		js.Set(k, complex(io[k], 0))
		js.Set(3+k, ic[k])
	}
	jb, err := MulVec(c, js)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range []float64{1.360652, -0.229514, 0.590166, 2.590166, 0.344262, -0.754096} { // p. 280
		if math.Abs(real(jb.At(i))-w) > 1e-5 {
			t.Errorf("J^b[%d] = %.6f, Happ %.6f", i, real(jb.At(i)), w)
		}
	}
	// Eqs. (3.13), (3.24): I^b = C I^s with I^s = (I°, 0); i^b = J^b − I^b
	// (p. 283), whose open-path parts A^o i^b vanish.
	isv := NewVec[Orthogonal](6)
	for k := 0; k < 3; k++ {
		isv.Set(k, complex(io[k], 0))
	}
	ib, err := MulVec(c, isv)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range []float64{1.360652, -0.229514, 1.590166, 1.590166, 1.344262, 0.245904} {
		got := real(jb.At(i) - ib.At(i))
		if math.Abs(got-w) > 1e-5 {
			t.Errorf("i^b[%d] = %.6f, Happ %.6f", i, got, w)
		}
	}
	for r := 0; r < 3; r++ {
		var s float64
		for k := 0; k < 6; k++ {
			s += asb[r][k] * real(jb.At(k)-ib.At(k))
		}
		if math.Abs(s) > 1e-12 {
			t.Errorf("open-path part %d of i^s = %g, want 0", r, s)
		}
	}
	// Kron's power invariant: the input power e·J is the same in the
	// primitive and the orthogonal frame (Happ Ch. IV §5).
	var pb, ps float64
	for i := 0; i < 6; i++ {
		pb += eb[i] * real(jb.At(i))
		ps += real(es.At(i)) * real(js.At(i))
	}
	if math.Abs(pb-ps) > 1e-12 {
		t.Errorf("power e·J: primitive %g, orthogonal %g", pb, ps)
	}
}

// toC converts a real matrix to complex.
func toC(a [][]float64) [][]complex128 {
	out := make([][]complex128, len(a))
	for i, r := range a {
		out[i] = make([]complex128, len(r))
		for j, v := range r {
			out[i][j] = complex(v, 0)
		}
	}
	return out
}

func near(t *testing.T, what string, got complex128, want, tol float64) {
	t.Helper()
	if math.Abs(real(got)-want) > tol || math.Abs(imag(got)) > tol {
		t.Errorf("%s = %.7f, Happ %.6f", what, real(got), want)
	}
}

// Happ Table 12.1 (p. 237) run on the Appendix I network: the twelve-step
// factorized algorithm must give the closed-path currents of eq. (3.76),
// and the branch voltages V_b and open-path voltages E_o that Happ prints
// on p. 284 (Appendix I §2, from the orthogonal equations directly).
func TestHappTable12_1(t *testing.T) {
	zbb := [][]float64{
		{1, 0, 0, 0, 0, -2}, {0, 2, 1, 0, 0, 0}, {0, 1, 2, 0, 0, 0},
		{0, 0, 0, 4, -2, 0}, {0, 0, 0, -2, 2, 0}, {-2, 0, 0, 0, 0, 5},
	}
	co := [][]float64{{0, 0, 0}, {0, 0, 0}, {0, -1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}}
	cc := [][]float64{{1, 1, 0}, {1, 0, -1}, {0, 1, 1}, {0, 1, 1}, {0, 1, 0}, {0, 0, 1}}
	eb := []complex128{1, 2, 3, 1, 4, 2}
	io := []complex128{1, 1, 1}
	s, err := PiecewiseZ(toC(zbb), toC(co), toC(cc), eb, io)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range []float64{1, 82, 15} {
		if cmplx.Abs(s.Ic[i]-complex(w/61, 0)) > 1e-14 {
			t.Errorf("i^c[%d] = %v, exact %g/61", i, s.Ic[i], w)
		}
	}
	for i, w := range []float64{2.868844, 0.131138, 0.950818, 9.672140, -4.491808, -6.491784} { // p. 284
		near(t, fmt.Sprintf("V_b[%d]", i), s.Vb[i], w, 5e-5)
	}
	for i, w := range []float64{8.491808, 10.540966, 8.672140} { // p. 284
		near(t, fmt.Sprintf("E_o[%d]", i), s.Eo[i], w, 5e-5)
	}
	for i, w := range []float64{1.360652, -0.229514, 0.590166, 2.590166, 0.344262, -0.754096} { // p. 280
		near(t, fmt.Sprintf("J^b[%d]", i), s.Jb[i], w, 1e-5)
	}
	// Power invariant: e_b·J_b = e_c·i^c + e_o·I° (Happ Ch. IV §5), with
	// e_c = C_cᵀ e_b = (3, 9, 4) and e_o = C_oᵀ e_b = (−4, −5, 1).
	orth := 3*s.Ic[0] + 9*s.Ic[1] + 4*s.Ic[2] + (-4 - 5 + 1)
	if cmplx.Abs(s.Power-orth) > 1e-12 {
		t.Errorf("power: primitive %v, orthogonal %v", s.Power, orth)
	}
	// Kron's output invariant: the junction-pair power E_o·I° equals the
	// power the junction-pair sources deliver, (V_b − e_b)·I^b with
	// I^b = C_o I°.
	var eo, vb complex128
	ib := []complex128{0, 0, -1, 1, -1, -1}
	for k := range io {
		eo += s.Eo[k] * io[k]
	}
	for i := range ib {
		vb += (s.Vb[i] - eb[i]) * ib[i]
	}
	if cmplx.Abs(eo-vb) > 1e-12 {
		t.Errorf("E_o·I° = %v, (V_b − e_b)·I^b = %v", eo, vb)
	}
	// A non-square frame is a dimension error, a dependent closed path a
	// singular Z_cc.
	if _, err := PiecewiseZ(toC(zbb), toC(co), toC([][]float64{{1}, {1}, {0}, {0}, {0}, {0}}), eb, io); !errors.Is(err, ErrDimension) {
		t.Errorf("non-square frame: %v", err)
	}
	dep := [][]float64{{1, 1, 1}, {1, 1, 1}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}
	if _, err := PiecewiseZ(toC(zbb), toC(co), toC(dep), eb, io); !errors.Is(err, ErrSingular) {
		t.Errorf("dependent closed paths: %v", err)
	}
}

// happFig129 is the single-common-bus network of Happ Fig. 12.9 (p. 258)
// as its subdivisions' nodal admittances and links. Every branch is 1 Ω
// except the link from 2A to 2B, which is 0 (p. 261). The subdivision
// trees are the ones whose radial equivalents are Fig. 12.13's Z_TT
// (p. 262): A is the chain G–3A–1A–2A, B the star 2B–{G, 1B, 3B}, C the
// chain G–2C–1C. Bus order: A (1A, 2A, 3A), B (1B, 2B, 3B), C (1C, 2C).
// Link orientations follow Fig. 12.14(b): e_c' is V1B − V1A, V2A − V2B,
// V2C − V3B, V3A − V1C.
func happFig129() ([][][]complex128, []Link) {
	ya := [][]float64{{2, -1, -1}, {-1, 1, 0}, {-1, 0, 2}}
	yb := [][]float64{{1, -1, 0}, {-1, 3, -1}, {0, -1, 1}}
	yc := [][]float64{{1, -1}, {-1, 2}}
	links := []Link{
		{From: BusRef{1, 0}, To: BusRef{0, 0}, Z: 1}, // c1: 1B → 1A
		{From: BusRef{0, 1}, To: BusRef{1, 1}, Z: 0}, // c2: 2A → 2B, Z_b2 = 0
		{From: BusRef{2, 1}, To: BusRef{1, 2}, Z: 1}, // c3: 2C → 3B
		{From: BusRef{0, 2}, To: BusRef{2, 0}, Z: 1}, // c4: 3A → 1C
	}
	return [][][]complex128{toC(ya), toC(yb), toC(yc)}, links
}

// Happ's two-level piecewise solution, Figs. 12.12–12.14 (pp. 262–263):
// building blocks (Z_TT, Z_cc⁻¹), every step of the algorithm for
// I^T = 0.5 at all eight buses, and the final bus voltages.
func TestHappSingleCommonBus(t *testing.T) {
	y, links := happFig129()
	tl, err := NewTwoLevel(y, links)
	if err != nil {
		t.Fatal(err)
	}
	zTT := [][][]float64{ // Fig. 12.13
		{{2, 2, 1}, {2, 3, 1}, {1, 1, 1}},
		{{2, 1, 1}, {1, 1, 1}, {1, 1, 2}},
		{{2, 1}, {1, 1}},
	}
	for k, zk := range zTT {
		got := tl.RadialZ(k)
		for i := range zk {
			for j, w := range zk[i] {
				near(t, fmt.Sprintf("Z_TT[%d][%d][%d]", k, i, j), got[i][j], w, 1e-14)
			}
		}
	}
	zccInv := [][]float64{ // Fig. 12.13, Z_cc⁻¹ (printed to six digits)
		{0.370370, 0.259259, 0.037037, 0.037037},
		{0.259259, 0.481481, -0.074074, -0.074074},
		{0.037037, -0.074074, 0.303704, 0.103704},
		{0.037037, -0.074074, 0.103704, 0.303704},
	}
	zcc := tl.Zcc()
	for j := range zccInv {
		e := make([]complex128, 4)
		e[j] = 1
		col, err := solveSystem(zcc, e)
		if err != nil {
			t.Fatal(err)
		}
		for i := range col {
			near(t, fmt.Sprintf("Z_cc⁻¹[%d][%d]", i, j), col[i], zccInv[i][j], 1e-6)
		}
	}
	inj := [][]complex128{{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5}}
	s, err := tl.Solve(inj)
	if err != nil {
		t.Fatal(err)
	}
	flat := func(v [][]complex128) []complex128 {
		var out []complex128
		for _, x := range v {
			out = append(out, x...)
		}
		return out
	}
	check := func(what string, got []complex128, want []float64) {
		t.Helper()
		for i, w := range want {
			near(t, fmt.Sprintf("%s[%d]", what, i), got[i], w, 1e-6)
		}
	}
	// Fig. 12.14 (a)–(f); Happ truncates to six digits.
	check("V_T⁰", flat(s.V0), []float64{2.5, 3.0, 1.5, 2.0, 1.5, 2.0, 1.5, 1.0})
	check("e_c'", s.Ecp, []float64{-0.5, 1.5, -1.0, 0.0})
	check("i^c", s.Ic, []float64{0.166666, 0.666666, -0.433333, -0.233333})
	check("I^T'", flat(s.Inj), []float64{0.166666, -0.666666, 0.233333, -0.166666, 0.666666, -0.433333, -0.233333, 0.433333})
	check("V_T¹", flat(s.V1), []float64{-0.766666, -1.433333, -0.266666, -0.1, 0.066666, -0.366666, -0.033333, 0.2})
	check("V_T", flat(s.V), []float64{1.733333, 1.566666, 1.233333, 1.9, 1.566666, 1.633333, 1.466666, 1.2})
	// Exact values: i^c = (1/6, 2/3, −13/30, −7/30).
	for i, w := range []float64{1.0 / 6, 2.0 / 3, -13.0 / 30, -7.0 / 30} {
		if cmplx.Abs(s.Ic[i]-complex(w, 0)) > 1e-14 {
			t.Errorf("i^c[%d] = %.16f, exact %.16f", i, real(s.Ic[i]), w)
		}
	}
	if cmplx.Abs(s.Supplied-s.Absorbed) > 1e-13*cmplx.Abs(s.Supplied) {
		t.Errorf("power supplied %v, absorbed %v", s.Supplied, s.Absorbed)
	}
}

// The same network as a netlist: Happ's two-level solution must equal the
// MNA operating point and the Aitchison-form tearing of diakoptics.go
// (cutting links b1, b3, b4 as RL1, RL3, RL4; b2 is the 0 V source VL2 and stays), to 1e−12.
func TestHappTwinMNA(t *testing.T) {
	src := `happ fig 12.9
RA1 n3a 0 1
RA2 n1a n3a 1
RA3 n2a n1a 1
RB1 n2b 0 1
RB2 n1b n2b 1
RB3 n3b n2b 1
RC1 n2c 0 1
RC2 n1c n2c 1
RL1 n1b n1a 1
VL2 n2a n2b 0
RL3 n2c n3b 1
RL4 n3a n1c 1
`
	buses := []string{"n1a", "n2a", "n3a", "n1b", "n2b", "n3b", "n1c", "n2c"}
	var b strings.Builder
	b.WriteString(src)
	for _, n := range buses {
		fmt.Fprintf(&b, "I%s 0 %s 0.5\n", n, n)
	}
	b.WriteString(".end\n")
	c := parse(t, b.String())
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	torn, err := TornOP(c, Tearing{Torn: []string{"RL1", "RL3", "RL4"}}, nil)
	if err != nil {
		t.Fatal(err)
	}
	y, links := happFig129()
	tl, err := NewTwoLevel(y, links)
	if err != nil {
		t.Fatal(err)
	}
	s, err := tl.Solve([][]complex128{{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5}})
	if err != nil {
		t.Fatal(err)
	}
	i := 0
	for k := range s.V {
		for _, v := range s.V[k] {
			mna, ok1 := op.Voltage(buses[i])
			tv, ok2 := torn.Voltage(buses[i])
			if !ok1 || !ok2 {
				t.Fatalf("node %s missing", buses[i])
			}
			if math.Abs(real(v)-mna) > 1e-12 || math.Abs(tv-mna) > 1e-12 {
				t.Errorf("%s: Happ %.15f, MNA %.15f, torn %.15f", buses[i], real(v), mna, tv)
			}
			i++
		}
	}
}

// A mutation inside one subdivision refactors only that block; the result
// equals a fresh build. A subdivision with a bus that cannot reach the
// common bus is reported as singular.
func TestHappSubdivisionUpdate(t *testing.T) {
	y, links := happFig129()
	tl, err := NewTwoLevel(y, links)
	if err != nil {
		t.Fatal(err)
	}
	inj := [][]complex128{{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5}}
	// Subdivision C at 1 kHz: its branch to G becomes 1 Ω ∥ 10 µF.
	yc := toC([][]float64{{1, -1}, {-1, 2}})
	yc[1][1] += complex(0, 2*math.Pi*1e3*10e-6)
	if err := tl.SetSubdivision(2, yc); err != nil {
		t.Fatal(err)
	}
	got, err := tl.Solve(inj)
	if err != nil {
		t.Fatal(err)
	}
	y[2] = yc
	fresh, err := NewTwoLevel(y, links)
	if err != nil {
		t.Fatal(err)
	}
	want, err := fresh.Solve(inj)
	if err != nil {
		t.Fatal(err)
	}
	for k := range got.V {
		for i := range got.V[k] {
			if cmplx.Abs(got.V[k][i]-want.V[k][i]) > 1e-14 {
				t.Errorf("V[%d][%d]: updated %v, fresh %v", k, i, got.V[k][i], want.V[k][i])
			}
		}
	}
	// Remove C's branch to G: 1C and 2C float.
	if err := tl.SetSubdivision(2, toC([][]float64{{1, -1}, {-1, 1}})); !errors.Is(err, ErrSingular) || !strings.Contains(err.Error(), "subdivision 2") {
		t.Errorf("floating subdivision: %v", err)
	}
	if _, err := NewTwoLevel(y, []Link{{From: BusRef{0, 0}, To: BusRef{0, 1}, Z: 1}}); !errors.Is(err, ErrDimension) {
		t.Errorf("link inside a subdivision: %v", err)
	}
}

// appendixI returns the Appendix I primitive network (eq. 1.3, p. 276).
func appendixI() (zbb [][]complex128, eb []complex128) {
	return toC([][]float64{
		{1, 0, 0, 0, 0, -2}, {0, 2, 1, 0, 0, 0}, {0, 1, 2, 0, 0, 0},
		{0, 0, 0, 4, -2, 0}, {0, 0, 0, -2, 2, 0}, {-2, 0, 0, 0, 0, 5},
	}), []complex128{1, 2, 3, 1, 4, 2}
}

// Happ Appendix I §3 (pp. 289–295): the same network with the open paths
// of Fig. I.3, which do not all end at one datum. Table I.2 (p. 291) gives
// C^b_s and A^s_b. The bus currents of Fig. I.1 become I° = (1, 1, 2)
// (p. 290, the change of open-path frame of Ch. VIII). The branch
// currents are frame-independent, so J^b is the same as in §2; the
// open-path voltages become E_o = (8.491808, 1.868826, 8.672140) (eq.
// 3.33, p. 293), and the junction-pair power E_o·I° is invariant.
func TestHappAppendixI3(t *testing.T) {
	zbb, eb := appendixI()
	cbs := [][]float64{
		{0, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 0, -1}, {0, -1, 0, 0, 1, 1},
		{0, -1, 1, 0, 1, 1}, {-1, 0, 0, 0, 1, 0}, {0, -1, 0, 0, 0, 1},
	}
	asb := [][]float64{
		{0, 0, 1, 0, -1, -1}, {1, -1, -1, 0, 0, 0}, {0, 0, -1, 1, 0, 0},
		{1, 0, -1, 0, 0, 1}, {0, 0, 1, 0, 0, -1}, {1, -1, -1, 0, 0, 1},
	}
	for i := 0; i < 6; i++ {
		for j := 0; j < 6; j++ {
			var s float64
			for k := 0; k < 6; k++ {
				s += asb[i][k] * cbs[k][j]
			}
			if want := map[bool]float64{true: 1, false: 0}[i == j]; s != want {
				t.Fatalf("Table I.2: (A·C)[%d][%d] = %g", i, j, s)
			}
		}
	}
	co, cc := make([][]float64, 6), make([][]float64, 6)
	for i, r := range cbs {
		co[i], cc[i] = r[:3], r[3:]
	}
	io := []complex128{1, 1, 2}
	s, err := PiecewiseZ(zbb, toC(co), toC(cc), eb, io)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range []float64{1, 82, 15} {
		if cmplx.Abs(s.Ic[i]-complex(w/61, 0)) > 1e-14 {
			t.Errorf("i^c[%d] = %v, exact %g/61", i, s.Ic[i], w)
		}
	}
	for i, w := range []float64{1.360652, -0.229514, 0.590166, 2.590166, 0.344262, -0.754096} { // eq. 2.9, p. 290
		near(t, fmt.Sprintf("J^b[%d]", i), s.Jb[i], w, 1e-5)
	}
	for i, w := range []float64{8.491808, 1.868826, 8.672140} { // eq. 3.33, p. 293
		near(t, fmt.Sprintf("E_o[%d]", i), s.Eo[i], w, 5e-5)
	}
	// I^b = C_o I° is the same branch-current excitation as in §2
	// (eq. 3.16, p. 290: (0, 0, −1, 1, −1, −1)).
	for i, w := range []float64{0, 0, -1, 1, -1, -1} {
		var v float64
		for k := range io {
			v += co[i][k] * real(io[k])
		}
		if v != w {
			t.Errorf("I^b[%d] = %g, Happ %g", i, v, w)
		}
	}
	// The same E_o by the nodal route with closed paths ignored (Ch. XI
	// §3), from the junction-pair rows of Table I.2's A^s_b alone.
	ybb := inverse(t, zbb)
	eo, ebv, jb, err := OpenPathsY(ybb, toC(asb[:3]), eb, toC([][]float64{{0, 0, -1, 1, -1, -1}})[0])
	if err != nil {
		t.Fatal(err)
	}
	for k := range eo {
		if cmplx.Abs(eo[k]-s.Eo[k]) > 1e-12 {
			t.Errorf("E_o[%d]: Ch. XI nodal %v, Table 12.1 %v", k, eo[k], s.Eo[k])
		}
	}
	for i, w := range []float64{1.868844, -1.868862, -2.049182, 8.672140, -8.491808, -8.491784} { // eq. 3.41, p. 294
		near(t, fmt.Sprintf("E_b[%d]", i), ebv[i], w, 5e-5)
		if cmplx.Abs(jb[i]-s.Jb[i]) > 1e-12 {
			t.Errorf("J^b[%d]: nodal %v, Table 12.1 %v", i, jb[i], s.Jb[i])
		}
	}
	// Junction-pair power is invariant under the change of open paths:
	// E_o·I° here equals (8.491808, 10.540966, 8.672140)·(1, 1, 1) of §2.
	var pNew complex128
	for k := range io {
		pNew += s.Eo[k] * io[k]
	}
	if math.Abs(real(pNew)-(8.491808+10.540966+8.672140)) > 1e-4 {
		t.Errorf("E_o·I° = %v, §2 frame %g", pNew, 8.491808+10.540966+8.672140)
	}
}

// inverse forms Y_bb = Z_bb⁻¹ column by column (test helper only; the
// algorithms never invert, but Ch. XI §3 starts from the primitive
// admittance).
func inverse(t *testing.T, z [][]complex128) [][]complex128 {
	t.Helper()
	n := len(z)
	y := make([][]complex128, n)
	for i := range y {
		y[i] = make([]complex128, n)
	}
	for j := 0; j < n; j++ {
		e := make([]complex128, n)
		e[j] = 1
		x, err := solveSystem(z, e)
		if err != nil {
			t.Fatal(err)
		}
		for i := range x {
			y[i][j] = x[i]
		}
	}
	return y
}

// Happ Ch. XI on the Appendix I network: ignoring the open paths (mesh
// solution, §2) and ignoring the closed paths (nodal solution, §3) both
// give the Table 12.1 branch currents; and eliminating instead of ignoring
// gives the Schur-complement identities (Z_cc)⁻¹ = Y^cc − Y^co (Y^oo)⁻¹
// Y^oc (eq. 11.10) and (Y^oo)⁻¹ = Z_oo − Z_oc (Z_cc)⁻¹ Z_co (eq. 11.28),
// with Y^ss = (Z_ss)⁻¹. Ignoring inactive axes is how kronspice skips
// unmodeled ports (CLAUDE.md §3a.1 item 4; Kron 1939 p. 429).
func TestHappChapterXI(t *testing.T) {
	zbb, eb := appendixI()
	cbs := [][]float64{
		{0, 0, 0, 1, 1, 0}, {0, 0, 0, 1, 0, -1}, {0, -1, 0, 0, 1, 1},
		{0, 0, 1, 0, 1, 1}, {-1, 0, 0, 0, 1, 0}, {0, -1, 0, 0, 0, 1},
	}
	asb := [][]float64{
		{0, 0, 1, 0, -1, -1}, {1, -1, -1, 0, 0, 0}, {-1, 1, 0, 1, 0, 0},
		{1, 0, -1, 0, 0, 1}, {0, 0, 1, 0, 0, -1}, {1, -1, -1, 0, 0, 1},
	}
	co, cc := make([][]float64, 6), make([][]float64, 6)
	for i, r := range cbs {
		co[i], cc[i] = r[:3], r[3:]
	}
	io := []complex128{1, 1, 1}
	ref, err := PiecewiseZ(zbb, toC(co), toC(cc), eb, io)
	if err != nil {
		t.Fatal(err)
	}
	ib := make([]complex128, 6)
	for i := range ib {
		for k := range io {
			ib[i] += complex(co[i][k], 0) * io[k]
		}
	}
	jm, em, err := ClosedPathsZ(zbb, toC(cc), eb, ib)
	if err != nil {
		t.Fatal(err)
	}
	eo, en, jn, err := OpenPathsY(inverse(t, zbb), toC(asb[:3]), eb, ib)
	if err != nil {
		t.Fatal(err)
	}
	for i := range jm {
		if cmplx.Abs(jm[i]-ref.Jb[i]) > 1e-12 || cmplx.Abs(jn[i]-ref.Jb[i]) > 1e-12 {
			t.Errorf("J^b[%d]: mesh %v, nodal %v, Table 12.1 %v", i, jm[i], jn[i], ref.Jb[i])
		}
		if cmplx.Abs(em[i]-en[i]) > 1e-12 {
			t.Errorf("E_b[%d]: mesh %v, nodal %v", i, em[i], en[i])
		}
	}
	for k := range eo {
		if cmplx.Abs(eo[k]-ref.Eo[k]) > 1e-12 {
			t.Errorf("E_o[%d]: nodal %v, Table 12.1 %v", k, eo[k], ref.Eo[k])
		}
	}
	// Z_ss = Cᵀ Z C, then Y^ss = (Z_ss)⁻¹ and the two Schur identities.
	c := toC(cbs)
	zss := make([][]complex128, 6)
	for i := range zss {
		zss[i] = make([]complex128, 6)
		for j := range zss[i] {
			for a := 0; a < 6; a++ {
				for b := 0; b < 6; b++ {
					zss[i][j] += c[a][i] * zbb[a][b] * c[b][j]
				}
			}
		}
	}
	yss := inverse(t, zss)
	block := func(m [][]complex128, r0, c0 int) [][]complex128 {
		out := make([][]complex128, 3)
		for i := range out {
			out[i] = append([]complex128(nil), m[r0+i][c0:c0+3]...)
		}
		return out
	}
	mul := func(a, b [][]complex128) [][]complex128 {
		out := make([][]complex128, 3)
		for i := range out {
			out[i] = make([]complex128, 3)
			for j := range out[i] {
				for k := 0; k < 3; k++ {
					out[i][j] += a[i][k] * b[k][j]
				}
			}
		}
		return out
	}
	sub := func(a, b [][]complex128) [][]complex128 {
		out := make([][]complex128, 3)
		for i := range out {
			out[i] = make([]complex128, 3)
			for j := range out[i] {
				out[i][j] = a[i][j] - b[i][j]
			}
		}
		return out
	}
	isIdentity := func(what string, m [][]complex128) {
		t.Helper()
		for i := range m {
			for j := range m[i] {
				want := complex(0, 0)
				if i == j {
					want = 1
				}
				if cmplx.Abs(m[i][j]-want) > 1e-12 {
					t.Errorf("%s: [%d][%d] = %v", what, i, j, m[i][j])
				}
			}
		}
	}
	zoo, zoc, zco, zcc := block(zss, 0, 0), block(zss, 0, 3), block(zss, 3, 0), block(zss, 3, 3)
	yoo, yoc, yco, ycc := block(yss, 0, 0), block(yss, 0, 3), block(yss, 3, 0), block(yss, 3, 3)
	isIdentity("Z_cc·[Y^cc − Y^co (Y^oo)⁻¹ Y^oc] (eq. 11.10)", mul(zcc, sub(ycc, mul(yco, mul(inverse(t, yoo), yoc)))))
	isIdentity("Y^oo·[Z_oo − Z_oc (Z_cc)⁻¹ Z_co] (eq. 11.28)", mul(yoo, sub(zoo, mul(zoc, mul(inverse(t, zcc), zco)))))
}
