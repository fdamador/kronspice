// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"sort"
	"strings"
	"sync"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// tearable lists a circuit's R, C and L elements that carry no chord or
// radiation edge.
func tearable(t *testing.T, c *circuit.Circuit, dc bool) []string {
	t.Helper()
	var out []string
	for _, e := range c.Elements {
		switch e.Kind {
		case circuit.Resistor, circuit.Capacitor:
		case circuit.Inductor:
			if dc {
				continue
			}
		default:
			continue
		}
		if _, err := Tear(c, Tearing{Torn: []string{e.Name}}, nil); err == nil {
			out = append(out, e.Name)
		}
	}
	return out
}

// cuts returns the tearings tried on a circuit: every single tearable
// branch, every second one, and all of them.
func cuts(names []string) [][]string {
	var out [][]string
	for _, n := range names {
		out = append(out, []string{n})
	}
	var alt []string
	for i := 0; i < len(names); i += 2 {
		alt = append(alt, names[i])
	}
	if len(alt) > 1 {
		out = append(out, alt)
	}
	if len(names) > 2 {
		out = append(out, names)
	}
	return out
}

// CLAUDE.md §3a.1 torn-solution test: for every linear golden circuit and
// a set of random ones, under many cuts, the torn operating point and AC
// response equal MNA to 1e−9 relative, the frame invariants agree with the
// untorn solve (power is invariant under tearing, Roth 1959 Sec. 5), the
// counts obey m = b − N + S after tearing, and a cut that leaves a
// singular block is reported, never solved.
func TestTornTwin(t *testing.T) {
	solved, split, singular, worst := 0, 0, 0, 0.0
	for name, c := range linearCircuits(t, 60) {
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			continue
		}
		var mnaI []float64
		for _, e := range op.Elements {
			mnaI = append(mnaI, e.Current)
		}
		untorn, err := Build(c)
		if err != nil {
			t.Fatal(err)
		}
		before := untorn.Counts()
		ref, err := untorn.Solve(NodeForm, 0, true)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for _, cut := range cuts(tearable(t, c, true)) {
			d, err := Tear(c, Tearing{Torn: cut}, nil)
			if err != nil {
				t.Fatalf("%s %v: %v", name, cut, err)
			}
			after := d.Counts()
			if after.Meshes != after.Branches-after.Junctions+after.SubNetworks ||
				after.Meshes != before.Meshes-len(cut)+(after.SubNetworks-before.SubNetworks) {
				t.Errorf("%s %v: counts %+v after tearing %+v", name, cut, after, before)
			}
			sol, st, err := d.Solve(0, true)
			if errors.Is(err, ErrSingular) {
				singular++
				if !strings.Contains(err.Error(), "torn piece") && !strings.Contains(err.Error(), "nearly floating") {
					t.Errorf("%s %v: singular cut not explained: %v", name, cut, err)
				}
				continue
			}
			if err != nil {
				t.Errorf("%s %v: %v", name, cut, err)
				continue
			}
			solved++
			if st.Pieces > 1 {
				split++
			}
			r, err := d.n.opResult(sol)
			if err != nil {
				t.Fatal(err)
			}
			w, at := worstRel(real2c(r.NodeVoltages), real2c(op.NodeVoltages))
			worst = math.Max(worst, w)
			if w > twinTol {
				t.Errorf("%s %v OP: node %s off by %.3g", name, cut, op.NodeNames[at], w)
			}
			var kI []float64
			for _, e := range r.Elements {
				kI = append(kI, e.Current)
			}
			if w, at := worstRel(real2c(kI), real2c(mnaI)); w > twinTol {
				t.Errorf("%s %v OP: current of %s off by %.3g", name, cut, op.Elements[at].Name, w)
			}
			scale := 1e-300
			for k := range sol.E.v {
				scale += cmplx.Abs(sol.E.v[k] * sol.I.v[k])
			}
			if e := cmplx.Abs(sol.PrimitivePower-ref.PrimitivePower) / scale; e > twinTol {
				t.Errorf("%s %v: impressed power %v torn, %v untorn", name, cut, sol.PrimitivePower, ref.PrimitivePower)
			}
			if e := cmplx.Abs(sol.Tellegen) / scale; e > twinTol {
				t.Errorf("%s %v: Tellegen sum %v", name, cut, sol.Tellegen)
			}
		}
		// AC at three frequencies, inductors tearable too.
		mna, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 1, FStart: 10, FStop: 1e7}, analysis.Options{})
		if err != nil {
			continue
		}
		freqs := mna.Freq
		for _, cut := range cuts(tearable(t, c, false)) {
			r, err := TornAC(c, freqs, Tearing{Torn: cut}, nil)
			if errors.Is(err, ErrSingular) {
				singular++
				continue
			}
			if err != nil {
				t.Errorf("%s %v AC: %v", name, cut, err)
				continue
			}
			solved++
			for k := range freqs {
				var got, want []complex128
				for n := 1; n < len(mna.V); n++ {
					got, want = append(got, r.V[n][k]), append(want, mna.V[n][k])
				}
				w, at := worstRel(got, want)
				worst = math.Max(worst, w)
				if w > twinTol {
					t.Errorf("%s %v AC at %g Hz: node %s off by %.3g", name, cut, freqs[k], mna.NodeNames[at+1], w)
					break
				}
			}
		}
	}
	if solved < 200 || split < 30 {
		t.Errorf("only %d torn solves (%d into several pieces)", solved, split)
	}
	t.Logf("%d torn solves (%d into several pieces), %d singular cuts reported; worst %.2g", solved, split, singular, worst)
}

const tornLadder = `torn ladder
V1 in 0 DC 1 AC 1
R1 in a 1k
C1 a 0 1n
R2 a b 1k
C2 b 0 1n
RT b c 2k
R3 c d 1k
C3 c 0 1n
C4 d 0 2n
R4 d 0 10k
.end
`

// Two pieces joined by RT: the interconnection has order 1 and its
// unknown is RT's current.
func TestTornPieces(t *testing.T) {
	c := parse(t, tornLadder)
	d, err := Tear(c, Tearing{Torn: []string{"RT"}}, nil)
	if err != nil {
		t.Fatal(err)
	}
	sol, st, err := d.Solve(2*math.Pi*1e5, false)
	if err != nil {
		t.Fatal(err)
	}
	sort.Ints(st.Sizes)
	if st.Pieces != 2 || st.Interconnection != 1 || fmt.Sprint(st.Sizes) != "[2 4]" {
		t.Errorf("stats %+v", st)
	}
	ref, err := d.n.Solve(NodeForm, 2*math.Pi*1e5, false)
	if err != nil {
		t.Fatal(err)
	}
	k, _ := d.n.branchOf("RT")
	if e := cmplx.Abs(sol.I.At(k)-ref.I.At(k)) / cmplx.Abs(ref.I.At(k)); e > 1e-12 {
		t.Errorf("RT current %v torn, %v untorn", sol.I.At(k), ref.I.At(k))
	}
}

func TestTornSingularPiece(t *testing.T) {
	c := parse(t, `floating piece
V1 1 0 1
R1 1 2 1k
R2 2 0 1k
R3 2 3 1k
R4 3 4 1k
C4 4 0 1n
.end
`)
	// At the operating point C4 is open, so the piece {3, 4} left by
	// tearing R3 has no reference; in AC C4 gives it one.
	if _, err := TornAC(c, []float64{1e3}, Tearing{Torn: []string{"R3"}}, nil); err != nil {
		t.Fatal(err)
	}
	_, err := TornOP(c, Tearing{Torn: []string{"R3"}}, nil)
	if err == nil {
		t.Fatal("expected a singular piece")
	}
	if !errors.Is(err, ErrSingular) || !strings.Contains(err.Error(), "Lemma 2") || !strings.Contains(err.Error(), "{3 4}") {
		t.Errorf("error %v", err)
	}
}

func TestTornErrors(t *testing.T) {
	c := parse(t, `errs
V1 1 0 DC 1 AC 1
L1 1 2 1m
L2 2 0 1m
K1 L1 L2 0.5
R1 2 3 1k
R2 3 0 1k
.end
`)
	for _, tc := range []struct{ torn, want string }{
		{"V1", "only R, C and L"},
		{"L1", "chord"},
		{"R9", "no element"},
	} {
		if _, err := Tear(c, Tearing{Torn: []string{tc.torn}}, nil); err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("tear %s: %v", tc.torn, err)
		}
	}
	if _, err := Tear(c, Tearing{Torn: []string{"R1", "r1"}}, nil); err == nil {
		t.Error("tearing twice accepted")
	}
	c2 := parse(t, "l\nV1 1 0 DC 1 AC 1\nL1 1 2 1m\nR1 2 0 1k\n.end\n")
	if _, err := TornOP(c2, Tearing{Torn: []string{"L1"}}, nil); err == nil || !strings.Contains(err.Error(), "short circuit") {
		t.Errorf("torn inductor at DC: %v", err)
	}
}

// The cache: a repeated solve hits every piece and is bit-identical; a
// change inside one piece refactors only that piece; concurrent solves
// sharing one cache agree (run with -race).
func TestTornCache(t *testing.T) {
	cache := NewTornCache(0)
	c := parse(t, tornLadder)
	w := 2 * math.Pi * 3e4
	solve := func(c *circuit.Circuit) (*Solution, TearStats) {
		t.Helper()
		d, err := Tear(c, Tearing{Torn: []string{"RT"}}, cache)
		if err != nil {
			t.Fatal(err)
		}
		s, st, err := d.Solve(w, false)
		if err != nil {
			t.Fatal(err)
		}
		return s, st
	}
	s1, st1 := solve(c)
	s2, st2 := solve(c)
	if st1.Misses != 2 || st2.Hits != 2 || st2.Misses != 0 {
		t.Errorf("first %+v, second %+v", st1, st2)
	}
	for k := range s1.E.v {
		if s1.E.v[k] != s2.E.v[k] || s1.I.v[k] != s2.I.v[k] {
			t.Fatalf("cached solve differs at branch %d", k)
		}
	}
	nc := parse(t, strings.Replace(tornLadder, "C3 c 0 1n", "C3 c 0 1.5n", 1))
	_, st3 := solve(nc)
	if st3.Hits != 1 || st3.Misses != 1 {
		t.Errorf("after a change in one piece: %+v", st3)
	}
	nocache, _, err := func() (*Solution, TearStats, error) {
		d, _ := Tear(nc, Tearing{Torn: []string{"RT"}}, nil)
		return d.Solve(w, false)
	}()
	if err != nil {
		t.Fatal(err)
	}
	s3, _ := solve(nc)
	for k := range s3.E.v {
		if s3.E.v[k] != nocache.E.v[k] {
			t.Fatalf("cached and uncached solves differ at branch %d", k)
		}
	}
	var wg sync.WaitGroup
	for g := 0; g < 8; g++ {
		wg.Add(1)
		go func(g int) {
			defer wg.Done()
			src := strings.Replace(tornLadder, "R4 d 0 10k", fmt.Sprintf("R4 d 0 %dk", 10+g%3), 1)
			cc := parse(t, src)
			for i := 0; i < 20; i++ {
				d, _ := Tear(cc, Tearing{Torn: []string{"RT"}}, cache)
				s, _, err := d.Solve(w, false)
				if err != nil {
					t.Error(err)
					return
				}
				ref, _ := d.n.Solve(NodeForm, w, false)
				if wr, _ := worstRel(s.E.v, ref.E.v); wr > twinTol {
					t.Errorf("goroutine %d: off by %g", g, wr)
				}
			}
		}(g)
	}
	wg.Wait()
	if h, m, n := cache.Stats(); h == 0 || m == 0 || n == 0 {
		t.Errorf("cache stats %d %d %d", h, m, n)
	}
}

// Sparse assembly of Y′ equals the dense tensor transform Aᵀ·Y·A (Kron
// 1939 Ch. XIV p. 360) on networks without apparent axes.
func TestNodeAssemblyEqualsTransform(t *testing.T) {
	c := parse(t, `rc
I1 0 1 DC 1 AC 1
R1 1 2 1k
C1 2 0 1n
R2 2 3 2k
G1 3 0 1 2 1m
R3 3 0 3k
C2 1 3 5p
.end
`)
	n, err := Build(c)
	if err != nil {
		t.Fatal(err)
	}
	rel := n.relations(2*math.Pi*1e5, false)
	ns, err := n.assembleNode(rel, nil)
	if err != nil {
		t.Fatal(err)
	}
	if ns.na != 0 {
		t.Fatalf("%d apparent axes", ns.na)
	}
	sys := ns.dense()
	b := len(n.branches)
	Y := NewMat[Primitive, Primitive](b, b)
	for k := 0; k < b; k++ {
		for _, tm := range ns.yrow[k] {
			Y.Set(k, tm.k, tm.c)
		}
	}
	A := n.g.incidence()
	Yp, err := Transform(A, Y)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < ns.nv; i++ {
		for j := 0; j < ns.nv; j++ {
			if d := cmplx.Abs(Yp.At(i, j) - sys[i][j]); d > 1e-15*(1+cmplx.Abs(sys[i][j])) {
				t.Errorf("Y′[%d][%d] = %v, assembled %v", i, j, Yp.At(i, j), sys[i][j])
			}
		}
	}
}
