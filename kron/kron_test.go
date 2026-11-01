// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

func parse(t testing.TB, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// randomCircuit is a random connected linear circuit of n nodes: a ring
// of resistors through ground (so it is connected), random extra R, C, L
// branches (the first two inductors coupled by a K), independent sources
// and controlled sources of all four kinds.
func randomCircuit(r *rand.Rand, nodes int) string {
	var b strings.Builder
	b.WriteString("random\n")
	node := func() int { return r.IntN(nodes + 1) }
	pair := func() (int, int) {
		a := node()
		c := node()
		for c == a {
			c = node()
		}
		return a, c
	}
	val := func(lo, hi float64) float64 { return math.Pow(10, lo+(hi-lo)*r.Float64()) }
	for i := 0; i < nodes; i++ {
		fmt.Fprintf(&b, "Rr%d %d %d %g\n", i, i, i+1, val(1, 5))
	}
	fmt.Fprintf(&b, "Rr%d %d 0 %g\n", nodes, nodes, val(1, 5))
	var inductors []string
	for i := 0; i < nodes; i++ {
		a, c := pair()
		switch r.IntN(3) {
		case 0:
			fmt.Fprintf(&b, "R%d %d %d %g\n", i, a, c, val(0, 6))
		case 1:
			fmt.Fprintf(&b, "C%d %d %d %g\n", i, a, c, val(-12, -6))
		case 2:
			// An inductor in series with a resistor, so DC never shorts a
			// source and inductor loops stay resistive.
			fmt.Fprintf(&b, "L%d %d m%d %g\nRl%d m%d %d %g\n", i, a, i, val(-9, -3), i, i, c, val(0, 3))
			inductors = append(inductors, fmt.Sprintf("L%d", i))
		}
	}
	// Mutual inductance (a chord in Kron's frame) between two of them.
	if len(inductors) >= 2 {
		fmt.Fprintf(&b, "K1 %s %s %g\n", inductors[0], inductors[1], 0.05+0.9*r.Float64())
	}
	a, c := pair()
	fmt.Fprintf(&b, "V1 %d %d DC %g AC %g %g\nRv1 %d %d %g\n", a, nodes+1, val(-1, 1), val(-1, 0), 360*r.Float64(), nodes+1, c, val(1, 3))
	a, c = pair()
	fmt.Fprintf(&b, "I1 %d %d DC %g AC %g\n", a, c, val(-4, -2), val(-4, -2))
	a, c = pair()
	cp, cn := pair()
	fmt.Fprintf(&b, "G1 %d %d %d %d %g\n", a, c, cp, cn, val(-5, -3))
	a, c = pair()
	cp, cn = pair()
	fmt.Fprintf(&b, "E1 %d e1 %d %d %g\nRe1 e1 %d %g\n", a, cp, cn, val(-1, 0.5), c, val(2, 4))
	a, c = pair()
	fmt.Fprintf(&b, "F1 %d %d V1 %g\n", a, c, val(-1, 0.5))
	a, c = pair()
	fmt.Fprintf(&b, "H1 %d h1 V1 %g\nRh1 h1 %d %g\n", a, val(1, 3), c, val(2, 4))
	b.WriteString(".end\n")
	return b.String()
}

// Random linear circuits with every supported element: mesh, node-pair
// and MNA agree to 1e−9 at DC and at three frequencies; Kirchhoff's laws
// and Tellegen's theorem hold for the Kron solutions; the frame power
// e·i equals e′·i′.
func TestTwinRandom(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	checked, singular := 0, 0
	worst := 0.0
	for trial := 0; trial < 300; trial++ {
		src := randomCircuit(r, 3+r.IntN(8))
		c := parse(t, src)
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			singular++ // e.g. a controlled-source loop with no solution: skip
			continue
		}
		var mnaI []float64
		for _, e := range op.Elements {
			mnaI = append(mnaI, e.Current)
		}
		freqs := []float64{10, 1e4, 1e7}
		mna, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: 10, FStop: 10}, analysis.Options{})
		_ = mna
		for _, form := range []Form{MeshForm, NodeForm} {
			k, err := OP(c, form)
			if err != nil {
				t.Fatalf("trial %d %s: %v\n%s", trial, form, err, src)
			}
			w1, _ := worstRel(real2c(k.NodeVoltages), real2c(op.NodeVoltages))
			var kI []float64
			for _, e := range k.Elements {
				kI = append(kI, e.Current)
			}
			w2, at := worstRel(real2c(kI), real2c(mnaI))
			worst = math.Max(worst, math.Max(w1, w2))
			if w1 > twinTol || w2 > twinTol {
				t.Fatalf("trial %d %s OP: voltages %.3g, currents %.3g (%s)\n%s", trial, form, w1, w2, op.Elements[at].Name, src)
			}
			for _, f := range freqs {
				spec := analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}
				m, err := analysis.AC(c, spec, analysis.Options{})
				if err != nil {
					continue
				}
				ka, err := AC(c, m.Freq, form)
				if err != nil {
					t.Fatalf("trial %d %s AC %g: %v\n%s", trial, form, f, err, src)
				}
				var got, want []complex128
				for n := 1; n < len(m.V); n++ {
					got, want = append(got, ka.V[n][0]), append(want, m.V[n][0])
				}
				for e := range m.I {
					got, want = append(got, ka.I[e][0]), append(want, m.I[e][0])
				}
				w, _ := worstRel(got[:len(m.V)-1], want[:len(m.V)-1])
				wi, _ := worstRel(got[len(m.V)-1:], want[len(m.V)-1:])
				worst = math.Max(worst, math.Max(w, wi))
				if w > twinTol || wi > twinTol {
					t.Fatalf("trial %d %s AC %g Hz: voltages %.3g, currents %.3g\n%s", trial, form, f, w, wi, src)
				}
			}
		}
		// Laws and invariants on the Kron solutions themselves.
		n, err := Build(c)
		if err != nil {
			t.Fatal(err)
		}
		for _, form := range []Form{MeshForm, NodeForm} {
			for _, omega := range []float64{0, 2 * math.Pi * 1e4} {
				sol, err := n.Solve(form, omega, omega == 0)
				if err != nil {
					t.Fatal(err)
				}
				scale := 0.0
				for k := 0; k < sol.E.Len(); k++ {
					scale += cmplx.Abs(sol.E.At(k) * sol.I.At(k))
				}
				if cmplx.Abs(sol.Tellegen) > 1e-12*scale+1e-300 {
					t.Fatalf("trial %d %s: Tellegen Σ E·i = %g against %g", trial, form, sol.Tellegen, scale)
				}
				if d := cmplx.Abs(sol.PrimitivePower - sol.FramePower); d > 1e-12*math.Max(cmplx.Abs(sol.PrimitivePower), scale)+1e-300 {
					t.Fatalf("trial %d %s: e·i = %g, e′·i′ = %g", trial, form, sol.PrimitivePower, sol.FramePower)
				}
				// KCL at every node and KVL around every mesh.
				A, C := n.g.incidence(), n.g.loops()
				kcl, _ := TMulVec(A, sol.I)
				kvl, _ := TMulVec(C, sol.E)
				if w := maxAbs(kcl.v); w > 1e-12*maxAbs(sol.I.v)*float64(len(n.branches)) {
					t.Fatalf("trial %d %s: KCL residual %g", trial, form, w)
				}
				if w := maxAbs(kvl.v); w > 1e-12*maxAbs(sol.E.v)*float64(len(n.branches)) {
					t.Fatalf("trial %d %s: KVL residual %g", trial, form, w)
				}
			}
		}
		checked++
	}
	t.Logf("%d random circuits checked (%d without an MNA solution skipped); worst twin error %.3g", checked, singular, worst)
	if checked < 250 {
		t.Errorf("only %d circuits checked", checked)
	}
}

// The connection tensors: Aᵀ·C = 0 (every mesh satisfies KCL), C has full
// column rank (its link rows are the identity), and the counts satisfy
// m = b − N + S, p = N − S (Kron 1939 p. 75), also for a network of
// several sub-networks.
func TestConnectionTensors(t *testing.T) {
	r := rand.New(rand.NewPCG(3, 4))
	for trial := 0; trial < 200; trial++ {
		n := 2 + r.IntN(12)
		b := r.IntN(25)
		from, to := make([]int, b), make([]int, b)
		for k := range from {
			from[k], to[k] = r.IntN(n), r.IntN(n)
		}
		g := newGraph(n, from, to)
		A, C := g.incidence(), g.loops()
		AC, err := TMul(A, C)
		if err != nil {
			t.Fatal(err)
		}
		for _, v := range AC.a {
			if v != 0 {
				t.Fatalf("trial %d: Aᵀ·C ≠ 0", trial)
			}
		}
		_, m := C.Dims()
		S := g.subNetworks()
		if m != b-n+S {
			t.Fatalf("trial %d: m = %d, b − N + S = %d", trial, m, b-n+S)
		}
		trees := 0
		for _, in := range g.tree {
			if in {
				trees++
			}
		}
		if trees != n-S {
			t.Fatalf("trial %d: %d tree branches, N − S = %d", trial, trees, n-S)
		}
		col := 0
		for k := 0; k < b; k++ {
			if g.tree[k] {
				continue
			}
			for j := 0; j < m; j++ {
				want := complex128(0)
				if j == col {
					want = 1
				}
				if C.At(k, j) != want {
					t.Fatalf("trial %d: link rows of C are not the identity", trial)
				}
			}
			col++
		}
	}
	// Counts of a real network: a bridge (6 branches, 4 nodes) has 3 meshes.
	c := parse(t, "bridge\nV1 1 0 10\nR1 1 2 1k\nR2 1 3 2k\nR3 2 0 3k\nR4 3 0 4k\nR5 2 3 5k\n.op\n")
	nw, err := Build(c)
	if err != nil {
		t.Fatal(err)
	}
	if got := nw.Counts(); got != (Counts{Branches: 6, Junctions: 4, SubNetworks: 1, Meshes: 3, JunctionPairs: 3}) {
		t.Fatalf("counts %+v", got)
	}
}

// Kron's transformation law: Z′ = Cₜ·Z·C for a diagonal Z is the sum, over
// branches, of each branch impedance times the outer product of its row of
// C (meshes sharing a branch couple through it).
func TestTransformLaw(t *testing.T) {
	c := parse(t, "bridge\nV1 1 0 10\nR1 1 2 1k\nR2 1 3 2k\nR3 2 0 3k\nR4 3 0 4k\nR5 2 3 5k\n.op\n")
	nw, _ := Build(c)
	C := nw.g.loops()
	b, m := C.Dims()
	Z := NewMat[Primitive, Primitive](b, b)
	for k, br := range nw.branches {
		Z.Set(k, k, complex(br.value, 0))
	}
	Zp, err := Transform(C, Z)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			want := complex128(0)
			for k := 0; k < b; k++ {
				want += C.At(k, i) * Z.At(k, k) * C.At(k, j)
			}
			if Zp.At(i, j) != want || Zp.At(i, j) != Zp.At(j, i) {
				t.Fatalf("Z′[%d][%d] = %v, want %v (symmetric)", i, j, Zp.At(i, j), want)
			}
		}
	}
	if _, err := Mul(Z, Z); err != nil {
		t.Fatal(err)
	}
	if _, err := MulVec(C, NewVec[Mesh](m+1)); !errors.Is(err, ErrDimension) {
		t.Fatalf("dimension mismatch not reported: %v", err)
	}
}

// Broken circuits return errors, never panic: a nonlinear element in the
// mesh frame (the node-pair frame solves it, see nonlinear.go), a loop of
// voltage sources, a unity-gain controlled source across its own control.
func TestKronErrors(t *testing.T) {
	cases := []struct {
		src  string
		want error
	}{
		{"d\nV1 1 0 1\nR1 1 2 1k\nD1 2 0 dm\n.model dm D\n.op\n", ErrNonlinear},
		{"vloop\nV1 1 0 1\nV2 1 0 2\nR1 1 0 1k\n.op\n", ErrSingular},
		// A unity-gain VCVS across its own control pair: V(1) = V(1).
		{"eloop\nE1 1 0 1 0 1\nR1 1 0 1k\n.op\n", ErrSingular},
	}
	for k, tc := range cases {
		c := parse(t, tc.src)
		for _, form := range []Form{MeshForm, NodeForm} {
			_, err := OP(c, form)
			if k == 0 && form == NodeForm {
				if err != nil {
					t.Errorf("diode in the node-pair frame: %v", err)
				}
				continue
			}
			if !errors.Is(err, tc.want) {
				t.Errorf("%s / %s: %v, want %v", strings.SplitN(tc.src, "\n", 2)[0], form, err, tc.want)
			}
		}
	}
}
