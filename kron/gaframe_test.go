// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/solver"
)

// passiveCircuit is a random connected network of resistors, capacitors,
// inductors (each in series with a resistor) and independent sources.
func passiveCircuit(r *rand.Rand, nodes int) string {
	var b strings.Builder
	b.WriteString("passive\n")
	node := func() int { return r.IntN(nodes + 1) }
	pair := func() (int, int) {
		a, c := node(), node()
		for c == a {
			c = node()
		}
		return a, c
	}
	val := func(lo, hi float64) float64 { return math.Pow(10, lo+(hi-lo)*r.Float64()) }
	for i := 0; i < nodes; i++ {
		fmt.Fprintf(&b, "Rr%d %d %d %g\n", i, i, i+1, val(1, 5))
	}
	for i := 0; i < nodes+2; i++ {
		a, c := pair()
		switch r.IntN(3) {
		case 0:
			fmt.Fprintf(&b, "R%d %d %d %g\n", i, a, c, val(0, 6))
		case 1:
			fmt.Fprintf(&b, "C%d %d %d %g\n", i, a, c, val(-12, -6))
		case 2:
			fmt.Fprintf(&b, "L%d %d m%d %g\nRl%d m%d %d %g\n", i, a, i, val(-9, -3), i, i, c, val(0, 3))
		}
	}
	a, c := pair()
	fmt.Fprintf(&b, "V1 %d v1 DC 1\nRv1 v1 %d 50\n", a, c)
	a, c = pair()
	fmt.Fprintf(&b, "I1 %d %d DC 1m\n", a, c)
	b.WriteString(".end\n")
	return b.String()
}

func relErr(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(math.Max(math.Abs(a), math.Abs(b)), 1e-300)
}

// The blade of the mesh axes squares to det(Cᵀ·Z·C) and the blade of the
// node-pair axes to det(Aᵀ·Y·A), three ways: the explicit outer product in
// the branch algebra, the product of the rejections, and Gaussian
// elimination on Kron's transformed tensor (CLAUDE.md §3b rule 1, M9e).
func TestGABladeIsDetZPrime(t *testing.T) {
	r := rand.New(rand.NewPCG(11, 12))
	checked := 0
	for trial := 0; trial < 300; trial++ {
		c := parse(t, passiveCircuit(r, 2+r.IntN(4)))
		nw, err := Build(c)
		if err != nil {
			t.Fatal(err)
		}
		mesh, node, err := nw.GAFrames()
		if err != nil {
			t.Fatal(err)
		}
		for _, f := range []*GAFrame{mesh, node} {
			det := GramDet(f.Gram())
			_, rej, ok := f.Independence()
			if !ok {
				t.Fatal("passive frame has a negative weight")
			}
			if relErr(rej, det) > 1e-9 {
				t.Fatalf("%s frame: rejections give %.12g, det(Gram) %.12g", f.Kind, rej, det)
			}
			if len(f.Branches) <= maxExplicitBranches {
				bl, alg, err := f.Blade()
				if err != nil {
					t.Fatal(err)
				}
				b2 := 1.0
				if alg != nil {
					b2 = alg.Norm2(bl)
				}
				if relErr(b2, det) > 1e-9 {
					t.Fatalf("%s frame: ⟨B·B̃⟩ = %.12g, det(Gram) = %.12g", f.Kind, b2, det)
				}
				if g := bl.Grades(); len(f.Axes) > 0 && det != 0 && (len(g) != 1 || g[0] != len(f.Axes)) {
					t.Fatalf("%s blade has grades %v, want %d", f.Kind, g, len(f.Axes))
				}
				checked++
			}
		}
	}
	if checked < 400 {
		t.Errorf("only %d explicit blades checked", checked)
	}
}

// Kron's duality between the two frames for a network of resistors: the
// mesh and node blades satisfy det(Cᵀ·R·C) = det(Aᵀ·G·A)·Π R_k (both are
// sums over spanning trees: Kirchhoff's matrix-tree theorem), and for unit
// resistors det(Aᵀ·A) counts the spanning trees — 16 for the complete
// graph on four nodes (Cayley, n^(n−2)).
func TestGAFrameDuality(t *testing.T) {
	r := rand.New(rand.NewPCG(13, 14))
	for trial := 0; trial < 100; trial++ {
		nodes := 2 + r.IntN(5)
		var b strings.Builder
		b.WriteString("resistors\n")
		k := 0
		for i := 0; i < nodes; i++ {
			fmt.Fprintf(&b, "R%d %d %d %g\n", k, i, i+1, math.Pow(10, 3*r.Float64()))
			k++
		}
		for j := 0; j < nodes+2; j++ {
			a, c := r.IntN(nodes+1), r.IntN(nodes+1)
			if a == c {
				continue
			}
			fmt.Fprintf(&b, "R%d %d %d %g\n", k, a, c, math.Pow(10, 3*r.Float64()))
			k++
		}
		b.WriteString(".end\n")
		nw, err := Build(parse(t, b.String()))
		if err != nil {
			t.Fatal(err)
		}
		mesh, node, _ := nw.GAFrames()
		prodR := 1.0
		for _, z := range mesh.Metric {
			prodR *= z
		}
		dm, dn := GramDet(mesh.Gram()), GramDet(node.Gram())
		if relErr(dm, dn*prodR) > 1e-9 {
			t.Fatalf("det Z′ = %.12g, det Y′·ΠR = %.12g", dm, dn*prodR)
		}
	}
	k4 := parse(t, "K4\nR1 0 1 1\nR2 0 2 1\nR3 0 3 1\nR4 1 2 1\nR5 1 3 1\nR6 2 3 1\n.end\n")
	nw, _ := Build(k4)
	_, node, _ := nw.GAFrames()
	bl, alg, err := node.Blade()
	if err != nil {
		t.Fatal(err)
	}
	if got := alg.Norm2(bl); math.Abs(got-16) > 1e-12 {
		t.Errorf("K4 node blade² = %g, want 16 spanning trees", got)
	}
}

// The linear golden circuits' frames: blade² = det(Gram), and every
// golden circuit (linear or not) passes the precheck — it must never
// reject a circuit the simulator solves.
func TestGAFrameGoldenCircuits(t *testing.T) {
	files, _ := filepath.Glob("../testdata/netlists/*.cir")
	linear := 0
	for _, f := range files {
		c, _ := load(t, f)
		if _, err := Precheck(c, PrecheckOptions{DC: true}); err != nil {
			t.Errorf("%s: precheck rejected a golden circuit: %v", filepath.Base(f), err)
		}
		nw, err := Build(c)
		if err != nil {
			continue
		}
		mesh, node, err := nw.GAFrames()
		if err != nil {
			continue
		}
		linear++
		for _, fr := range []*GAFrame{mesh, node} {
			det := GramDet(fr.Gram())
			if _, rej, _ := fr.Independence(); relErr(rej, det) > 1e-9 {
				t.Errorf("%s %s: rejections %.10g, det %.10g", filepath.Base(f), fr.Kind, rej, det)
			}
		}
	}
	if linear < 8 {
		t.Errorf("only %d linear golden circuits had GA frames", linear)
	}
}

// Degenerate topologies are rejected before any matrix is built, with the
// loop named; every rejected circuit is also singular for MNA (no false
// positives), checked on random circuits seeded with source loops.
func TestPrecheck(t *testing.T) {
	cases := []struct {
		src, want string
		dc        bool
	}{
		{"V loop\nV1 1 0 1\nV2 1 0 2\nR1 1 0 1k\n.end\n", "V2, V1", false},
		{"E loop\nV1 1 0 1\nE1 2 0 1 0 2\nV2 2 1 1\nR1 2 0 1k\n.end\n", "V2", false},
		{"V-L loop\nV1 1 0 1\nL1 1 0 1m\nR1 1 0 1k\n.end\n", "L1, V1", true},
	}
	for _, tc := range cases {
		c := parse(t, tc.src)
		_, err := Precheck(c, PrecheckOptions{DC: tc.dc})
		if !errors.Is(err, ErrSingular) || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%q: got %v, want ErrSingular naming %s", strings.SplitN(tc.src, "\n", 2)[0], err, tc.want)
		}
		if _, opErr := analysis.OP(c, analysis.Options{}); opErr == nil {
			t.Errorf("%q: MNA solved a circuit the precheck rejects", tc.src)
		}
	}
	// The V–L loop passes without the DC option: it is fine in AC and UIC
	// transient analyses.
	if _, err := Precheck(parse(t, cases[2].src), PrecheckOptions{}); err != nil {
		t.Errorf("V–L loop rejected without DC: %v", err)
	}
	// Ill-conditioning is reported, not rejected.
	rep, err := Precheck(parse(t, "ill\nV1 1 0 1\nR1 1 2 1e-6\nR2 2 0 1e9\nR3 2 3 1\nR4 3 0 1\n.end\n"), PrecheckOptions{DC: true, Conditioning: true})
	if err != nil || !rep.Linear || math.Min(rep.MeshMin, rep.NodeMin) > 1e-3 {
		t.Errorf("ill-conditioned network: %+v, %v", rep, err)
	}
	// Random circuits with extra ideal sources: every rejection is singular
	// for MNA.
	r := rand.New(rand.NewPCG(15, 16))
	rejected := 0
	for trial := 0; trial < 300; trial++ {
		src := passiveCircuit(r, 2+r.IntN(4))
		nodes := 6
		for j := 0; j < 2; j++ {
			a, c := r.IntN(nodes), r.IntN(nodes)
			if a == c {
				continue
			}
			src = strings.Replace(src, ".end\n", fmt.Sprintf("Vx%d %d %d DC %g\n.end\n", j, a, c, r.Float64()), 1)
		}
		d := parseMaybe(src)
		if d == nil {
			continue
		}
		_, err := Precheck(d, PrecheckOptions{DC: true})
		if err == nil {
			continue
		}
		rejected++
		if _, opErr := analysis.OP(d, analysis.Options{}); !errors.Is(opErr, solver.ErrSingularMatrix) {
			t.Fatalf("precheck rejected (%v) but MNA gave %v for\n%s", err, opErr, src)
		}
	}
	if rejected < 10 {
		t.Errorf("only %d random circuits rejected; the test is not exercising the precheck", rejected)
	}
}

// BenchmarkPrecheckVsLU compares the GA rejection pass with a dense LU of
// the same frame's Gram matrix on a 12×12 resistor grid.
func BenchmarkPrecheckVsLU(b *testing.B) {
	var sb strings.Builder
	sb.WriteString("grid\n")
	const n = 12
	id := func(i, j int) int { return i*n + j + 1 }
	k := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if j+1 < n {
				fmt.Fprintf(&sb, "R%d %d %d 1k\n", k, id(i, j), id(i, j+1))
				k++
			}
			if i+1 < n {
				fmt.Fprintf(&sb, "R%d %d %d 1k\n", k, id(i, j), id(i+1, j))
				k++
			}
		}
	}
	fmt.Fprintf(&sb, "Rg 1 0 1\nV1 %d 0 1\n.end\n", id(n-1, n-1))
	c := parse(b, sb.String())
	nw, _ := Build(c)
	_, node, _ := nw.GAFrames()
	b.Run("rejection", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			node.Independence()
		}
	})
	b.Run("gram+LU", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			GramDet(node.Gram())
		}
	})
	b.Run("precheck", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err := Precheck(c, PrecheckOptions{DC: true}); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("precheck+conditioning", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err := Precheck(c, PrecheckOptions{DC: true, Conditioning: true}); err != nil {
				b.Fatal(err)
			}
		}
	})
	b.Run("MNA-OP", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			if _, err := analysis.OP(c, analysis.Options{}); err != nil {
				b.Fatal(err)
			}
		}
	})
}

// parseMaybe parses a netlist, nil when it does not make a valid circuit.
func parseMaybe(src string) *circuit.Circuit {
	d, err := netlist.ParseString(src)
	if err != nil {
		return nil
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil
	}
	return c
}
