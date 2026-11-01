// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// ladder is an n-section RC ladder driven by a voltage source.
func ladder(n int) string {
	var b strings.Builder
	b.WriteString("ladder\nV1 1 0 DC 1 AC 1\n")
	for i := 1; i <= n; i++ {
		fmt.Fprintf(&b, "R%d %d %d 1k\nC%d %d 0 1n\n", i, i, i+1, i, i+1)
	}
	b.WriteString(".end\n")
	return b.String()
}

func BenchmarkOP(b *testing.B) {
	for _, n := range []int{10, 40} {
		d, _ := netlist.ParseString(ladder(n))
		c, _ := circuit.FromDeck(d)
		for _, form := range []Form{MeshForm, NodeForm} {
			b.Run(fmt.Sprintf("%s/%d", form, n), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					if _, err := OP(c, form); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
		b.Run(fmt.Sprintf("mna/%d", n), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				if _, err := analysis.OP(c, analysis.Options{}); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// One AC point of a 40-section ladder: untorn node form, torn at the
// middle resistor with a cold cache (every piece factored), and torn with
// a warm cache (both pieces reused, as for the unchanged pieces of a GA
// candidate).
func BenchmarkTornAC(b *testing.B) {
	d, _ := netlist.ParseString(ladder(40))
	c, _ := circuit.FromDeck(d)
	w := 2 * math.Pi * 1e5
	b.Run("untorn", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			n, _ := Build(c)
			if _, err := n.Solve(NodeForm, w, false); err != nil {
				b.Fatal(err)
			}
		}
	})
	for _, warm := range []bool{false, true} {
		name := "torn/cold"
		if warm {
			name = "torn/warm"
		}
		b.Run(name, func(b *testing.B) {
			cache := NewTornCache(0)
			for i := 0; i < b.N; i++ {
				if !warm {
					cache = nil
				}
				t, err := Tear(c, Tearing{Torn: []string{"R20"}}, cache)
				if err != nil {
					b.Fatal(err)
				}
				if _, _, err := t.Solve(w, false); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// Kron's frames are built for whatever the GA produces: any circuit that
// passes validation is solved or rejected with an error, never a panic.
func FuzzKron(f *testing.F) {
	f.Add(ladder(3))
	f.Add("b\nV1 1 0 10\nR1 1 2 1k\nR2 1 3 2k\nR3 2 0 3k\nR4 3 0 4k\nR5 2 3 5k\n.op\n")
	f.Add("k\nV1 1 0 DC 1 AC 1\nR1 1 2 10\nL1 2 0 1m\nL2 3 0 4m\nK1 L1 L2 0.9\nR2 3 0 100\n.ac dec 3 1 1meg\n")
	f.Add("c\nV1 1 0 DC 1 AC 1\nE1 2 0 1 0 2\nR1 2 3 1k\nG1 3 0 2 0 1m\nF1 3 0 V1 0.5\nH1 4 0 V1 100\nR2 4 0 1k\nL1 3 0 1m\nC1 2 0 1n\n.ac dec 3 1 1meg\n")
	f.Fuzz(func(t *testing.T, src string) {
		d, err := netlist.ParseString(src)
		if err != nil {
			return
		}
		c, err := circuit.FromDeck(d)
		if err != nil {
			return
		}
		for _, form := range []Form{MeshForm, NodeForm} {
			_, _ = OP(c, form)
			_, _ = AC(c, []float64{1, 1e6}, form)
		}
	})
}

// A rectifier driving a fixed RC ladder (diodeLadder, jacobian_test.go):
// the transient by MNA and with every Newton Jacobian torn at RT.
//
// "grid" is the same rectifier driving a 12×12 RC mesh through RT.
func BenchmarkJacobianTran(b *testing.B) {
	for _, plant := range []string{"ladder", "grid"} {
		src := diodeLadder
		if plant == "grid" {
			src = rectifierGrid(12)
		}
		b.Run(plant, func(b *testing.B) { benchJacobianTran(b, src) })
	}
}

// rectifierGrid is a half-wave rectifier driving an n×n RC mesh.
func rectifierGrid(n int) string {
	var s strings.Builder
	s.WriteString("rectifier into a grid\nV1 in 0 SIN(0 5 1k)\nR1 in a 100\nD1 a b dmod\nC1 b 0 1u\nRL b 0 10k\nRT b g0_0 1k\n")
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Fprintf(&s, "CG%d_%d g%d_%d 0 10n\n", i, j, i, j)
			if i+1 < n {
				fmt.Fprintf(&s, "RV%d_%d g%d_%d g%d_%d 100\n", i, j, i, j, i+1, j)
			}
			if j+1 < n {
				fmt.Fprintf(&s, "RH%d_%d g%d_%d g%d_%d 100\n", i, j, i, j, i, j+1)
			}
		}
	}
	fmt.Fprintf(&s, "RG g%d_%d 0 1k\n.model dmod D(IS=1e-14 N=1.05 RS=0.5 CJO=2p)\n.end\n", n-1, n-1)
	return s.String()
}

func benchJacobianTran(b *testing.B, src string) {
	d, _ := netlist.ParseString(src)
	c, _ := circuit.FromDeck(d)
	spec := analysis.TranSpec{Step: 10e-6, Stop: 3e-3}
	b.Run("mna", func(b *testing.B) {
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			if _, err := analysis.Tran(c, spec, analysis.TranOptions{}); err != nil {
				b.Fatal(err)
			}
		}
	})
	for _, tol := range []float64{0, 0.01} {
		b.Run(fmt.Sprintf("torn/reuse=%g", tol), func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				jt, _ := NewJacobianTearing(c, Tearing{Torn: []string{"RT"}})
				jt.ReuseTol = tol
				if _, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{Solver: jt}}); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
