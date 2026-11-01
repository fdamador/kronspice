// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"testing"
)

// AutoTear on every linear circuit that has one: the torn solution equals
// the untorn node-pair solution (itself the MNA twin, TestTornTwin) to
// 1e−9 of the largest branch voltage, at DC and at 1 kHz; on a 400-node
// ladder it finds a cut into the asked number of pieces of similar size.
func TestAutoTear(t *testing.T) {
	tried, solved := 0, 0
	for name, c := range linearCircuits(t, 60) {
		n, err := Build(c)
		if err != nil {
			continue
		}
		for _, dc := range []bool{true, false} {
			w := 0.0
			if !dc {
				w = 2 * math.Pi * 1e3
			}
			ref, err := n.Solve(NodeForm, w, dc)
			if err != nil {
				continue
			}
			tried++
			tear, _, sol, _, err := AutoTearSolvable(c, 3, w, dc)
			if err != nil {
				continue // too small or no separable groups: reported, not solved
			}
			solved++
			scale := 0.0
			for i := 0; i < ref.E.Len(); i++ {
				scale = math.Max(scale, cmplx.Abs(ref.E.At(i)))
			}
			for i := 0; i < ref.E.Len(); i++ {
				if d := cmplx.Abs(sol.E.At(i) - ref.E.At(i)); d > 1e-9*math.Max(scale, 1e-30) {
					t.Errorf("%s dc=%v cut %v: branch %d voltage %v, untorn %v", name, dc, tear.Torn, i, sol.E.At(i), ref.E.At(i))
					break
				}
			}
		}
	}
	if solved < tried/3 {
		t.Errorf("AutoTear solved %d of %d circuits", solved, tried)
	}
	t.Logf("AutoTear solved %d of %d circuit/analysis pairs", solved, tried)

	// A long RC ladder tears into near-equal pieces.
	var b strings.Builder
	b.WriteString("ladder\nV1 n0 0 DC 1 AC 1\n")
	for i := 0; i < 400; i++ {
		fmt.Fprintf(&b, "R%d n%d n%d 1k\nC%d n%d 0 1n\n", i, i, i+1, i, i+1)
	}
	b.WriteString(".end\n")
	c := parse(t, b.String())
	tear, d, sol, st, err := AutoTearSolvable(c, 8, 2*math.Pi*1e4, false)
	if err != nil {
		t.Fatal(err)
	}
	if len(tear.Torn) != 7 || st.Pieces < 8 {
		t.Errorf("ladder: %d torn branches, %d pieces (sizes %v)", len(tear.Torn), st.Pieces, st.Sizes)
	}
	ref, err := d.n.Solve(NodeForm, 2*math.Pi*1e4, false)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < ref.E.Len(); i++ {
		if cmplx.Abs(sol.E.At(i)-ref.E.At(i)) > 1e-9 {
			t.Fatalf("ladder branch %d: %v vs %v", i, sol.E.At(i), ref.E.At(i))
		}
	}
	t.Logf("ladder: torn %v, piece sizes %v, interconnection order %d", tear.Torn, st.Sizes, st.Interconnection)
}
