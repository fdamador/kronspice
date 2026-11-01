// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// Stamps of the four controlled sources (Vlach & Singhal, Table 4.3), with
// output on unknown 0 (node 1), control on unknown 1 (node 2), and for F/H
// a controlling branch at unknown 2.
func TestControlledStamps(t *testing.T) {
	cases := []struct {
		kind circuit.Kind
		want map[[2]int]float64
	}{
		{circuit.VCCS, map[[2]int]float64{{0, 1}: 2}},
		{circuit.CCCS, map[[2]int]float64{{0, 2}: 2}},
		{circuit.VCVS, map[[2]int]float64{{0, 3}: 1, {3, 0}: 1, {3, 1}: -2}},
		{circuit.CCVS, map[[2]int]float64{{0, 3}: 1, {3, 0}: 1, {3, 2}: -2}},
	}
	for _, c := range cases {
		e := circuit.Element{Name: "X", Kind: c.kind, Nodes: [2]int{1, 0}, Value: 2,
			CtrlNodes: [2]int{2, 0}, CtrlSource: "vc"}
		d := NewControlled(e)
		ctx := NewSetupContext(3)
		ctx.AllocBranch() // unknown 2: the controlling source's branch
		ctx.RegisterBranch("vc", 2)
		if err := d.Setup(ctx); err != nil {
			t.Fatal(err)
		}
		if err := d.Link(ctx); err != nil {
			t.Fatal(err)
		}
		n := ctx.Size()
		m := &denseStub{n: n, a: make([]float64, n*n)}
		rhs := make([]float64, n)
		if err := d.Load(&LoadContext{M: m, RHS: rhs}); err != nil {
			t.Fatal(err)
		}
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if got, want := m.a[i*n+j], c.want[[2]int{i, j}]; got != want {
					t.Errorf("%v: A[%d][%d] = %g, want %g", c.kind, i, j, got, want)
				}
			}
		}
		for _, v := range rhs {
			if v != 0 {
				t.Errorf("%v: rhs %v", c.kind, rhs)
			}
		}
		// All delivered power is model-supplied.
		op := OperatingPoint{Voltage: 2, Current: -3, PowerAbsorbed: -6}
		if got := d.ModelSuppliedPower(op); got != 6 {
			t.Errorf("%v: model-supplied %g", c.kind, got)
		}
	}
}
