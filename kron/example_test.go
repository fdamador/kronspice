// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron_test

import (
	"fmt"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
)

// A Wheatstone bridge solved in Kron's mesh frame and node-pair frame.
func Example() {
	d, _ := netlist.ParseString("bridge\nV1 1 0 10\nR1 1 2 1k\nR2 1 3 2k\nR3 2 0 3k\nR4 3 0 4k\nR5 2 3 5k\n.op\n")
	c, _ := circuit.FromDeck(d)
	n, _ := kron.Build(c)
	fmt.Printf("%+v\n", n.Counts())
	for _, form := range []kron.Form{kron.MeshForm, kron.NodeForm} {
		r, _ := kron.OP(c, form)
		fmt.Printf("%s: V(2) = %.6f V, V(3) = %.6f V\n", form, r.NodeVoltages[2], r.NodeVoltages[3])
	}
	// Output:
	// {Branches:6 Junctions:4 SubNetworks:1 Meshes:3 JunctionPairs:3}
	// mesh: V(2) = 7.411765 V, V(3) = 6.823529 V
	// node: V(2) = 7.411765 V, V(3) = 6.823529 V
}
