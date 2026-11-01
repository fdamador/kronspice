// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis_test

import (
	"fmt"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

func ExampleOP() {
	c, err := circuit.NewBuilder("divider").
		VoltageSource("V1", "in", "0", 10).
		Resistor("R1", "in", "out", 1e3).
		Resistor("R2", "out", "0", 1e3).
		Build()
	if err != nil {
		fmt.Println(err)
		return
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		fmt.Println(err)
		return
	}
	v, _ := r.Voltage("out")
	fmt.Printf("v(out) = %g V\n", v)
	fmt.Printf("inflow %g W, dissipated %g W\n", r.Ledger.Inflow, r.Ledger.Dissipation)
	// Output:
	// v(out) = 5 V
	// inflow 0.05 W, dissipated 0.05 W
}
