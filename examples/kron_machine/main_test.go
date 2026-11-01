// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
)

// The motor reaches nearly synchronous speed in both of Kron's frames,
// which agree; the pumped tank grows at Mathieu's rate.
func TestExample(t *testing.T) {
	s1 := Primitive(0.12, 0.12, 0.11)
	s2, err := kron.RotatingFrame(s1, kron.RotorFrame(4, [][2]int{{2, 3}}, 1))
	if err != nil {
		t.Fatal(err)
	}
	var final [2]float64
	for k, spec := range []*circuit.MachineSpec{s1, s2} {
		c, err := Motor(spec)
		if err != nil {
			t.Fatal(err)
		}
		wf, err := analysis.Tran(c, analysis.TranSpec{Step: 5e-4, Stop: 1.5, UIC: true}, analysis.TranOptions{Options: analysis.Options{RelTol: 1e-5}})
		if err != nil {
			t.Fatal(err)
		}
		w, _ := wf.NodeVoltage("w")
		final[k] = w[len(w)-1]
	}
	if final[0] < 300 || final[0] > 314.16 || math.Abs(final[0]-final[1]) > 1e-3*final[0] {
		t.Errorf("final speeds %g (frame I), %g (frame II)", final[0], final[1])
	}
	w0 := 2 * math.Pi * 1e3
	for _, damp := range []float64{0, 0.5} {
		rate, _, _ := parametric(1e-3, 0.05, 1e3, damp)
		if want := 0.05 * w0 / 4 * (1 - damp); math.Abs(rate-want) > 0.01*want {
			t.Errorf("damping %g: growth %g /s, Mathieu %g /s", damp, rate, want)
		}
	}
}
