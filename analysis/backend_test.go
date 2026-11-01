// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/solver"
)

// rcMesh builds a w×h mesh of resistors with a capacitor from every node to
// ground, a pulsed source (with AC 1 V) driving one corner through 1 kΩ, a
// DC supply on the far corner, and diode clamps along one edge, so OP and
// transient are nonlinear. Values are random but reproducible.
func rcMesh(w, h int, seed int64) (*circuit.Circuit, error) {
	rng := rand.New(rand.NewSource(seed))
	b := circuit.NewBuilder(fmt.Sprintf("rc mesh %dx%d", w, h))
	b.Model("dclamp", "d", map[string]float64{"is": 1e-14, "n": 1.0})
	node := func(x, y int) string { return fmt.Sprintf("n%d_%d", x, y) }
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			if x+1 < w {
				b.Resistor(fmt.Sprintf("Rx%d_%d", x, y), node(x, y), node(x+1, y), 100+900*rng.Float64())
			}
			if y+1 < h {
				b.Resistor(fmt.Sprintf("Ry%d_%d", x, y), node(x, y), node(x, y+1), 100+900*rng.Float64())
			}
			b.Capacitor(fmt.Sprintf("C%d_%d", x, y), node(x, y), "0", 1e-9*(1+rng.Float64()))
		}
	}
	for x := 0; x < w; x++ {
		b.Diode(fmt.Sprintf("D%d", x), node(x, h-1), "0", "dclamp")
	}
	b.Add(circuit.Element{Name: "Vin", Kind: circuit.VoltageSource, Value: 0, HasDC: true, ACMag: 1,
		Wave: circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 2, 1e-7, 1e-8, 1e-8, 5e-7, 2e-6}}}, "in", "0")
	b.Resistor("Rin", "in", node(0, 0), 1000)
	b.VoltageSource("Vdd", "vdd", "0", 1.5)
	b.Resistor("Rdd", "vdd", node(w-1, h-1), 500)
	return b.Build()
}

func maxRel(a, b []float64) float64 {
	worst, scale := 0.0, 1e-30
	for i := range a {
		worst = math.Max(worst, math.Abs(a[i]-b[i]))
		scale = math.Max(scale, math.Abs(b[i]))
	}
	return worst / scale
}

// Dense and sparse backends agree to 1e-9 relative on OP, AC and transient
// of a 100-node nonlinear mesh, and Auto picks sparse at that size.
func TestBackendsAgreeOnLargeCircuit(t *testing.T) {
	c, err := rcMesh(10, 10, 1)
	if err != nil {
		t.Fatal(err)
	}
	s, err := New(c)
	if err != nil {
		t.Fatal(err)
	}
	if s.Backend() != solver.SparseLU {
		t.Fatalf("auto backend at %d unknowns: %v", len(s.rhs), s.Backend())
	}
	ops := map[solver.Backend][]float64{}
	acs := map[solver.Backend][]complex128{}
	trs := map[solver.Backend][]float64{}
	times := map[solver.Backend][]float64{}
	for _, be := range []solver.Backend{solver.DenseLU, solver.SparseLU} {
		r, err := OP(c, Options{Backend: be})
		if err != nil {
			t.Fatalf("%v OP: %v", be, err)
		}
		ops[be] = r.NodeVoltages
		ac, err := AC(c, ACSpec{Variation: "dec", Points: 10, FStart: 1e3, FStop: 1e9}, Options{Backend: be})
		if err != nil {
			t.Fatalf("%v AC: %v", be, err)
		}
		v, _ := ac.NodeVoltage("n9_9")
		acs[be] = v
		wf, err := Tran(c, TranSpec{Step: 1e-8, Stop: 2e-6}, TranOptions{Options: Options{Backend: be}})
		if err != nil {
			t.Fatalf("%v tran: %v", be, err)
		}
		v2, _ := wf.NodeVoltage("n5_5")
		trs[be], times[be] = v2, wf.Time
		if wf.Ledger.MaxRelativeResidual > 1e-9 {
			t.Errorf("%v: ledger residual %g", be, wf.Ledger.MaxRelativeResidual)
		}
	}
	d, sp := solver.DenseLU, solver.SparseLU
	if r := maxRel(ops[sp], ops[d]); r > 1e-9 {
		t.Errorf("OP differs by %g", r)
	}
	worst, scale := 0.0, 0.0
	for i := range acs[d] {
		worst = math.Max(worst, cmplx.Abs(acs[sp][i]-acs[d][i]))
		scale = math.Max(scale, cmplx.Abs(acs[d][i]))
	}
	if worst > 1e-9*scale {
		t.Errorf("AC differs by %g", worst/scale)
	}
	// Round-off can in principle change an adaptive step; require the same
	// time grid here and say so if it ever differs.
	if len(times[sp]) != len(times[d]) || maxRel(times[sp], times[d]) > 1e-12 {
		t.Fatalf("time grids differ: %d vs %d points", len(times[sp]), len(times[d]))
	}
	if r := maxRel(trs[sp], trs[d]); r > 1e-9 {
		t.Errorf("transient differs by %g", r)
	}
}

// Small circuits stay on the dense backend under Auto.
func TestAutoPicksDenseForSmallCircuits(t *testing.T) {
	c, err := rcMesh(3, 3, 1)
	if err != nil {
		t.Fatal(err)
	}
	s, err := New(c)
	if err != nil {
		t.Fatal(err)
	}
	if s.Backend() != solver.DenseLU {
		t.Fatalf("auto backend at %d unknowns: %v", len(s.rhs), s.Backend())
	}
}

// BenchmarkMeshTran measures a whole transient (Newton, LTE control,
// ledger) on meshes of growing size with each backend.
func BenchmarkMeshTran(b *testing.B) {
	for _, w := range []int{4, 8, 16, 24} {
		c, err := rcMesh(w, w, 1)
		if err != nil {
			b.Fatal(err)
		}
		for _, be := range []solver.Backend{solver.DenseLU, solver.SparseLU} {
			b.Run(fmt.Sprintf("%s/nodes=%d", be, w*w), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					if _, err := Tran(c, TranSpec{Step: 2e-8, Stop: 1e-6}, TranOptions{Options: Options{Backend: be}}); err != nil {
						b.Fatal(err)
					}
				}
			})
		}
	}
}
