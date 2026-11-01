// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Every golden circuit — linear and nonlinear — solved in Kron's node-pair
// frame with state-dependent primitive tensors gives the MNA operating
// point, and where the netlist has an .ac line, the MNA small-signal
// response.
func TestNonlinearKronTwin(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	nonlinear, worstOP, worstAC := 0, 0.0, 0.0
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".cir")
		c, d := load(t, f)
		if _, err := Build(c); err != nil {
			nonlinear++
		}
		// MNA run with tight Newton tolerances, so both solutions are
		// exact to round-off and can be compared at 1e−9.
		mna, err := analysis.OP(c, analysis.Options{RelTol: 1e-12, VnTol: 1e-15, AbsTol: 1e-18})
		if err != nil {
			t.Fatalf("%s MNA: %v", name, err)
		}
		k, err := NonlinearOP(c, NonlinearOptions{})
		if err != nil {
			t.Errorf("%s Kron: %v", name, err)
			continue
		}
		// Voltages against the larger of their scale and 1 V (a circuit
		// biased at 0 V has femtovolt round-off for a scale).
		w := 0.0
		at := 0
		scale := math.Max(maxAbs(real2c(mna.NodeVoltages)), 1)
		for i := range mna.NodeVoltages {
			if e := math.Abs(k.NodeVoltages[i]-mna.NodeVoltages[i]) / scale; e > w {
				w, at = e, i
			}
		}
		worstOP = math.Max(worstOP, w)
		if w > 1e-9 {
			t.Errorf("%s OP: node %s off by %.3g", name, mna.NodeNames[at], w)
		}
		for _, a := range d.Analyses {
			if a.Kind != "ac" {
				continue
			}
			ref, err := analysis.AC(c, analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, analysis.Options{})
			if err != nil {
				t.Fatal(err)
			}
			got, err := NonlinearAC(c, ref.Freq, NonlinearOptions{})
			if err != nil {
				t.Errorf("%s AC: %v", name, err)
				continue
			}
			for fi := range ref.Freq {
				var g, w0 []complex128
				for n := 1; n < len(ref.V); n++ {
					g, w0 = append(g, got.V[n][fi]), append(w0, ref.V[n][fi])
				}
				w, at := worstRel(g, w0)
				worstAC = math.Max(worstAC, w)
				if w > 1e-9 {
					t.Errorf("%s AC at %g Hz: node %s off by %.3g", name, ref.Freq[fi], ref.NodeNames[at+1], w)
					break
				}
			}
		}
	}
	if nonlinear < 15 {
		t.Errorf("only %d nonlinear golden circuits", nonlinear)
	}
	t.Logf("%d golden circuits (%d nonlinear): worst OP difference %.2g, worst AC %.2g", len(files), nonlinear, worstOP, worstAC)
}

// Kron's transformation of the primitive device tensors, Y′ = Cᵀ·Y·C with
// Y the direct sum of every device's Y_d and C their incidence, done as
// typed tensors, equals the scattered assembly used by the solve.
func TestNonlinearTransformLaw(t *testing.T) {
	for _, f := range []string{"bjt_ce_op.cir", "mos1_cs_ac.cir", "diode_fwd.cir"} {
		c, _ := load(t, filepath.Join("..", "testdata", "netlists", f))
		prims, x, err := NonlinearPrimitives(c, NonlinearOptions{})
		if err != nil {
			t.Fatal(err)
		}
		n := len(x)
		P := 0
		for _, p := range prims {
			P += len(p.Axes)
		}
		Y := NewMat[Primitive, Primitive](P, P)
		C := NewMat[Primitive, Node](P, n)
		off := 0
		for _, p := range prims {
			for a, ga := range p.Axes {
				C.Set(off+a, ga, 1)
				for b := range p.Axes {
					Y.Set(off+a, off+b, complex(p.Y[a][b], 0))
				}
			}
			off += len(p.Axes)
		}
		Yp, err := Transform(C, Y)
		if err != nil {
			t.Fatal(err)
		}
		S, _ := connect(n, prims)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				if d := cmplx.Abs(Yp.At(i, j) - complex(S[i][j], 0)); d > 1e-12*(1+math.Abs(S[i][j])) {
					t.Errorf("%s: Y′[%d][%d] = %v, assembled %g", f, i, j, Yp.At(i, j), S[i][j])
				}
			}
		}
		t.Logf("%s: %d devices, %d primitive axes → %d network axes", f, len(prims), P, n)
	}
}

// A machine is one more primitive network in Kron's frame: the torque
// operating point of analysis.TestMachineTorqueOP, solved here.
func TestNonlinearKronMachine(t *testing.T) {
	s := primitiveMachine(10e-3, 12e-3, 8e-3, 1)
	s.Theta0 = 0.7
	c, err := circuit.NewBuilder("machine OP").Domain(exchange.Mechanical, "w").
		MachineModel("mm", s).
		CurrentSource("I1", "0", "a", 2).CurrentSource("I3", "0", "c", 3).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}, {"c", "0"}, {"d", "0"}}, "w", "0").
		Resistor("RB", "b", "0", 1e9).Resistor("RD", "d", "0", 1e9).Resistor("RW", "w", "0", 4).Build()
	if err != nil {
		t.Fatal(err)
	}
	ref, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	k, err := NonlinearOP(c, NonlinearOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if w, _ := worstRel(real2c(k.NodeVoltages), real2c(ref.NodeVoltages)); w > 1e-9 {
		t.Errorf("machine OP off by %g", w)
	}
	T := -8e-3 * math.Sin(0.7) * 2 * 3 // −M·sin θ·I_α·I_a
	if w := k.NodeVoltages[c.NodeIndex("w")]; math.Abs(w-4*T) > 1e-9*math.Abs(4*T) {
		t.Errorf("ω = %g, want %g", w, 4*T)
	}
}
