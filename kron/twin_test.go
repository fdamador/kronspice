// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// twinTol is the twin-test tolerance (CLAUDE.md §3a rule 4): 1e−9
// relative, each value against the largest value of its kind in the same
// solution (a node voltage of 1e−15 V next to 5 V is compared at 5e−9 V).
const twinTol = 1e-9

func load(t testing.TB, path string) (*circuit.Circuit, *netlist.Deck) {
	t.Helper()
	src, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	d, err := netlist.ParseString(string(src))
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c, d
}

// close reports the worst error of got against want relative to the
// largest |want| (and the index where it occurs).
func worstRel(got, want []complex128) (float64, int) {
	scale := maxAbs(want)
	if scale == 0 {
		scale = 1
	}
	w, at := 0.0, -1
	for i := range want {
		if e := cmplx.Abs(got[i]-want[i]) / scale; e > w {
			w, at = e, i
		}
	}
	return w, at
}

func real2c(x []float64) []complex128 {
	out := make([]complex128, len(x))
	for i, v := range x {
		out[i] = complex(v, 0)
	}
	return out
}

// Every linear golden circuit: the operating point (and the AC sweep when
// the netlist has one) from the mesh form and from the node-pair form
// equals the MNA result to 1e−9 relative — node voltages and element
// currents.
func TestTwinGolden(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	linear, worst := 0, 0.0
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".cir")
		c, d := load(t, f)
		if _, err := Build(c); errors.Is(err, ErrNonlinear) {
			continue
		} else if err != nil {
			t.Errorf("%s: %v", name, err)
			continue
		}
		linear++
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			t.Errorf("%s: MNA: %v", name, err)
			continue
		}
		var mnaI []float64
		for _, e := range op.Elements {
			mnaI = append(mnaI, e.Current)
		}
		for _, form := range []Form{MeshForm, NodeForm} {
			r, err := OP(c, form)
			if err != nil {
				t.Errorf("%s/%s: %v", name, form, err)
				continue
			}
			w, at := worstRel(real2c(r.NodeVoltages), real2c(op.NodeVoltages))
			worst = math.Max(worst, w)
			if w > twinTol {
				t.Errorf("%s/%s OP: node %s off by %.3g relative", name, form, op.NodeNames[at], w)
			}
			var kI []float64
			for _, e := range r.Elements {
				kI = append(kI, e.Current)
			}
			w, at = worstRel(real2c(kI), real2c(mnaI))
			worst = math.Max(worst, w)
			if w > twinTol {
				t.Errorf("%s/%s OP: current of %s off by %.3g relative", name, form, op.Elements[at].Name, w)
			}
		}
		for _, a := range d.Analyses {
			if a.Kind != "ac" {
				continue
			}
			spec := analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}
			mna, err := analysis.AC(c, spec, analysis.Options{})
			if err != nil {
				t.Errorf("%s: MNA AC: %v", name, err)
				continue
			}
			for _, form := range []Form{MeshForm, NodeForm} {
				r, err := AC(c, mna.Freq, form)
				if err != nil {
					t.Errorf("%s/%s AC: %v", name, form, err)
					continue
				}
				for k := range mna.Freq {
					var got, want []complex128
					for n := 1; n < len(mna.V); n++ {
						got, want = append(got, r.V[n][k]), append(want, mna.V[n][k])
					}
					w, at := worstRel(got, want)
					worst = math.Max(worst, w)
					if w > twinTol {
						t.Errorf("%s/%s AC at %g Hz: node %s off by %.3g", name, form, mna.Freq[k], mna.NodeNames[at+1], w)
						break
					}
				}
			}
		}
	}
	t.Logf("%d linear golden circuits, mesh and node forms; worst difference from MNA %.3g", linear, worst)
	if linear < 10 {
		t.Errorf("only %d linear golden circuits", linear)
	}
}
