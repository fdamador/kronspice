// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// M9g: the open-system ledger balances on every golden circuit in every
// analysis its netlist requests, plus the operating point (CLAUDE.md
// §3c; back-filled over M2–M9), and Tellegen's sum Σ v·i vanishes at the
// operating point of every golden circuit made of two-terminal elements.
func TestGoldenLedgers(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	worst := map[string]float64{}
	count := map[string]int{}
	note := func(kind string, rr float64) {
		count[kind]++
		worst[kind] = math.Max(worst[kind], rr)
	}
	tellegen := 0
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".cir")
		src, err := os.ReadFile(f)
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
		if c.Validate() == nil {
			r, err := OP(c, Options{})
			if err != nil {
				t.Fatalf("%s: %v", name, err)
			}
			note("op", r.Ledger.RelativeResidual())
			if r.Ledger.Scale > 0 {
				note("op raw", math.Abs(r.Ledger.Residual)/r.Ledger.Scale)
			}
			cur := make([]float64, len(c.Elements))
			two := true
			for k, e := range c.Elements {
				if len(e.Extra) > 0 || e.Kind == circuit.Mutual {
					two = false
					break
				}
				el, _ := r.Element(e.Name)
				cur[k] = el.Current
			}
			if two {
				sum, scale, err := em.Tellegen(c, r.NodeVoltages, cur)
				if err != nil || math.Abs(sum) > math.Max(1e-12*scale, output.LedgerAbsFloor) {
					t.Errorf("%s: Tellegen Σ v·i = %g of %g (%v)", name, sum, scale, err)
				}
				tellegen++
			}
		}
		for _, a := range d.Analyses {
			switch a.Kind {
			case "dc":
				var sw []DCSweep
				for _, x := range a.Sweeps {
					sw = append(sw, DCSweep{Source: x.Source, Start: x.Start, Stop: x.Stop, Step: x.Step})
				}
				r, err := DC(c, sw, Options{})
				if err != nil {
					t.Fatalf("%s: %v", name, err)
				}
				note("dc", r.MaxRelativeResidual)
			case "ac":
				r, err := AC(c, ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, Options{})
				if err != nil {
					t.Fatalf("%s: %v", name, err)
				}
				note("ac", r.Ledger.MaxRelativeResidual)
			case "tran":
				wf, err := Tran(c, TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC}, TranOptions{})
				if err != nil {
					t.Fatalf("%s: %v", name, err)
				}
				note("tran", wf.Ledger.MaxRelativeResidual)
				raw := 0.0
				for k, res := range wf.Ledger.Residual {
					sc := wf.Ledger.Inflow[k] + wf.Ledger.Outflow[k] + wf.Ledger.Dissipation[k] + math.Abs(wf.Ledger.StorageRate[k]) + math.Abs(wf.Ledger.ModelSupplied[k])
					if sc > 0 {
						raw = math.Max(raw, math.Abs(res)/sc)
					}
				}
				note("tran raw", raw)
			}
		}
	}
	// "raw" is |residual| / scale point by point before the round-off
	// floor is applied: logged, not limited (a point where every power is
	// itself round-off can show a large raw ratio).
	limit := map[string]float64{"op": 1e-12, "dc": 1e-9, "ac": 1e-9, "tran": 1e-9, "op raw": math.Inf(1), "tran raw": math.Inf(1)}
	for kind, w := range worst {
		if w > limit[kind] {
			t.Errorf("%s: worst ledger residual %g > %g", kind, w, limit[kind])
		}
		t.Logf("%-4s %2d runs, worst relative residual %.2g", kind, count[kind], w)
	}
	t.Logf("Tellegen checked on %d two-terminal circuits", tellegen)
	if count["op"] < 30 || count["tran"] < 10 || count["ac"] < 4 || tellegen < 10 {
		t.Errorf("coverage: %v, Tellegen %d", count, tellegen)
	}
}
