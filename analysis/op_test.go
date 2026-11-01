// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/internal/golden"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/solver"
)

func near(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestDivider(t *testing.T) {
	c, err := circuit.NewBuilder("div").
		VoltageSource("V1", "in", "0", 10).
		Resistor("R1", "in", "out", 1e3).
		Resistor("R2", "out", "0", 1e3).
		Build()
	if err != nil {
		t.Fatal(err)
	}
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("out"); !near(v, 5, 1e-12) {
		t.Fatalf("v(out) = %g", v)
	}
	v1, _ := r.Element("V1")
	if !near(v1.Current, -5e-3, 1e-15) {
		t.Fatalf("i(V1) = %g, want -5mA (current into + terminal)", v1.Current)
	}
	l := r.Ledger
	if !near(l.Inflow, 0.05, 1e-15) || !near(l.Dissipation, 0.05, 1e-15) || l.RelativeResidual() > 1e-12 {
		t.Fatalf("ledger = %+v", l)
	}
	// Every element appears with every channel in the ledger.
	if len(l.Rows) != 3*exchange.NumChannels {
		t.Fatalf("ledger rows = %d", len(l.Rows))
	}
}

func TestCurrentSourceDirection(t *testing.T) {
	// SPICE: "I1 0 a 2m" drives 2 mA from node 0 through the source into node a.
	c, _ := circuit.NewBuilder("i").CurrentSource("I1", "0", "a", 2e-3).Resistor("R1", "a", "0", 1e3).Build()
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("a"); !near(v, 2, 1e-12) {
		t.Fatalf("v(a) = %g, want 2", v)
	}
	if !near(r.Ledger.Inflow, 4e-3, 1e-15) {
		t.Fatalf("inflow = %g", r.Ledger.Inflow)
	}
}

func TestAbsorbingSourceIsOutflow(t *testing.T) {
	// V2 (5 V) is charged by V1 (9 V) through R1: V2 absorbs power, which
	// leaves the circuit through its electrical port.
	c, _ := circuit.NewBuilder("charge").
		VoltageSource("V1", "a", "0", 9).
		VoltageSource("V2", "b", "0", 5).
		Resistor("R1", "a", "b", 100).Build()
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	l := r.Ledger
	i := 4.0 / 100
	if !near(l.Inflow, 9*i, 1e-12) || !near(l.Outflow, 5*i, 1e-12) || !near(l.Dissipation, 4*i, 1e-12) {
		t.Fatalf("ledger = %+v", l)
	}
}

func TestSingularVoltageLoop(t *testing.T) {
	c, err := circuit.NewBuilder("loop").
		VoltageSource("V1", "a", "0", 1).
		VoltageSource("V2", "a", "0", 2).
		Resistor("R1", "a", "0", 1).Build()
	if err != nil {
		t.Fatal(err)
	}
	if _, err := OP(c, Options{}); !errors.Is(err, solver.ErrSingularMatrix) {
		t.Fatalf("err = %v, want ErrSingularMatrix", err)
	}
}

func TestFloatingNodeRejected(t *testing.T) {
	c := &circuit.Circuit{
		NodeNames: []string{"0", "a", "b"},
		Elements: []circuit.Element{
			{Name: "R1", Kind: circuit.Resistor, Nodes: [2]int{1, 0}, Value: 1},
			{Name: "I1", Kind: circuit.CurrentSource, Nodes: [2]int{2, 0}, Value: 1},
		},
	}
	if _, err := OP(c, Options{}); !errors.Is(err, circuit.ErrFloatingNode) {
		t.Fatalf("err = %v, want ErrFloatingNode", err)
	}
}

// TestGolden runs every netlist in testdata/netlists and compares with
// testdata/golden (CLAUDE.md §8 rules 4–5): name,value files for .op and
// time tables for .tran.
func TestGolden(t *testing.T) {
	files, err := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	if err != nil || len(files) == 0 {
		t.Fatalf("no netlists: %v", err)
	}
	for _, f := range files {
		for _, be := range []solver.Backend{solver.DenseLU, solver.SparseLU} {
			name := strings.TrimSuffix(filepath.Base(f), ".cir")
			gpath := filepath.Join("..", "testdata", "golden", name+".csv")
			t.Run(name+"/"+be.String(), func(t *testing.T) { checkGolden(t, f, gpath, be) })
		}
	}
}

// TestNgspiceTwin compares every netlist that has an ngspice reference in
// testdata/golden/ngspice (written by internal/golden/ref/ngspice_all.sh
// from the same netlist, run through ngspice itself) with the same
// tolerances as TestGolden. This is the SPICE3-compatibility check of
// CLAUDE.md §1 goal 4: testdata/golden holds independent references
// (exact solutions, Python re-implementations of the model equations),
// this directory holds what ngspice actually computes. Netlists without
// an ngspice file are skipped, so the test passes before the references
// are generated.
func TestNgspiceTwin(t *testing.T) {
	refs, err := filepath.Glob(filepath.Join("..", "testdata", "golden", "ngspice", "*.csv"))
	if err != nil {
		t.Fatal(err)
	}
	if len(refs) == 0 {
		t.Skip("no ngspice references yet: run internal/golden/ref/ngspice_all.sh in the tools container")
	}
	for _, g := range refs {
		name := strings.TrimSuffix(filepath.Base(g), ".csv")
		f := filepath.Join("..", "testdata", "netlists", name+".cir")
		if _, err := os.Stat(f); err != nil {
			t.Errorf("%s: ngspice reference without a netlist", name)
			continue
		}
		t.Run(name, func(t *testing.T) { checkGolden(t, f, g, solver.Auto) })
	}
}

// checkGolden runs netlist f with backend be and compares it with the
// reference in gpath: a name,value file for .op, or a table for .dc, .ac
// or .tran (CLAUDE.md §8 rules 4–5).
func checkGolden(t *testing.T, f, gpath string, be solver.Backend) {
	t.Helper()
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
	isTable, err := golden.IsTable(gpath)
	if err != nil {
		t.Fatal(err)
	}
	if isTable {
		for _, a := range d.Analyses {
			if a.Kind != "dc" {
				continue
			}
			var sw []DCSweep
			for _, x := range a.Sweeps {
				sw = append(sw, DCSweep{Source: x.Source, Start: x.Start, Stop: x.Stop, Step: x.Step})
			}
			r, err := DC(c, sw, Options{Backend: be})
			if err != nil {
				t.Fatal(err)
			}
			want, err := golden.LoadTable(gpath)
			if err != nil {
				t.Fatal(err)
			}
			get := func(n string) ([]float64, bool) {
				if strings.HasPrefix(n, "v(") {
					return r.NodeVoltage(strings.TrimSuffix(strings.TrimPrefix(n, "v("), ")"))
				}
				name := strings.TrimSuffix(strings.TrimPrefix(n, "i("), ")")
				for i, e := range r.ElementNames {
					if strings.EqualFold(e, name) {
						return r.I[i], true
					}
				}
				return nil, false
			}
			for _, m := range golden.CompareSweep(r.Values[0], get, want) {
				t.Error(m)
			}
			if r.MaxRelativeResidual > 1e-9 {
				t.Errorf("ledger residual %g", r.MaxRelativeResidual)
			}
			return
		}
		for _, a := range d.Analyses {
			if a.Kind != "ac" {
				continue
			}
			r, err := AC(c, ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, Options{Backend: be})
			if err != nil {
				t.Fatal(err)
			}
			want, err := golden.LoadTable(gpath)
			if err != nil {
				t.Fatal(err)
			}
			for _, m := range golden.CompareAC(r.Freq, r.NodeVoltage, want) {
				t.Error(m)
			}
			return
		}
		var spec *netlist.Analysis
		for i := range d.Analyses {
			if d.Analyses[i].Kind == "tran" {
				spec = &d.Analyses[i]
			}
		}
		if spec == nil {
			t.Fatal("golden table but no .tran in netlist")
		}
		wf, err := Tran(c, TranSpec{Step: spec.Step, Stop: spec.Stop, Start: spec.Start, MaxStep: spec.MaxStep, UIC: spec.UIC}, TranOptions{Options: Options{Backend: be}})
		if err != nil {
			t.Fatal(err)
		}
		want, err := golden.LoadTable(gpath)
		if err != nil {
			t.Fatal(err)
		}
		get := func(n string) ([]float64, bool) {
			if strings.HasPrefix(n, "v(") {
				return wf.NodeVoltage(strings.TrimSuffix(strings.TrimPrefix(n, "v("), ")"))
			}
			return wf.ElementCurrent(strings.TrimSuffix(strings.TrimPrefix(n, "i("), ")"))
		}
		compare := golden.CompareTable
		if strings.Contains(filepath.ToSlash(gpath), "/golden/ngspice/") {
			// Dense ngspice reference: check the simulation at its own
			// steps (golden.CompareTableAtSim).
			compare = golden.CompareTableAtSim
		}
		for _, m := range compare(wf.Time, get, want) {
			t.Error(m)
		}
		if wf.Ledger.MaxRelativeResidual > 1e-9 {
			t.Errorf("ledger residual %g", wf.Ledger.MaxRelativeResidual)
		}
		return
	}
	r, err := OP(c, Options{Backend: be})
	if err != nil {
		t.Fatal(err)
	}
	want, err := golden.Load(gpath)
	if err != nil {
		t.Fatal(err)
	}
	for _, m := range golden.Compare(golden.Values(r), want) {
		t.Error(m)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-12 {
		t.Errorf("ledger residual %g", rr)
	}
}

func TestParallelSimsAreIndependent(t *testing.T) {
	errs := make(chan error, 16)
	for g := 0; g < 16; g++ {
		go func(g int) {
			v := float64(g + 1)
			c, err := circuit.NewBuilder("p").VoltageSource("V1", "a", "0", v).
				Resistor("R1", "a", "b", 1e3).Resistor("R2", "b", "0", 1e3).Build()
			if err != nil {
				errs <- err
				return
			}
			s, err := New(c)
			if err != nil {
				errs <- err
				return
			}
			for k := 0; k < 200; k++ {
				r, err := s.OP(Options{})
				if err != nil {
					errs <- err
					return
				}
				if got, _ := r.Voltage("b"); got != v/2 {
					errs <- fmt.Errorf("goroutine %d: v(b) = %g", g, got)
					return
				}
			}
			errs <- nil
		}(g)
	}
	for i := 0; i < 16; i++ {
		if err := <-errs; err != nil {
			t.Fatal(err)
		}
	}
}

func ladder(n int) *circuit.Circuit {
	b := circuit.NewBuilder("ladder").VoltageSource("V1", "n0", "0", 1)
	for i := 0; i < n; i++ {
		b.Resistor(fmt.Sprintf("Rs%d", i), fmt.Sprintf("n%d", i), fmt.Sprintf("n%d", i+1), 1e3)
		b.Resistor(fmt.Sprintf("Rp%d", i), fmt.Sprintf("n%d", i+1), "0", 2e3)
	}
	c, err := b.Build()
	if err != nil {
		panic(err)
	}
	return c
}

func BenchmarkOPLadder30(b *testing.B) {
	s, err := New(ladder(30))
	if err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := s.OP(Options{}); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkSolveLadder30(b *testing.B) {
	s, err := New(ladder(30))
	if err != nil {
		b.Fatal(err)
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if err := s.solve(); err != nil {
			b.Fatal(err)
		}
	}
}
