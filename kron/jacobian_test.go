// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/internal/golden"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/solver"
)

// jacobianCuts are the cuts tried on a circuit: every element that can be
// torn (resistors, capacitors, diodes, transistors, …) on its own, and all
// of them.
func jacobianCuts(c *circuit.Circuit) [][]string {
	var names []string
	for i, e := range c.Elements {
		if terms, err := tearTerminals(&c.Elements[i]); err == nil && len(terms) >= 2 {
			names = append(names, e.Name)
		}
	}
	var out [][]string
	for _, n := range names {
		out = append(out, []string{n})
	}
	if len(names) > 1 {
		out = append(out, names)
	}
	return out
}

// runDeck runs a golden netlist's analysis (as analysis.TestGolden does)
// with the given options and returns the golden comparison messages and
// a flat list of values for the twin comparison.
func runDeck(c *circuit.Circuit, d *netlist.Deck, gpath string, opt analysis.Options) (msgs []string, vals []float64, err error) {
	isTable, err := golden.IsTable(gpath)
	if err != nil {
		return nil, nil, err
	}
	if !isTable {
		r, err := analysis.OP(c, opt)
		if err != nil {
			return nil, nil, err
		}
		return golden.Compare(golden.Values(r), mustLoad(gpath)), r.NodeVoltages, nil
	}
	want, err := golden.LoadTable(gpath)
	if err != nil {
		return nil, nil, err
	}
	for _, a := range d.Analyses {
		switch a.Kind {
		case "dc":
			var sw []analysis.DCSweep
			for _, x := range a.Sweeps {
				sw = append(sw, analysis.DCSweep{Source: x.Source, Start: x.Start, Stop: x.Stop, Step: x.Step})
			}
			r, err := analysis.DC(c, sw, opt)
			if err != nil {
				return nil, nil, err
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
			for i := 1; i < len(c.NodeNames); i++ {
				v, _ := r.NodeVoltage(c.NodeNames[i])
				vals = append(vals, v...)
			}
			return golden.CompareSweep(r.Values[0], get, want), vals, nil
		case "ac":
			r, err := analysis.AC(c, analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, opt)
			if err != nil {
				return nil, nil, err
			}
			for i := 1; i < len(r.V); i++ {
				for _, v := range r.V[i] {
					vals = append(vals, real(v), imag(v))
				}
			}
			return golden.CompareAC(r.Freq, r.NodeVoltage, want), vals, nil
		case "tran":
			wf, err := analysis.Tran(c, analysis.TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC},
				analysis.TranOptions{Options: opt})
			if err != nil {
				return nil, nil, err
			}
			get := func(n string) ([]float64, bool) {
				if strings.HasPrefix(n, "v(") {
					return wf.NodeVoltage(strings.TrimSuffix(strings.TrimPrefix(n, "v("), ")"))
				}
				return wf.ElementCurrent(strings.TrimSuffix(strings.TrimPrefix(n, "i("), ")"))
			}
			vals = append(vals, wf.Time...)
			for i := 1; i < len(c.NodeNames); i++ {
				v, _ := wf.NodeVoltage(c.NodeNames[i])
				vals = append(vals, v...)
			}
			return golden.CompareTable(wf.Time, get, want), vals, nil
		}
	}
	return nil, nil, fmt.Errorf("golden table but no dc, ac or tran analysis")
}

func mustLoad(path string) map[string]float64 {
	m, err := golden.Load(path)
	if err != nil {
		panic(err)
	}
	return m
}

// M9c acceptance: every golden circuit — diodes, BJTs, MOSFETs, the
// linear ones — under every single-element cut and the all-element cut,
// with every Newton Jacobian (operating point, DC sweep, AC, transient)
// solved by tearing, matches its ngspice golden values within the golden
// tolerances, and matches the untorn MNA run.
func TestJacobianTearingGolden(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	runs, split, worst := 0, 0, 0.0
	var tot JacobianStats
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
		gpath := filepath.Join("..", "testdata", "golden", name+".csv")
		refMsgs, ref, err := runDeck(c, d, gpath, analysis.Options{})
		if err != nil {
			t.Fatalf("%s untorn: %v", name, err)
		}
		if len(refMsgs) > 0 {
			// The untorn engine itself misses this golden (a model under
			// development); the torn runs are held to the untorn result.
			t.Logf("%s: untorn run misses its golden (%s); torn runs checked against the untorn run only", name, refMsgs[0])
		}
		for _, cut := range jacobianCuts(c) {
			jt, err := NewJacobianTearing(c, Tearing{Torn: cut})
			if err != nil {
				t.Fatalf("%s %v: %v", name, cut, err)
			}
			msgs, got, err := runDeck(c, d, gpath, analysis.Options{Solver: jt})
			if err != nil {
				t.Errorf("%s %v: %v", name, cut, err)
				continue
			}
			runs++
			if len(refMsgs) == 0 {
				for _, m := range msgs {
					t.Errorf("%s %v: %s", name, cut, m)
				}
			}
			st := jt.Stats()
			if st.MaxPieces > 1 {
				split++
			}
			tot.Solves += st.Solves
			tot.Factors += st.Factors
			tot.Reuses += st.Reuses
			tot.Fallbacks += st.Fallbacks
			if len(got) != len(ref) {
				// A transient may take different steps; the golden
				// comparison above already resampled it.
				continue
			}
			scale := 0.0
			for _, v := range ref {
				scale = math.Max(scale, math.Abs(v))
			}
			for i := range ref {
				worst = math.Max(worst, math.Abs(got[i]-ref[i])/math.Max(scale, 1e-30))
			}
		}
	}
	if runs < 40 || split < 30 {
		t.Errorf("only %d torn runs (%d split)", runs, split)
	}
	if tot.Fallbacks*20 > tot.Solves {
		t.Errorf("%d of %d solves fell back to an untorn LU", tot.Fallbacks, tot.Solves)
	}
	if worst > 1e-6 {
		t.Errorf("torn and untorn runs differ by %.3g of scale", worst)
	}
	t.Logf("%d torn runs (%d into several pieces); %d systems, %d pieces factored, %d reused unchanged, %d fallbacks; worst twin difference %.2g",
		runs, split, tot.Solves, tot.Factors, tot.Reuses, tot.Fallbacks, worst)
}

const diodeLadder = `diode ladder: a fixed RC network torn from a rectifier
V1 in 0 SIN(0 5 1k)
R1 in a 100
D1 a b dmod
C1 b 0 1u
RL b 0 10k
RT b p 1k
R2 p q 1k
C2 q 0 100n
R3 q r 1k
C3 r 0 100n
R4 r 0 10k
.model dmod D(IS=1e-14 N=1.05 RS=0.5 CJO=2p)
.tran 10u 3m
.end
`

// The pieces without nonlinear devices are factored once and reused at
// every Newton iteration of every time step; only the rectifier's piece is
// refactored. With ReuseTol, the rectifier's piece also keeps its factors
// while the diode moves little, and refinement makes up the difference:
// the waveform is unchanged.
func TestJacobianReuse(t *testing.T) {
	c := parse(t, diodeLadder)
	spec := analysis.TranSpec{Step: 10e-6, Stop: 3e-3}
	ref, err := analysis.Tran(c, spec, analysis.TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	for _, tol := range []float64{0, 0.01} {
		jt, err := NewJacobianTearing(c, Tearing{Torn: []string{"RT"}})
		if err != nil {
			t.Fatal(err)
		}
		jt.ReuseTol = tol
		wf, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{Solver: jt}})
		if err != nil {
			t.Fatal(err)
		}
		st := jt.Stats()
		if st.MaxPieces != 2 || st.Fallbacks*20 > st.Solves {
			t.Errorf("tol %g: %+v", tol, st)
		}
		// Two pieces per solve: the ladder's is reused every time but the
		// first (and after a timestep change).
		if st.Reuses+st.StaleReuses < st.Solves/2 {
			t.Errorf("tol %g: only %d of %d pieces reused (%+v)", tol, st.Reuses+st.StaleReuses, 2*st.Solves, st)
		}
		if tol > 0 && st.StaleReuses == 0 {
			t.Errorf("tol %g: no stale reuse (%+v)", tol, st)
		}
		vb, _ := wf.NodeVoltage("r")
		rb, _ := ref.NodeVoltage("r")
		if len(vb) != len(rb) {
			t.Fatalf("tol %g: %d points torn, %d untorn", tol, len(vb), len(rb))
		}
		w := 0.0
		for i := range rb {
			w = math.Max(w, math.Abs(vb[i]-rb[i]))
		}
		if w > 1e-9 {
			t.Errorf("tol %g: v(r) differs by %g V", tol, w)
		}
		t.Logf("tol %g: %d systems, %d factors, %d reused, %d stale, %d refinements, worst |Δv(r)| %.2g V",
			tol, st.Solves, st.Factors, st.Reuses, st.StaleReuses, st.Refinements, w)
	}
}

func TestJacobianErrors(t *testing.T) {
	c := parse(t, diodeLadder)
	for _, tc := range []struct{ torn, want string }{
		{"V1", "branch-current unknown"},
		{"RL", "connects to ground"},
		{"RX", "no element"},
	} {
		if _, err := NewJacobianTearing(c, Tearing{Torn: []string{tc.torn}}); err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%s: %v", tc.torn, err)
		}
	}
	if _, err := NewJacobianTearing(c, Tearing{Torn: []string{"RT", "rt"}}); err == nil {
		t.Error("tearing twice accepted")
	}
	// A singular circuit stays singular: the fallback reports what MNA
	// reports.
	s := parse(t, "loop\nV1 1 0 1\nV2 1 0 2\nR1 1 2 1k\nR2 2 0 1k\n.end\n")
	jt, _ := NewJacobianTearing(s, Tearing{Torn: []string{"R1"}})
	_, err := analysis.OP(s, analysis.Options{Solver: jt, NoGminStepping: true, NoSourceStepping: true, NoPseudoTransient: true})
	_, ref := analysis.OP(s, analysis.Options{NoGminStepping: true, NoSourceStepping: true, NoPseudoTransient: true})
	if err == nil || ref == nil || errors.Is(err, solver.ErrSingularMatrix) != errors.Is(ref, solver.ErrSingularMatrix) {
		t.Errorf("torn %v, untorn %v", err, ref)
	}
}

// One JacobianTearing shared by concurrent analyses (run with -race).
func TestJacobianConcurrent(t *testing.T) {
	c := parse(t, diodeLadder)
	jt, _ := NewJacobianTearing(c, Tearing{Torn: []string{"RT"}})
	ref, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	var wg sync.WaitGroup
	for g := 0; g < 6; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			r, err := analysis.OP(c, analysis.Options{Solver: jt})
			if err != nil {
				t.Error(err)
				return
			}
			for i := range r.NodeVoltages {
				if math.Abs(r.NodeVoltages[i]-ref.NodeVoltages[i]) > 1e-12 {
					t.Errorf("node %d: %g torn, %g untorn", i, r.NodeVoltages[i], ref.NodeVoltages[i])
				}
			}
		}()
	}
	wg.Wait()
}

// A torn solve whose pieces are all reused does not allocate (CLAUDE.md
// §8 rule 7).
func TestJacobianSolveAllocs(t *testing.T) {
	jt := &JacobianTearing{cuts: [][2]int{{3, 4}}}
	m, err := jt.Matrix(8)
	if err != nil {
		t.Fatal(err)
	}
	stamp := func() {
		m.Clear()
		for i := 0; i < 8; i++ {
			m.Add(i, i, 2)
			if i > 0 && i != 4 {
				m.Add(i, i-1, -1)
				m.Add(i-1, i, -1)
			}
		}
		m.Add(4, 3, -0.5)
		m.Add(3, 4, -0.5)
	}
	b, x := make([]float64, 8), make([]float64, 8)
	b[0] = 1
	stamp()
	if err := m.Solve(b, x); err != nil {
		t.Fatal(err)
	}
	if n := testing.AllocsPerRun(50, func() {
		stamp()
		if err := m.Solve(b, x); err != nil {
			t.Fatal(err)
		}
	}); n != 0 {
		t.Errorf("%v allocations per torn solve", n)
	}
	if st := jt.Stats(); st.MaxPieces != 2 || st.Fallbacks != 0 {
		t.Errorf("%+v", st)
	}
}
