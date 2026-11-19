// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// An amplifier grown from an embryo that holds only a 12 V rail, a source
// behind 1 kΩ and a 10 kΩ load (examples/amplifier_topology): parts R, C
// and a 2N3904, target 20 dB at 1 kHz. The written best netlist is
// re-simulated from scratch and must give 20 ± 0.5 dB with a transistor in
// its forward-active region and a closed AC ledger.
func TestAmplifierTopologyExample(t *testing.T) {
	if testing.Short() || raceOn {
		t.Skip("about 50 s: 200 generations of 1000")
	}
	checkGrownAmplifier(t, "config.json")
}

// The same with Koza's developmental encoding and MAP-Elites over four
// descriptors (transistor count, gain error, bias error, program size).
func TestAmplifierTopologyMapElites(t *testing.T) {
	if testing.Short() || raceOn {
		t.Skip("about 65 s: 200 batches of 1000")
	}
	checkGrownAmplifier(t, "config_map_elites.json")
}

// The same with Kron's connection-tensor encoding: a fixed library of six
// R, three C and two Q slots, evolved as their interconnection.
func TestAmplifierTopologyKron(t *testing.T) {
	if testing.Short() || raceOn {
		t.Skip("about 3 min: 200 batches of 1000")
	}
	checkGrownAmplifier(t, "config_kron.json")
}

// With ECE2 coils allowed (κ_a = 0.75, register entry H48) and no bias
// objective, evolution drops the transistor: a coil and a resistor reach
// 20 dB, all of the gain's power booked as model-supplied by the
// hypothetical scalar term. Beside it, the standard prediction: the same
// circuit with the ECE term removed (Maxwell) is a passive divider and
// cannot exceed 0 dB.
func TestAmplifierTopologyECE2(t *testing.T) {
	if testing.Short() || raceOn {
		t.Skip("about 100 s: 200 generations of 1000")
	}
	c, db := grownAmplifier(t, "config_ece2.json")
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 10, FStart: 1e3, FStop: 1e3}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	supplied := 0.0
	for _, p := range r.Ledger.ModelSupplied {
		supplied += p
	}
	if !(supplied > 0) {
		t.Errorf("gain %.2f dB with no model-supplied power booked", db)
	}
	m := c.Clone()
	m.ECE = nil
	for i := range m.Elements {
		delete(m.Elements[i].Params, "kappa_a")
		delete(m.Elements[i].Params, "kappa_b")
	}
	rm, err := analysis.AC(m, analysis.ACSpec{Variation: "dec", Points: 10, FStart: 1e3, FStop: 1e3}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(rm, "out", "0")
	dbm := 20 * math.Log10(cmplx.Abs(h[0]))
	if dbm > 0 {
		t.Errorf("standard (κ = 0) gain %.2f dB > 0 dB", dbm)
	}
	t.Logf("ECE2: %.3f dB with %.3g W model-supplied; same circuit under Maxwell: %.2f dB", db, supplied, dbm)
}

func checkGrownAmplifier(t *testing.T, config string) {
	t.Helper()
	c, _ := grownAmplifier(t, config)
	bias, err := ActiveBias(2, 1e-4, 0.1).Eval(&Eval{Ctx: context.Background(), Circuit: c})
	if err != nil || bias != 0 {
		t.Errorf("no transistor in its active region (error %g, %v)", bias, err)
	}
}

// grownAmplifier runs a config as cmd/evolve does, re-simulates the
// written best netlist and checks 20 ± 0.5 dB at 1 kHz and the AC ledger.
func grownAmplifier(t *testing.T, config string) (*circuit.Circuit, float64) {
	t.Helper()
	ex, err := LoadExperiment(filepath.Join("..", "examples", "amplifier_topology", config))
	if err != nil {
		t.Fatal(err)
	}
	ex.Options.Timeout = 30 * time.Second // generous: no candidate may time out (§7)
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 11, 18, 23, 0, 0, 0, time.UTC))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, nil, rd.Observe)
	if err != nil {
		t.Fatal(err)
	}
	rd.Close()
	for _, h := range st.History {
		if h.TimedOut != 0 {
			t.Fatalf("generation %d: %d candidates timed out; the result depends on machine speed", h.Generation, h.TimedOut)
		}
	}
	src, err := os.ReadFile(filepath.Join(rd.Path, "best.cir"))
	if err != nil {
		t.Fatal(err)
	}
	d, err := netlist.ParseString(string(src))
	if err != nil {
		t.Fatalf("%v\n%s", err, src)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 10, FStart: 1e3, FStop: 1e3}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(r, "out", "0")
	db := 20 * math.Log10(cmplx.Abs(h[0]))
	if math.Abs(db-20) > 0.5 {
		t.Errorf("gain at 1 kHz %.3f dB, want 20 ± 0.5", db)
	}
	if r.Ledger.MaxRelativeResidual > 1e-9 {
		t.Errorf("AC ledger residual %g", r.Ledger.MaxRelativeResidual)
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	var lines []string
	for _, e := range c.Elements[len(ex.Problem.Template.Elements):] {
		line := e.Name + " " + c.NodeNames[e.Nodes[0]] + " " + c.NodeNames[e.Nodes[1]]
		if e.Kind == circuit.BJT {
			q, _ := op.Element(e.Name)
			line += " b=" + c.NodeNames[e.Extra[0]]
			line += " (V_CE " + fmtG(q.Voltage) + " V, I_C " + fmtG(q.Current*1e3) + " mA)"
		} else if e.Params["kappa_a"] != 0 || e.Params["kappa_b"] != 0 {
			line += " " + fmtG(e.Value) + " (ECE κ_a " + fmtG(e.Params["kappa_a"]) + ", κ_b " + fmtG(e.Params["kappa_b"]) + ")"
		} else {
			line += " " + fmtG(e.Value)
		}
		lines = append(lines, line)
	}
	t.Logf("best after %d generations (%d evaluations): %.3f dB at 1 kHz, %d parts:\n  %s",
		st.Generation, st.Evaluations, db, len(lines), strings.Join(lines, "\n  "))
	return c, db
}

func fmtG(v float64) string { return strconv.FormatFloat(v, 'g', 4, 64) }
