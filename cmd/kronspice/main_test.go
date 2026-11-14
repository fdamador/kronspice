// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"bytes"
	"encoding/json"
	"io"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

const divider = `Divider
V1 in 0 DC 10
R1 in out 1k
R2 out 0 1k
.op
.end
`

const rcTran = `RC step
V1 in 0 PULSE(0 1 0 1n 1n 1 2)
R1 in out 1k
C1 out 0 1u
.tran 10u 5m
.end
`

const rcAC = `RC low-pass
V1 in 0 AC 1
R1 in out 1k
C1 out 0 159.15n
.ac dec 2 10 100k
.end
`

const dcSweep = `Sweep
V1 in 0 DC 0
R1 in out 1k
R2 out 0 1k
.dc V1 0 2 1
.end
`

func runJSON(t *testing.T, netlist string) (*output.Document, error) {
	t.Helper()
	var buf bytes.Buffer
	err := run(strings.NewReader(netlist), &buf, options{json: true, backend: solver.Auto})
	var doc output.Document
	if derr := json.Unmarshal(buf.Bytes(), &doc); derr != nil {
		t.Fatalf("output is not JSON: %v\n%s", derr, buf.String())
	}
	return &doc, err
}

func TestJSONOperatingPoint(t *testing.T) {
	doc, err := runJSON(t, divider)
	if err != nil || !doc.OK || doc.Error != nil {
		t.Fatalf("run: err %v, ok %v, error %+v", err, doc.OK, doc.Error)
	}
	if doc.Tool != "kronspice" || doc.Title != "Divider" || len(doc.Analyses) != 1 || doc.Analyses[0].Kind != "op" {
		t.Fatalf("document header: %+v", doc)
	}
	op := doc.Analyses[0].OP
	if op == nil {
		t.Fatal("no op section")
	}
	if v := op.Nodes["out"]; math.Abs(v-5) > 1e-9 {
		t.Errorf("v(out) = %g, want 5", v)
	}
	if _, ground := op.Nodes["0"]; ground {
		t.Error("ground listed among nodes")
	}
	if len(op.Elements) != 3 {
		t.Errorf("elements: %d, want 3", len(op.Elements))
	}
	l := op.Ledger
	if l == nil {
		t.Fatal("no ledger")
	}
	if math.Abs(l.Inflow-0.05) > 1e-9 || math.Abs(l.Dissipation-0.05) > 1e-9 || l.RelativeResidual > 1e-9 {
		t.Errorf("ledger: inflow %g dissipation %g relative residual %g", l.Inflow, l.Dissipation, l.RelativeResidual)
	}
	for _, row := range l.Rows {
		if row.Status == "not modeled" {
			t.Errorf("unmodeled row listed: %+v", row)
		}
	}
}

func TestJSONTransient(t *testing.T) {
	doc, err := runJSON(t, rcTran)
	if err != nil || !doc.OK {
		t.Fatalf("run: %v %+v", err, doc.Error)
	}
	tr := doc.Analyses[0].Tran
	if tr == nil || len(tr.Time) == 0 {
		t.Fatal("no transient data")
	}
	if n := len(tr.Time); n != 501 {
		t.Errorf("resampled points: %d, want 501 on a 10 µs grid to 5 ms", n)
	}
	vout := tr.V["out"]
	if len(vout) != len(tr.Time) {
		t.Fatalf("v(out) has %d points for %d times", len(vout), len(tr.Time))
	}
	// 5 ms is five time constants: 1 − e⁻⁵ = 0.9933.
	if got := vout[len(vout)-1]; math.Abs(got-0.9933) > 2e-3 {
		t.Errorf("v(out) at 5 ms = %g, want 0.9933", got)
	}
	if len(tr.Ledger.StoredEnergy) != len(tr.Time) || tr.Ledger.MaxRelativeResidual > 1e-6 {
		t.Errorf("ledger: %d stored-energy points, max relative residual %g", len(tr.Ledger.StoredEnergy), tr.Ledger.MaxRelativeResidual)
	}
	if tr.Stats.Accepted == 0 {
		t.Error("no accepted steps reported")
	}
}

func TestJSONACAndDC(t *testing.T) {
	doc, err := runJSON(t, rcAC)
	if err != nil || !doc.OK {
		t.Fatalf("ac run: %v %+v", err, doc.Error)
	}
	ac := doc.Analyses[0].AC
	if ac == nil || len(ac.Freq) != 9 {
		t.Fatalf("ac: %+v", ac)
	}
	// Corner at 1 kHz: |H| = 1/√2, phase −45°.
	k := 4 // 10, 31.6, 100, 316, 1000 Hz
	if math.Abs(ac.Freq[k]-1000) > 1e-6 {
		t.Fatalf("freq[4] = %g", ac.Freq[k])
	}
	s := ac.V["out"]
	if math.Abs(s.Mag[k]-1/math.Sqrt2) > 1e-4 || math.Abs(s.PhaseDeg[k]+45) > 0.01 {
		t.Errorf("v(out) at 1 kHz: %g ∠ %g°", s.Mag[k], s.PhaseDeg[k])
	}
	if len(ac.Ledger.Dissipation) != len(ac.Freq) {
		t.Errorf("ac ledger has %d points", len(ac.Ledger.Dissipation))
	}

	doc, err = runJSON(t, dcSweep)
	if err != nil || !doc.OK {
		t.Fatalf("dc run: %v %+v", err, doc.Error)
	}
	dc := doc.Analyses[0].DC
	if dc == nil || len(dc.Sources) != 1 || len(dc.Values[0]) != 3 {
		t.Fatalf("dc: %+v", dc)
	}
	if got := dc.V["out"]; len(got) != 3 || math.Abs(got[2]-1) > 1e-9 {
		t.Errorf("v(out) over the sweep: %v, want [0 0.5 1]", got)
	}
}

func TestJSONErrors(t *testing.T) {
	cases := []struct{ name, netlist, class string }{
		{"syntax", "Bad\nR1 in\n.op\n.end\n", "syntax"},
		{"floating", "Floating\nV1 in 0 DC 1\nR1 in 0 1k\nR2 x y 1k\n.op\n.end\n", "invalid"},
		{"no ground", "No ground\nV1 a b DC 1\nR1 a b 1k\n.op\n.end\n", "invalid"},
	}
	for _, c := range cases {
		doc, err := runJSON(t, c.netlist)
		if err == nil || doc.OK || doc.Error == nil {
			t.Errorf("%s: expected a failure, got err %v doc %+v", c.name, err, doc)
			continue
		}
		if doc.Error.Class != c.class {
			t.Errorf("%s: class %q, want %q (%s)", c.name, doc.Error.Class, c.class, doc.Error.Message)
		}
		if doc.Analyses == nil {
			t.Errorf("%s: analyses should be an empty list, not null", c.name)
		}
	}
}

func TestTextReportFromStdin(t *testing.T) {
	var buf bytes.Buffer
	if err := run(strings.NewReader(divider), &buf, options{backend: solver.Auto}); err != nil {
		t.Fatal(err)
	}
	out := buf.String()
	for _, want := range []string{"Divider", "v(out)", "Energy ledger", "residual"} {
		if !strings.Contains(out, want) {
			t.Errorf("text report lacks %q:\n%s", want, out)
		}
	}
	if strings.Contains(out, `"ok"`) {
		t.Error("text mode wrote JSON")
	}
}

const negTank = `Negative-resistance tank
V1 in 0 AC 1
R1 in a 50
Rn a b -20
L1 b c 1u
C1 c 0 1n
Rleak c 0 1e9
.ac dec 10 1meg 30meg
.end
`

// -smith measures the port with the circuit's own AC source set to zero,
// keeps negative resistance (|Γ| > 1), and writes the extended chart.
func TestSmithFlag(t *testing.T) {
	dir := t.TempDir()
	opt := options{json: true, backend: solver.Auto, smithPlus: "a", smithMinus: "0", z0: 50,
		smithSVG: dir + "/chart.svg", smithCSV: dir + "/smith.csv"}
	var buf bytes.Buffer
	if err := run(strings.NewReader(negTank), &buf, opt); err != nil {
		t.Fatal(err)
	}
	var doc output.Document
	if err := json.Unmarshal(buf.Bytes(), &doc); err != nil {
		t.Fatal(err)
	}
	sm := doc.Analyses[0].Smith
	if sm == nil || sm.Plus != "a" || sm.Z0 != 50 {
		t.Fatalf("smith doc %+v", sm)
	}
	// Port a sees Rn + L1 + (C1 ∥ Rleak) to ground, and R1 to the shorted
	// source: Re Z of the series branch is −20 Ω near resonance, so the
	// parallel 50 ∥ (−20 + jX) is active there.
	if sm.MaxGamma <= 1 || len(sm.ActiveBands) == 0 {
		t.Fatalf("max |Γ| %g, bands %v", sm.MaxGamma, sm.ActiveBands)
	}
	for _, f := range []string{"chart.svg", "smith.csv"} {
		if b, err := os.ReadFile(dir + "/" + f); err != nil || len(b) < 100 {
			t.Fatalf("%s: %v (%d bytes)", f, err, len(b))
		}
	}
	// Text report.
	buf.Reset()
	opt.json, opt.smithSVG, opt.smithCSV = false, "", ""
	if err := run(strings.NewReader(negTank), &buf, opt); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(buf.String(), "Smith chart: port a–0") || !strings.Contains(buf.String(), "Negative resistance") {
		t.Fatalf("report:\n%s", buf.String())
	}
	// No .ac analysis, or an unknown port: an error, not a silent skip.
	if err := run(strings.NewReader(divider), io.Discard, options{backend: solver.Auto, smithPlus: "out", smithMinus: "0", z0: 50}); err == nil {
		t.Fatal("want an error without .ac")
	}
	if _, err := runJSONOpt(t, rcAC, options{json: true, backend: solver.Auto, smithPlus: "nope", smithMinus: "0", z0: 50}); err == nil {
		t.Fatal("want an error for an unknown node")
	}
}

func runJSONOpt(t *testing.T, netlist string, opt options) (*output.Document, error) {
	t.Helper()
	var buf bytes.Buffer
	err := run(strings.NewReader(netlist), &buf, opt)
	var doc output.Document
	if derr := json.Unmarshal(buf.Bytes(), &doc); derr != nil {
		t.Fatalf("output is not JSON: %v\n%s", derr, buf.String())
	}
	if err != nil && (doc.Error == nil || doc.Error.Class != "invalid") {
		t.Fatalf("error class %+v", doc.Error)
	}
	return &doc, err
}
