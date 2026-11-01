// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// The M7 acceptance test, run as `cmd/evolve` runs it: the 2N3904
// common-emitter example evolves, from E-series values only, to the spec
// 30 ± 0.5 dB mid-band gain, −3 dB edges at 100 Hz and 200 kHz within 5 %
// and V_CE ≥ 3 V, checked on the written best netlist re-simulated from
// scratch. The tolerances are those of the 5 % parts the design is built
// from: with E24/E12 values in ~4 % steps and five coupled targets the
// weighted error settles near 0.02, the floor of the grid (single-part
// moves from the GA's best do not improve it), so the spec, not the
// fitness value, is what the test holds. Over seeds 1–12, 10 meet it; the
// other two miss the gain by under 0.1 dB beyond 0.5 dB (doc/M7.md).
func TestAmplifierExample(t *testing.T) {
	if testing.Short() {
		t.Skip("about 25 s: the full 300-generation run")
	}
	ex, err := LoadExperiment(filepath.Join("..", "examples", "amplifier_gain", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 10, 24, 22, 0, 0, 0, time.UTC))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, nil, rd.Observe)
	if err != nil {
		t.Fatal(err)
	}
	rd.Close()
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
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 20, FStart: 10, FStop: 1e8}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(r, "o", "0")
	g, _ := analysis.GainAt(r.Freq, h, 1e4)
	lo, hi, err := analysis.Cutoff3dB(r.Freq, h)
	if err != nil {
		t.Fatal(err)
	}
	if db := 20 * math.Log10(g); math.Abs(db-30) > 0.5 {
		t.Errorf("gain at 10 kHz %.3f dB, want 30 ± 0.5", db)
	}
	if math.Abs(hi-2e5)/2e5 > 0.05 || math.Abs(lo-100)/100 > 0.05 {
		t.Errorf("band edges %.4g Hz and %.4g Hz, want 100 Hz and 200 kHz within 5 %%", lo, hi)
	}
	if r.Ledger.MaxRelativeResidual > 1e-9 {
		t.Errorf("AC ledger residual %g", r.Ledger.MaxRelativeResidual)
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	q, _ := op.Element("Q1")
	if q.Voltage < 3 {
		t.Errorf("V_CE %g V, want ≥ 3 V", q.Voltage)
	}
	t.Logf("best after %d generations: gain %.3f dB, edges %.4g Hz / %.4g Hz, V_CE %s, I_C %s",
		st.Generation, 20*math.Log10(g), lo, hi, output.Eng(q.Voltage, "V"), output.Eng(q.Current, "A"))
}
