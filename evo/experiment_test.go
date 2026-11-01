// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"encoding/csv"
	"errors"
	"math"
	"os"
	"path/filepath"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// The M5 acceptance test, run exactly as `cmd/evolve` runs it: the RC
// low-pass example reaches its target cutoff within 1 %, every output file
// is written, and the best netlist re-simulates to the same cutoff.
func TestRCLowpassExample(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "rc_lowpass_values", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Date(2026, 10, 24, 18, 0, 0, 0, time.UTC))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, nil, rd.Observe)
	if err != nil {
		t.Fatal(err)
	}
	rd.Close()
	r, c := st.Best.Values[0], st.Best.Values[1]
	if fc := 1 / (2 * math.Pi * r * c); math.Abs(fc-1000) > 10 {
		t.Fatalf("fc = %g Hz", fc)
	}
	for _, f := range []string{"config.json", "stats.csv", "best.cir", "checkpoint.json", "best/gen0000.cir"} {
		if _, err := os.Stat(filepath.Join(rd.Path, f)); err != nil {
			t.Errorf("missing %s: %v", f, err)
		}
	}
	// stats.csv: header + one row per generation.
	f, _ := os.Open(filepath.Join(rd.Path, "stats.csv"))
	rows, err := csv.NewReader(f).ReadAll()
	f.Close()
	if err != nil || len(rows) != len(st.History)+1 {
		t.Fatalf("stats.csv: %d rows, %v", len(rows), err)
	}
	// The best netlist is plain SPICE: re-parse it and measure again.
	src, _ := os.ReadFile(filepath.Join(rd.Path, "best.cir"))
	d, err := netlist.ParseString(string(src))
	if err != nil {
		t.Fatalf("%v\n%s", err, src)
	}
	ck, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	res, err := analysis.AC(ck, analysis.ACSpec{Variation: "dec", Points: 100, FStart: 10, FStop: 1e5}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	h, _ := analysis.Transfer(res, "out", "0")
	_, fhi, _ := analysis.Cutoff3dB(res.Freq, h)
	if math.Abs(fhi-1000)/1000-st.Best.Fitness > 1e-12 {
		t.Fatalf("re-simulated cutoff %g Hz, fitness %g", fhi, st.Best.Fitness)
	}
}

// A run stopped part-way resumes from its run directory and appends to its
// statistics instead of starting over.
func TestResumeFromRunDir(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "negative_resistance_gain", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	ex.Options.HasTarget = false
	ex.Options.Generations = 6
	full, err := Run(context.Background(), ex.Problem, ex.Options, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	rd, err := NewRunDir(t.TempDir(), ex, time.Unix(0, 0))
	if err != nil {
		t.Fatal(err)
	}
	stop := errors.New("stop")
	_, err = Run(context.Background(), ex.Problem, ex.Options, nil, func(st *State) error {
		if err := rd.Observe(st); err != nil {
			return err
		}
		if st.Generation == 2 {
			return stop
		}
		return nil
	})
	rd.Close()
	if !errors.Is(err, stop) {
		t.Fatal(err)
	}
	cp, err := LoadCheckpoint(filepath.Join(rd.Path, "checkpoint.json"))
	if err != nil {
		t.Fatal(err)
	}
	rd2, err := OpenRunDir(rd.Path, ex, cp)
	if err != nil {
		t.Fatal(err)
	}
	st, err := Run(context.Background(), ex.Problem, ex.Options, cp, rd2.Observe)
	rd2.Close()
	if err != nil {
		t.Fatal(err)
	}
	if st.Best.Fitness != full.Best.Fitness || len(st.History) != len(full.History) {
		t.Fatalf("resumed best %g (%d gens), uninterrupted %g (%d gens)", st.Best.Fitness, len(st.History),
			full.Best.Fitness, len(full.History))
	}
	f, _ := os.Open(filepath.Join(rd.Path, "stats.csv"))
	rows, _ := csv.NewReader(f).ReadAll()
	f.Close()
	if len(rows) != 1+len(full.History) {
		t.Fatalf("stats.csv has %d rows, want %d", len(rows), 1+len(full.History))
	}
}

func TestExamplesLoad(t *testing.T) {
	dirs, _ := filepath.Glob(filepath.Join("..", "examples", "*", "config.json"))
	if len(dirs) < 3 {
		t.Fatalf("found %d example configs", len(dirs))
	}
	for _, p := range dirs {
		if _, err := LoadExperiment(p); err != nil {
			t.Errorf("%s: %v", p, err)
		}
	}
}

func TestConfigErrors(t *testing.T) {
	nl := "x\nV1 in 0 AC 1\nR1 in out 1k\nC1 out 0 1n\n.op\n"
	bad := []Config{
		{Genes: []GeneConfig{{Element: "R1", Min: 1, Max: 2}}},                                                               // no objectives
		{Genes: []GeneConfig{{Element: "R1", Min: 1, Max: 2}}, Objectives: []ObjectiveConfig{{Kind: "cutoff", TargetHz: 1}}}, // no .ac
		{Genes: []GeneConfig{{Element: "R1", Min: 1, Max: 2}}, Objectives: []ObjectiveConfig{{Kind: "nope"}}},                // unknown kind
		{Genes: []GeneConfig{{Element: "R1", Min: 1, Max: 2, Sign: "sideways"}}, Objectives: []ObjectiveConfig{{Kind: "op_voltage"}}},
		{Genes: []GeneConfig{{Element: "C1", Min: 1, Max: 2, Sign: "either"}}, Objectives: []ObjectiveConfig{{Kind: "op_voltage"}}},
		{Genes: []GeneConfig{{Element: "R1", Min: 1, Max: 2}}, Objectives: []ObjectiveConfig{{Kind: "op_voltage"}}, GA: GAConfig{Population: 1}},
	}
	for i, c := range bad {
		if _, err := NewExperiment(c, nl); err == nil {
			t.Errorf("config %d accepted", i)
		}
	}
}

func TestControlLines(t *testing.T) {
	got := controlLines("t\n.options reltol=1e-4\n.op\n.tran 1u\n+ 1m\n.ac dec 10 1 1k\n.end\n")
	want := []string{".op", ".tran 1u 1m", ".ac dec 10 1 1k"}
	if len(got) != len(want) {
		t.Fatalf("got %q", got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("got %q, want %q", got, want)
		}
	}
}
