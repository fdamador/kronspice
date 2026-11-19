// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func TestMapElitesCells(t *testing.T) {
	m := &MapElites{Descriptors: []Descriptor{
		DescParts(0, 10, 5),       // bins of 2 parts
		DescTerm("e", 0, 0, 1, 4), // bins of 0.25
	}}
	p := Problem{}
	at := func(parts int, term float64) Individual {
		return Individual{Values: make([]float64, parts), Terms: []float64{term}, Fitness: term}
	}
	for _, c := range []struct {
		ind  Individual
		want int
		ok   bool
	}{
		{at(0, 0), 0, true},
		{at(3, 0.3), 1 + 1*5, true},
		{at(99, 7), 4 + 3*5, true}, // beyond the range: end bins
		{Individual{Failure: "singular"}, 0, false},
		{Individual{Values: []float64{1}}, 0, false}, // no term: declined
	} {
		got, ok := m.cell(p, c.ind)
		if ok != c.ok || (ok && got != c.want) {
			t.Errorf("%+v: cell %d %v, want %d %v", c.ind, got, ok, c.want, c.ok)
		}
	}
	// A cell keeps its best; a tie keeps the earlier candidate.
	arch := map[int]Individual{}
	a, b, d := at(1, 0.1), at(1, 0.05), at(1, 0.05)
	a.Error, b.Error, d.Error = "a", "b", "d"
	m.insert(p, arch, []Individual{a, b, d})
	if len(arch) != 1 || arch[0].Error != "b" {
		t.Errorf("archive %+v", arch)
	}
	if err := (&MapElites{}).validate(1); !errors.Is(err, ErrBadOptions) {
		t.Errorf("no descriptors: %v", err)
	}
	bad := &MapElites{Descriptors: []Descriptor{DescParts(1, 1, 3)}}
	if err := bad.validate(1); !errors.Is(err, ErrBadOptions) {
		t.Errorf("empty range: %v", err)
	}
}

func TestMapElitesCount(t *testing.T) {
	top, emb := ampTopology(), ampEmbryo(t)
	p := Problem{Topology: top, Template: emb}
	g := Genome{Parts: []Part{{Kind: 'Q', A: 2, B: 0, C: 1}, {Kind: 'Q', A: 3, B: 2, C: 1}, {Kind: 'R', A: 3, B: 1, X: 4}}}
	v, ok := DescCount(circuit.BJT, 0, 4, 4).Of(p, Individual{Genome: g})
	if !ok || v != 2 {
		t.Errorf("BJT count %g %v, want 2", v, ok)
	}
}

// A MAP-Elites run on the filter keeps a growing archive of distinct
// cells, never loses a cell's quality, and resumes bit-identically.
func TestMapElitesRun(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "filter_topology", "config.json"))
	if err != nil {
		t.Fatal(err)
	}
	opt := ex.Options
	opt.Population, opt.Generations = 60, 8
	opt.Algorithm = MAPElites
	opt.MapElites = &MapElites{Descriptors: []Descriptor{
		DescParts(0, 20, 10),
		DescCount(circuit.Inductor, 0, 6, 6),
		DescTerm("mask", 0, 0, 2, 8),
	}}
	var cells []int
	var best []float64
	full, err := Run(context.Background(), ex.Problem, opt, nil, func(st *State) error {
		cells = append(cells, st.History[len(st.History)-1].Cells)
		best = append(best, st.Best.Fitness)
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	for i := 1; i < len(cells); i++ {
		if cells[i] < cells[i-1] || best[i] > best[i-1] {
			t.Fatalf("archive shrank or got worse: cells %v best %v", cells, best)
		}
	}
	if cells[len(cells)-1] < 10 {
		t.Errorf("only %d cells filled: %v", cells[len(cells)-1], cells)
	}
	seen := map[int]bool{}
	for _, ind := range full.Population {
		k, ok := opt.MapElites.cell(ex.Problem, ind)
		if !ok || seen[k] {
			t.Fatalf("archive member out of place or duplicated cell %d", k)
		}
		seen[k] = true
	}
	stop := errors.New("stop")
	path := filepath.Join(t.TempDir(), "cp.json")
	_, err = Run(context.Background(), ex.Problem, opt, nil, func(st *State) error {
		if st.Generation == 3 {
			if err := SaveCheckpoint(path, st); err != nil {
				return err
			}
			return stop
		}
		return nil
	})
	if !errors.Is(err, stop) {
		t.Fatal(err)
	}
	cp, err := LoadCheckpoint(path)
	if err != nil {
		t.Fatal(err)
	}
	resumed, err := Run(context.Background(), ex.Problem, opt, cp, nil)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(full.Population, resumed.Population) || !reflect.DeepEqual(full.Best, resumed.Best) {
		t.Fatal("resumed MAP-Elites run differs")
	}
	t.Logf("cells per generation %v, best %v", cells, best)
	// Local search is refused with MAP-Elites.
	opt.LocalSearch = &LocalSearch{}
	if _, err := Run(context.Background(), ex.Problem, opt, nil, nil); !errors.Is(err, ErrBadOptions) {
		t.Errorf("local search with MAP-Elites: %v", err)
	}
}

// ga.precheck changes no result: the developmental filter run with
// kron.Precheck on equals the run without it, candidate for candidate.
func TestPrecheckConfigSameResult(t *testing.T) {
	ex, err := LoadExperiment(filepath.Join("..", "examples", "filter_topology", "config_developmental.json"))
	if err != nil {
		t.Fatal(err)
	}
	ex.Options.Population, ex.Options.Generations = 60, 10
	plain, err := Run(context.Background(), ex.Problem, ex.Options, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	ex.Config.GA.Precheck = &PrecheckConfig{DC: true}
	raw := ex.Config
	ex2, err := NewExperiment(raw, "Low-pass embryo\nVIN vs 0 DC 0 AC 1\nRS vs in 1k\nRL out 0 1k\n.ac dec 20 10 100k\n")
	if err != nil {
		t.Fatal(err)
	}
	if ex2.Problem.Precheck == nil || !ex2.Problem.Precheck.DC {
		t.Fatal("precheck not set from the config")
	}
	ex2.Options.Population, ex2.Options.Generations = 60, 10
	checked, err := Run(context.Background(), ex2.Problem, ex2.Options, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	for i := range plain.Population {
		a, b := plain.Population[i], checked.Population[i]
		if a.Fitness != b.Fitness && !(a.Failure != "" && b.Failure != "") {
			t.Fatalf("candidate %d: %g (%s) without precheck, %g (%s) with", i, a.Fitness, a.Failure, b.Fitness, b.Failure)
		}
	}
	if plain.Best.Fitness != checked.Best.Fitness {
		t.Errorf("best %g without, %g with precheck", plain.Best.Fitness, checked.Best.Fitness)
	}
}
