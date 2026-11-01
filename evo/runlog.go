// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"
)

// RunDir writes an experiment's outputs (CLAUDE.md §4, experiment outputs):
//
//	runs/<timestamp>-<name>/config.json       copy of the configuration
//	                        stats.csv         one row per generation
//	                        best/genNNNN.cir  best netlist of each generation
//	                        best.cir          best netlist of the run
//	                        checkpoint.json   state after the latest generation
//	                        front.csv         current Pareto front (NSGA-II)
//	                        front/NNN.cir     final Pareto front netlists (NSGA-II)
//
// Best netlists are plain SPICE that ngspice, LTspice and QSpice can open.
type RunDir struct {
	Path string
	ex   *Experiment
	csv  *csv.Writer
	f    *os.File
	last int // generations already written to stats.csv
}

var unsafeName = regexp.MustCompile(`[^A-Za-z0-9._-]+`)

// NewRunDir creates base/<timestamp>-<name>/ and writes the config copy.
func NewRunDir(base string, ex *Experiment, now time.Time) (*RunDir, error) {
	name := unsafeName.ReplaceAllString(ex.Config.Name, "_")
	dir := filepath.Join(base, now.Format("20060102-150405")+"-"+name)
	return openRunDir(dir, ex, false)
}

// OpenRunDir reopens an existing run directory to resume it; statistics
// are appended after the generations already recorded.
func OpenRunDir(dir string, ex *Experiment, resumed *State) (*RunDir, error) {
	rd, err := openRunDir(dir, ex, true)
	if err != nil {
		return nil, err
	}
	rd.last = len(resumed.History)
	return rd, nil
}

func openRunDir(dir string, ex *Experiment, resume bool) (*RunDir, error) {
	if err := os.MkdirAll(filepath.Join(dir, "best"), 0o755); err != nil {
		return nil, err
	}
	rd := &RunDir{Path: dir, ex: ex}
	flags := os.O_CREATE | os.O_WRONLY | os.O_TRUNC
	if resume {
		flags = os.O_CREATE | os.O_WRONLY | os.O_APPEND
	}
	if !resume && ex.Raw != nil {
		if err := os.WriteFile(filepath.Join(dir, "config.json"), ex.Raw, 0o644); err != nil {
			return nil, err
		}
	}
	f, err := os.OpenFile(filepath.Join(dir, "stats.csv"), flags, 0o644)
	if err != nil {
		return nil, err
	}
	rd.f, rd.csv = f, csv.NewWriter(f)
	if !resume {
		hdr := []string{"generation", "best", "mean", "worst", "diversity_decades", "failed", "evaluations", "best_model_supplied_w"}
		if ex.Options.Algorithm == NSGA2 {
			hdr = append(hdr, "front_size")
		}
		for _, g := range ex.Problem.Space {
			hdr = append(hdr, "best_"+g.Label())
		}
		rd.csv.Write(hdr)
	}
	return rd, nil
}

// Observe is an Observer that records each generation.
func (rd *RunDir) Observe(st *State) error {
	for ; rd.last < len(st.History); rd.last++ {
		h := st.History[rd.last]
		ff := func(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }
		row := []string{strconv.Itoa(h.Generation), ff(h.Best), ff(h.Mean), ff(h.Worst), ff(h.Diversity),
			strconv.Itoa(h.Failed), strconv.Itoa(h.Evaluations), ff(h.BestSupplied)}
		if rd.ex.Options.Algorithm == NSGA2 {
			row = append(row, strconv.Itoa(h.FrontSize))
		}
		for _, v := range st.Population[h.BestIndex].Values {
			row = append(row, ff(v))
		}
		if err := rd.csv.Write(row); err != nil {
			return err
		}
		ind := st.Population[h.BestIndex]
		if err := rd.writeNetlist(filepath.Join(rd.Path, "best", fmt.Sprintf("gen%04d.cir", h.Generation)), ind,
			fmt.Sprintf("best of generation %d", h.Generation)); err != nil {
			return err
		}
	}
	rd.csv.Flush()
	if err := rd.csv.Error(); err != nil {
		return err
	}
	if err := rd.writeNetlist(filepath.Join(rd.Path, "best.cir"), st.Best, "best of run"); err != nil {
		return err
	}
	if rd.ex.Options.Algorithm == NSGA2 {
		if err := rd.writeFront(st); err != nil {
			return err
		}
	}
	return SaveCheckpoint(filepath.Join(rd.Path, "checkpoint.json"), st)
}

func (rd *RunDir) writeNetlist(path string, ind Individual, what string) error {
	c, _, err := rd.ex.Problem.circuitFor(ind.Genome)
	if err != nil {
		return nil // an invalid best has no netlist; stats.csv records it
	}
	title := fmt.Sprintf("%s: %s, fitness %g", rd.ex.Config.Name, what, ind.Fitness)
	if ind.ModelSupplied > 0 {
		title += fmt.Sprintf(", model-supplied %g W", ind.ModelSupplied)
	}
	c.Title = strings.ReplaceAll(title, "\n", " ")
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := c.WriteNetlist(f, rd.ex.Control...); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// writeFront rewrites front.csv (one row per Pareto-front member: its
// objective errors and gene values, sorted by the first objective) and,
// once the run is done, one netlist per member under front/.
func (rd *RunDir) writeFront(st *State) error {
	front := sortedFront(st.Front)
	f, err := os.Create(filepath.Join(rd.Path, "front.csv"))
	if err != nil {
		return err
	}
	w := csv.NewWriter(f)
	hdr := []string{"member", "fitness", "model_supplied_w"}
	for _, ob := range rd.ex.Problem.Objectives {
		hdr = append(hdr, ob.Name)
	}
	if rd.ex.Problem.Topology == nil {
		for _, g := range rd.ex.Problem.Space {
			hdr = append(hdr, g.Label())
		}
	}
	w.Write(hdr)
	ff := func(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }
	for i, ind := range front {
		row := []string{strconv.Itoa(i), ff(ind.Fitness), ff(ind.ModelSupplied)}
		for _, v := range ind.Terms {
			row = append(row, ff(v))
		}
		if rd.ex.Problem.Topology == nil {
			for _, v := range ind.Values {
				row = append(row, ff(v))
			}
		}
		w.Write(row)
	}
	w.Flush()
	if err := w.Error(); err != nil {
		f.Close()
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}
	if !st.Done {
		return nil
	}
	dir := filepath.Join(rd.Path, "front")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	for i, ind := range front {
		if err := rd.writeNetlist(filepath.Join(dir, fmt.Sprintf("%03d.cir", i)), ind,
			fmt.Sprintf("Pareto front member %d", i)); err != nil {
			return err
		}
	}
	return nil
}

// sortedFront is a copy of a front sorted by its first objective (stable,
// so ties keep population order): the order of front.csv and front/*.cir.
func sortedFront(f []Individual) []Individual {
	out := append([]Individual(nil), f...)
	sort.SliceStable(out, func(a, b int) bool { return out[a].Terms[0] < out[b].Terms[0] })
	return out
}

// Close flushes and closes the statistics file.
func (rd *RunDir) Close() error {
	rd.csv.Flush()
	return rd.f.Close()
}
