// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command evolve runs a genetic-algorithm experiment from a JSON config.
//
//	go run ./cmd/evolve -config examples/rc_lowpass_values/config.json
//	go run ./cmd/evolve -config cfg.json -resume runs/20261024-183000-x
//
// Outputs go to runs/<timestamp>-<name>/ (config copy, stats.csv, best
// netlist per generation, best.cir, checkpoint.json; with "algorithm":
// "nsga2" also front.csv and front/*.cir, the final Pareto front).
package main

import (
	"context"
	"flag"
	"fmt"
	"os"
	"os/signal"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/fdamador/kronspice/evo"
)

// version is set at release build time (-ldflags "-X main.version=v1.2.3").
var version = "dev"

func main() {
	showVersion := flag.Bool("version", false, "print the version and exit")
	cfgPath := flag.String("config", "", "experiment config (JSON)")
	out := flag.String("out", "", "base directory for run outputs (default: the config's out_dir, or runs)")
	resume := flag.String("resume", "", "run directory to resume from its checkpoint.json")
	workers := flag.Int("workers", 0, "worker goroutines (default: config, or all CPUs)")
	quiet := flag.Bool("quiet", false, "print only the final result")
	flag.Parse()
	if *showVersion {
		fmt.Println("evolve", version)
		return
	}
	if *cfgPath == "" {
		fmt.Fprintln(os.Stderr, "evolve: -config is required")
		flag.Usage()
		os.Exit(2)
	}
	if err := run(*cfgPath, *out, *resume, *workers, *quiet); err != nil {
		fmt.Fprintln(os.Stderr, "evolve:", err)
		os.Exit(1)
	}
}

func run(cfgPath, out, resume string, workers int, quiet bool) error {
	ex, err := evo.LoadExperiment(cfgPath)
	if err != nil {
		return err
	}
	if workers > 0 {
		ex.Options.Workers = workers
	}
	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt)
	defer stop()

	var st *evo.State
	var rd *evo.RunDir
	if resume != "" {
		if st, err = evo.LoadCheckpoint(filepath.Join(resume, "checkpoint.json")); err != nil {
			return err
		}
		if rd, err = evo.OpenRunDir(resume, ex, st); err != nil {
			return err
		}
	} else {
		base := out
		if base == "" {
			base = ex.Config.OutDir
		}
		if base == "" {
			base = "runs"
		}
		if rd, err = evo.NewRunDir(base, ex, time.Now()); err != nil {
			return err
		}
	}
	defer rd.Close()

	labels := make([]string, len(ex.Problem.Space))
	for i, g := range ex.Problem.Space {
		labels[i] = g.Label()
	}
	observe := func(s *evo.State) error {
		if err := rd.Observe(s); err != nil {
			return err
		}
		if !quiet && len(s.History) > 0 {
			h := s.History[len(s.History)-1]
			unit := "dec"
			if ex.Problem.Topology != nil {
				unit = "parts"
			}
			front := ""
			if ex.Options.Algorithm == evo.NSGA2 {
				front = fmt.Sprintf("  front %d", h.FrontSize)
			}
			fmt.Printf("gen %4d  best %-12.6g mean %-12.6g failed %3d  diversity %.3g %s%s\n",
				h.Generation, h.Best, h.Mean, h.Failed, h.Diversity, unit, front)
		}
		return nil
	}
	st, err = evo.Run(ctx, ex.Problem, ex.Options, st, observe)
	if st != nil {
		b := st.Best
		fmt.Printf("\nbest fitness %g after generation %d (%d evaluations)\n", b.Fitness, st.Generation, st.Evaluations)
		if ex.Problem.Topology != nil {
			if c, _, err := ex.Problem.Circuit(b.Genome); err == nil {
				n := len(ex.Problem.Template.Elements)
				fmt.Printf("  %d evolved parts:\n", len(c.Elements)-n)
				for _, e := range c.Elements[n:] {
					fmt.Printf("  %-8s %-6s %-6s %g\n", e.Name, c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]], e.Value)
				}
			}
		} else {
			for i, v := range b.Values {
				fmt.Printf("  %-12s %g\n", labels[i], v)
			}
		}
		for i, ob := range ex.Problem.Objectives {
			if i < len(b.Terms) {
				fmt.Printf("  objective %-28s error %g (weight %g)\n", ob.Name, b.Terms[i], ob.Weight)
			}
		}
		if b.Failure != "" {
			fmt.Printf("  failed: %s: %s\n", b.Failure, b.Error)
		}
		if b.ModelSupplied > 0 {
			fmt.Printf("  model-supplied power %g W (idealized active elements: a real circuit needs a supply for it)\n", b.ModelSupplied)
		}
		if ex.Options.Algorithm == evo.NSGA2 {
			fmt.Printf("\nPareto front: %d non-dominated designs (front.csv, front/*.cir)\n", len(st.Front))
			fmt.Printf("  %-6s", "member")
			for _, ob := range ex.Problem.Objectives {
				fmt.Printf(" %14.14s", ob.Name)
			}
			fmt.Println()
			front := append([]evo.Individual(nil), st.Front...)
			sort.SliceStable(front, func(a, b int) bool { return front[a].Terms[0] < front[b].Terms[0] })
			for i, ind := range front {
				fmt.Printf("  %-6d", i)
				for _, v := range ind.Terms {
					fmt.Printf(" %14.6g", v)
				}
				fmt.Println()
			}
		}
		fmt.Printf("outputs in %s\n", strings.TrimSuffix(rd.Path, string(filepath.Separator)))
	}
	return err
}
