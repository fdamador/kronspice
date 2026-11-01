// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kronspice simulates a SPICE netlist.
//
//	kronspice [-o out.csv] [-all-channels] circuit.cir
//
// It runs the analyses the netlist requests (.op, .dc, .ac, .tran; .op if
// none) and
// prints node voltages, element operating points and the energy ledger. With
// -o it writes the results as CSV; for a transient, the waveforms are
// resampled on the .tran print step.
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

// version is set at release build time (-ldflags "-X main.version=v1.2.3").
var version = "dev"

func main() {
	showVersion := flag.Bool("version", false, "print the version and exit")
	out := flag.String("o", "", "write results as CSV to this file (the last analysis in the netlist)")
	all := flag.Bool("all-channels", false, "list every external channel of every element, including unmodeled ones")
	backend := flag.String("solver", "auto", "linear solver: auto (dense below 30 unknowns, sparse above), dense or sparse")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "usage: kronspice [flags] circuit.cir\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if *showVersion {
		fmt.Println("kronspice", version)
		return
	}
	if flag.NArg() != 1 {
		flag.Usage()
		os.Exit(2)
	}
	be, ok := map[string]solver.Backend{"auto": solver.Auto, "dense": solver.DenseLU, "sparse": solver.SparseLU}[*backend]
	if !ok {
		fmt.Fprintf(os.Stderr, "kronspice: unknown -solver %q (auto, dense or sparse)\n", *backend)
		os.Exit(2)
	}
	if err := run(flag.Arg(0), *out, *all, be); err != nil {
		fmt.Fprintln(os.Stderr, "kronspice:", err)
		os.Exit(1)
	}
}

func run(path, csvPath string, all bool, be solver.Backend) error {
	f, err := os.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()
	deck, err := netlist.Parse(f)
	if err != nil {
		return err
	}
	c, err := circuit.FromDeck(deck)
	if err != nil {
		return err
	}
	analyses := deck.Analyses
	if len(analyses) == 0 {
		fmt.Fprintln(os.Stderr, "kronspice: no analysis requested; running .op")
		analyses = []netlist.Analysis{{Kind: "op"}}
	}
	var csvWrite func(*os.File) error
	for i, a := range analyses {
		if i > 0 {
			fmt.Println()
		}
		switch a.Kind {
		case "op":
			r, err := analysis.OP(c, analysis.Options{Backend: be})
			if err != nil {
				return err
			}
			if err := output.WriteReport(os.Stdout, r, output.ReportOptions{AllChannels: all}); err != nil {
				return err
			}
			csvWrite = func(f *os.File) error { return output.WriteCSV(f, r) }
		case "tran":
			spec := analysis.TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC}
			wf, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{Backend: be}})
			if err != nil {
				return err
			}
			if err := output.WriteTranReport(os.Stdout, wf); err != nil {
				return err
			}
			step := a.Step
			csvWrite = func(f *os.File) error { return output.WriteTranCSV(f, wf, step) }
		case "ac":
			r, err := analysis.AC(c, analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, analysis.Options{Backend: be})
			if err != nil {
				return err
			}
			if err := output.WriteACReport(os.Stdout, r); err != nil {
				return err
			}
			csvWrite = func(f *os.File) error { return output.WriteACCSV(f, r) }
		case "dc":
			var sweeps []analysis.DCSweep
			for _, sw := range a.Sweeps {
				sweeps = append(sweeps, analysis.DCSweep{Source: sw.Source, Start: sw.Start, Stop: sw.Stop, Step: sw.Step})
			}
			r, err := analysis.DC(c, sweeps, analysis.Options{Backend: be})
			if err != nil {
				return err
			}
			if err := output.WriteDCReport(os.Stdout, r); err != nil {
				return err
			}
			csvWrite = func(f *os.File) error { return output.WriteDCCSV(f, r) }
		}
	}
	if csvPath != "" && csvWrite != nil {
		cf, err := os.Create(csvPath)
		if err != nil {
			return err
		}
		if err := csvWrite(cf); err != nil {
			cf.Close()
			return err
		}
		return cf.Close()
	}
	return nil
}
