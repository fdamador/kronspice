// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kronspice simulates a SPICE netlist.
//
//	kronspice [-o out.csv] [-json] [-all-channels] [-solver auto|dense|sparse] [circuit.cir]
//
// It runs the analyses the netlist requests (.op, .dc, .ac, .tran; .op if
// none) and prints node voltages, element operating points and the energy
// ledger. The netlist is read from the named file, or from standard input
// when no file is given or the name is "-". With -o it writes the last
// analysis as CSV; for a transient, the waveforms are resampled on the
// .tran print step. With -json it writes one JSON document to standard
// output instead of the text report (output.Document): every analysis in
// order, and on failure the error with its class, so that scripts and AI
// assistants can drive the simulator without parsing prose.
package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/solver"
)

// version is set at release build time (-ldflags "-X main.version=v1.2.3").
var version = "dev"

// options are the command's flags.
type options struct {
	csvPath     string
	allChannels bool
	json        bool
	backend     solver.Backend
}

func main() {
	showVersion := flag.Bool("version", false, "print the version and exit")
	out := flag.String("o", "", "write results as CSV to this file (the last analysis in the netlist)")
	all := flag.Bool("all-channels", false, "list every external channel of every element, including unmodeled ones")
	jsonOut := flag.Bool("json", false, "write a JSON document to stdout instead of the text report")
	backend := flag.String("solver", "auto", "linear solver: auto (dense below 30 unknowns, sparse above), dense or sparse")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "usage: kronspice [flags] [circuit.cir | -]\n  (no file or \"-\": read the netlist from standard input)\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if *showVersion {
		fmt.Println("kronspice", version)
		return
	}
	if flag.NArg() > 1 {
		flag.Usage()
		os.Exit(2)
	}
	be, ok := map[string]solver.Backend{"auto": solver.Auto, "dense": solver.DenseLU, "sparse": solver.SparseLU}[*backend]
	if !ok {
		fmt.Fprintf(os.Stderr, "kronspice: unknown -solver %q (auto, dense or sparse)\n", *backend)
		os.Exit(2)
	}
	in := io.Reader(os.Stdin)
	if path := flag.Arg(0); path != "" && path != "-" {
		f, err := os.Open(path)
		if err != nil {
			fmt.Fprintln(os.Stderr, "kronspice:", err)
			os.Exit(1)
		}
		defer f.Close()
		in = f
	}
	opt := options{csvPath: *out, allChannels: *all, json: *jsonOut, backend: be}
	if err := run(in, os.Stdout, opt); err != nil {
		if !opt.json {
			fmt.Fprintln(os.Stderr, "kronspice:", err)
		}
		os.Exit(1)
	}
}

// run reads a netlist, runs its analyses and writes the report (text or
// JSON) to w. In JSON mode an error is reported inside the document as well
// as returned, and analyses completed before it are kept.
func run(in io.Reader, w io.Writer, opt options) error {
	doc := &output.Document{Tool: "kronspice", Version: version, Analyses: []output.AnalysisDoc{}}
	err := simulate(in, w, opt, doc)
	if opt.json {
		doc.OK = err == nil
		if err != nil {
			doc.Error = &output.ErrorInfo{Class: errorClass(err), Message: err.Error()}
		}
		if werr := output.WriteJSON(w, doc); werr != nil {
			return werr
		}
	}
	return err
}

func simulate(in io.Reader, w io.Writer, opt options, doc *output.Document) error {
	deck, err := netlist.Parse(in)
	if err != nil {
		return err
	}
	doc.Title = deck.Title
	c, err := circuit.FromDeck(deck)
	if err != nil {
		return err
	}
	analyses := deck.Analyses
	if len(analyses) == 0 {
		fmt.Fprintln(os.Stderr, "kronspice: no analysis requested; running .op")
		analyses = []netlist.Analysis{{Kind: "op"}}
	}
	text := w
	if opt.json {
		text = io.Discard
	}
	var csvWrite func(io.Writer) error
	for i, a := range analyses {
		if i > 0 {
			fmt.Fprintln(text)
		}
		switch a.Kind {
		case "op":
			r, err := analysis.OP(c, analysis.Options{Backend: opt.backend})
			if err != nil {
				return err
			}
			if err := output.WriteReport(text, r, output.ReportOptions{AllChannels: opt.allChannels}); err != nil {
				return err
			}
			doc.Analyses = append(doc.Analyses, output.AnalysisDoc{Kind: "op", OP: output.OPDocument(r)})
			csvWrite = func(f io.Writer) error { return output.WriteCSV(f, r) }
		case "tran":
			spec := analysis.TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC}
			wf, err := analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{Backend: opt.backend}})
			if err != nil {
				return err
			}
			if err := output.WriteTranReport(text, wf); err != nil {
				return err
			}
			step := a.Step
			doc.Analyses = append(doc.Analyses, output.AnalysisDoc{Kind: "tran", Tran: output.TranDocument(wf, step)})
			csvWrite = func(f io.Writer) error { return output.WriteTranCSV(f, wf, step) }
		case "ac":
			r, err := analysis.AC(c, analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}, analysis.Options{Backend: opt.backend})
			if err != nil {
				return err
			}
			if err := output.WriteACReport(text, r); err != nil {
				return err
			}
			doc.Analyses = append(doc.Analyses, output.AnalysisDoc{Kind: "ac", AC: output.ACDocument(r)})
			csvWrite = func(f io.Writer) error { return output.WriteACCSV(f, r) }
		case "dc":
			var sweeps []analysis.DCSweep
			for _, sw := range a.Sweeps {
				sweeps = append(sweeps, analysis.DCSweep{Source: sw.Source, Start: sw.Start, Stop: sw.Stop, Step: sw.Step})
			}
			r, err := analysis.DC(c, sweeps, analysis.Options{Backend: opt.backend})
			if err != nil {
				return err
			}
			if err := output.WriteDCReport(text, r); err != nil {
				return err
			}
			doc.Analyses = append(doc.Analyses, output.AnalysisDoc{Kind: "dc", DC: output.DCDocument(r)})
			csvWrite = func(f io.Writer) error { return output.WriteDCCSV(f, r) }
		}
	}
	if opt.csvPath != "" && csvWrite != nil {
		cf, err := os.Create(opt.csvPath)
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

// errorClass names an error's class with the same words evo.Penalty uses
// for its penalty classes (CLAUDE.md §4), plus "syntax" for netlist errors,
// so a caller can branch on the class instead of the message.
func errorClass(err error) string {
	switch {
	case errors.Is(err, netlist.ErrSyntax), errors.Is(err, netlist.ErrBadValue):
		return "syntax"
	case errors.Is(err, context.DeadlineExceeded), errors.Is(err, context.Canceled):
		return "timeout"
	case errors.Is(err, em.ErrLedgerResidual):
		return "ledger"
	case errors.Is(err, analysis.ErrNoConvergence), errors.Is(err, analysis.ErrTimestepTooSmall), errors.Is(err, analysis.ErrTooManySteps):
		return "no-convergence"
	case errors.Is(err, solver.ErrSingularMatrix), errors.Is(err, solver.ErrNonFinite):
		return "singular"
	case errors.Is(err, circuit.ErrNoGround), errors.Is(err, circuit.ErrFloatingNode), errors.Is(err, circuit.ErrBadValue),
		errors.Is(err, circuit.ErrBadModel), errors.Is(err, circuit.ErrSelfConnected), errors.Is(err, circuit.ErrDuplicateName),
		errors.Is(err, circuit.ErrUnknownKind), errors.Is(err, circuit.ErrEmptyCircuit), errors.Is(err, circuit.ErrTooLarge),
		errors.Is(err, circuit.ErrBadWaveform), errors.Is(err, em.ErrNotQuasiStatic), errors.Is(err, em.ErrOutsideModelValidity),
		errors.Is(err, analysis.ErrBadACSpec), errors.Is(err, analysis.ErrBadDCSpec), errors.Is(err, analysis.ErrBadTranSpec),
		errors.Is(err, analysis.ErrUnknownSource):
		return "invalid"
	}
	return "error"
}
