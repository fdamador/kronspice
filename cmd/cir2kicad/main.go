// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command cir2kicad converts SPICE netlists to KiCad 8 schematics
// (.kicad_sch next to each .cir, or in -o). Connectivity is carried by net
// labels and checked by reading the schematic back (kicad.Read). With
// -results (default) each netlist also gets a results sheet: the MNA
// solution and Kron's automatically torn solution side by side
// (NAME_results.md, NAME_results.csv).
//
//	go run ./cmd/cir2kicad examples/water_vibration/netlists/*.cir
package main

import (
	"bytes"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/fdamador/kronspice/kicad"
	"github.com/fdamador/kronspice/report"
)

var version = "dev"

func main() {
	out := flag.String("o", "", "output directory (default: next to each netlist)")
	results := flag.Bool("results", true, "also write NAME_results.md and NAME_results.csv: MNA and Kron torn solutions")
	pieces := flag.Int("pieces", 4, "pieces asked of the automatic Kron tearing for the results sheet")
	freq := flag.Float64("freq", 0, "AC frequency for the results sheet [Hz] (0: middle of the netlist's .ac sweep, or DC)")
	showVersion := flag.Bool("version", false, "print the version and exit")
	flag.Usage = func() {
		fmt.Fprintf(flag.CommandLine.Output(), "usage: cir2kicad [-o dir] netlist.cir...\n")
		flag.PrintDefaults()
	}
	flag.Parse()
	if *showVersion {
		fmt.Println("cir2kicad", version)
		return
	}
	if flag.NArg() == 0 {
		flag.Usage()
		os.Exit(2)
	}
	failed := false
	for _, f := range flag.Args() {
		if err := convert(f, *out, *results, report.Options{Pieces: *pieces, Freq: *freq}); err != nil {
			fmt.Fprintf(os.Stderr, "%s: %v\n", f, err)
			failed = true
		}
	}
	if failed {
		os.Exit(1)
	}
}

func convert(path, dir string, results bool, ro report.Options) error {
	src, err := os.ReadFile(path)
	if err != nil {
		return err
	}
	base := strings.TrimSuffix(filepath.Base(path), filepath.Ext(path))
	var b bytes.Buffer
	if err := kicad.Write(&b, string(src), kicad.Options{Project: base}); err != nil {
		return err
	}
	if _, err := kicad.Read(b.String()); err != nil {
		return fmt.Errorf("schematic check: %v", err)
	}
	if dir == "" {
		dir = filepath.Dir(path)
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	dst := filepath.Join(dir, base+".kicad_sch")
	if err := os.WriteFile(dst, b.Bytes(), 0o644); err != nil {
		return err
	}
	fmt.Println(dst)
	if !results {
		return nil
	}
	ro.Source = filepath.Base(path)
	sh, err := report.Compute(string(src), ro)
	if err != nil {
		return fmt.Errorf("results sheet: %v", err)
	}
	var md, csv bytes.Buffer
	sh.Markdown(&md, 40)
	sh.CSV(&csv)
	for name, b := range map[string][]byte{base + "_results.md": md.Bytes(), base + "_results.csv": csv.Bytes()} {
		p := filepath.Join(dir, name)
		if err := os.WriteFile(p, b, 0o644); err != nil {
			return err
		}
		fmt.Println(p)
	}
	return nil
}
