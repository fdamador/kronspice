// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command bearden_meg: the Motionless Electromagnetic Generator (Patrick,
// Bearden et al., US 6,362,718, 2002; register entry H10) and a two-magnet
// "knife edge" relative, as magnetic circuits in kronspice's magnetic
// domain with the full energy ledger. The claim under test is a
// coefficient of performance above 1 (load power over drive power, the
// magnets supplying the difference); the standard prediction is below 1,
// with the magnets supplying no net power over whole cycles.
//
// Each case is also written as an ngspice twin (ngspice/*.cir, the same
// circuit in ngspice syntax) so the two simulators can be read against
// each other; see ngspice_run.sh and TestNgspiceTwin.
//
//	go run ./examples/bearden_meg            # sweeps and the table
//	go run ./examples/bearden_meg -write     # regenerate netlists/ and ngspice/
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
)

// Sweep is one named group of runs.
type Sweep struct {
	Name  string
	Specs []Spec
}

// Sweeps returns the runs of the study.
func Sweeps() []Sweep {
	var out []Sweep
	add := func(name string, ss ...Spec) { out = append(out, Sweep{name, ss}) }
	var ss []Spec
	for _, v := range []float64{0.5, 2, 5, 20, 50, 100} {
		ss = append(ss, Spec{Kind: Knife, V: v})
	}
	add("knife: drive 0.5-100 V", ss...)
	ss = nil
	for _, r := range []float64{1, 10, 100, 300, 1e3, 3e3, 1e4, 1e5} {
		ss = append(ss, Spec{Kind: Knife, RL: r})
	}
	add("knife: load 1 ohm-100 kohm", ss...)
	ss = nil
	for _, f := range []float64{1e3, 5e3, 10e3, 50e3} {
		ss = append(ss, Spec{Kind: Knife, F: f})
	}
	add("knife: frequency 1-50 kHz", ss...)
	ss = nil
	for _, m := range []float64{0.9, 0.95, 0.99, 1, 1.01, 1.05, 1.1} {
		ss = append(ss, Spec{Kind: Knife, Mismatch: m})
	}
	add("knife: magnet mismatch", ss...)
	ss = nil
	for _, g := range []float64{0.02, 0.05, 0.1, 0.15, 0.3} {
		ss = append(ss, Spec{Kind: Knife, Gap: g})
	}
	add("knife: magnet spacing 2-30 cm", ss...)
	ss = nil
	for _, l := range []float64{1e-8, 3e-8, 1e-7, 3e-7, 1e-6} {
		ss = append(ss, Spec{Kind: Knife, NullLeak: l})
	}
	add("knife: air leakage at the null", ss...)
	ss = nil
	for _, r := range []float64{1e3, 1e4} {
		for _, c := range []float64{30e-9, 60e-9, 80e-9, 100e-9, 150e-9, 300e-9} {
			ss = append(ss, Spec{Kind: Knife, RL: r, Cres: c})
		}
	}
	add("knife: series-resonant drive", ss...)
	ss = nil
	for _, v := range []float64{1, 2, 5, 10, 20} {
		ss = append(ss, Spec{Kind: MEG, V: v})
	}
	add("meg: drive 1-20 V", ss...)
	ss = nil
	for _, r := range []float64{10, 100, 1e3, 1e4} {
		ss = append(ss, Spec{Kind: MEG, RL: r})
	}
	add("meg: load 10 ohm-10 kohm", ss...)
	ss = nil
	for _, f := range []float64{1e3, 5e3, 10e3, 25e3} {
		ss = append(ss, Spec{Kind: MEG, F: f})
	}
	add("meg: frequency 1-25 kHz", ss...)
	return out
}

// Twins are the cases written as kronspice and ngspice netlists and
// compared between the two simulators.
func Twins() []Spec {
	return []Spec{
		{Name: "knife_baseline", Kind: Knife},
		{Name: "knife_20V", Kind: Knife, V: 20},
		{Name: "knife_resonant", Kind: Knife, RL: 1e3, Cres: 80e-9},
		{Name: "knife_leak", Kind: Knife, NullLeak: 1e-7},
		{Name: "meg_baseline", Kind: MEG},
		{Name: "meg_10ohm", Kind: MEG, RL: 10},
	}
}

func main() {
	write := flag.Bool("write", false, "regenerate netlists/ and ngspice/ next to this file")
	dir := flag.String("dir", "examples/bearden_meg", "example folder (for -write and results.csv)")
	flag.Parse()
	if *write {
		if err := writeNetlists(*dir); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		return
	}
	if err := run(*dir); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func writeNetlists(dir string) error {
	for _, sub := range []string{"netlists", "ngspice"} {
		if err := os.MkdirAll(filepath.Join(dir, sub), 0o755); err != nil {
			return err
		}
	}
	for _, s := range Twins() {
		if err := os.WriteFile(filepath.Join(dir, "netlists", s.Name+".cir"), []byte(s.Kron()), 0o644); err != nil {
			return err
		}
		if err := os.WriteFile(filepath.Join(dir, "ngspice", s.Name+".cir"), []byte(s.Ngspice(s.Name)), 0o644); err != nil {
			return err
		}
		fmt.Println("wrote", s.Name)
	}
	return nil
}

func run(dir string) error {
	rows := [][]string{{"sweep", "kind", "V", "f_Hz", "RL_ohm", "mismatch", "Cres_F", "P_drive_W", "P_magnets_W", "P_load_W", "COP", "ledger_rel_residual"}}
	best, failed := 0.0, 0
	fmt.Println("Bearden MEG and two-magnet knife edge, kronspice magnetic domain, full ledger")
	fmt.Println("COP = load power / drive power over whole cycles; claim > 1, standard prediction < 1")
	fmt.Println()
	for _, sw := range Sweeps() {
		fmt.Println(sw.Name)
		for _, s := range sw.Specs {
			r, _, err := s.Run()
			if err != nil {
				fmt.Printf("  %-58s FAILED: %v\n", s.title(), err)
				failed++
				continue
			}
			d := s.withDefaults()
			fmt.Printf("  %-58s drive %10.4g W  magnets %+9.2g W  load %10.4g W  COP %.4f  residual %.1e\n",
				s.title(), r.Pdrive, r.Pmagnet, r.Pload, r.COP, r.RelResidual(s))
			if r.COP > best {
				best = r.COP
			}
			g := func(x float64) string { return strconv.FormatFloat(x, 'g', 6, 64) }
			rows = append(rows, []string{sw.Name, d.Kind.String(), g(d.V), g(d.F), g(d.RL), g(d.Mismatch), g(d.Cres),
				g(r.Pdrive), g(r.Pmagnet), g(r.Pload), g(r.COP), g(r.RelResidual(s))})
		}
		fmt.Println()
	}
	fmt.Printf("highest COP over %d runs: %.4f (%d runs did not converge)\n", len(rows)-1, best, failed)
	f, err := os.Create(filepath.Join(dir, "results.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	if err := w.WriteAll(rows); err != nil {
		return err
	}
	fmt.Println("wrote", filepath.Join(dir, "results.csv"))
	return nil
}
