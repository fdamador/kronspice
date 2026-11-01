// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package export writes the three files that go with every example
// circuit: the SPICE netlist (NAME.cir), a KiCad 8 schematic
// (NAME.kicad_sch, connectivity checked by reading it back) and the results
// sheet with the MNA and Kron torn solutions (NAME_results.md and .csv).
// The Go code that builds the circuit stays the source; these files are
// generated from it.
//
// Examples call it from an export test, run with
//
//	KRONSPICE_EXPORT=1 go test ./examples/... -run TestExport
//
// Without KRONSPICE_EXPORT the same test only checks, in memory, that the
// netlist parses back, the schematic reads back and the sheet computes.
package export

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kicad"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/report"
)

// Files is what Source produces, in memory.
type Files struct {
	Netlist, Schematic, Markdown, CSV []byte
	Sheet                             *report.Sheet // nil if the sheet failed
	SheetErr                          error
}

// Source builds the files for a netlist text. A results sheet that cannot
// be computed (for example an operating point that does not converge) is
// written as a sheet stating the error; the netlist and schematic are
// still produced.
func Source(name, src string, opt report.Options) (*Files, error) {
	if _, err := netlist.ParseString(src); err != nil {
		return nil, fmt.Errorf("export %s: netlist does not parse: %w", name, err)
	}
	f := &Files{Netlist: []byte(src)}
	var b bytes.Buffer
	if err := kicad.Write(&b, src, kicad.Options{Project: name}); err != nil {
		return nil, fmt.Errorf("export %s: schematic: %w", name, err)
	}
	if _, err := kicad.Read(b.String()); err != nil {
		return nil, fmt.Errorf("export %s: schematic check: %w", name, err)
	}
	f.Schematic = b.Bytes()
	if opt.Source == "" {
		opt.Source = name + ".cir"
	}
	sh, err := report.Compute(src, opt)
	var md, csv bytes.Buffer
	if err != nil {
		f.SheetErr = err
		fmt.Fprintf(&md, "# Results: %s\n\nNetlist `%s`. The results sheet could not be computed: %v\n", name, opt.Source, err)
	} else {
		f.Sheet = sh
		sh.Markdown(&md, opt.MaxRows)
		sh.CSV(&csv)
	}
	f.Markdown, f.CSV = md.Bytes(), csv.Bytes()
	return f, nil
}

// Circuit writes c as a netlist with the given control cards (".ac …",
// ".tran …"), then proceeds as Source.
func Circuit(name string, c *circuit.Circuit, opt report.Options, control ...string) (*Files, error) {
	var b bytes.Buffer
	if err := c.WriteNetlist(&b, control...); err != nil {
		return nil, fmt.Errorf("export %s: %w", name, err)
	}
	return Source(name, b.String(), opt)
}

// Write puts the files in dir.
func (f *Files) Write(dir, name string) ([]string, error) {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return nil, err
	}
	var out []string
	for _, x := range []struct {
		suffix string
		b      []byte
	}{{".cir", f.Netlist}, {".kicad_sch", f.Schematic}, {"_results.md", f.Markdown}, {"_results.csv", f.CSV}} {
		if x.b == nil {
			continue
		}
		p := filepath.Join(dir, name+x.suffix)
		if err := os.WriteFile(p, x.b, 0o644); err != nil {
			return out, err
		}
		out = append(out, p)
	}
	return out, nil
}

// Item is one circuit of an example.
type Item struct {
	Name    string
	Build   func() (*circuit.Circuit, error) // or Src
	Src     string
	Control []string // control cards for Build's netlist
	Opt     report.Options
	// Tol is the largest accepted |torn − MNA| / max|MNA| (default 1e-9,
	// CLAUDE.md §3a rule 4). Raise it only for an ill-conditioned circuit,
	// with the reason at the call site.
	Tol float64
}

// Test exports every item to dir ("export" inside the example) when
// KRONSPICE_EXPORT is set, and otherwise only checks them.
func Test(t *testing.T, dir string, items ...Item) {
	t.Helper()
	write := os.Getenv("KRONSPICE_EXPORT") != ""
	for _, it := range items {
		var f *Files
		var err error
		if it.Build != nil {
			c, berr := it.Build()
			if berr != nil {
				t.Errorf("%s: build: %v", it.Name, berr)
				continue
			}
			f, err = Circuit(it.Name, c, it.Opt, it.Control...)
		} else {
			f, err = Source(it.Name, it.Src, it.Opt)
		}
		if err != nil {
			t.Errorf("%v", err)
			continue
		}
		if f.SheetErr != nil {
			t.Logf("%s: results sheet: %v", it.Name, f.SheetErr)
		} else if f.Sheet.KronNote != "" {
			t.Logf("%s: %s", it.Name, f.Sheet.KronNote)
		} else if f.Sheet.MaxDiff > tolOf(it) {
			t.Errorf("%s: Kron torn solution differs from MNA by %.3g, tolerance %.3g (CLAUDE.md §3a rule 4)", it.Name, f.Sheet.MaxDiff, tolOf(it))
		}
		if write {
			ps, err := f.Write(dir, it.Name)
			if err != nil {
				t.Errorf("%s: %v", it.Name, err)
			}
			for _, p := range ps {
				t.Log(p)
			}
		}
	}
}

func tolOf(it Item) float64 {
	if it.Tol > 0 {
		return it.Tol
	}
	return 1e-9
}
