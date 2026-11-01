// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command regen regenerates golden reference files by running ngspice in
// batch mode on every netlist in a directory.
//
//	go run ./internal/golden/cmd/regen -netlists testdata/netlists -out testdata/golden
//
// It requires ngspice on PATH. Each netlist must contain a .op line; regen
// appends a .control block that prints all operating-point values.
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
)

func main() {
	dir := flag.String("netlists", "testdata/netlists", "directory of .cir files")
	out := flag.String("out", "testdata/golden", "directory for golden .csv files")
	bin := flag.String("ngspice", "ngspice", "ngspice executable")
	flag.Parse()

	if _, err := exec.LookPath(*bin); err != nil {
		fmt.Fprintf(os.Stderr, "regen: %s not found on PATH: %v\n", *bin, err)
		os.Exit(1)
	}
	files, err := filepath.Glob(filepath.Join(*dir, "*.cir"))
	if err != nil || len(files) == 0 {
		fmt.Fprintf(os.Stderr, "regen: no netlists in %s\n", *dir)
		os.Exit(1)
	}
	failed := false
	for _, f := range files {
		vals, err := run(*bin, f)
		if err != nil {
			fmt.Fprintf(os.Stderr, "regen: %s: %v\n", f, err)
			failed = true
			continue
		}
		name := strings.TrimSuffix(filepath.Base(f), ".cir") + ".csv"
		if err := write(filepath.Join(*out, name), vals); err != nil {
			fmt.Fprintf(os.Stderr, "regen: %v\n", err)
			failed = true
			continue
		}
		fmt.Printf("%s: %d values\n", name, len(vals))
	}
	if failed {
		os.Exit(1)
	}
}

// line matches ngspice "print all" output, e.g. "v(out) = 5.000000e+00" or
// "v1#branch = -5.00000e-03".
var line = regexp.MustCompile(`^\s*([A-Za-z0-9_#().:+-]+)\s*=\s*([-+0-9.eE]+)\s*$`)

func run(bin, path string) (map[string]string, error) {
	src, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	// Insert a control block before .end.
	text := string(src)
	ctl := ".control\nop\nprint all\n.endc\n"
	if i := strings.LastIndex(strings.ToLower(text), ".end"); i >= 0 {
		text = text[:i] + ctl + text[i:]
	} else {
		text += ctl + ".end\n"
	}
	tmp, err := os.CreateTemp("", "regen-*.cir")
	if err != nil {
		return nil, err
	}
	defer os.Remove(tmp.Name())
	if _, err := tmp.WriteString(text); err != nil {
		return nil, err
	}
	tmp.Close()

	cmd := exec.Command(bin, "-b", tmp.Name())
	var stdout bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stdout
	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("ngspice failed: %v\n%s", err, stdout.String())
	}
	vals := map[string]string{}
	sc := bufio.NewScanner(&stdout)
	for sc.Scan() {
		m := line.FindStringSubmatch(sc.Text())
		if m == nil {
			continue
		}
		key := strings.ToLower(m[1])
		if strings.HasSuffix(key, "#branch") {
			key = "i(" + strings.TrimSuffix(key, "#branch") + ")"
		} else if !strings.HasPrefix(key, "v(") {
			key = "v(" + key + ")"
		}
		vals[key] = m[2]
	}
	if len(vals) == 0 {
		return nil, fmt.Errorf("no values parsed from ngspice output:\n%s", stdout.String())
	}
	return vals, nil
}

func write(path string, vals map[string]string) error {
	keys := make([]string, 0, len(vals))
	for k := range vals {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	var b strings.Builder
	b.WriteString("# reference: ngspice batch run (internal/golden/cmd/regen)\nname,value\n")
	for _, k := range keys {
		fmt.Fprintf(&b, "%s,%s\n", k, vals[k])
	}
	return os.WriteFile(path, []byte(b.String()), 0o644)
}
