// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
	"math"
	"os"
	"strings"
)

func create(path string) *os.File {
	f, err := os.Create(path)
	if err != nil {
		log.Fatal(err)
	}
	return f
}

// writePaschen writes breakdown voltage against pressure × gap for every
// gas, with the 500 psi points for the Part 2 gaps marked.
func writePaschen(path string) {
	f := create(path)
	defer f.Close()
	fmt.Fprint(f, "pd_torr_cm")
	for _, g := range Gases {
		fmt.Fprintf(f, ",%s", strings.ToLower(strings.Fields(g.Name)[0]))
	}
	fmt.Fprintln(f)
	for e := -1.0; e <= 3.5; e += 0.025 {
		pd := math.Pow(10, e)
		fmt.Fprintf(f, "%.6g", pd)
		for _, g := range Gases {
			v := g.Breakdown(torr, pd/100)
			if math.IsInf(v, 0) {
				fmt.Fprint(f, ",")
				continue
			}
			fmt.Fprintf(f, ",%.6g", v)
		}
		fmt.Fprintln(f)
	}
}

// writeMarx writes the fired generator's output voltage against time.
func writeMarx(path string, m Marx) {
	w, err := run(m.FireNetlist(), 0)
	if err != nil {
		log.Fatal(err)
	}
	out, _ := w.NodeVoltage("out")
	f := create(path)
	defer f.Close()
	fmt.Fprintln(f, "t_s,v_out")
	for k := range w.Time {
		fmt.Fprintf(f, "%.6g,%.6g\n", w.Time[k], out[k])
	}
}

// writePump writes the pump's varactor and output voltages against time,
// thinned to at most about 20,000 rows (the full run has 400,000 points).
func writePump(path string, p Pump) {
	w, err := run(p.Netlist(), 0)
	if err != nil {
		log.Fatal(err)
	}
	v, _ := w.NodeVoltage("v")
	out, _ := w.NodeVoltage("out")
	f := create(path)
	defer f.Close()
	fmt.Fprintln(f, "t_s,v_varactor,v_out")
	step := len(w.Time)/20000 + 1
	for k := 0; k < len(w.Time); k += step {
		fmt.Fprintf(f, "%.6g,%.6g,%.6g\n", w.Time[k], v[k], out[k])
	}
}
