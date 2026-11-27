// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command switch_relaxation runs two relaxation oscillators built from one
// parameterized subcircuit (.param, .subckt with params:, a hysteretic S
// switch) and compares each period with the exact piecewise-exponential
// value.
//
//	go run ./examples/switch_relaxation
package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// Cell is one oscillator: R charges C from vcc until the switch turns on at
// hi = VT + VH; RON then discharges it to lo = VT − VH.
type Cell struct {
	Node             string
	R, C, Ron, Vcc   float64
	Lo, Hi           float64
	Exact, Simulated float64 // period [s]
	Cycles           int
}

// ExactPeriod is the charge time RC·ln((vcc − lo)/(vcc − hi)) plus the
// discharge time through R ∥ RON towards vcc·RON/(R + RON).
func (c Cell) ExactPeriod() float64 {
	charge := c.R * c.C * math.Log((c.Vcc-c.Lo)/(c.Vcc-c.Hi))
	vinf := c.Vcc * c.Ron / (c.R + c.Ron)
	tau := c.C * c.R * c.Ron / (c.R + c.Ron)
	return charge + tau*math.Log((c.Hi-vinf)/(c.Lo-vinf))
}

// Result holds both cells and the run's ledger mismatch.
type Result struct {
	Cells    []Cell
	Mismatch float64 // integration mismatch relative to the energy in
}

// Run simulates the netlist in dir.
func Run(dir string) (*Result, error) {
	d, err := netlist.ParseFile(filepath.Join(dir, "relax.cir"))
	if err != nil {
		return nil, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil, err
	}
	a := d.Analyses[0]
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: a.Step, Stop: a.Stop, UIC: a.UIC}, analysis.TranOptions{})
	if err != nil {
		return nil, err
	}
	p := d.Params
	cells := []Cell{
		{Node: "a", R: p["rc"], C: p["cc"], Ron: 1e3, Vcc: p["vcc"], Lo: 2, Hi: 3},
		{Node: "b", R: 2 * p["rc"], C: p["cc"], Ron: 1e3, Vcc: p["vcc"], Lo: 1.5, Hi: 3.5},
	}
	for i := range cells {
		v, ok := wf.NodeVoltage(cells[i].Node)
		if !ok {
			return nil, fmt.Errorf("no node %s", cells[i].Node)
		}
		// Turn-on instants are the local maxima of the capacitor voltage.
		var on []float64
		for k := 1; k+1 < len(v); k++ {
			if v[k] > v[k-1] && v[k] >= v[k+1] && v[k] > cells[i].Hi-0.1 {
				on = append(on, wf.Time[k])
			}
		}
		cells[i].Exact = cells[i].ExactPeriod()
		if len(on) >= 2 {
			cells[i].Cycles = len(on) - 1
			cells[i].Simulated = (on[len(on)-1] - on[0]) / float64(len(on)-1)
		}
	}
	r := &Result{Cells: cells}
	if l := wf.Ledger; l.EnergyIn > 0 {
		r.Mismatch = math.Abs(l.IntegrationMismatch) / l.EnergyIn
	}
	return r, nil
}

func main() {
	_, file, _, _ := runtime.Caller(0)
	r, err := Run(filepath.Dir(file))
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	fmt.Println("Relaxation oscillators from one subcircuit (S switch, .param, .subckt params:)")
	fmt.Println("cell  R       window      exact period  simulated     error     cycles")
	for _, c := range r.Cells {
		fmt.Printf("%-4s  %-6.0f  %.1f–%.1f V   %.5f ms    %.5f ms   %+.2e  %d\n", c.Node, c.R, c.Lo, c.Hi,
			c.Exact*1e3, c.Simulated*1e3, c.Simulated/c.Exact-1, c.Cycles)
	}
	fmt.Printf("ledger integration mismatch %.2e of the energy in\n", r.Mismatch)
}
