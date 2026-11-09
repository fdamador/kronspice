// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// -viz dir: the numbers behind viz/render.py's figure.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

func writeViz(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	create := func(name string) (*os.File, error) { return os.Create(filepath.Join(dir, name)) }

	f, err := create("elements.csv")
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "element,burn_kJ_per_g,vaporize_kJ_per_g,return")
	for _, e := range Elements {
		fmt.Fprintf(f, "%s,%g,%g,%g\n", e.Name, e.OxidationHeat()/1e3, e.VaporizationCost()/1e3, e.Return())
	}
	f.Close()

	f, err = create("burntime.csv")
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "T_K,r_nm,t_s")
	for T := 1000.0; T <= 3500; T += 25 {
		for _, r := range []float64{5, 10, 25, 50, 100} {
			fmt.Fprintf(f, "%g,%g,%g\n", T, r, BurnTime(r*1e-9, T, 0.21))
		}
	}
	f.Close()

	// the welder's last three cycles
	w, err := TuneWeld(70, 22.5, 60, 120)
	if err != nil {
		return err
	}
	b := circuit.NewBuilder("AC stick welder and arc")
	b.VoltageWave("VW", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 70 * math.Sqrt2, 60}})
	b.Inductor("LW", "s", "a", w.L)
	b.NonlinearCC("ARC", "a", "0", arcCurve(22.5))
	ckt, err := b.BuildTransient()
	if err != nil {
		return err
	}
	T := 1 / 60.0
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: T / 400, Stop: 20 * T, MaxStep: T / 400}, analysis.TranOptions{})
	if err != nil {
		return err
	}
	ia, _ := wf.ElementCurrent("Vs_ARC")
	va, _ := wf.NodeVoltage("a")
	vs, _ := wf.NodeVoltage("s")
	f, err = create("weld.csv")
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "t_ms,v_source,v_arc,i_arc")
	for k, t := range wf.Time {
		if t >= 17*T {
			fmt.Fprintf(f, "%g,%g,%g,%g\n", (t-17*T)*1e3, vs[k], va[k], ia[k])
		}
	}
	f.Close()

	// line power against the ball's resistance
	f, err = create("line.csv")
	if err != nil {
		return err
	}
	fmt.Fprintln(f, "a_m,R_ohm,P_W")
	for _, a := range []float64{0.15, 2.5} {
		l := Line{Vphase: 35e3 / math.Sqrt(3), F: 50, A: a, D: 20}
		l.Cm = 4 * math.Pi * eps0 * a * a / 20
		l.Cb = 4 * math.Pi * eps0 * a
		for e := 0.0; e <= 12; e += 0.1 {
			p, err := LinePower(l, math.Pow(10, e))
			if err != nil {
				return err
			}
			fmt.Fprintf(f, "%g,%g,%g\n", a, math.Pow(10, e), p)
		}
	}
	return f.Close()
}
