// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render_part56.py (Parts 5 and 6), as the run computed it.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
)

func writeViz56(dir string, p5a Part5a, p5b []FieldResult, p6 Part6) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	e := DefaultERing()
	T := 2 * math.Pi / e.Omega0()
	// 5a: stored energy and port powers for chosen drives, time in periods.
	for i, r := range p5a.Set.Runs {
		var tp, w, pin, pref []float64
		pk := 0.0
		for _, x := range r.PInc {
			pk = math.Max(pk, x)
		}
		for k := range r.Time {
			tp = append(tp, r.Time[k]/T)
			w = append(w, r.Stored[k]/r.EInc)
			pin = append(pin, r.PInc[k]/pk)
			pref = append(pref, r.PRef[k]/pk)
		}
		if err := writeCSV(dir, fmt.Sprintf("p5a_run%d.csv", i), "time_periods,stored_over_einc,pinc_over_peak,pref_over_peak", tp, w, pin, pref); err != nil {
			return err
		}
	}
	names := ""
	for i, r := range p5a.Set.Runs {
		names += fmt.Sprintf("%d,%s,%.6f,%.6f\n", i, r.Name, r.Absorbed, r.PeakStored/r.EInc)
	}
	if err := os.WriteFile(filepath.Join(dir, "p5a_runs.csv"), []byte(names), 0o644); err != nil {
		return err
	}
	for _, tr := range []TrainRun{p5a.Locked, p5a.Anti} {
		var c, in, out, heat, peak []float64
		for i, cy := range tr.Cycles[:len(tr.Cycles)-1] {
			c, in, out, heat, peak = append(c, float64(i+1)), append(in, cy.EInc), append(out, cy.EOut), append(heat, cy.Heat), append(peak, cy.Peak)
		}
		name := "p5a_train_locked.csv"
		if tr.Name != p5a.Locked.Name {
			name = "p5a_train_anti.csv"
		}
		if err := writeCSV(dir, name, "cycle,e_in,e_out,heat,peak", c, in, out, heat, peak); err != nil {
			return err
		}
	}
	// 5b: per grid, the implosion and burst histories; the plane of the finest.
	for _, r := range p5b {
		col := func(ss []FieldSample, f func(FieldSample) float64) []float64 {
			var o []float64
			for _, s := range ss {
				o = append(o, f(s))
			}
			return o
		}
		for _, run := range []struct {
			name string
			ss   []FieldSample
		}{{"implosion", r.Implosion}, {"burst", r.Burst}} {
			if err := writeCSV(dir, fmt.Sprintf("p5b_%s_c%g.csv", run.name, r.Cells), "time_periods,inside,leak,surface_net,exact_inc",
				col(run.ss, func(s FieldSample) float64 { return s.T }),
				col(run.ss, func(s FieldSample) float64 { return s.Inside / r.EInc }),
				col(run.ss, func(s FieldSample) float64 { return (s.Outside + s.Absorbed) / r.EInc }),
				col(run.ss, func(s FieldSample) float64 { return (s.SurfIn - s.SurfOut) / r.EInc }),
				col(run.ss, func(s FieldSample) float64 { return s.ExactInc / r.EInc })); err != nil {
				return err
			}
		}
	}
	if len(p5b) > 0 {
		r := p5b[len(p5b)-1]
		var x, z, ey, bx, bz, ring []float64
		for _, p := range r.Plane {
			x, z, ey, bx, bz = append(x, p.X), append(z, p.Z), append(ey, p.Ey), append(bx, p.Bx), append(bz, p.Bz)
			ring = append(ring, map[bool]float64{false: 0, true: 1}[p.InRing])
		}
		if err := writeCSV(dir, "p5b_plane.csv", "x_m,z_m,ey_peak,bx_quarter,bz_quarter,in_ring", x, z, ey, bx, bz, ring); err != nil {
			return err
		}
		var c, f0, q, ab []float64
		for _, r := range p5b {
			c, f0, q, ab = append(c, r.Cells), append(f0, r.F0), append(q, r.Q), append(ab, r.Peak/(r.Peak+r.LeakAtPeak))
		}
		if err := writeCSV(dir, "p5b_summary.csv", "cells,f0_hz,q,absorbed", c, f0, q, ab); err != nil {
			return err
		}
	}
	// 6: patterns of the three types (plain ring), beats, angles, three states, implosion.
	for _, x := range p6.Persist {
		if x.Ring != "plain" {
			continue
		}
		var j []float64
		for k := range x.Snap[0] {
			j = append(j, float64(k))
		}
		cols := [][]float64{j}
		hdr := "node"
		for s := range x.Snap {
			cols = append(cols, x.Snap[s])
			hdr += fmt.Sprintf(",u_t%d", s)
		}
		if err := writeCSV(dir, "p6_type_"+x.Type+".csv", hdr, cols...); err != nil {
			return err
		}
	}
	for i, b := range append(append(append([]Beat{}, p6.Beats...), p6.Thetas...), p6.Crossed) {
		var t []float64
		for _, s := range b.T {
			t = append(t, s*b.DfExact)
		}
		if err := writeCSV(dir, fmt.Sprintf("p6_beat%d.csv", i), "time_beats,p", t, b.P); err != nil {
			return err
		}
	}
	var d, dfe, dff, dfo []float64
	for _, b := range p6.Beats {
		d, dfe, dff, dfo = append(d, b.Delta), append(dfe, b.DfExact), append(dff, b.DfFit), append(dfo, b.DfFirst)
	}
	if err := writeCSV(dir, "p6_df.csv", "delta,df_exact,df_kronspice,df_first_order", d, dfe, dff, dfo); err != nil {
		return err
	}
	th := p6.Three
	var tt []float64
	for _, s := range th.T {
		tt = append(tt, s*th.Periods/th.T[len(th.T)-1])
	}
	if err := writeCSV(dir, "p6_three.csv", "time_periods,p1,p2,p3,pred1,pred2,pred3,full1,full2,full3", tt, th.P[0], th.P[1], th.P[2],
		th.Pred[0], th.Pred[1], th.Pred[2], th.PredFull[0], th.PredFull[1], th.PredFull[2]); err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(dir, "p6_three_meta.csv"), []byte(fmt.Sprintf("%d,%d,%d,%d\n", th.From, p6.Cluster.Assign[0], p6.Cluster.Assign[1], p6.Cluster.Assign[2])), 0o644); err != nil {
		return err
	}
	for i, ex := range p6.Explosions {
		if err := writeCSV(dir, fmt.Sprintf("p6_explode%d.csv", i), "time_periods,energy,p_conj", ex.T, ex.E, ex.PConj); err != nil {
			return err
		}
	}
	for i, rb := range append(append([]Rebuild{}, p6.Rebuilds...), p6.FRebuild, p6.FGrowth) {
		if err := writeCSV(dir, fmt.Sprintf("p6_rebuild%d.csv", i), "time_periods,energy,p_conj", rb.T, rb.E, rb.PConj); err != nil {
			return err
		}
	}
	return writeCSV(dir, "p6_fexplode.csv", "time_periods,energy,p_conj", p6.FExplode.T, p6.FExplode.E, p6.FExplode.PConj)
}
