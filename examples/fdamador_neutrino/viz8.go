// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render_part8.py (triplet.png) and viz/plot_triplet.py
// (structure_m+1.png, structure_m0.png, structure_m-1.png), as the run
// computed it.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// stateTag names the files of the triplet states +1, 0, −1.
var stateTag = [3]string{"p1", "0", "m1"}

func writeViz8(dir string, p Part8) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	var meta strings.Builder
	kv := func(k string, v float64) { fmt.Fprintf(&meta, "%s,%.9g\n", k, v) }
	meta.WriteString("key,value\n")
	ti := p.Indep
	for a := 0; a < 3; a++ {
		kv(fmt.Sprintf("gram_triplet_%d", a), ti.GramTriplet[a])
		kv(fmt.Sprintf("gram_types_%d", a), math.Abs(ti.GramTypes[a]))
		kv(fmt.Sprintf("omega_line_%d", a), ti.Omega[a])
	}
	for _, x := range []struct {
		name string
		c    TripletCl
	}{{"recip", p.Recip.Best}, {"cap", p.RecipCap.Best}, {"rot", p.Rot.Best}} {
		c := x.c
		kv(x.name+"_th12", c.Th12)
		kv(x.name+"_th23", c.Th23)
		kv(x.name+"_th13", c.Th13)
		kv(x.name+"_ratio2", c.Ratio2)
		kv(x.name+"_dcp", c.DeltaCP)
		kv(x.name+"_leak", c.Leak)
		kv(x.name+"_cost", c.Cost2)
		kv(x.name+"_nu", c.Knobs.Nu)
	}
	m := NuFIT6()
	kv("nufit_th12", m.Th12)
	kv("nufit_s12", m.S12)
	kv("nufit_th23", m.Th23)
	kv("nufit_s23", m.S23)
	kv("nufit_th13", m.Th13)
	kv("nufit_s13", m.S13)
	kv("nufit_ratio", m.Ratio)
	kv("nufit_dcp", m.DeltaCP)
	kv("torque_dj", p.Torque.DeltaJ)
	kv("torque_int", p.Torque.TorqueAt)
	kv("torque_df", p.Torque.Df)
	kv("om_ratio", p.OmRatio)
	kv("w0", p.OmRatioW[0])
	kv("w1", p.OmRatioW[1])
	kv("burst", 1881.0)
	for i, rm := range p.RotMasses {
		kv(fmt.Sprintf("rm%d_m", i), float64(rm.M))
		kv(fmt.Sprintf("rm%d_a", i), rm.A)
		kv(fmt.Sprintf("rm%d_omega", i), rm.Omega)
		kv(fmt.Sprintf("rm%d_rho2", i), rm.Rho2)
	}
	if err := os.WriteFile(filepath.Join(dir, "p8_meta.csv"), []byte(meta.String()), 0o644); err != nil {
		return err
	}
	// 8b kronspice run (time in periods of the mean cluster ω)
	r := p.RecipRun
	c := p.Recip.Best
	wref := (c.Omega[0] + c.Omega[1] + c.Omega[2]) / 3
	tp := make([]float64, len(r.T))
	for i, t := range r.T {
		tp[i] = t * wref / (2 * math.Pi)
	}
	if err := writeCSV(dir, "p8_triplet_run.csv", "t_periods,pe,pmu,ptau,pe_full,pmu_full,ptau_full,pe_3,pmu_3,ptau_3",
		tp, r.P[0], r.P[1], r.P[2], r.PFull[0], r.PFull[1], r.PFull[2], r.P3[0], r.P3[1], r.P3[2]); err != nil {
		return err
	}
	// 8c torque
	tq := p.Torque
	tb := make([]float64, len(tq.T))
	for i, t := range tq.T {
		tb[i] = t * tq.Df
	}
	if err := writeCSV(dir, "p8_torque.csv", "t_beats,j,torque_int,u", tb, tq.Jw, tq.TorqueInt, tq.U); err != nil {
		return err
	}
	// field (every grid's summary; the finest grid's fields for the figures)
	if len(p.Field) == 0 {
		return nil
	}
	var fm strings.Builder
	fm.WriteString("cells,state,m,f,jzA,jzM,ub_min,ub_max,mz_uB,mperp_uB,p_uBc,rho2_a2\n")
	e := DefaultERing()
	for _, res := range p.Field {
		for q, s := range res.States {
			w := 2 * math.Pi * s.F
			ax, pp := s.MomentPerQuantum()
			fmt.Fprintf(&fm, "%g,%s,%d,%.9g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", res.Cells, stateTag[q], s.M, s.F, s.JzAMean*w, s.JzMMean*w,
				s.UBMin, s.UBMax, ax/MuB, pp/MuB, s.DipolePerQuantum()*C0/MuB, s.Rho2/(e.A*e.A))
		}
	}
	if err := os.WriteFile(filepath.Join(dir, "p8_field_meta.csv"), []byte(fm.String()), 0o644); err != nil {
		return err
	}
	res := p.Field[len(p.Field)-1]
	for q, s := range res.States {
		w := 2 * math.Pi * s.F
		var th, jzA, jzM, ue, ub, mx, my, mz, px, py, pz []float64
		for k := range s.Theta {
			th = append(th, s.Theta[k])
			jzA, jzM = append(jzA, s.JzA[k]*w), append(jzM, s.JzM[k]*w)
			ue, ub = append(ue, s.UE[k]), append(ub, s.UB[k])
			mx, my, mz = append(mx, s.Mt[k][0]), append(my, s.Mt[k][1]), append(mz, s.Mt[k][2])
			px, py, pz = append(px, s.Pt[k][0]), append(py, s.Pt[k][1]), append(pz, s.Pt[k][2])
		}
		if err := writeCSV(dir, "p8_cycle_"+stateTag[q]+".csv", "theta,ue,ub,jzA_w,jzM_w,mx,my,mz,px,py,pz", th, ue, ub, jzA, jzM, mx, my, mz, px, py, pz); err != nil {
			return err
		}
		// planes at θ = π/4 (t = T/8): z = 0 and y = 0
		const th8 = math.Pi / 4
		L := 1.8 * e.A
		np := int(math.Round(2*L/res.Dx)) + 1
		var xy, xz [][3]float64
		for i := 0; i < np; i++ {
			for j := 0; j < np; j++ {
				u, v := -L+float64(i)*res.Dx, -L+float64(j)*res.Dx
				xy = append(xy, [3]float64{u, v, 0})
				xz = append(xz, [3]float64{u, 0, v})
			}
		}
		for _, pl := range []struct {
			name string
			pts  [][3]float64
		}{{"plane", xy}, {"mer", xz}} {
			es, bs, err := res.Sample(q, th8, pl.pts)
			if err != nil {
				return err
			}
			if err := writePts(dir, "p8_"+pl.name+"_"+stateTag[q]+".csv", pl.pts, es, bs); err != nil {
				return err
			}
		}
		// 3-D arrows: E on the ring's centre circle, B on circles beside the tube and on the axis
		var pts [][3]float64
		for k := 0; k < 16; k++ {
			ph := 2 * math.Pi * float64(k) / 16
			pts = append(pts, [3]float64{e.A * math.Cos(ph), e.A * math.Sin(ph), 0})
		}
		for _, rz := range [][2]float64{{e.A - 1.7*e.R, 0}, {e.A + 1.7*e.R, 0}, {e.A, 1.7 * e.R}, {e.A, -1.7 * e.R}} {
			for k := 0; k < 12; k++ {
				ph := 2 * math.Pi * (float64(k) + 0.5) / 12
				pts = append(pts, [3]float64{rz[0] * math.Cos(ph), rz[0] * math.Sin(ph), rz[1]})
			}
		}
		for k := -2; k <= 2; k++ {
			pts = append(pts, [3]float64{0, 0, float64(k) * 0.5 * e.A})
		}
		es, bs, err := res.Sample(q, th8, pts)
		if err != nil {
			return err
		}
		if err := writePts(dir, "p8_pts_"+stateTag[q]+".csv", pts, es, bs); err != nil {
			return err
		}
	}
	return nil
}

func writePts(dir, name string, pts, es, bs [][3]float64) error {
	n := len(pts)
	cols := make([][]float64, 9)
	for i := range cols {
		cols[i] = make([]float64, n)
	}
	for k := 0; k < n; k++ {
		for c := 0; c < 3; c++ {
			cols[c][k], cols[3+c][k], cols[6+c][k] = pts[k][c], es[k][c], bs[k][c]
		}
	}
	return writeCSV(dir, name, "x,y,z,ex,ey,ez,bx,by,bz", cols...)
}
