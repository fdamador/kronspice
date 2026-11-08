// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render_part9.py (pair.png), viz/plot_pair_structure.py
// (structure_pair.png) and viz/plot_pair_cycle.py (pair_cycle.gif), as the
// run computed it.

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
)

func writeViz9(dir string, p Part9) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	write := func(name string, b *strings.Builder) error {
		return os.WriteFile(filepath.Join(dir, name), []byte(b.String()), 0o644)
	}
	e := DefaultERing()
	var b strings.Builder
	b.WriteString("geom,d_over_a,kd,k,rm_r11,rm_dipole,f_in,q_in,qrad_in,f_anti,q_anti,qrad_anti,q_in_kron,q_anti_kron\n")
	for _, r := range p.Pairs {
		fmt.Fprintf(&b, "%s,%.6g,%.6g,%.8g,%.8g,%.8g,%.9g,%.6g,%.6g,%.9g,%.6g,%.6g,%.6g,%.6g\n", map[Geom]string{Coaxial: "coaxial", Coplanar: "side"}[r.G],
			r.D/e.A, r.KD, r.ML, r.RmR, r.RmDip, r.In.F, r.In.Q, r.In.QRad, r.Anti.F, r.Anti.Q, r.Anti.QRad, r.QInK, r.QAntiK)
	}
	if err := write("p9_pairs.csv", &b); err != nil {
		return err
	}
	// a finer coaxial curve from the matrix alone (no kronspice), for the chart
	b.Reset()
	b.WriteString("geom,d_over_a,f_in,q_in,qrad_in,f_anti,q_anti,qrad_anti,k,f_single\n")
	for _, g := range []Geom{Coaxial, Coplanar} {
		lo := 0.5
		if g == Coplanar {
			lo = 2.5
		}
		for dd := lo; dd <= 12.0001; dd *= 1.06 {
			a := NewArray(e, g, dd*e.A, 2)
			in, an := a.PairModes()
			m, _ := a.Coupling(0, 1, e.Omega0())
			fmt.Fprintf(&b, "%s,%.6g,%.9g,%.6g,%.6g,%.9g,%.6g,%.6g,%.8g,%.9g\n", map[Geom]string{Coaxial: "coaxial", Coplanar: "side"}[g],
				dd, in.F, in.Q, in.QRad, an.F, an.Q, an.QRad, m/a.LSelf(e.Omega0()), p.Single.F)
		}
	}
	if err := write("p9_curve.csv", &b); err != nil {
		return err
	}
	for k, bt := range p.Beats {
		b.Reset()
		b.WriteString("t_periods,e0,e1,u,i0,i1,v0,v1\n")
		a := NewArray(e, Coaxial, bt.D, 2)
		r := bt.Run
		T := 1 / ((bt.FIn + bt.FAnti) / 2)
		for t := range r.Time {
			fmt.Fprintf(&b, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", r.Time[t]/T, r.ringEnergy(a, 0, t), r.ringEnergy(a, 1, t), r.U[t],
				r.I[0][t], r.I[1][t], r.VC[0][t], r.VC[1][t])
		}
		if err := write(fmt.Sprintf("p9_beat%d.csv", k), &b); err != nil {
			return err
		}
	}
	b.Reset()
	b.WriteString("geom,d_over_a,dmdd,f_in_a,f_in_k,f_anti_a,f_anti_k,f_rand_early,s_early,f_rand_late,s_late\n")
	for i, f := range p.Forces {
		fmt.Fprintf(&b, "%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", map[Geom]string{Coaxial: "coaxial", Coplanar: "side"}[p.ForceGeoms[i]],
			f.D/e.A, f.DMdd, f.FInA, f.FInK, f.FAntiA, f.FAntiK, f.FRandEarly, f.FRandEarlySigma, f.FRandLate, f.FRandLateSigma)
	}
	if err := write("p9_force.csv", &b); err != nil {
		return err
	}
	// a finer analytic force curve (coaxial and side by side)
	b.Reset()
	b.WriteString("geom,d_over_a,f_in,f_anti\n")
	for _, g := range []Geom{Coaxial, Coplanar} {
		lo := 0.6
		if g == Coplanar {
			lo = 2.5
		}
		for dd := lo; dd <= 8.0001; dd *= 1.06 {
			a := NewArray(e, g, dd*e.A, 2)
			in, an := a.PairModes()
			var fi, fa float64
			for _, x := range []struct {
				m  AMode
				sg float64
				o  *float64
			}{{in, 1, &fi}, {an, -1, &fa}} {
				m, _ := a.Coupling(0, 1, x.m.W)
				*x.o = x.sg * dMdd(e, g, dd*e.A, x.m.W) / (2 * (a.LSelf(x.m.W) + x.sg*m)) * e.A
			}
			fmt.Fprintf(&b, "%s,%.6g,%.6g,%.6g\n", map[Geom]string{Coaxial: "coaxial", Coplanar: "side"}[g], dd, fi, fa)
		}
	}
	if err := write("p9_forcecurve.csv", &b); err != nil {
		return err
	}
	b.Reset()
	b.WriteString("n,f_bright,q_bright,f_dark,q_dark,qrad_dark,q_dark_kron\n")
	for _, r := range p.Chain {
		fmt.Fprintf(&b, "%d,%.9g,%.6g,%.9g,%.6g,%.6g,%.6g\n", r.N, r.Bright.F, r.Bright.Q, r.Dark.F, r.Dark.Q, r.DarkRad.Q, r.QDarkK)
	}
	if err := write("p9_chain.csv", &b); err != nil {
		return err
	}
	// mode vectors of the N = 8 chain's darkest and brightest modes (for the chart)
	b.Reset()
	b.WriteString("n,j,dark,bright\n")
	for _, r := range p.Chain {
		if r.N != 8 && r.N != 16 {
			continue
		}
		for j := range r.Dark.V {
			fmt.Fprintf(&b, "%d,%d,%.6g,%.6g\n", r.N, j, real(r.Dark.V[j]), real(r.Bright.V[j]))
		}
	}
	if err := write("p9_chainvec.csv", &b); err != nil {
		return err
	}
	if len(p.Field) > 0 {
		b.Reset()
		b.WriteString("label,cells,d_over_a,nodes,f,q,q_energy,residual,seconds,circuit_f,circuit_q\n")
		for i, r := range p.Field {
			fmt.Fprintf(&b, "%s,%g,%.6g,%d,%.9g,%.6g,%.6g,%.3g,%.4g,%.9g,%.9g\n", p.FieldLabel[i], r.Cells, r.D/e.A, r.Nodes, r.F, r.Q, r.QEnergy,
				r.Residual, r.Seconds, p.FieldPredF[i], p.FieldPred[i])
		}
		if err := write("p9_field.csv", &b); err != nil {
			return err
		}
		b.Reset()
		b.WriteString("label,t_periods,amp\n")
		n := 0
		for i, r := range p.Field {
			for k := range r.AT {
				fmt.Fprintf(&b, "%s,%.5g,%.6g\n", p.FieldLabel[i]+fmt.Sprintf(" r/%g", r.Cells), r.AT[k], r.Amp[k])
				n++
			}
		}
		// rows read back with -pairfieldcsv carry no traces: keep the file the field run wrote
		if n > 0 {
			if err := write("p9_fieldamp.csv", &b); err != nil {
				return err
			}
		}
	}
	b.Reset()
	b.WriteString("key,value\n")
	kv := func(k string, v float64) { fmt.Fprintf(&b, "%s,%.9g\n", k, v) }
	kv("f_single", p.Single.F)
	kv("q_single", p.Single.Q)
	kv("qrad_single", p.Single.QRad)
	kv("a", e.A)
	kv("r", e.R)
	kv("chain_d_over_a", p.ChainD/e.A)
	kv("chain_exp", p.ChainExp)
	kv("chain_coef", p.ChainCoef)
	kv("n_needed", p.NNeeded)
	kv("q_needed", DefaultBench().Run().QNeeded)
	kv("sigma_weak", p.Nu.SigmaWeak)
	kv("sigma_model", p.Nu.SigmaModel)
	kv("shift_2a", p.Nu.Shift2a)
	kv("contact_v", p.Nu.ContactV)
	for k, bt := range p.Beats {
		kv(fmt.Sprintf("beat%d_d_over_a", k), bt.D/e.A)
		kv(fmt.Sprintf("beat%d_df", k), bt.DfMatrix)
		kv(fmt.Sprintf("beat%d_fmean", k), (bt.FIn+bt.FAnti)/2)
		kv(fmt.Sprintf("beat%d_qin", k), bt.QIn)
		kv(fmt.Sprintf("beat%d_qanti", k), bt.QAnti)
	}
	return write("p9_meta.csv", &b)
}
