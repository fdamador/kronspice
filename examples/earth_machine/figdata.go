// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Figure data: CSV files under viz/data, written by -data, read by
// viz/plot_layers.py. Every number comes from the runs above.

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"sort"

	"github.com/fdamador/kronspice/analysis"
)

func writeCSV(path, header string, rows [][]float64) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, header)
	for _, r := range rows {
		for i, v := range r {
			if i > 0 {
				fmt.Fprint(f, ",")
			}
			fmt.Fprintf(f, "%.10g", v)
		}
		fmt.Fprintln(f)
	}
	return nil
}

func writeFigData(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	p := func(name string) string { return filepath.Join(dir, name) }

	// Layer 2: profiles, field and charge per shell; switch-on.
	a := StandardAtmosphere()
	s0, H := a.Calibrate()
	two, _ := a.TwoScale(1500, 5000)
	var rows [][]float64
	for z := 0.0; z <= 60e3; z += 250 {
		s1 := s0 * math.Exp(z/H)
		rows = append(rows, []float64{z, s1, two(z), a.J / s1, a.J / two(z)})
	}
	if err := writeCSV(p("l2_profile.csv"), "z,sigma_exp,sigma_two,E_exp,E_two", rows); err != nil {
		return err
	}
	// Space-charge per shell junction (single exponential column).
	sh := a.Shells()
	rows = nil
	I := a.GlobalCurrent()
	cum := 0.0
	for k := 1; k < len(sh); k++ {
		q := sh[k-1].C*I*sh[k-1].R - sh[k].C*I*sh[k].R
		cum += q
		rows = append(rows, []float64{sh[k].Z0, q, cum})
	}
	if err := writeCSV(p("l2_charge.csv"), "z,q,cum", rows); err != nil {
		return err
	}
	t, vi, ex, err := a.SwitchOn(3000, 2)
	if err != nil {
		return err
	}
	rows = nil
	for k := range t {
		rows = append(rows, []float64{t[k], vi[k], ex[k]})
	}
	if err := writeCSV(p("l2_switchon.csv"), "t,vi_circuit,vi_exact", rows); err != nil {
		return err
	}

	// Layer 3: constant-height and two-height spectra.
	lossy := FitLossy(MeasuredSchumann[0], 5)
	sp, err := lossy.Sweep(2, 45, 861)
	if err != nil {
		return err
	}
	th, _ := fitTwoHeight(5)
	sp2, err := th.Spectrum(181, 3, 42, 0.1)
	if err != nil {
		return err
	}
	rows = nil
	for k := range sp.F {
		rows = append(rows, []float64{sp.F[k], sp.Pole[k], sp.Equator[k]})
	}
	if err := writeCSV(p("l3_spectra.csv"), "f,const_pole,const_equator", rows); err != nil {
		return err
	}
	if err := writeTwoHeightData(dir, sp2); err != nil {
		return err
	}
	ideal := IdealCavity()
	rows = nil
	for _, nb := range []int{31, 61, 91, 181, 361} {
		c := ideal
		c.Bands = nb
		m := c.Modes(5)
		r := []float64{float64(nb)}
		for n := 1; n <= 5; n++ {
			r = append(r, m[n-1]/LegendreFreq(n)-1)
		}
		rows = append(rows, r)
	}
	if err := writeCSV(p("l3_convergence.csv"), "bands,e1,e2,e3,e4,e5", rows); err != nil {
		return err
	}
	rows = nil
	for n := 1; n <= 5; n++ {
		fc := lossy.Modes(5)[n-1]
		w0 := 2 * math.Pi * fc
		fd := math.Sqrt(w0*w0-lossy.Alpha*lossy.Alpha) / (2 * math.Pi)
		f2, q2 := th.Mode(n)
		rows = append(rows, []float64{float64(n), LegendreFreq(n), fd, w0 / (2 * lossy.Alpha), f2, q2, MeasuredSchumann[n-1]})
	}
	if err := writeCSV(p("l3_modes.csv"), "n,legendre,const_f,const_q,two_f,two_q,measured", rows); err != nil {
		return err
	}

	// Layer 1 statistics: interval survival functions.
	r := Rikitake().StatRun(20000, 0.001, 50, 1e-2)
	iv := append([]float64(nil), r.Intervals...)
	sort.Float64s(iv)
	rows = nil
	for k, v := range iv {
		rows = append(rows, []float64{v, 1 - float64(k)/float64(len(iv))})
	}
	if err := writeCSV(p("l1_survival.csv"), "interval,survival", rows); err != nil {
		return err
	}
	rows = nil
	for _, mu := range []float64{0.5, 1, 2} {
		for _, aa := range []float64{0.5, 1, 2, 4} {
			s := withMuA(mu, aa).StatRun(20000, 0.002, 200, 1e-2)
			m, cv, mx := meanCV(s.Intervals)
			rows = append(rows, []float64{mu, aa, float64(len(s.Rev)), m, cv, mx / m, precursor(s)})
		}
	}
	if err := writeCSV(p("l1_mua.csv"), "mu,A,reversals,mean,cv,maxovermean,precursor", rows); err != nil {
		return err
	}
	_, per := Equilibrium(1)
	wn := 2 * math.Pi / per
	rows = nil
	for x := 0.2; x <= 3.01; x += 0.1 {
		q := Rikitake()
		q.GMod, q.Wd = 0.3, x*wn
		s := q.StatRun(10000, 0.002, 200, 1e-2)
		rows = append(rows, []float64{x, s.RMS, float64(len(s.Rev)) / 9800})
	}
	if err := writeCSV(p("l1_drive.csv"), "wd_over_wn,rms,rate", rows); err != nil {
		return err
	}
	rows = nil
	for _, k := range []float64{-0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.75} {
		q := Rikitake()
		q.K = k
		s := q.StatRun(20000, 0.002, 200, 1e-2)
		rows = append(rows, []float64{k, float64(len(s.Rev)), s.RMS})
	}
	if err := writeCSV(p("l1_feedback.csv"), "K,reversals,rms", rows); err != nil {
		return err
	}
	if err := writeStackedData(dir); err != nil {
		return err
	}

	// Links: the unified network and the dipole → V_I chain.
	t2, _ := fitTwoHeight(5)
	_, hm := t2.Heights(complex(8, 0))
	u := Unified{Bands: 61, Layers: 40, Top: 100e3, Sigma: two, Hm: real(hm), RLoss: imag(hm) / real(hm), FRef: 8}
	c, err := u.Circuit(I, 1)
	if err != nil {
		return err
	}
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 391, FStart: 3, FStop: 42}, analysis.Options{})
	if err != nil {
		return err
	}
	vp, _ := ac.NodeVoltage("t0")
	ve, _ := ac.NodeVoltage(fmt.Sprintf("t%d", u.Bands/2))
	rows = nil
	for k := range ac.Freq {
		rows = append(rows, []float64{ac.Freq[k], cmplx.Abs(vp[k]), cmplx.Abs(ve[k])})
	}
	if err := writeCSV(p("link_unified.csv"), "f,pole,equator", rows); err != nil {
		return err
	}
	tr := Rikitake().RunRK4(400, 0.05)
	var rms float64
	k0 := len(tr.T) / 4
	for _, x := range tr.I1[k0:] {
		rms += x * x
	}
	rms = math.Sqrt(rms / float64(len(tr.I1)-k0))
	rows = nil
	for k := k0; k < len(tr.T); k++ {
		m := math.Min(1, math.Abs(tr.I1[k])/rms)
		q := 1 + (2.0-1)*(1-m)
		rows = append(rows, []float64{tr.T[k], tr.I1[k], a.VI / math.Sqrt(q)})
	}
	return writeCSV(p("link_dipole.csv"), "t,I1,VI", rows)
}

func writeTwoHeightData(dir string, sp Spectrum) error {
	var rows [][]float64
	for k := range sp.F {
		rows = append(rows, []float64{sp.F[k], sp.Pole[k], sp.Equator[k]})
	}
	return writeCSV(filepath.Join(dir, "l3_two.csv"), "f,pole,equator", rows)
}

// writeStackedData: side-by-side and stacked (gap 0.53 r_o and 3 r_o)
// traces for fig10 (viz/plot_stacked.py).
func writeStackedData(dir string) error {
	side := Rikitake().RunRK4(400, 0.05)
	near := StackedRikitake(StackedGapHomopolar).RunRK4(400, 0.05)
	far := StackedRikitake(3).RunRK4(400, 0.05)
	var rows [][]float64
	for k := range side.T {
		rows = append(rows, []float64{side.T[k], side.I1[k], near.I1[k], far.I1[k],
			side.W1[k] - side.W2[k], near.W1[k] - near.W2[k], far.W1[k] - far.W2[k]})
	}
	return writeCSV(filepath.Join(dir, "stacked_traces.csv"), "t,I1_side,I1_near,I1_far,A_side,A_near,A_far", rows)
}
