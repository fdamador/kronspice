// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// normalCircuit is an ordinary two-winding transformer with the given
// inductances and coupling, in the same rig as the spring transformer: the
// comparison case for the turns ratio.
func normalCircuit(r Rig, lp, ls, k, rp, rs, rload float64, ac bool) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("normal transformer")
	src := circuit.Element{Name: "VSRC", Kind: circuit.VoltageSource}
	if ac {
		src.ACMag = 1
	} else {
		src.Wave = circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, r.Vpk, r.Freq}}
	}
	b.Add(src, "s", "0")
	b.Resistor("RSRC", "s", "a", r.Rsrc)
	b.Resistor("RWIN", "a", "a1", rp)
	b.Inductor("LIN", "a1", "0", lp)
	b.Resistor("RWOUT", "o", "o1", rs)
	b.Inductor("LOUT", "o1", "0", ls)
	b.Resistor("RLOAD", "o", "0", rload)
	b.Mutual("K0", "LIN", "LOUT", k)
	return b.BuildTransient()
}

func writeCSV(path, head string, rows [][]float64) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, head)
	for _, row := range rows {
		for i, v := range row {
			if i > 0 {
				fmt.Fprint(f, ",")
			}
			fmt.Fprintf(f, "%.9g", v)
		}
		fmt.Fprintln(f)
	}
	return nil
}

// WriteData writes the report's data to dir:
//
//	wave.csv   — one design run at 1 kHz: input and output voltage and
//	             current, volts per winding turn of each winding, dΦ/dt in
//	             the passes, and the normal 100:100 transformer beside it;
//	freq.csv   — AC sweep 1 Hz–1 MHz, loaded: |v_s/v_p|, phase, |Z_in|,
//	             efficiency, per core;
//	ratio.csv  — ratio against the passes enclosed and the core turns,
//	             spring and normal;
//	passes.csv — dΦ/dt per pass (relative) for the plain and design springs.
func WriteData(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	r := DefaultRig()
	s := DefaultSpring(Materials[2])
	pri, sec := designWindings(s.N)
	axis := Winding{Name: "axis", Azimuth: -math.Pi / 2, First: 0, Passes: s.N, Turns: 1}

	// ---- wave.csv ----
	ws := []Winding{pri, sec, axis}
	_, fs, err := s.Inductance(ws)
	if err != nil {
		return err
	}
	cr, err := s.RunCircuit(r, pri, sec)
	if err != nil {
		return err
	}
	c, _, L, err := s.build(r, pri, sec, axis, cr.Rload, false)
	if err != nil {
		return err
	}
	per := 1 / r.Freq
	w, err := analysis.Tran(c, analysis.TranSpec{Step: per / 400, Stop: 20 * per, MaxStep: per / 400}, analysis.TranOptions{})
	if err != nil {
		return err
	}
	vin, _ := w.NodeVoltage("a")
	vout, _ := w.NodeVoltage("o")
	vax, _ := w.NodeVoltage("x")
	iin, _ := w.ElementCurrent("LIN")
	iout, _ := w.ElementCurrent("LOUT")
	iax, _ := w.ElementCurrent("LAX")
	segIn, _ := s.segments(pri)
	segOut, _ := s.segments(sec)
	// pass flux of pass k at the primary's azimuth, from the three winding currents
	passFlux := func(i int, k int) float64 {
		j := segIn[k]
		return fs[0].Seg[j]*iin[i] + fs[1].Seg[j]*iout[i] + fs[2].Seg[j]*iax[i]
	}
	_ = segOut
	// normal 100:100 transformer: same primary, same coupling, ratio 1
	kps := L[0][1] / math.Sqrt(L[0][0]*L[1][1])
	nc, err := normalCircuit(r, L[0][0], L[0][0], kps, s.WindingResistance(pri, r.CuDia), s.WindingResistance(pri, r.CuDia), cr.Rload*25, false)
	if err != nil {
		return err
	}
	wn, err := analysis.Tran(nc, analysis.TranSpec{Step: per / 400, Stop: 20 * per, MaxStep: per / 400}, analysis.TranOptions{})
	if err != nil {
		return err
	}
	nin, _ := wn.NodeVoltage("a")
	nout, _ := wn.NodeVoltage("o")
	var rows [][]float64
	t0 := 18 * per
	for i := 1; i < len(w.Time); i++ {
		if w.Time[i] < t0 {
			continue
		}
		dt := w.Time[i] - w.Time[i-1]
		d0 := (passFlux(i, 0) - passFlux(i-1, 0)) / dt
		d5 := (passFlux(i, 5) - passFlux(i-1, 5)) / dt
		d9 := (passFlux(i, 9) - passFlux(i-1, 9)) / dt
		// normal transformer at the same instant (nearest sample)
		jn := nearest(wn.Time, w.Time[i])
		rows = append(rows, []float64{(w.Time[i] - t0) * 1e3, vin[i], vout[i], vax[i], iin[i], iout[i],
			vin[i] / float64(pri.Turns), vout[i] / float64(sec.Turns), d0, d5, d9, nin[jn], nout[jn], nin[jn] / float64(pri.Turns), nout[jn] / float64(pri.Turns)})
	}
	if err := writeCSV(filepath.Join(dir, "wave.csv"),
		"t_ms,v_in,v_out,v_axis,i_in,i_out,vpt_in,vpt_out,dphi_pass1,dphi_pass6,dphi_pass10,normal_v_in,normal_v_out,normal_vpt_in,normal_vpt_out", rows); err != nil {
		return err
	}

	// ---- freq.csv ----
	rows = nil
	spec := analysis.ACSpec{Variation: "dec", Points: 40, FStart: 1, FStop: 1e6}
	var cols [][][]float64
	for _, m := range Materials {
		sm := DefaultSpring(m)
		crm, err := sm.RunCircuit(r, pri, sec)
		if err != nil {
			return err
		}
		c, _, _, err := sm.build(r, pri, sec, axis, crm.Rload, true)
		if err != nil {
			return err
		}
		res, err := analysis.AC(c, spec, analysis.Options{})
		if err != nil {
			return err
		}
		va, _ := res.NodeVoltage("a")
		vo, _ := res.NodeVoltage("o")
		vs, _ := res.NodeVoltage("s")
		col := make([][]float64, len(res.Freq))
		for i := range res.Freq {
			g := vo[i] / va[i]
			iIn := (vs[i] - va[i]) / complex(r.Rsrc, 0)
			zin := va[i] / iIn
			pin := real(va[i] * cmplx.Conj(iIn))
			pout := real(vo[i]*cmplx.Conj(vo[i])) / crm.Rload
			col[i] = []float64{res.Freq[i], cmplx.Abs(g), cmplx.Phase(g) * 180 / math.Pi, cmplx.Abs(zin), pout / pin}
		}
		cols = append(cols, col)
	}
	for i := range cols[0] {
		row := []float64{cols[0][i][0]}
		for _, col := range cols {
			row = append(row, col[i][1:]...)
		}
		rows = append(rows, row)
	}
	head := "f"
	for _, tag := range []string{"galfenol", "steel", "nano"} {
		head += fmt.Sprintf(",ratio_%s,phase_%s,zin_%s,eta_%s", tag, tag, tag, tag)
	}
	if err := writeCSV(filepath.Join(dir, "freq.csv"), head, rows); err != nil {
		return err
	}

	// ---- ratio.csv ----
	rows = nil
	for ns := 1; ns <= 10; ns++ {
		q := sec
		q.First, q.Passes = 0, ns
		rr, err := s.Ratio(pri, q)
		if err != nil {
			return err
		}
		rows = append(rows, []float64{0, float64(ns), rr, IdealRatio(pri, q), float64(q.Turns) / float64(pri.Turns)})
	}
	for _, n := range []int{1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30, 40} {
		sn := DefaultSpring(Materials[2])
		sn.N = n
		p, q := designWindings(n)
		if n < 2 {
			q.First, q.Passes = 0, 1
		}
		rr, err := sn.Ratio(p, q)
		if err != nil {
			return err
		}
		rows = append(rows, []float64{1, float64(n), rr, IdealRatio(p, q), float64(q.Turns) / float64(p.Turns)})
	}
	if err := writeCSV(filepath.Join(dir, "ratio.csv"), "sweep,x,ratio,ideal,normal", rows); err != nil {
		return err
	}

	// ---- passes.csv ----
	rows = nil
	plain := DefaultSpring(Materials[2])
	plain.A, plain.P, plain.RetArea = 1e-3, 3e-3, 1
	steel := DefaultSpring(Materials[1])
	gal := DefaultSpring(Materials[0])
	var pf [][]float64
	for _, sp := range []Spring{gal, steel, plain, s} {
		f, err := sp.unitField(pri)
		if err != nil {
			return err
		}
		v := sp.PassFlux(f)
		mx := 0.0
		for _, x := range v {
			mx = math.Max(mx, x)
		}
		for i := range v {
			v[i] /= mx
		}
		pf = append(pf, v)
	}
	for k := 0; k < s.N; k++ {
		rows = append(rows, []float64{float64(k + 1), pf[0][k], pf[1][k], pf[2][k], pf[3][k]})
	}
	if err := writeCSV(filepath.Join(dir, "passes.csv"), "pass,galfenol_design,steel_design,nano_plain,nano_design", rows); err != nil {
		return err
	}
	return writeLedger(dir)
}

func nearest(ts []float64, t float64) int {
	lo, hi := 0, len(ts)-1
	for hi-lo > 1 {
		m := (lo + hi) / 2
		if ts[m] < t {
			lo = m
		} else {
			hi = m
		}
	}
	if math.Abs(ts[lo]-t) < math.Abs(ts[hi]-t) {
		return lo
	}
	return hi
}

// LedgerRow is one core-turn count's energy ledger over a 20-cycle run.
type LedgerRow struct {
	N                             int
	Ratio                         float64 // AC open-circuit v_s/v_p
	DPhi                          float64 // peak dΦ/dt per pass [Wb/s]
	Bpk                           float64 // peak core flux density [T]
	EIn, EBack, ELoad, ECu, ESrcR float64 // gross in, returned to the source, load, copper, source R [J]
	EStore, Mismatch              float64 // [J]
	Pin, Pout                     float64 // steady state, last 5 cycles [W]
}

// LedgerSweep runs the nanocrystalline spring with N core turns, the
// primary round all N passes and the secondary round 2 (100 turns each),
// at the design drive (20 V peak, 1 kHz) into one fixed load (the N = 10
// design's 632 Ω), and books the whole run's ledger.
func LedgerSweep(ns []int) ([]LedgerRow, error) {
	r := DefaultRig()
	s0 := DefaultSpring(Materials[2])
	p0, q0 := designWindings(s0.N)
	c0, err := s0.RunCircuit(r, p0, q0)
	if err != nil {
		return nil, err
	}
	r.Load = c0.Rload
	var out []LedgerRow
	for _, n := range ns {
		s := DefaultSpring(Materials[2])
		s.N = n
		p, q := designWindings(n)
		cr, err := s.RunCircuit(r, p, q)
		if err != nil {
			return nil, fmt.Errorf("N = %d: %w", n, err)
		}
		out = append(out, LedgerRow{N: n, Ratio: cr.RatioAC,
			DPhi: math.Sqrt2 * cr.VinRMS / float64(p.Turns*p.Passes), Bpk: cr.Bpk,
			EIn: cr.EIn, EBack: cr.EOut, ELoad: cr.ELoad, ECu: cr.EDis - cr.ELoad - cr.ESrcR, ESrcR: cr.ESrcR,
			EStore: cr.EStore, Mismatch: cr.Mismatch, Pin: cr.Pin, Pout: cr.Pout})
	}
	return out, nil
}

// writeLedger writes ledger.csv for the report.
func writeLedger(dir string) error {
	rows, err := LedgerSweep([]int{2, 4, 6, 10, 16, 20, 30, 40})
	if err != nil {
		return err
	}
	var out [][]float64
	for _, l := range rows {
		out = append(out, []float64{float64(l.N), l.Ratio, l.DPhi, l.Bpk, l.EIn, l.EBack, l.EIn - l.EBack, l.ELoad, l.ECu + l.ESrcR, l.EStore, l.Mismatch, l.Pin, l.Pout})
	}
	return writeCSV(filepath.Join(dir, "ledger.csv"), "N,ratio,dphi_pk,b_pk,e_in,e_back,e_net,e_load,e_copper,e_stored,residual,p_in,p_out", out)
}
