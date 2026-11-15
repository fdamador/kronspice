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
)

// Tags name the cases in the report's CSV columns, in Cases order.
var Tags = []string{"plain", "ece_t", "ece2_t", "ece_s", "ece2_s"}

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

// interp is y(t) by linear interpolation on the sorted series ts, ys.
func interp(ts, ys []float64, t float64) float64 {
	n := len(ts)
	if t <= ts[0] {
		return ys[0]
	}
	if t >= ts[n-1] {
		return ys[n-1]
	}
	lo, hi := 0, n-1
	for hi-lo > 1 {
		m := (lo + hi) / 2
		if ts[m] <= t {
			lo = m
		} else {
			hi = m
		}
	}
	a := (t - ts[lo]) / (ts[hi] - ts[lo])
	return ys[lo] + a*(ys[hi]-ys[lo])
}

// cumulative is the running trapezoidal integral of p over ts.
func cumulative(ts, p []float64) []float64 {
	e := make([]float64, len(ts))
	for k := 1; k < len(ts); k++ {
		e[k] = e[k-1] + 0.5*(p[k]+p[k-1])*(ts[k]-ts[k-1])
	}
	return e
}

// WriteData writes the report's data to dir:
//
//	ac.csv      — the AC sweep: |I|, phase [deg] and model-supplied power [W]
//	              per case;
//	tran.csv    — the 0.3 s transient on a 50 µs grid: source current and the
//	              running energies (source in, heat, model-supplied, stored) per
//	              case;
//	summary.csv — one row per case: the numbers of run_output.txt.
func WriteData(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	var acRows [][]float64
	var sum [][]float64
	const dt, stop = 50e-6, 0.3
	nt := int(math.Round(stop/dt)) + 1
	tr := make([][]float64, nt)
	for k := range tr {
		tr[k] = []float64{float64(k) * dt}
	}
	acHead, trHead := "f", "t"
	for n, k := range Cases {
		tag := Tags[n]
		acHead += fmt.Sprintf(",mag_%s,phase_%s,pmodel_%s", tag, tag, tag)
		trHead += fmt.Sprintf(",i_%s,ein_%s,eheat_%s,emodel_%s,stored_%s", tag, tag, tag, tag, tag)

		f, i, r, err := sweep(k.Netlist())
		if err != nil {
			return fmt.Errorf("%s: %w", k.Name, err)
		}
		if acRows == nil {
			acRows = make([][]float64, len(f))
			for j := range f {
				acRows[j] = []float64{f[j]}
			}
		}
		for j := range f {
			// V1's current runs n+ → n− through the source; the loop current is −I.
			acRows[j] = append(acRows[j], cmplx.Abs(i[j]), cmplx.Phase(-i[j])*180/math.Pi, r.Ledger.ModelSupplied[j])
		}

		ck, err := build(k.Netlist())
		if err != nil {
			return err
		}
		wf, err := analysis.Tran(ck, analysis.TranSpec{Step: 20e-6, Stop: stop, MaxStep: 20e-6}, analysis.TranOptions{})
		if err != nil {
			return fmt.Errorf("%s transient: %w", k.Name, err)
		}
		iv, ok := wf.ElementCurrent("V1")
		if !ok {
			return fmt.Errorf("%s: no transient current for V1", k.Name)
		}
		l := wf.Ledger
		ein, eheat, emod := cumulative(wf.Time, l.Inflow), cumulative(wf.Time, l.Dissipation), cumulative(wf.Time, l.ModelSupplied)
		for j := range tr {
			t := tr[j][0]
			tr[j] = append(tr[j], -interp(wf.Time, iv, t), interp(wf.Time, ein, t), interp(wf.Time, eheat, t),
				interp(wf.Time, emod, t), interp(wf.Time, l.StoredEnergy, t))
		}

		res, err := Run(k)
		if err != nil {
			return err
		}
		sum = append(sum, []float64{float64(n), k.S, res.FPeak, res.FWant, res.Q, res.QWant, res.IPeak, res.TwinDiff,
			res.PSource, res.PHeat, res.PModel, res.EIn, res.EOut, res.EHeat, res.EModel, res.EStored, res.Mismatch})
	}
	if err := writeCSV(filepath.Join(dir, "ac.csv"), acHead, acRows); err != nil {
		return err
	}
	if err := writeCSV(filepath.Join(dir, "tran.csv"), trHead, tr); err != nil {
		return err
	}
	return writeCSV(filepath.Join(dir, "summary.csv"),
		"case,s,f_peak,f_want,q,q_want,i_peak,twin_diff,p_source,p_heat,p_model,e_in,e_out,e_heat,e_model,e_stored,mismatch", sum)
}
