// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"encoding/csv"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"
)

// writeCSV writes columns (equal length) with a header.
func writeCSV(path string, head []string, cols ...[]float64) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	if err := w.Write(head); err != nil {
		return err
	}
	for j := range cols[0] {
		row := make([]string, len(cols))
		for c := range cols {
			row[c] = strconv.FormatFloat(cols[c][j], 'g', 9, 64)
		}
		if err := w.Write(row); err != nil {
			return err
		}
	}
	w.Flush()
	return w.Error()
}

// writeBenchData writes the bench waveforms and spectra for the figures.
func writeBenchData(dir string, cs []BenchCase) error {
	var head []string
	cols := [][]float64{cs[0].Run.T}
	head = append(head, "t")
	for k, c := range cs {
		head = append(head, "case"+strconv.Itoa(k))
		cols = append(cols, c.Run.VLoad)
	}
	if err := writeCSV(filepath.Join(dir, "bench_wave.csv"), head, cols...); err != nil {
		return err
	}
	b := cs[0].B
	n := 4 * b.MP
	f := make([]float64, n)
	sp := make([][]float64, len(cs))
	for k := 1; k <= n; k++ {
		f[k-1] = float64(k) * b.F0
	}
	for c := range cs {
		sp[c] = make([]float64, n)
		for k := 1; k <= n; k++ {
			sp[c][k-1] = 2 * cmplx.Abs(cs[c].Run.Load[k])
		}
	}
	hs := []string{"f"}
	for k := range cs {
		hs = append(hs, "case"+strconv.Itoa(k))
	}
	return writeCSV(filepath.Join(dir, "bench_spectrum.csv"), hs, append([][]float64{f}, sp...)...)
}

// writeTeslaData writes each winding's waveforms (sparks + signal, the
// signal alone, the mixing part) and the mixing spectrum near fs.
func writeTeslaData(dir string, rs []TCResult) error {
	for _, r := range rs {
		tag := "ordinary"
		if r.Coil.Bifilar {
			tag = "bifilar"
		}
		if err := writeCSV(filepath.Join(dir, "tesla_"+tag+"_wave.csv"), []string{"t", "both", "sparks", "signal", "mix"},
			r.Both.T, r.Both.Top, r.Sparks.Top, r.Signal.Top, r.Mix); err != nil {
			return err
		}
		f0 := 1 / r.Coil.Window
		ks := int(math.Round(r.Coil.Fs / f0))
		var f, mix, sig []float64
		for k := ks - 120; k <= ks+120; k++ {
			f = append(f, float64(k)*f0)
			mix = append(mix, 2*cmplx.Abs(dftReal(r.Mix, k)))
			sig = append(sig, 2*cmplx.Abs(dftReal(r.Signal.Top, k)))
		}
		if err := writeCSV(filepath.Join(dir, "tesla_"+tag+"_spectrum.csv"), []string{"f", "mix", "signal"}, f, mix, sig); err != nil {
			return err
		}
	}
	return nil
}
