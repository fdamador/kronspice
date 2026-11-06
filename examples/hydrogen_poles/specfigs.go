// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
)

// SpectrumFigures writes Figs. S1–S3.
func SpectrumFigures(dir string, sp Spectrum) error {
	// S1: the two 2p networks' responses near resonance
	f := NewFig(960, 460, "Fig. S1 — The fine structure as two kronspice resonances",
		"Im(1/V) at the drive of Kron's 2p network (h = 0.05 a0) with the spin–orbit term for j = ½ and j = 3⁄2; each crosses zero at its level")
	const h, rmax = 0.05, 100
	lad12 := NewLadderV(Potential(1, 0.5, true), h, rmax)
	lad32 := NewLadderV(Potential(1, 1.5, true), h, rmax)
	lad0 := NewLadderV(Potential(1, 0.5, false), h, rmax)
	c := sp.Levels[1].Fine // 2p without spin–orbit at h/2
	toMHz := MuRedH * HaHz / 1e6
	var es []float64
	for i := 0; i <= 120; i++ {
		es = append(es, c+(-12000+24000*float64(i)/120)/toMHz)
	}
	ys := make([][]float64, 3)
	errs := make([]error, 3)
	parallel(3, func(i int) {
		ys[i], errs[i] = []Ladder{lad12, lad32, lad0}[i].Response(es, 4)
	})
	for _, e := range errs {
		if e != nil {
			return e
		}
	}
	mx := 0.0
	for _, y := range ys {
		for _, v := range y {
			mx = math.Max(mx, math.Abs(v))
		}
	}
	xs := make([]float64, len(es))
	for i := range es {
		xs[i] = (es[i] - c) * toMHz / 1000
	}
	ax := f.Axes(90, 80, 820, 280, -12, 12, -1.05, 1.05, false, false)
	ax.Grid(Ticks(-12, 12, 3), Ticks(-1, 1, 0.5), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "frequency offset from the 2p level without spin–orbit (GHz)", "Im(1/V), scaled")
	ax.Line([]float64{-12, 12}, []float64{0, 0}, colInk2, 1)
	cols := []string{colBlue, colOrange, colMuted}
	for i, y := range ys {
		s := make([]float64, len(y))
		for k := range y {
			s[k] = y[k] / mx
		}
		ax.Line(xs, s, cols[i], 2)
	}
	for i, l := range []Level{sp.Levels[2], sp.Levels[3]} {
		x := (l.Fine - c) * toMHz / 1000
		ax.Dot(x, 0, 5, cols[i])
	}
	f.Legend(90, 420, []string{"2p½ network", "2p3⁄2 network", "no spin–orbit"}, cols)
	f.Text(90, 446, fmt.Sprintf("The two zero crossings are %.1f MHz apart at this grid; Richardson over h = 0.1 and 0.05 a0 gives %.1f MHz, measured 10 969.1 MHz.",
		(sp.Levels[3].Fine-sp.Levels[2].Fine)*toMHz, sp.FS/1e6), "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figS1_fine_structure.svg")); err != nil {
		return err
	}

	// S2: every line against measurement, across the spectrum
	f = NewFig(1000, 520, "Fig. S2 — Hydrogen's lines found by kronspice and the model, against measurement",
		"difference from the measured value in parts per million, placed by frequency; the bands show where each line is observed")
	ax = f.Axes(90, 90, 860, 300, 1e8, 1e16, -40, 80, true, false)
	bands := []struct {
		lo, hi float64
		name   string
	}{{1e8, 3e11, "radio / microwave"}, {3e11, 4e14, "infrared"}, {4e14, 7.9e14, "visible"}, {7.9e14, 1e16, "ultraviolet"}}
	for i, b := range bands {
		fill := colGrid
		if i%2 == 1 {
			fill = colSurface
		}
		f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" fill-opacity="0.6"/>`, ax.PX(b.lo), ax.Y0, ax.PX(b.hi)-ax.PX(b.lo), ax.H, fill))
		f.Text((ax.PX(b.lo)+ax.PX(b.hi))/2, ax.Y0+16, b.name, "middle", colInk2, 11, 400)
	}
	ax.Grid(Decades(1e8, 1e16, 1), Ticks(-40, 80, 20), func(v float64) string { return Pow10Label(v) + " Hz" },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "frequency (log scale)", "found − measured (ppm)")
	ax.Line([]float64{1e8, 1e16}, []float64{0, 0}, colInk2, 1)
	labels := []string{"Lyα", "Hα", "Paα", "Brα", "2p fine structure", "21 cm", "ESR (1 T)"}
	for i, l := range sp.Lines {
		ppm := 1e6 * (l.Hz/l.Measured - 1)
		col := colBlue
		if l.How != "kronspice resonance" {
			col = colOrange
		}
		ax.Dot(l.Hz, ppm, 6, col)
		dy := -12.0
		if i%2 == 1 {
			dy = 20
		}
		f.Text(ax.PX(l.Hz), ax.PY(ppm)+dy, labels[i], "middle", colInk, 12, 500)
	}
	f.Legend(90, 460, []string{"kronspice circuit resonance", "this example's model (Part 1)"}, []string{colBlue, colOrange})
	f.Text(90, 488, "All within 60 ppm. The remainders are physics the models leave out: the Dirac level shifts and the Lamb shift for the optical lines; Breit, recoil", "start", colInk2, 12, 400)
	f.Text(90, 506, "and proton-structure terms for 21 cm; the bound-state g-factor correction (−α²/3) for ESR. The ESR point is the frequency per tesla, plotted at 1 T.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figS2_lines.svg")); err != nil {
		return err
	}

	// S3: H₂, HD, D₂ vibration circuits
	f = NewFig(960, 460, "Fig. S3 — Molecular vibration as an LC resonance: H₂, HD and D₂",
		"masses as capacitors, the Heitler–London spring as an inductor (damped to Q = 200 for display); dashed: measured ω_e")
	var ws []float64
	for i := 0; i <= 600; i++ {
		ws = append(ws, 2000+3000*float64(i)/600)
	}
	masses := [3][2]float64{{MpMe, MpMe}, {MpMe, MdMe}, {MdMe, MdMe}}
	resp := make([][]float64, 3)
	parallel(3, func(i int) {
		resp[i], errs[i] = VibResponse(sp.K, masses[i][0], masses[i][1], ws)
	})
	for _, e := range errs[:3] {
		if e != nil {
			return e
		}
	}
	ax = f.Axes(90, 80, 820, 280, 2000, 5000, 0, 1.05, false, false)
	ax.Grid(Ticks(2000, 5000, 500), Ticks(0, 1, 0.25), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "wavenumber (cm⁻¹)", "|V| at the second mass, scaled")
	vcols := []string{colBlue, colAqua, colOrange}
	names := []string{"H₂", "HD", "D₂"}
	var leg []string
	for i := range resp {
		mx := 0.0
		for _, v := range resp[i] {
			mx = math.Max(mx, v)
		}
		s := make([]float64, len(resp[i]))
		for k := range s {
			s[k] = resp[i][k] / mx
		}
		ax.Line(ws, s, vcols[i], 2)
		ax.VLine(measWe[i], vcols[i], 1.2, true)
		leg = append(leg, fmt.Sprintf("%s: %.0f cm⁻¹ (measured %.0f)", names[i], sp.Vib[i], measWe[i]))
	}
	f.Legend(90, 420, leg, vcols)
	f.Text(90, 446, "Each peak sits 13.5 % below the measured line: the Heitler–London spring is too soft. The isotope ratios, set by the masses alone, agree to 0.07 %.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "figS3_vibration.svg"))
}
