// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
	"strings"
)

// Figures writes the field-line and scan figures.
func Figures(dir string, u, f4 Bar) error {
	if err := figFieldLines(dir, u, f4); err != nil {
		return err
	}
	return figScan(dir, u, f4)
}

// contour returns the segments of the level set ψ = c on a grid (marching
// squares, linear interpolation along cell edges).
func contour(xs, ys []float64, v [][]float64, c float64) [][4]float64 {
	var seg [][4]float64
	interp := func(x0, y0, v0, x1, y1, v1 float64) (float64, float64) {
		t := (c - v0) / (v1 - v0)
		return x0 + t*(x1-x0), y0 + t*(y1-y0)
	}
	for i := 0; i+1 < len(xs); i++ {
		for j := 0; j+1 < len(ys); j++ {
			// corners: a (i, j), b (i+1, j), c2 (i+1, j+1), d (i, j+1)
			type pt struct{ x, y, v float64 }
			p := [4]pt{{xs[i], ys[j], v[i][j]}, {xs[i+1], ys[j], v[i+1][j]}, {xs[i+1], ys[j+1], v[i+1][j+1]}, {xs[i], ys[j+1], v[i][j+1]}}
			var hits [][2]float64
			for e := 0; e < 4; e++ {
				a, b := p[e], p[(e+1)%4]
				if (a.v-c)*(b.v-c) < 0 {
					x, y := interp(a.x, a.y, a.v, b.x, b.y, b.v)
					hits = append(hits, [2]float64{x, y})
				}
			}
			for k := 0; k+1 < len(hits); k += 2 {
				seg = append(seg, [4]float64{hits[k][0], hits[k][1], hits[k+1][0], hits[k+1][1]})
			}
		}
	}
	return seg
}

// figFieldLines: the field lines (contours of ψ = ρA_φ) of the uniform bar
// and of the four-segment bar, in a plane through the axis, axis drawn
// horizontally as in Russell's drawing.
func figFieldLines(dir string, u, f4 Bar) error {
	f := NewFig(1000, 930, "Fig. M1 — Field lines of a bar magnet: one pair against four pairs of poles",
		"contours of the flux function ψ = ρA_φ in a plane through the axis (axis horizontal, as in Russell's drawing); the bar is the grey outline")
	panel := func(y0 float64, b Bar, title string) {
		zext, rext := 1.0*b.L, 0.42*b.L
		ax := f.Axes(70, y0, 860, 860*rext*2/(2*zext), -zext, zext, -rext, rext, false, false)
		sheet := b.Sheet(400)
		nz, nr := 220, 54
		zs := make([]float64, nz+1)
		rs := make([]float64, nr+1)
		for i := range zs {
			zs[i] = -zext + 2*zext*float64(i)/float64(nz)
		}
		for j := range rs {
			rs[j] = rext * float64(j) / float64(nr)
		}
		psi := make([][]float64, len(zs))
		peak := 0.0
		for i, z := range zs {
			psi[i] = make([]float64, len(rs))
			for j, r := range rs {
				_, _, p := b.Field(r, z, sheet)
				psi[i][j] = p
				peak = math.Max(peak, math.Abs(p))
			}
		}
		f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="%s" stroke-width="1"/>`, ax.X0, ax.Y0, ax.W, ax.H, colGrid))
		var sb strings.Builder
		for _, frac := range []float64{0.004, 0.012, 0.03, 0.06, 0.1, 0.16, 0.24, 0.34, 0.46, 0.6, 0.75, 0.9} {
			for _, sg := range []float64{1, -1} {
				for _, s := range contour(zs, rs, psi, sg*frac*peak) {
					for _, m := range []float64{1, -1} {
						fmt.Fprintf(&sb, "M%.1f %.1fL%.1f %.1f", ax.PX(s[0]), ax.PY(m*s[1]), ax.PX(s[2]), ax.PY(m*s[3]))
					}
				}
			}
		}
		f.Raw(fmt.Sprintf(`<path d="%s" fill="none" stroke="%s" stroke-width="1.2" stroke-linecap="round"/>`, sb.String(), colBlue))
		// the bar and its segments
		f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" fill-opacity="0.18" stroke="%s" stroke-width="1.5"/>`,
			ax.PX(-b.L/2), ax.PY(b.A), ax.PX(b.L/2)-ax.PX(-b.L/2), ax.PY(-b.A)-ax.PY(b.A), colMuted, colInk2))
		for _, p := range b.Poles() {
			col, lab := colOrange, "N"
			if p.Q < 0 {
				col, lab = colInk, "S"
			}
			x := ax.PX(p.Z)
			f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="5" fill="%s" stroke="%s" stroke-width="2"/>`, x, ax.PY(0), col, colSurface))
			f.Text(x, ax.PY(-b.A)+16, lab, "middle", colInk, 12, 600)
		}
		f.Text(70, y0-10, title, "start", colInk, 14, 600)
	}
	panel(100, u, "Uniform bar: magnetic charge on the two end faces — one pair; the outside lines run end to end")
	panel(530, f4, "Four segments S N | N S | S N | N S: the lines bow from pole to pole along the axis — a chain of lenses")
	f.Text(70, 916, "Field lines are contours of the flux function ψ; dots mark the magnetic charge −∇·M (orange N, black S). 400 sheet loops per bar.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "figM1_fieldlines.svg"))
}

// figScan: the bench test, B_z along a line beside the bar.
func figScan(dir string, u, f4 Bar) error {
	f := NewFig(960, 480, "Fig. M2 — The bench test: B_z along a line 3 mm beside the bar",
		"a Hall probe run parallel to the axis; sign changes count the poles it passes")
	gap := u.A + 3e-3
	zs, bu := u.Scan(gap, 601, u.Sheet(2000))
	_, b4 := f4.Scan(gap, 601, f4.Sheet(2000))
	mt := func(v []float64) []float64 {
		o := make([]float64, len(v))
		for i := range v {
			o[i] = v[i] * 1e3
		}
		return o
	}
	zmm := make([]float64, len(zs))
	for i := range zs {
		zmm[i] = zs[i] * 1e3
	}
	lo, hi := -250.0, 250.0
	a := f.Axes(90, 80, 830, 300, zmm[0], zmm[len(zmm)-1], lo, hi, false, false)
	a.Grid(Ticks(-100, 100, 25), Ticks(lo, hi, 100), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "position along the bar, z (mm)", "B_z (mT)")
	f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" fill-opacity="0.15"/>`, a.PX(-u.L/2*1e3), a.Y0, a.PX(u.L/2*1e3)-a.PX(-u.L/2*1e3), a.H, colMuted))
	f.Text(a.PX(0), a.Y0+16, "the bar", "middle", colInk2, 12, 400)
	a.Line([]float64{zmm[0], zmm[len(zmm)-1]}, []float64{0, 0}, colInk2, 1)
	a.Line(zmm, mt(bu), colBlue, 2)
	a.Line(zmm, mt(b4), colOrange, 2)
	f.Legend(90, 440, []string{fmt.Sprintf("uniform bar: %d sign changes", SignChanges(bu)), fmt.Sprintf("four segments: %d sign changes", SignChanges(b4))}, []string{colBlue, colOrange})
	f.Text(90, 466, "One sign change per pole site: a uniform bar gives 2; four pairs end to end give 5 (8 if separated). More than 2 on an ordinary bar supports Russell.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "figM2_scan.svg"))
}
