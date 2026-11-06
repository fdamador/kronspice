// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
)

// Figure 8: hydrogen under pressure. Two stacked charts on their own axes:
// the levels against the wall radius, and Lyman α against the pressure.
func figPressure(dir string, _, _ Spectrum) error {
	f := NewFig(960, 820, "Fig. 8 — High pressure, short waves: hydrogen compressed by a wall in Kron's network",
		"the network's outer short circuit moved in to radius r_c; Michels, de Boer & Bijl (1937)")
	var rcs, e1, e2s, e2p, pr, ly []float64
	for rc := 1.5; rc <= 10.0001; rc += 0.1 {
		c := Compress(rc)
		rcs, e1, e2s, e2p = append(rcs, rc), append(e1, c.E1s), append(e2s, c.E2s), append(e2p, c.E2p)
		if c.PressureGPa > 1e-3 {
			pr, ly = append(pr, c.PressureGPa), append(ly, c.LymanNm)
		}
	}
	a := f.Axes(90, 90, 830, 300, 1.5, 10, -0.6, 1.6, false, false) // curves above 1.6 Ha are clipped
	a.Grid(Ticks(2, 10, 1), Ticks(-0.6, 1.6, 0.4), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%.1f", v) }, "wall radius r_c (a0)", "energy (Ha)")
	a.Line([]float64{1.5, 10}, []float64{0, 0}, colMuted, 1)
	f.Text(a.PX(9.9), a.PY(0)-6, "zero: unbound above", "end", colInk2, 12, 400)
	a.Line(rcs, e1, colBlue, 2)
	a.Line(rcs, e2s, colOrange, 2)
	a.Line(rcs, e2p, colAqua, 2)
	f.Text(a.PX(9.8), a.PY(-0.5)-8, "1s", "end", colInk, 13, 600)
	f.Text(a.PX(5.1), a.PY(Compress(5.1).E2s)-8, "2s", "start", colInk, 13, 600)
	f.Text(a.PX(4.3), a.PY(Compress(4.3).E2p)+18, "2p", "start", colInk, 13, 600)
	for _, x := range []float64{1.835, 2} {
		a.VLine(x, colInk2, 1, true)
	}
	f.Text(a.PX(1.835)-4, a.PY(-0.5), "1.835: ionized", "end", colInk2, 12, 400)
	f.Text(a.PX(2)+6, a.PY(-0.5)+4, "2.0: E_1s = −1/8 exactly", "start", colInk2, 12, 400)
	f.Legend(90, 448, []string{"1s", "2s", "2p"}, []string{colBlue, colOrange, colAqua})

	b := f.Axes(90, 510, 830, 220, 1e-3, 3e3, 0, 130, true, false)
	b.Grid(Decades(1e-3, 1e3, 1), Ticks(0, 120, 20), Pow10Label, func(v float64) string { return fmt.Sprintf("%g", v) },
		"pressure on the wall (GPa)", "Lyman α (nm)")
	b.Line(pr, ly, colBlue, 2)
	b.VLine(364, colInk, 1.5, true)
	f.Text(b.PX(364)-6, b.Y0+14, "Earth's centre, 364 GPa", "end", colInk, 12, 400)
	b.Line([]float64{1e-3, 3e3}, []float64{121.57, 121.57}, colMuted, 1)
	f.Text(b.PX(1.2e-3), b.PY(121.57)+16, "free atom 121.57 nm", "start", colInk2, 12, 400)
	f.Text(90, 790, "The atom's pressure on the wall is P = −dE_1s/dV. Squeezed, every level rises and every line shortens; the 2s = 2p of the free atom splits.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig8_pressure.svg"))
}

// Figure 9: one sphere (1s) and two spheres meeting at a plane (2p_z).
func figOrbitals(dir string, _, _ Spectrum) error {
	f := NewFig(960, 540, "Fig. 9 — One sphere and two: the 1s and 2p_z states of hydrogen",
		"electron density in a plane through the axis, each scaled to its own peak; Kron's network gives both levels")
	ramp := []string{"#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6"}
	panel := func(x0 float64, title string, dens func(x, z float64) float64, ext float64) {
		ax := f.Axes(x0, 90, 360, 360, -ext, ext, -ext, ext, false, false)
		peak := 0.0
		const n = 72
		cell := 2 * ext / n
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				peak = math.Max(peak, dens(-ext+(float64(i)+0.5)*cell, -ext+(float64(j)+0.5)*cell))
			}
		}
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				x, z := -ext+float64(i)*cell, -ext+float64(j)*cell
				d := dens(x+cell/2, z+cell/2) / peak
				k := int(math.Sqrt(d) * float64(len(ramp)+1))
				if k == 0 {
					continue
				}
				if k > len(ramp) {
					k = len(ramp)
				}
				px, py := ax.PX(x), ax.PY(z+cell)
				f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill="%s"/>`, px, py, ax.PX(x+cell)-px+0.3, ax.PY(z)-py+0.3, ramp[k-1]))
			}
		}
		ax.Grid(Ticks(-ext, ext, ext/2), Ticks(-ext, ext, ext/2), func(v float64) string { return fmt.Sprintf("%g", v) },
			func(v float64) string { return fmt.Sprintf("%g", v) }, "x (a0)", "z (a0)")
		f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="4" fill="%s" stroke="%s" stroke-width="2"/>`, ax.PX(0), ax.PY(0), colInk, colSurface))
		f.Text(x0+180, 80, title, "middle", colInk, 14, 600)
	}
	panel(90, "1s: E = −1/2 Ha, one sphere", func(x, z float64) float64 { return math.Exp(-2 * math.Hypot(x, z)) }, 4)
	panel(560, "2p_z: E = −1/8 Ha, two lobes, nodal plane z = 0", func(x, z float64) float64 {
		r := math.Hypot(x, z)
		return z * z * math.Exp(-r)
	}, 10)
	f.Text(90, 528, "Russell's two spheres meeting at a central plane, with a point between them, have the shape of 2p_z: the electron is never on the plane z = 0.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig9_orbitals.svg"))
}
