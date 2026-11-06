// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
	"strings"
)

// contour returns the segments of the level set v = c on a grid (marching
// squares, linear interpolation along cell edges).
func contour(xs, ys []float64, v [][]float64, c float64) [][4]float64 {
	var seg [][4]float64
	for i := 0; i+1 < len(xs); i++ {
		for j := 0; j+1 < len(ys); j++ {
			p := [4][3]float64{{xs[i], ys[j], v[i][j]}, {xs[i+1], ys[j], v[i+1][j]}, {xs[i+1], ys[j+1], v[i+1][j+1]}, {xs[i], ys[j+1], v[i][j+1]}}
			var hits [][2]float64
			for e := 0; e < 4; e++ {
				a, b := p[e], p[(e+1)%4]
				if (a[2]-c)*(b[2]-c) < 0 {
					t := (c - a[2]) / (b[2] - a[2])
					hits = append(hits, [2]float64{a[0] + t*(b[0]-a[0]), a[1] + t*(b[1]-a[1])})
				}
			}
			for k := 0; k+1 < len(hits); k += 2 {
				seg = append(seg, [4]float64{hits[k][0], hits[k][1], hits[k+1][0], hits[k+1][1]})
			}
		}
	}
	return seg
}

// meridian draws one field-line panel in the (x, z) plane through the axis
// (z up, x = ±ρ): a density shading, flux-function contours, and a ring of
// radius ringR coloured by the sign of B_r (orange N, blue S).
type meridianSpec struct {
	X0, Y0, Side float64
	Ext          float64
	Title, Sub   string
	Density      func(rho, z float64) float64
	Flux         func(rho, z float64) float64
	Br           func(r, th float64) float64
	RingR        float64
	Dashed       float64      // radius of a dashed reference circle (0 = none)
	Marks        [][3]float64 // (z, moment sign, size) arrows on the axis
}

func meridian(f *Fig, s meridianSpec) {
	sc := s.Side / (2 * s.Ext)
	cx, cy := s.X0+s.Side/2, s.Y0+s.Side/2
	px := func(x, z float64) (float64, float64) { return cx + sc*x, cy - sc*z }
	f.Text(s.X0, s.Y0-24, s.Title, "start", colInk, 14, 600)
	if s.Sub != "" {
		f.Text(s.X0, s.Y0-8, s.Sub, "start", colInk2, 12, 400)
	}
	// density shading
	const nd = 64
	dpeak := 0.0
	dens := make([][]float64, nd)
	for i := range dens {
		dens[i] = make([]float64, nd)
		for j := range dens[i] {
			x := -s.Ext + (float64(i)+0.5)*2*s.Ext/nd
			z := -s.Ext + (float64(j)+0.5)*2*s.Ext/nd
			if s.Density != nil {
				dens[i][j] = s.Density(math.Abs(x), z)
			}
			dpeak = math.Max(dpeak, dens[i][j])
		}
	}
	var sb strings.Builder
	cell := s.Side / nd
	for i := range dens {
		for j := range dens[i] {
			if dpeak == 0 {
				continue
			}
			t := math.Pow(dens[i][j]/dpeak, 0.5)
			if t < 0.02 {
				continue
			}
			fmt.Fprintf(&sb, `<rect x="%.1f" y="%.1f" width="%.2f" height="%.2f" fill-opacity="%.3f"/>`, s.X0+float64(i)*cell, s.Y0+float64(nd-1-j)*cell, cell+0.3, cell+0.3, 0.55*t)
		}
	}
	f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s"/>`, s.X0, s.Y0, s.Side, s.Side, colSurface, colGrid))
	f.Raw(fmt.Sprintf(`<g fill="%s">%s</g>`, colAqua, sb.String()))
	// flux contours
	const ng = 161
	xs, zs := make([]float64, ng), make([]float64, ng)
	for i := range xs {
		xs[i] = s.Ext * float64(i) / float64(ng-1)
		zs[i] = -s.Ext + 2*s.Ext*float64(i)/float64(ng-1)
	}
	psi := make([][]float64, ng)
	parallel(ng, func(i int) {
		psi[i] = make([]float64, ng)
		for j := range zs {
			psi[i][j] = s.Flux(xs[i], zs[j])
		}
	})
	peak := 0.0
	for i := range psi {
		for j := range psi[i] {
			if xs[i] > 0.04*s.Ext && math.Hypot(xs[i], zs[j]) > 0.06*s.Ext {
				peak = math.Max(peak, math.Abs(psi[i][j]))
			}
		}
	}
	sb.Reset()
	for _, frac := range []float64{0.002, 0.006, 0.015, 0.03, 0.055, 0.09, 0.14, 0.21, 0.3, 0.42, 0.57, 0.75, 0.95} {
		for _, sg := range []float64{1, -1} {
			for _, sgm := range contour(xs, zs, psi, sg*frac*peak) {
				for _, m := range []float64{1, -1} {
					x1, y1 := px(m*sgm[0], sgm[1])
					x2, y2 := px(m*sgm[2], sgm[3])
					fmt.Fprintf(&sb, "M%.1f %.1fL%.1f %.1f", x1, y1, x2, y2)
				}
			}
		}
	}
	f.Raw(fmt.Sprintf(`<path d="%s" fill="none" stroke="%s" stroke-width="1" stroke-opacity="0.85" stroke-linecap="round"/>`, sb.String(), colInk))
	// ring coloured by the sign of B_r
	if s.RingR > 0 {
		const n = 180
		for k := 0; k < n; k++ {
			a0 := 2 * math.Pi * float64(k) / n
			a1 := 2 * math.Pi * float64(k+1) / n
			th := math.Abs(math.Mod(0.5*(a0+a1)+math.Pi, 2*math.Pi) - math.Pi) // polar angle from +z
			col := colOrange
			if s.Br(s.RingR, th) < 0 {
				col = colBlue
			}
			x1, y1 := px(s.RingR*math.Sin(a0), s.RingR*math.Cos(a0))
			x2, y2 := px(s.RingR*math.Sin(a1), s.RingR*math.Cos(a1))
			f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="5"/>`, x1, y1, x2, y2, col))
		}
	}
	if s.Dashed > 0 && s.Dashed < s.Ext {
		f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="%.1f" fill="none" stroke="%s" stroke-width="1.2" stroke-dasharray="5 4"/>`, cx, cy, s.Dashed*sc, colInk2))
	}
	for _, m := range s.Marks {
		x, y := px(0, m[0])
		l := m[2] * sc
		xo := x - 0.32*l*boolf(l > 30) // electron arrows beside the axis, proton arrows on it
		for _, st := range []struct {
			col string
			w   float64
		}{{colSurface, 6}, {colOrange, 2.6}} {
			f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%g" stroke-linecap="round"/>`, xo, y+m[1]*l/2, xo, y-m[1]*l/2, st.col, st.w))
		}
		arrowHead(f, xo, y-m[1]*l/2, -m[1]*math.Pi/2, 8, colOrange)
	}
	// axis label
	f.Text(s.X0+s.Side-6, s.Y0+s.Side-8, fmt.Sprintf("±%.0f a0", s.Ext), "end", colInk2, 11, 400)
}

func arrowHead(f *Fig, x, y, ang, s float64, col string) {
	p := func(da float64) string { return fmt.Sprintf("%.1f,%.1f", x-s*math.Cos(ang+da), y-s*math.Sin(ang+da)) }
	f.Raw(fmt.Sprintf(`<polygon points="%.1f,%.1f %s %s" fill="%s"/>`, x+0.3*s*math.Cos(ang), y+0.3*s*math.Sin(ang), p(0.45), p(-0.45), col))
}

// AtomFigures writes Figs. H1–H4 and the sphere animation.
func AtomFigures(dir string, as []AtomResult) error {
	// H1: field lines of the four levels
	f := NewFig(1180, 1210, "Fig. H1 — The magnetic field of one hydrogen atom, level by level",
		"field lines (contours of ρA_φ) in a plane through the spin axis; shading = electron density |ψ|²; ring at 8 a0 coloured by the sign of B_r (orange N, blue S)")
	subs := []string{
		"one pair: B_r is N over one cap, S over the other",
		"j = ½: no octupole (2j < 3) — one pair",
		"j = 3⁄2, m = 3⁄2: spin octupole present, but the dipole leads everywhere",
		"j = 3⁄2, m = ½: four zones N S N S round the atom out to 15.6 a0 (dashed)",
	}
	for i, a := range as {
		a := a
		x0 := 40 + float64(i%2)*570
		y0 := 120 + float64(i/2)*540
		meridian(f, meridianSpec{X0: x0, Y0: y0, Side: 500, Ext: 18, Title: a.F.S.Name, Sub: subs[i],
			Density: func(rho, z float64) float64 { return a.F.S.Density(math.Hypot(rho, z), math.Atan2(rho, z)) },
			Flux:    a.F.Flux, Br: a.F.Br, RingR: 8, Dashed: a.Rc * boolf(i == 3)})
	}
	f.Text(40, 1196, "Electron field from its orbital current and spin magnetization (Pauli spinors, g_s = 2.0023), exact Legendre-series Green's function; proton spin dipole included.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figH1_atom_fieldlines.svg")); err != nil {
		return err
	}

	// H2: B_r on spheres of growing radius
	f = NewFig(1100, 470, "Fig. H2 — B_r round the atom on spheres of growing radius (each curve scaled to its own peak)",
		"from the north pole (0°) to the south pole (180°); each sign change separates an N zone from an S zone")
	radii := []float64{3, 6, 10, 15, 24}
	cols := []string{colInk, colBlue, colAqua, colOrange, colMuted}
	for p, idx := range []int{3, 2} {
		a := as[idx]
		ax := f.Axes(80+float64(p)*520, 90, 440, 260, 0, 180, -1.1, 1.1, false, false)
		ax.Grid(Ticks(0, 180, 30), Ticks(-1, 1, 0.5), func(v float64) string { return fmt.Sprintf("%g°", v) },
			func(v float64) string { return fmt.Sprintf("%g", v) }, "polar angle θ", "B_r / max|B_r|")
		f.Text(80+float64(p)*520, 82, a.F.S.Name, "start", colInk, 14, 600)
		ax.Line([]float64{0, 180}, []float64{0, 0}, colInk2, 1)
		for k, r := range radii {
			th, v := make([]float64, 361), make([]float64, 361)
			mx := 0.0
			for i := range th {
				th[i] = 180 * float64(i) / 360
				v[i] = a.F.Br(r, math.Max(1e-6, th[i]*math.Pi/180))
				mx = math.Max(mx, math.Abs(v[i]))
			}
			for i := range v {
				v[i] /= mx
			}
			ax.Line(th, v, cols[k], 2)
		}
	}
	var names []string
	for _, r := range radii {
		names = append(names, fmt.Sprintf("r = %g a0", r))
	}
	f.Legend(80, 430, names, cols)
	f.Text(80, 456, "m = ½: three sign changes (four zones, N S N S) out to 15.6 a0, then one (two zones). m = 3⁄2: one sign change at every radius.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figH2_sphere_scans.svg")); err != nil {
		return err
	}

	// H3: zones against radius
	f = NewFig(960, 440, "Fig. H3 — How many N/S zones a sphere round the atom sees, against its radius",
		"two zones = one pair of poles (dipole); four zones = the axial octupole's bands")
	ax := f.Axes(90, 80, 820, 260, 0.5, 60, 0, 5, true, false)
	ax.Grid([]float64{0.5, 1, 2, 5, 10, 20, 50}, []float64{0, 1, 2, 3, 4, 5}, func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "sphere radius (a0, log scale)", "zones")
	scols := []string{colMuted, colAqua, colBlue, colOrange}
	names = nil
	for i, a := range as {
		y := make([]float64, len(a.Zones))
		for k := range y {
			y[k] = float64(a.Zones[k]) + 0.06*float64(i-1)
		}
		ax.Line(a.Radii, y, scols[i], 2.5)
		names = append(names, a.F.S.Name)
	}
	ax.VLine(as[3].Rc, colOrange, 1, true)
	f.Text(ax.PX(as[3].Rc)+6, 96, fmt.Sprintf("r_c = %.1f a0", as[3].Rc), "start", colOrange, 12, 500)
	f.Legend(90, 400, names, scols)
	f.Text(90, 426, "Only 2p3⁄2 with m = ±½ shows four zones outside the atom. The ground state, and every j = ½ level, shows one pair at every radius.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figH3_zones.svg")); err != nil {
		return err
	}

	// H4: the bar's bench scan on the atom
	f = NewFig(960, 440, "Fig. H4 — The octave bar's bench test on the atom: B_z along a line 5 a0 beside the spin axis",
		"the scan of examples/octave_magnet: a dipole gives two sign changes; more mean more pole sites along the line")
	ax = f.Axes(90, 80, 820, 260, -30, 30, -1.1, 1.1, false, false)
	ax.Grid(Ticks(-30, 30, 10), Ticks(-1, 1, 0.5), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "position along the axis z (a0)", "B_z / max|B_z|")
	ax.Line([]float64{-30, 30}, []float64{0, 0}, colInk2, 1)
	names = nil
	for i, a := range as {
		fa := a.F
		zs, b := LineScan(func(rho, z float64) float64 { _, bz := fa.Cyl(rho, z); return bz }, 5, 30, 601)
		mx := 0.0
		for _, v := range b {
			mx = math.Max(mx, math.Abs(v))
		}
		for k := range b {
			b[k] /= mx
		}
		ax.Line(zs, b, scols[i], 2)
		names = append(names, fmt.Sprintf("%s: %d", a.F.S.Name, SignChanges(b)))
	}
	f.Legend(90, 400, names, scols)
	f.Text(90, 426, "Sign changes in the legend. 2p3⁄2 m = ½ gives four at this distance: two N and two S sites along the line, the four-zone pattern seen from the side.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figH4_line_scan.svg")); err != nil {
		return err
	}
	return AtomGIF(filepath.Join(dir, "atom_spheres.gif"), []AtomField{as[1].F, as[2].F, as[3].F})
}

func boolf(b bool) float64 {
	if b {
		return 1
	}
	return 0
}

// MoleculeFigures writes Figs. M1–M4 and the bond animation.
func MoleculeFigures(dir string, ms []MolResult) error {
	// M1: Heitler–London energies
	f := NewFig(960, 460, "Fig. M1 — Why the two electron magnets pair off: the Heitler–London energies of H₂",
		"energy of two hydrogen atoms against their separation, relative to the separated atoms")
	ax := f.Axes(90, 80, 820, 280, 0.8, 6, -5, 6, false, false)
	ax.Grid(Ticks(1, 6, 1), Ticks(-4, 6, 2), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "proton separation R (a0)", "E − 2E_H (eV)")
	var rs, es, et []float64
	for i := 0; i <= 260; i++ {
		r := 0.8 + 5.2*float64(i)/260
		rs = append(rs, r)
		es = append(es, (HLEnergy(r, 1)+1)*Hartree)
		et = append(et, (HLEnergy(r, -1)+1)*Hartree)
	}
	ax.Line([]float64{0.8, 6}, []float64{0, 0}, colInk2, 1)
	ax.Line(rs, es, colBlue, 2.5)
	ax.Line(rs, et, colOrange, 2.5)
	re, de := HLMinimum()
	ax.Dot(re, -de*Hartree, 5, colBlue)
	ax.Dot(1.4011, -4.748, 5, colInk)
	f.Text(ax.PX(re)+8, ax.PY(-de*Hartree)+16, fmt.Sprintf("Heitler–London: %.2f a0, %.2f eV", re, de*Hartree), "start", colBlue, 12, 500)
	f.Text(ax.PX(1.4011)+8, ax.PY(-4.748)+4, "measured: 1.40 a0, 4.75 eV", "start", colInk, 12, 500)
	f.Legend(90, 420, []string{"singlet: electron spins opposite (↑↓) — binds", "triplet: spins parallel (↑↑) — repels at every R", "measured bond"}, []string{colBlue, colOrange, colInk})
	f.Text(90, 446, "The bond exists only with the electron spins opposite, so in the molecule the two electron magnets cancel.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figM1_energies.svg")); err != nil {
		return err
	}

	// M2: densities along the bond
	f = NewFig(960, 440, "Fig. M2 — Electron density and spin density along the bond axis (R = 1.40 a0)",
		"the spin density is what makes the electrons magnetic; in the bound (singlet) molecule it is zero everywhere")
	ax = f.Axes(90, 80, 820, 260, -4, 4, -0.05, 0.75, false, false)
	ax.Grid(Ticks(-4, 4, 1), Ticks(0, 0.7, 0.1), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%.1f", v) }, "position along the bond z (a0)", "per a0³")
	sing, trip := ms[0].F.M, ms[2].F.M
	var zs, d1, s1, s3 []float64
	for i := 0; i <= 400; i++ {
		z := -4 + 8*float64(i)/400
		zs = append(zs, z)
		d1 = append(d1, sing.Density(0, z))
		s1 = append(s1, sing.SpinDensity(0, z))
		s3 = append(s3, trip.SpinDensity(0, z))
	}
	ax.Line(zs, d1, colMuted, 2)
	ax.Line(zs, s3, colOrange, 2.5)
	ax.Line(zs, s1, colBlue, 2.5)
	for _, z := range []float64{-0.70055, 0.70055} {
		ax.VLine(z, colInk2, 1, true)
	}
	f.Legend(90, 400, []string{"electron density, singlet", "spin density, triplet (↑↑)", "spin density, singlet (↑↓) = 0"}, []string{colMuted, colOrange, colBlue})
	f.Text(90, 426, "Dashed lines: the protons. The triplet's spin density sits on both atoms (cusps at the protons); the singlet's vanishes identically.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figM2_spin_density.svg")); err != nil {
		return err
	}

	// M3: field lines of the four magnets
	f = NewFig(1180, 700, "Fig. M3 — The field of H₂'s four magnets, in a plane through the bond (bond vertical)",
		"field lines (contours of ρA_φ); shading = spin density; ring at 1 a0 coloured by the sign of B_r (orange N, blue S); arrows = the magnetic moments (orange)")
	titles := []string{"Singlet + ortho protons (bound)", "Triplet + para protons (unbound)", "Triplet + ortho: all four aligned (unbound)"}
	subs := []string{"electrons cancel; two proton dipoles of 0.0015 μB", "the two electron spins, −2.0 μB in all", "electron spins plus proton spins"}
	for k, idx := range []int{1, 2, 3} {
		m := ms[idx]
		mf := m.F
		var marks [][3]float64
		if mf.M.Triplet {
			marks = append(marks, [3]float64{-0.70055, -1, 0.9}, [3]float64{0.70055, -1, 0.9})
		}
		if mf.M.Ortho {
			marks = append(marks, [3]float64{-0.70055, 1, 0.45}, [3]float64{0.70055, 1, 0.45})
		}
		meridian(f, meridianSpec{X0: 40 + float64(k)*380, Y0: 130, Side: 360, Ext: 3.2, Title: titles[k], Sub: subs[k],
			Density: mf.M.SpinDensity, Flux: func(rho, z float64) float64 { _, _, p := mf.Field(rho, z); return p },
			Br: func(r, th float64) float64 {
				bp, bz, _ := mf.Field(r*math.Sin(th), r*math.Cos(th))
				return bp*math.Sin(th) + bz*math.Cos(th)
			}, RingR: 1.0, Marks: marks})
	}
	f.Text(40, 530, "Long arrows, beside the axis: electron spin moments (opposite to the spins, the electron being negative). Short arrows, on the axis: proton moments (drawn", "start", colInk2, 12, 400)
	f.Text(40, 548, "far larger than their 1/658 share). In the bound singlet only the protons remain, and the para form (proton spins opposite) has no field at all.", "start", colInk2, 12, 400)
	f.Text(40, 576, "Zones round the molecule (Fig. M4): ortho protons alone show four zones between 0.3 and 1.6 a0 — two dipoles end to end, each 658 times weaker than one", "start", colInk2, 12, 400)
	f.Text(40, 594, "electron's. With the triplet electrons the field is one pair at every radius except a thin shell round the protons (0.6–0.8 a0).", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figM3_h2_fieldlines.svg")); err != nil {
		return err
	}

	// M4: zones against radius
	f = NewFig(960, 440, "Fig. M4 — N/S zones on spheres round the H₂ molecule, against their radius",
		"centred between the protons (each at 0.70 a0 from the centre)")
	ax = f.Axes(90, 80, 820, 260, 0.2, 20, 0, 5, true, false)
	ax.Grid([]float64{0.2, 0.5, 1, 2, 5, 10, 20}, []float64{0, 1, 2, 3, 4, 5}, func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "sphere radius (a0, log scale)", "zones")
	mcols := []string{colBlue, colAqua, colOrange}
	var names []string
	names = nil
	for k, idx := range []int{1, 2, 3} {
		y := make([]float64, len(ms[idx].Zones))
		for i := range y {
			y[i] = float64(ms[idx].Zones[i]) + 0.06*float64(k-1)
		}
		ax.Line(ms[idx].Radii, y, mcols[k], 2.5)
		names = append(names, titles[k])
	}
	ax.VLine(0.70055, colInk2, 1, true)
	f.Text(ax.PX(0.70055)+6, 96, "protons", "start", colInk2, 12, 400)
	f.Legend(90, 400, names, mcols)
	f.Text(90, 426, "Singlet + para (bound, the commonest form at low temperature) has no field and no zones; it is left off the chart.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figM4_h2_zones.svg")); err != nil {
		return err
	}
	return BondGIF(filepath.Join(dir, "h2_bond.gif"))
}
