// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
	"sort"
)

// WriteFigures draws the report's figures into dir from the computed
// spectra (hydrogen network, light-atom network).
func WriteFigures(dir string, hyd, light Spectrum) error {
	steps := []func(string, Spectrum, Spectrum) error{figResponse, figSpectrum, figFrequencies, figHyperfine, figLoops, figScales, figAtom3D, figPressure, figOrbitals}
	for _, s := range steps {
		if err := s(dir, hyd, light); err != nil {
			return err
		}
	}
	return nil
}

func pct(v []float64, p float64) float64 {
	s := append([]float64(nil), v...)
	sort.Float64s(s)
	return s[int(p*float64(len(s)-1))]
}

// Figure 1: the frequency response of the whole network.
func figResponse(dir string, hyd, light Spectrum) error {
	f := NewFig(960, 560, "Fig. 1 — Frequency response of Kron's network: every pole is a level of the atom",
		"kronspice AC sweep of the driving-point impedance, plotted against x = 1/√(−2E); poles sit at x = n/Z")
	lo := math.Floor(math.Min(pct(hyd.LogZ, 0.01), pct(light.LogZ, 0.01)))
	hi := math.Ceil(math.Max(pct(hyd.LogZ, 0.995), pct(light.LogZ, 0.995)))
	a := f.Axes(90, 90, 830, 340, 0.75, 5.8, lo, hi, false, false)
	a.Grid(Ticks(1, 5.8, 1), Ticks(lo, hi, 1), func(v float64) string { return fmt.Sprintf("%g", v) },
		func(v float64) string { return Pow10Label(math.Pow(10, v)) }, "x = 1/√(−2E)   (E in hartree; a level n of charge product Z sits at x = n/Z)", "|Z_in| (Ω)")
	a.Line(hyd.X, hyd.LogZ, colBlue, 1.5)
	a.Line(light.X, light.LogZ, colOrange, 1.5)
	for k := 0; k < NLevels; k++ {
		for _, s := range []struct {
			sp  Spectrum
			col string
			dy  float64
		}{{hyd, colBlue, 14}, {light, colOrange, 30}} {
			x := 1 / math.Sqrt(-2*s.sp.Circuit[k])
			a.VLine(x, s.col, 1, true)
			f.Text(a.PX(x)+3, a.Y0+s.dy, fmt.Sprintf("%ds", k+1), "start", colInk, 12, 600)
		}
	}
	f.Legend(90, 512, []string{"hydrogen network (q = e, Z = 1)", "light atom (q = 0.910 e, Z = 0.910)"}, []string{colBlue, colOrange})
	f.Text(90, 540, fmt.Sprintf("Each network: %d (Z = 1) and %d (Z = 0.910) nodes, %d frequencies; the dashed lines are the poles pinned by bisection.", hyd.Nodes, light.Nodes, hyd.SweepPts), "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig1_response.svg"))
}

// Figure 2: the line spectrum from the network's poles against measurement.
func figSpectrum(dir string, hyd, light Spectrum) error {
	f := NewFig(960, 440, "Fig. 2 — The spectrum read from the poles, against measured hydrogen",
		"Lyman and Balmer lines from the levels of Fig. 1 (h → 0), vacuum wavelengths")
	a := f.Axes(230, 80, 700, 240, 80, 820, 0, 3, false, false)
	a.Grid(Ticks(100, 800, 100), nil, func(v float64) string { return fmt.Sprintf("%g", v) }, nil, "wavelength (nm)", "")
	rows := []struct {
		name  string
		col   string
		nm    func(Line) float64
		lines []Line
	}{
		{"measured hydrogen (NIST)", colInk, func(l Line) float64 { return l.MeasuredNm }, hyd.Lines()},
		{"network, q = e", colBlue, func(l Line) float64 { return l.Circuit }, hyd.Lines()},
		{"network, q = 0.910 e (H23)", colOrange, func(l Line) float64 { return l.Circuit }, light.Lines()},
	}
	for i, r := range rows {
		y0 := 2.85 - float64(i)
		f.Text(a.X0-14, a.PY(y0-0.35)+4, r.name, "end", colInk, 13, 400)
		for _, l := range r.lines {
			x := a.PX(r.nm(l))
			f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5" stroke-linecap="round"/>`, x, a.PY(y0-0.05), x, a.PY(y0-0.65), r.col))
		}
		for _, l := range r.lines {
			if l.Up-l.Lo == 1 || l.Lo == 2 {
				name := map[string]string{"Lyman α": "Lyα", "Balmer α (Hα)": "Hα", "Balmer β (Hβ)": "Hβ", "Balmer γ (Hγ)": "Hγ"}[l.Name]
				if name != "" {
					f.Text(a.PX(r.nm(l))+5, a.PY(y0-0.75)+4, name, "start", colInk2, 11, 400)
				}
			}
		}
	}
	f.Text(a.X0, 404, "Blue lies on black to the 0.01 nm shown (Bohr levels; the 10⁻⁵ left is fine structure and the Lamb shift).", "start", colInk2, 12, 400)
	f.Text(a.X0, 424, "Orange is the same pattern stretched by 1/Z² = 1.207: the light atom's lines are 21 % too long.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig2_spectrum.svg"))
}

// rowItem is one labelled value in a row-per-item dot plot.
type rowItem struct {
	name string
	v    float64
	col  string
	unit string
}

// rowDots draws one row per item on a shared logarithmic axis: the name on
// the left, a dot at the value, the value beside the dot.
func rowDots(f *Fig, items []rowItem, lo, hi float64, every int, xlabel string) {
	rowH := 30.0
	a := f.Axes(300, 80, 600, rowH*float64(len(items)), lo, hi, 0, 1, true, false)
	a.Grid(Decades(lo, hi, every), nil, Pow10Label, nil, xlabel, "")
	for i, it := range items {
		y := a.Y0 + rowH*(float64(i)+0.5)
		f.Text(a.X0-14, y+4, it.name, "end", colInk, 13, 400)
		x := a.PX(it.v)
		f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="6" fill="%s" stroke="%s" stroke-width="2"/>`, x, y, it.col, colSurface))
		anchor, dx := "start", 12.0
		if x > a.X0+a.W*0.75 {
			anchor, dx = "end", -12
		}
		f.Text(x+dx, y+4, fmt.Sprintf("%.3g %s", it.v, it.unit), anchor, colInk2, 12, 400)
	}
}

// Figure 3: every frequency in the model on one axis.
func figFrequencies(dir string, hyd, _ Spectrum) error {
	k := Knot(WvdM().R)
	items := []rowItem{
		{"Kron network poles (circuit units)", 1 / math.Sqrt(-0.5+NewLadderV0(1, 0.1)) / (2 * math.Pi), colMuted, "Hz"},
		{"21 cm line (measured)", HFMeasHz, colInk, "Hz"},
		{"Hα (measured)", C0 / 656.46e-9, colInk, "Hz"},
		{"1S–2S (measured)", OneSTwoSHz, colInk, "Hz"},
		{"electron ring, m_e c²/h (H23)", MeMeV * MeV / H, colOrange, "Hz"},
		{"knot squashing tone (H29)", k.ExciteMeV * MeV / H, colOrange, "Hz"},
		{"proton light ring, m_p c²/h (H31)", MpMeV * MeV / H, colBlue, "Hz"},
	}
	f := NewFig(960, 80+30*float64(len(items))+110, "Fig. 3 — Where the frequencies are: the lines are not the rings' own frequencies",
		"one logarithmic axis")
	rowDots(f, items, 1e-2, 1e24, 2, "frequency (Hz)")
	y := 80 + 30*float64(len(items)) + 70
	f.Text(40, y, "The atom's lines (black) are differences of Coulomb levels, 10⁵ below the electron ring's own frequency and 10⁸ below the proton's.", "start", colInk2, 12, 400)
	f.Text(40, y+18, "The network's poles are an analog: E = 1/ω² − V0 maps them onto the levels, so their frequencies carry the levels, not the light.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig3_frequencies.svg"))
}

func sign(v float64) float64 {
	if v < 0 {
		return -1
	}
	return 1
}

// NewLadderV0 is the shift V0 of a network (for the circuit-unit frequency).
func NewLadderV0(z, h float64) float64 { return z/h + 1 }

// Figure 4: the 21 cm line by model.
func figHyperfine(dir string, _, _ Spectrum) error {
	pr := Protons()
	type row struct {
		name string
		a    Atom
		col  string
	}
	rows := []row{
		{"q = 0.910 e + light ring", Atom{WvdM(), pr[0]}, colOrange},
		{"q = 0.910 e + Skyrme knot", Atom{WvdM(), pr[2]}, colOrange},
		{"q = 0.910 e + quark rings", Atom{WvdM(), pr[1]}, colOrange},
		{"q = e + light ring", Atom{WvdMUnit(), pr[0]}, colBlue},
		{"q = e + Skyrme knot", Atom{WvdMUnit(), pr[2]}, colBlue},
		{"q = e + quark rings", Atom{WvdMUnit(), pr[1]}, colBlue},
		{"measured e + measured p", Atom{Measured(), pr[3]}, colInk},
	}
	f := NewFig(960, 120+38*float64(len(rows))+70, "Fig. 4 — The 21 cm line from the current-loop circuit, model by model",
		"ν / measured 1420.406 MHz; kronspice loop interaction × 4 / h")
	a := f.Axes(260, 80, 640, 38*float64(len(rows)), 0, 1.1, 0, float64(len(rows)), false, false)
	a.Grid(Ticks(0, 1.1, 0.2), nil, func(v float64) string { return fmt.Sprintf("%.1f", v) }, nil, "ν / ν measured", "")
	a.VLine(1, colInk, 1.5, true)
	for i, r := range rows {
		hf, err := HyperfineCircuit(r.a, 8)
		if err != nil {
			return err
		}
		v := hf.Nu / HFMeasHz
		y := a.Y0 + 38*float64(i) + 9
		f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="20" rx="4" fill="%s"/>`, a.X0, y, a.PX(v)-a.X0, r.col))
		f.Text(a.X0-12, y+15, r.name, "end", colInk, 13, 400)
		f.Text(a.PX(v)+8, y+15, fmt.Sprintf("%.1f MHz", hf.Nu/1e6), "start", colInk2, 12, 400)
	}
	f.Legend(260, f.H-24, []string{"toroidal-photon charge 0.910 e", "charge set to e", "measured particles"}, []string{colOrange, colBlue, colInk})
	return f.Save(filepath.Join(dir, "fig4_hyperfine.svg"))
}

// Figure 5: the electron's spin cloud as current loops round the proton.
func figLoops(dir string, _, _ Spectrum) error {
	a0 := A0 * Recoil
	f := NewFig(960, 560, "Fig. 5 — The 21 cm circuit: the electron's spin cloud as coaxial current loops",
		"half-plane through the axis, equal scales; shading = 1s density |ψ|², dots = the 32 loops (area ∝ current)")
	ax := f.Axes(110, 80, 240, 384, 0, 1.6, -1.6, 1.6, false, false) // 150 px per a0 both ways
	ramp := []string{"#fcfcfb", "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5"}
	const cell = 0.04
	for rho := 0.0; rho < 1.6-1e-9; rho += cell {
		for z := -1.6; z < 1.6-1e-9; z += cell {
			r := math.Hypot(rho+cell/2, z+cell/2)
			k := int(math.Sqrt(math.Exp(-2*r)) * float64(len(ramp)))
			if k >= len(ramp) {
				k = len(ramp) - 1
			}
			if k == 0 {
				continue
			}
			x0, y0 := ax.PX(rho), ax.PY(z+cell)
			f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>`, x0, y0, ax.PX(rho+cell)-x0+0.4, ax.PY(z)-y0+0.4, ramp[k]))
		}
	}
	ax.Grid(Ticks(0, 1.6, 0.4), Ticks(-1.6, 1.6, 0.4), func(v float64) string { return fmt.Sprintf("%.1f", v) }, func(v float64) string { return fmt.Sprintf("%.1f", v) }, "ρ / a0", "z / a0")
	loops := SpinLoops(Measured().Moment(), a0, 8)
	imax := 0.0
	for _, l := range loops {
		imax = math.Max(imax, l.I)
	}
	for _, l := range loops {
		ax.Dot(l.Rho/a0, l.Z/a0, 3+6*math.Sqrt(l.I/imax), colOrange)
	}
	f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="5" fill="none" stroke="%s" stroke-width="2.5"/>`, ax.PX(0), ax.PY(0), colInk))
	f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>`, ax.PX(0)+5, ax.PY(0)+5, ax.PX(0.5), ax.PY(-1.35), colInk2))
	f.Text(ax.PX(0.52), ax.PY(-1.35)+4, "proton ring, 4×10⁻⁶ a0 (drawn enlarged)", "start", colInk, 12, 400)
	x := 470.0
	lines := []string{
		"A spin cloud is a magnetization M(r) = μ_e|ψ(r)|² ẑ.",
		"A magnetization is a current: J = ∇ × M = (2/a)·M·sin θ φ̂,",
		"so the cloud is a set of loops coaxial with the proton.",
		"",
		"Each loop and the proton ring is a kronspice inductor",
		"driven to its current, coupled to the ring by Maxwell's",
		"coaxial-circle mutual inductance.",
		"",
		"ΣM·I_p·I_k (solved currents) = μ_p·(2/3)μ0·M(0),",
		"Fermi's contact energy, to 1.3×10⁻⁵.",
		"The ledger's stored energy = ½ΣLI² + ΣMII to 4×10⁻¹⁶.",
		"",
		"Loops: 8 radial steps (exact weights for e^(−2r/a)),",
		"4 Gauss–Legendre angles (exact for sin²θ).",
	}
	for i, s := range lines {
		f.Text(x, 110+20*float64(i), s, "start", colInk, 13, 400)
	}
	return f.Save(filepath.Join(dir, "fig5_loops.svg"))
}

// Figure 6: sizes on one axis.
func figScales(dir string, _, _ Spectrum) error {
	wv := WvdM()
	items := []rowItem{
		{"electron size bound (g − 2)", 1e-3, colInk, "fm"},
		{"proton light ring (H31)", ProtonRingRadius() * 1e15, colBlue, "fm"},
		{"quark rings, rms (H31)", Protons()[1].R * 1e15, colBlue, "fm"},
		{"proton, measured charge radius", RpMeas * 1e15, colInk, "fm"},
		{"knot scale ℓ (H29)", Knot(wv.R).L * 1e15, colOrange, "fm"},
		{"electron ring λ_C/4π (H23)", wv.R * 1e15, colOrange, "fm"},
		{"Bohr radius a0", A0 * 1e15, colInk, "fm"},
	}
	f := NewFig(960, 80+30*float64(len(items))+90, "Fig. 6 — Sizes: both rings are far smaller than the atom", "one logarithmic axis")
	rowDots(f, items, 1e-4, 1e6, 1, "size (fm)")
	f.Text(40, 80+30*float64(len(items))+70, "The electron ring is 0.0036 a0 and the proton ring 4×10⁻⁶ a0: at the atom's scale both are points, so the atom is the Coulomb problem.", "start", colInk2, 12, 400)
	return f.Save(filepath.Join(dir, "fig6_scales.svg"))
}
