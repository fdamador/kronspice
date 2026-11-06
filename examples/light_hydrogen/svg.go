// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"html"
	"math"
	"os"
	"strings"
)

// A small SVG plotter for the report's figures (standard library only).
// Colours are the validated reference palette's light-mode steps; every
// figure carries its own light surface so it reads the same in light and
// dark viewers.
const (
	colSurface = "#fcfcfb"
	colInk     = "#0b0b0b"
	colInk2    = "#52514e"
	colGrid    = "#e4e3df"
	colBlue    = "#2a78d6" // categorical slot 1
	colOrange  = "#eb6834" // slot 2
	colAqua    = "#1baf7a" // slot 3
	colMuted   = "#9a988f"
)

// Fig is one SVG figure.
type Fig struct {
	W, H float64
	b    strings.Builder
}

// NewFig starts a figure with a title and an optional subtitle.
func NewFig(w, h float64, title, sub string) *Fig {
	f := &Fig{W: w, H: h}
	fmt.Fprintf(&f.b, `<svg xmlns="http://www.w3.org/2000/svg" width="%g" height="%g" viewBox="0 0 %g %g" font-family="system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif">`+"\n", w, h, w, h)
	fmt.Fprintf(&f.b, `<rect width="%g" height="%g" fill="%s"/>`+"\n", w, h, colSurface)
	f.Text(24, 30, title, "start", colInk, 17, 600)
	if sub != "" {
		f.Text(24, 52, sub, "start", colInk2, 13, 400)
	}
	return f
}

// Text writes a label at pixel coordinates.
func (f *Fig) Text(x, y float64, s, anchor, color string, size float64, weight int) {
	fmt.Fprintf(&f.b, `<text x="%.1f" y="%.1f" text-anchor="%s" fill="%s" font-size="%g" font-weight="%d">%s</text>`+"\n",
		x, y, anchor, color, size, weight, html.EscapeString(s))
}

// Raw appends SVG markup.
func (f *Fig) Raw(s string) { f.b.WriteString(s + "\n") }

// Save closes the figure and writes it.
func (f *Fig) Save(path string) error {
	f.b.WriteString("</svg>\n")
	return os.WriteFile(path, []byte(f.b.String()), 0o644)
}

// Axes is a plot area with linear or log scales.
type Axes struct {
	f                      *Fig
	X0, Y0, W, H           float64
	XMin, XMax, YMin, YMax float64
	XLog, YLog             bool
}

// Axes adds a plot area at pixel box (x0, y0, w, h).
func (f *Fig) Axes(x0, y0, w, h, xmin, xmax, ymin, ymax float64, xlog, ylog bool) *Axes {
	return &Axes{f, x0, y0, w, h, xmin, xmax, ymin, ymax, xlog, ylog}
}

func scale(v, lo, hi float64, log bool) float64 {
	if log {
		return (math.Log10(v) - math.Log10(lo)) / (math.Log10(hi) - math.Log10(lo))
	}
	return (v - lo) / (hi - lo)
}

// PX and PY map data to pixels.
func (a *Axes) PX(x float64) float64 { return a.X0 + a.W*scale(x, a.XMin, a.XMax, a.XLog) }
func (a *Axes) PY(y float64) float64 { return a.Y0 + a.H*(1-scale(y, a.YMin, a.YMax, a.YLog)) }

// Grid draws recessive grid lines, tick labels and axis titles. A nil tick
// list draws no grid on that axis.
func (a *Axes) Grid(xt, yt []float64, xf, yf func(float64) string, xlabel, ylabel string) {
	for _, x := range xt {
		p := a.PX(x)
		a.f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>`, p, a.Y0, p, a.Y0+a.H, colGrid))
		a.f.Text(p, a.Y0+a.H+18, xf(x), "middle", colInk2, 12, 400)
	}
	for _, y := range yt {
		p := a.PY(y)
		a.f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>`, a.X0, p, a.X0+a.W, p, colGrid))
		a.f.Text(a.X0-8, p+4, yf(y), "end", colInk2, 12, 400)
	}
	a.f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>`, a.X0, a.Y0+a.H, a.X0+a.W, a.Y0+a.H, colMuted))
	if xlabel != "" {
		a.f.Text(a.X0+a.W/2, a.Y0+a.H+42, xlabel, "middle", colInk2, 13, 400)
	}
	if ylabel != "" {
		x, y := a.X0-58, a.Y0+a.H/2
		a.f.Raw(fmt.Sprintf(`<text x="%.1f" y="%.1f" text-anchor="middle" fill="%s" font-size="13" transform="rotate(-90 %.1f %.1f)">%s</text>`,
			x, y, colInk2, x, y, html.EscapeString(ylabel)))
	}
}

// Line draws a 2 px polyline, clipped to the y range.
func (a *Axes) Line(xs, ys []float64, color string, width float64) {
	var pts []string
	flush := func() {
		if len(pts) > 1 {
			a.f.Raw(fmt.Sprintf(`<polyline fill="none" stroke="%s" stroke-width="%g" stroke-linejoin="round" points="%s"/>`, color, width, strings.Join(pts, " ")))
		}
		pts = pts[:0]
	}
	for i := range xs {
		y := math.Max(a.YMin, math.Min(a.YMax, ys[i]))
		if xs[i] < a.XMin || xs[i] > a.XMax || math.IsNaN(y) {
			flush()
			continue
		}
		pts = append(pts, fmt.Sprintf("%.1f,%.1f", a.PX(xs[i]), a.PY(y)))
	}
	flush()
}

// Dot draws a marker with a 2 px surface ring.
func (a *Axes) Dot(x, y, r float64, color string) {
	a.f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" stroke="%s" stroke-width="2"/>`, a.PX(x), a.PY(y), r, color, colSurface))
}

// VLine draws a vertical rule at x, optionally dashed.
func (a *Axes) VLine(x float64, color string, width float64, dash bool) {
	d := ""
	if dash {
		d = ` stroke-dasharray="4 4"`
	}
	p := a.PX(x)
	a.f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="%g"%s/>`, p, a.Y0, p, a.Y0+a.H, color, width, d))
}

// Legend draws a swatch row at pixel (x, y).
func (f *Fig) Legend(x, y float64, names, colors []string) {
	for i, n := range names {
		f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="14" height="4" rx="2" fill="%s"/>`, x, y-4, colors[i]))
		f.Text(x+20, y+1, n, "start", colInk, 13, 400)
		x += 28 + 7.2*float64(len([]rune(n)))
	}
}

// Ticks returns 1, 2, 5 × 10ⁿ ticks (log) or evenly spaced ones (linear).
func Ticks(lo, hi float64, step float64) []float64 {
	var t []float64
	for v := math.Ceil(lo/step) * step; v <= hi+1e-9*step; v += step {
		t = append(t, v)
	}
	return t
}

// Decades returns the powers of ten in [lo, hi], every k-th.
func Decades(lo, hi float64, k int) []float64 {
	var t []float64
	for e := int(math.Ceil(math.Log10(lo))); float64(e) <= math.Log10(hi)+1e-9; e += k {
		t = append(t, math.Pow(10, float64(e)))
	}
	return t
}

// Pow10Label writes 10ⁿ with a superscript exponent.
func Pow10Label(v float64) string {
	e := int(math.Round(math.Log10(v)))
	sup := map[rune]rune{'-': '⁻', '0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴', '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹'}
	var s strings.Builder
	s.WriteString("10")
	for _, r := range fmt.Sprint(e) {
		s.WriteRune(sup[r])
	}
	return s.String()
}
