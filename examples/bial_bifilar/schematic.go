// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"strings"
)

// svg accumulates a schematic drawing. Strokes and text use currentColor, so
// the drawing follows the page's colour scheme when embedded inline; the
// standalone file sets a white background and black ink.
type svg struct{ b strings.Builder }

func (s *svg) f(format string, a ...any) { fmt.Fprintf(&s.b, format+"\n", a...) }

func (s *svg) line(x1, y1, x2, y2 float64) {
	s.f(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>`, x1, y1, x2, y2)
}

func (s *svg) dot(x, y float64) {
	s.f(`<circle cx="%.1f" cy="%.1f" r="3.5" class="dot"/>`, x, y)
}

// label writes a component's name, value and an optional detail line
// beside (x, y).
func (s *svg) label(x, y float64, name, value, detail string) {
	dy := 0.0
	if detail != "" {
		dy = -7
	}
	s.f(`<text x="%.1f" y="%.1f" class="name">%s</text>`, x, y-4+dy, name)
	s.f(`<text x="%.1f" y="%.1f" class="value">%s</text>`, x, y+12+dy, value)
	if detail != "" {
		s.f(`<text x="%.1f" y="%.1f" class="detail">%s</text>`, x, y+26+dy, detail)
	}
}

const compLen = 56 // length of a drawn component along its wire [px]

// resistor draws a zigzag from y0 down to y0+compLen at x, with leads.
func (s *svg) resistor(x, y0 float64) {
	pts := []string{fmt.Sprintf("%.1f,%.1f", x, y0)}
	n := 6
	for i := 1; i <= n; i++ {
		dx := 9.0
		if i%2 == 0 {
			dx = -9
		}
		y := y0 + 8 + float64(i-1)*(compLen-16)/float64(n-1)
		pts = append(pts, fmt.Sprintf("%.1f,%.1f", x+dx, y))
	}
	pts = append(pts, fmt.Sprintf("%.1f,%.1f", x, y0+compLen))
	s.f(`<polyline points="%s" class="part"/>`, strings.Join(pts, " "))
}

func (s *svg) capacitor(x, y0 float64) {
	yc := y0 + compLen/2
	s.line(x, y0, x, yc-5)
	s.line(x, yc+5, x, y0+compLen)
	s.f(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="plate"/>`, x-16, yc-5, x+16, yc-5)
	s.f(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" class="plate"/>`, x-16, yc+5, x+16, yc+5)
}

// inductor draws four bumps to the right of the wire.
func (s *svg) inductor(x, y0 float64) {
	h := float64(compLen) / 4
	d := fmt.Sprintf("M %.1f %.1f", x, y0)
	for i := 0; i < 4; i++ {
		d += fmt.Sprintf(" a %.1f %.1f 0 0 1 0 %.1f", h/2, h/2, h)
	}
	s.f(`<path d="%s" class="part" fill="none"/>`, d)
}

// vsource draws a voltage source with + at the top.
func (s *svg) vsource(x, y0 float64) {
	yc := y0 + compLen/2
	s.line(x, y0, x, yc-18)
	s.line(x, yc+18, x, y0+compLen)
	s.f(`<circle cx="%.1f" cy="%.1f" r="18" class="part" fill="none"/>`, x, yc)
	s.f(`<text x="%.1f" y="%.1f" class="sign">+</text>`, x-4, yc-4)
	s.f(`<text x="%.1f" y="%.1f" class="sign">−</text>`, x-4, yc+13)
}

// isource draws a current source pushing current upwards.
func (s *svg) isource(x, y0 float64) {
	yc := y0 + compLen/2
	s.line(x, y0, x, yc-18)
	s.line(x, yc+18, x, y0+compLen)
	s.f(`<circle cx="%.1f" cy="%.1f" r="18" class="part" fill="none"/>`, x, yc)
	s.line(x, yc+10, x, yc-10)
	s.f(`<polyline points="%.1f,%.1f %.1f,%.1f %.1f,%.1f" class="part" fill="none"/>`, x-5, yc-4, x, yc-11, x+5, yc-4)
}

func (s *svg) ground(x, y float64) {
	s.line(x, y, x, y+12)
	for i, w := range []float64{18, 11, 4} {
		yy := y + 12 + float64(i)*5
		s.line(x-w, yy, x+w, yy)
	}
}

// column draws a series string of parts between the top rail and the
// bottom rail at x, spreading them evenly, with labels.
type part struct {
	kind                string // "R", "C", "L", "V", "I"
	name, value, detail string
}

func (s *svg) column(x, top, bottom float64, parts []part, labelDX float64) {
	gap := (bottom - top - float64(len(parts))*compLen) / float64(len(parts)+1)
	y := top
	for _, p := range parts {
		s.line(x, y, x, y+gap)
		y += gap
		switch p.kind {
		case "R":
			s.resistor(x, y)
		case "C":
			s.capacitor(x, y)
		case "L":
			s.inductor(x, y)
		case "V":
			s.vsource(x, y)
		case "I":
			s.isource(x, y)
		}
		s.label(x+labelDX, y+compLen/2, p.name, p.value, p.detail)
		y += compLen
	}
	s.line(x, y, x, bottom)
}

// schematicSVG draws the equivalent circuit of a measured design: the Bi–Al
// element and the coil in parallel between node a (top rail) and ground,
// with the AC probe. Values are the ones the simulation used.
func schematicSVG(r Result, standalone bool) string {
	d, el := r.Design, r.Element
	var s svg
	const w, h = 1000.0, 540.0
	const top, bottom = 100.0, 440.0
	style := `<style>
 .wires line, .wires polyline, .wires path, .wires circle { stroke: currentColor; stroke-width: 2; fill: none; }
 .wires .plate { stroke-width: 3.5; }
 .wires .dot { fill: currentColor; stroke: none; }
 text { fill: currentColor; font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; }
 .name { font-size: 14px; font-weight: 600; }
 .value { font-size: 12.5px; opacity: .85; }
 .detail { font-size: 11.5px; opacity: .6; }
 .sign { font-size: 13px; }
 .box { fill: none; stroke: currentColor; stroke-width: 1.2; stroke-dasharray: 6 5; opacity: .5; }
 .boxlabel { font-size: 14px; font-weight: 600; }
 .note { font-size: 12.5px; opacity: .75; }
 .node { font-size: 14px; font-style: italic; }
</style>`
	if standalone {
		style = strings.Replace(style, "<style>", "<style>\n svg { color: #111; background: #fff; }", 1)
	}
	s.f(`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" width="%.0f" height="%.0f" role="img" aria-label="Equivalent circuit">`, w, h, w, h)
	s.f("%s", style)
	if standalone {
		s.f(`<rect width="100%%" height="100%%" fill="#fff"/>`)
	}
	s.f(`<g class="wires">`)
	// Boxes: element 24..452, coil 492..812, probe at 880.
	const eL, eR, cL, cR = 24.0, 452.0, 492.0, 812.0
	xL, xR := 60.0, 900.0
	s.line(xL, top, xR, top)
	s.line(xL, bottom, xR, bottom)
	xg := (eR + cL) / 2
	s.ground(xg, bottom)
	s.dot(xg, bottom)
	s.f(`<text x="%.1f" y="%.1f" class="node">a</text>`, xL-18, top+5)
	s.f(`<text x="%.1f" y="%.1f" class="node">0</text>`, xg-4, bottom+48)

	f := eng
	// Bi–Al element.
	boxTitle := "Bi–Al element (dry): Seebeck EMF in series with the plates"
	switch d.Variant {
	case Wet:
		boxTitle = "Bi–Al element (wet): Al | electrolyte | Bi cell"
	case Thermopile:
		boxTitle = fmt.Sprintf("Bi–Al thermopile, %d pairs, charging the plates", d.Pairs)
	}
	s.f(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" class="box"/>`, eL, top-58, eR-eL, bottom-top+92)
	s.f(`<text x="%.1f" y="%.1f" class="boxlabel">%s</text>`, eL+12, top-36, boxTitle)
	switch d.Variant {
	case Dry, Thermopile:
		x := 150.0
		rDetail, vDetail := "Bi + Al plates", fmt.Sprintf("Seebeck, ΔT = %g K", d.A.DeltaT)
		if d.Variant == Thermopile {
			rDetail = fmt.Sprintf("%d legs + plates", 2*d.Pairs)
			vDetail = fmt.Sprintf("%d pairs, ΔT = %g K", d.Pairs, d.A.DeltaT)
		}
		s.column(x, top, bottom, []part{
			{"R", "Rint", f(el.Rint, "Ω"), rDetail},
			{"V", "Vbial", f(el.EMF, "V"), vDetail},
			{"C", "Cplate", f(el.Cplate, "F"), fmt.Sprintf("plates, ε = %.3g", d.DielectricEps)},
		}, 28)
		s.dot(x, top)
		s.dot(x, bottom)
		s.f(`<text x="%.1f" y="%.1f" class="note">energy in: thermal port (heat of the hot plates)</text>`, eL+12, bottom+26)
	case Wet:
		x := 140.0
		split := top + (bottom-top)*0.6
		s.column(x, top, split, []part{
			{"R", "Relec", f(el.Rint, "Ω"), "electrolyte"},
			{"V", "Vbial", f(el.EMF, "V"), "cell EMF"},
		}, 28)
		s.dot(x, split)
		xa, xb := x-50, x+50
		s.line(xa, split, xb, split)
		s.column(xa, split, bottom, []part{{"R", "Rct", f(el.Rct, "Ω"), ""}}, 22)
		s.column(xb, split, bottom, []part{{"C", "Cdl", f(el.Cdl, "F"), ""}}, 26)
		s.dot(xa, bottom)
		s.dot(xb, bottom)
		s.dot(x, top)
		xc := 330.0
		s.column(xc, top, bottom, []part{{"C", "Cgeo", f(el.Cgeo, "F"), "plates in water"}}, 26)
		s.dot(xc, top)
		s.dot(xc, bottom)
		s.f(`<text x="%.1f" y="%.1f" class="note">energy in: chemical port (aluminium plate, %.3g Ah)</text>`, eL+12, bottom+26, el.CapacityC/3600)
	}
	// Coil.
	kind := coilKind(d.Coil.Bifilar)
	s.f(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" class="box"/>`, cL, top-58, cR-cL, bottom-top+92)
	s.f(`<text x="%.1f" y="%.1f" class="boxlabel">%s coil, %d turns</text>`, cL+12, top-36, strings.ToUpper(kind[:1])+kind[1:], d.Coil.Turns)
	xs := 540.0
	s.column(xs, top, bottom, []part{{"C", "Cself", f(r.Cself, "F"), "turn-to-turn"}}, 26)
	s.dot(xs, top)
	s.dot(xs, bottom)
	xl := 690.0
	s.column(xl, top, bottom, []part{
		{"R", "Rcoil", f(r.Rcoil, "Ω"), "copper, at f0"},
		{"L", "Lcoil", f(r.L, "H"), ""},
	}, 26)
	s.dot(xl, top)
	s.dot(xl, bottom)
	s.f(`<text x="%.1f" y="%.1f" class="note">%s diameter, %s wire</text>`, cL+12, bottom+26, f(d.Coil.Diameter, "m"), f(d.Coil.WireDiameter, "m"))
	// Probe.
	xp := 870.0
	s.column(xp, top, bottom, []part{{"I", "Iprobe", "AC 1 A", "measures Z"}}, 26)
	s.f(`</g>`)
	s.f(`<text x="24" y="%.1f" class="note">Simulated: f0 = %s, Q = %.4g, stored at steady state %s (%s per m³)</text>`,
		h-22, eng(r.F0, "Hz"), r.Q, eng(r.StoredJ, "J"), eng(r.EnergyDensity, "J"))
	s.f(`</svg>`)
	return s.b.String()
}
