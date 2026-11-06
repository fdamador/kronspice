// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"path/filepath"
	"sort"
	"strings"
)

// Pole colours for the ring figures (validated palette, light mode).
var (
	rgbN    = [3]float64{235, 104, 52} // N: orange
	rgbS    = [3]float64{42, 120, 214} // S: blue
	rgbLine = [3]float64{38, 37, 35}   // field lines: ink
)

// Camera is an orthographic view: azimuth Az about z, elevation El above
// the mid-plane.
type Camera struct{ Az, El float64 }

// Project maps a world point to screen x (right), y (up) and depth (larger
// is farther from the camera).
func (c Camera) Project(p [3]float64) (x, y, d float64) {
	ca, sa := math.Cos(c.Az), math.Sin(c.Az)
	x1 := p[0]*ca + p[1]*sa
	y1 := -p[0]*sa + p[1]*ca
	ce, se := math.Cos(c.El), math.Sin(c.El)
	return x1, p[2]*ce + y1*se, y1*ce - p[2]*se
}

// Facing is n̂ · (direction to the camera).
func (c Camera) Facing(n [3]float64) float64 {
	ca, sa := math.Cos(c.Az), math.Sin(c.Az)
	y1 := -n[0]*sa + n[1]*ca
	return -y1*math.Cos(c.El) + n[2]*math.Sin(c.El)
}

// Prim is one drawable piece of the 3D scene: a face (polygon with an
// outward normal) or a field-line segment.
type Prim struct {
	P     [][3]float64
	N     [3]float64
	RGB   [3]float64
	Line  bool
	Width float64
}

// poleRGB is the colour of the ring's surface at azimuth phi in section k
// (on the outer face; the inner face has the opposite pole for radial
// sections): radial sections show their outer pole, tangential sections are
// N in the half toward the head of M and S in the other half.
func (r Ring) poleRGB(k int, phi float64, inner bool) [3]float64 {
	n := r.Sign[k] > 0
	if !r.Radial {
		p0, p1 := r.Span(k)
		n = (phi > 0.5*(p0+p1)) == (r.Sign[k] > 0)
		inner = false
	}
	if n != inner {
		return rgbN
	}
	return rgbS
}

func tint(c [3]float64, t float64) [3]float64 {
	return [3]float64{c[0] + (255-c[0])*t, c[1] + (255-c[1])*t, c[2] + (255-c[2])*t}
}

// Scene builds the ring's faces and its field lines.
func (r Ring) Scene(lines [][][3]float64) []Prim {
	var ps []Prim
	const ns = 8 // steps per section
	h := r.H / 2
	for k := 0; k < r.N(); k++ {
		p0, p1 := r.Span(k)
		for i := 0; i < ns; i++ {
			a := p0 + (p1-p0)*float64(i)/ns
			b := p0 + (p1-p0)*float64(i+1)/ns
			m := 0.5 * (a + b)
			at := func(rho, phi, z float64) [3]float64 { return [3]float64{rho * math.Cos(phi), rho * math.Sin(phi), z} }
			co, cin := r.poleRGB(k, m, false), r.poleRGB(k, m, true)
			for _, zz := range []float64{h, -h} {
				for j := 0; j < 2; j++ { // two radial strips so depth sorting stays local
					r0 := r.Ri + (r.Ro-r.Ri)*float64(j)/2
					r1 := r.Ri + (r.Ro-r.Ri)*float64(j+1)/2
					ps = append(ps, Prim{P: [][3]float64{at(r0, a, zz), at(r1, a, zz), at(r1, b, zz), at(r0, b, zz)},
						N: [3]float64{0, 0, math.Copysign(1, zz)}, RGB: tint(co, 0.45)})
				}
			}
			ps = append(ps, Prim{P: [][3]float64{at(r.Ro, a, -h), at(r.Ro, b, -h), at(r.Ro, b, h), at(r.Ro, a, h)},
				N: [3]float64{math.Cos(m), math.Sin(m), 0}, RGB: co})
			ps = append(ps, Prim{P: [][3]float64{at(r.Ri, a, -h), at(r.Ri, b, -h), at(r.Ri, b, h), at(r.Ri, a, h)},
				N: [3]float64{-math.Cos(m), -math.Sin(m), 0}, RGB: cin})
		}
	}
	for _, l := range lines {
		// every third traced point (0.9 mm steps) is enough to draw
		for i := 0; i+1 < len(l); i += 3 {
			j := min(i+3, len(l)-1)
			ps = append(ps, Prim{P: [][3]float64{l[i], l[j]}, RGB: rgbLine, Line: true, Width: 1.1})
		}
	}
	return ps
}

// painted is a projected primitive ready to draw, back to front.
type painted struct {
	XY    [][2]float64
	RGB   [3]float64
	Line  bool
	Width float64
	Depth float64
}

// Paint projects, culls back faces, shades and sorts the scene. scale and
// (cx, cy) map metres to pixels (y down).
func (c Camera) Paint(ps []Prim, scale, cx, cy float64) []painted {
	light := [3]float64{-0.35, -0.55, 0.76} // world direction toward the light
	var out []painted
	for _, p := range ps {
		var q painted
		q.Line, q.Width, q.RGB = p.Line, p.Width, p.RGB
		if !p.Line {
			if c.Facing(p.N) <= 1e-9 {
				continue
			}
			l := math.Max(0, p.N[0]*light[0]+p.N[1]*light[1]+p.N[2]*light[2])
			f := 0.72 + 0.28*l
			q.RGB = [3]float64{p.RGB[0] * f, p.RGB[1] * f, p.RGB[2] * f}
		}
		for _, v := range p.P {
			x, y, d := c.Project(v)
			q.XY = append(q.XY, [2]float64{cx + scale*x, cy - scale*y})
			q.Depth += d / float64(len(p.P))
		}
		if p.Line {
			q.Depth -= 1e-4 // a line lying on a face is drawn over it
		}
		out = append(out, q)
	}
	sort.SliceStable(out, func(i, j int) bool { return out[i].Depth > out[j].Depth })
	return out
}

func hex(c [3]float64) string {
	cl := func(v float64) int { return int(math.Max(0, math.Min(255, math.Round(v)))) }
	return fmt.Sprintf("#%02x%02x%02x", cl(c[0]), cl(c[1]), cl(c[2]))
}

// RingView is everything one 3D panel needs.
type RingView struct {
	Ring  Ring
	Lines [][][3]float64
	PoleN []float64 // azimuths of the N poles seen from outside
	PoleS []float64
}

// NewRingView traces the 3D field lines: from the outer N faces at three
// heights, four lines per pole at equal steps of flux (an even number, so
// none starts on a pole's centre line, which runs straight out to infinity).
func NewRingView(r Ring, cs []Charge) RingView {
	var seeds [][3]float64
	for _, z := range []float64{-0.35 * r.H, 0, 0.35 * r.H} {
		seeds = append(seeds, Seeds(cs, r.Ro+0.4e-3, z, 1, 4)...)
	}
	v := RingView{Ring: r, Lines: r.Lines(cs, seeds, 0.3e-3)}
	phi, br, _, _ := RingScan(cs, r.Ro+3e-3, 0, 720)
	v.PoleN = PeakAngles(phi, br)
	neg := make([]float64, len(br))
	for i := range br {
		neg[i] = -br[i]
	}
	v.PoleS = PeakAngles(phi, neg)
	return v
}

// svg3D draws one 3D panel into f.
func svg3D(f *Fig, v RingView, cam Camera, x0, y0, w, hgt float64) {
	scale := 0.27 * w / v.Ring.Ro
	cx, cy := x0+w/2, y0+hgt/2
	id := fmt.Sprintf("p3d%.0f", x0)
	f.Raw(fmt.Sprintf(`<clipPath id="%s"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath><g clip-path="url(#%s)">`, id, x0, y0+50, w, hgt-50, id))
	var sb strings.Builder
	for _, q := range cam.Paint(v.Ring.Scene(v.Lines), scale, cx, cy) {
		if q.Line {
			fmt.Fprintf(&sb, `<line class="fl" x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f"/>`+"\n", q.XY[0][0], q.XY[0][1], q.XY[1][0], q.XY[1][1])
			continue
		}
		var pts []string
		for _, p := range q.XY {
			pts = append(pts, fmt.Sprintf("%.1f,%.1f", p[0], p[1]))
		}
		fmt.Fprintf(&sb, `<polygon points="%s" fill="%[2]s" stroke="%[2]s" stroke-width="0.6" stroke-linejoin="round"/>`+"\n", strings.Join(pts, " "), hex(q.RGB))
	}
	f.Raw(sb.String() + "</g>")
	// pole labels on the visible outer side
	for _, set := range []struct {
		phis []float64
		lab  string
		rgb  [3]float64
	}{{v.PoleN, "N", rgbN}, {v.PoleS, "S", rgbS}} {
		for _, p := range set.phis {
			n := [3]float64{math.Cos(p), math.Sin(p), 0}
			if cam.Facing(n) < 0.2 {
				continue
			}
			rr := v.Ring.Ro + 0.5e-3
			x, y, _ := cam.Project([3]float64{rr * math.Cos(p), rr * math.Sin(p), 0})
			f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="8" fill="%s" fill-opacity="0.9" stroke="%s" stroke-width="1.5"/>`, cx+scale*x, cy-scale*y, colSurface, hex(set.rgb)))
			f.Text(cx+scale*x, cy-scale*y+4, set.lab, "middle", hex(set.rgb), 11, 700)
		}
	}
}

// RingFigures writes Figs. R1–R3 and the rotating GIFs.
func RingFigures(dir string, rad, tan, ctl Ring) error {
	cr, ct, cc := rad.Charges(0.5e-3, 1), tan.Charges(0.5e-3, 1), ctl.Charges(0.5e-3, 1)
	vr, vt := NewRingView(rad, cr), NewRingView(tan, ct)

	// Fig. R1: 3D
	f := NewFig(1200, 660, "Fig. R1 — A ring magnet of 12 sections, each magnetized opposite to its neighbours, and its field outside",
		"field lines traced in 3D from the outer N faces (four per pole at three heights, at equal steps of flux) until they return into an S face; orange N, blue S")
	f.Raw(fmt.Sprintf(`<style>.fl{stroke:%s;stroke-opacity:0.8;stroke-width:1.1;stroke-linecap:round}</style>`, hex(rgbLine)))
	cam := Camera{Az: -0.9, El: 0.5}
	svg3D(f, vr, cam, 20, 80, 580, 520)
	svg3D(f, vt, cam, 600, 80, 580, 520)
	f.Text(310, 96, "Radial: outward and inward in turn", "middle", colInk, 15, 600)
	f.Text(310, 116, "N and S alternate face by face; poles at the section centres", "middle", colInk2, 12, 400)
	f.Text(890, 96, "Along the ring: like poles meet at every joint", "middle", colInk, 15, 600)
	f.Text(890, 116, "the octave bar's S N | N S closed into a ring; poles at the joints", "middle", colInk2, 12, 400)
	f.Text(24, 630, "Ring 50 mm across, 30 mm bore, 10 mm tall, Br 1.25 T (declared). Each line leaves an N face, arcs over the gap, and enters the neighbouring S face:", "start", colInk2, 12, 400)
	f.Text(24, 648, "twelve arches round the outside, bulging above and below the mid-plane. Field from the surface charge σ = M·n̂ (checked against a point-dipole twin).", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figR1_ring3d.svg")); err != nil {
		return err
	}

	// Fig. R2: mid-plane
	f = NewFig(1200, 520, "Fig. R2 — The field in the ring's mid-plane (top view)",
		"by symmetry B lies in the mid-plane, so these lines are exact; lines from the outer N faces, arrows N → S; arrows in the sections are M")
	for i, c := range []struct {
		r     Ring
		cs    []Charge
		title string
		k     int
	}{{rad, cr, "Radial, alternating: 12 poles outside", 8}, {tan, ct, "Along the ring, alternating: 12 poles at the joints", 8}, {ctl, cc, "Control, all outward: no poles round the side", 72}} {
		midPlane(f, c.r, c.cs, 20+float64(i)*393, 80, 380, c.title, c.k)
	}
	f.Text(24, 506, "Outside lines leave the plane only in the control, where the flux of the all-N outer face returns over the top and bottom into the all-S bore.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figR2_midplane.svg")); err != nil {
		return err
	}

	// Fig. R3: the scan
	f = NewFig(1000, 500, "Fig. R3 — The bench test: B_r round the ring, 3 mm outside, in the mid-plane",
		"a Hall probe carried round the ring (or the ring turned under it); sign changes count the poles")
	var series [][]float64
	var phis []float64
	for _, cs := range [][]Charge{cr, ct, cc} {
		phi, br, _, _ := RingScan(cs, rad.Ro+3e-3, 0, 720)
		phis = phi
		mt := make([]float64, len(br))
		for i := range br {
			mt[i] = br[i] * 1e3
		}
		series = append(series, mt)
	}
	deg := make([]float64, len(phis))
	for i := range phis {
		deg[i] = phis[i] * 180 / math.Pi
	}
	a := f.Axes(90, 80, 860, 320, 0, 360, -400, 400, false, false)
	a.Grid(Ticks(0, 360, 30), Ticks(-400, 400, 200), func(v float64) string { return fmt.Sprintf("%g°", v) },
		func(v float64) string { return fmt.Sprintf("%g", v) }, "angle round the ring (section joints every 30°)", "B_r (mT)")
	a.Line([]float64{0, 360}, []float64{0, 0}, colInk2, 1)
	cols := []string{colBlue, colOrange, colMuted}
	for i, s := range series {
		a.Line(deg, s, cols[i], 2)
	}
	f.Legend(90, 462, []string{
		fmt.Sprintf("radial, alternating: %d sign changes", CyclicSignChanges(series[0])),
		fmt.Sprintf("along the ring, alternating: %d", CyclicSignChanges(series[1])),
		fmt.Sprintf("all outward: %d", CyclicSignChanges(series[2]))}, cols)
	f.Text(90, 488, "Twelve sign changes for both alternating rings, peaks shifted by half a section (15°): radial poles sit mid-section, along-the-ring poles at the joints.", "start", colInk2, 12, 400)
	if err := f.Save(filepath.Join(dir, "figR3_scan.svg")); err != nil {
		return err
	}

	// rotating views
	if err := RingGIF(filepath.Join(dir, "ring3d_radial.gif"), vr); err != nil {
		return err
	}
	return RingGIF(filepath.Join(dir, "ring3d_tangential.gif"), vt)
}

// midPlane draws one top-view panel: the sections with their M arrows and
// the field lines in z = 0, from the outer N faces.
func midPlane(f *Fig, r Ring, cs []Charge, x0, y0, side float64, title string, k int) {
	ext := 2.2 * r.Ro
	scale := side / (2 * ext)
	cx, cy := x0+side/2, y0+side/2+16
	px := func(p [3]float64) (float64, float64) { return cx + scale*p[0], cy - scale*p[1] }
	f.Text(x0+side/2, y0+4, title, "middle", colInk, 13, 600)
	f.Raw(fmt.Sprintf(`<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="%s"/>`, cx-side/2, cy-side/2, side, side, colGrid))
	f.Raw(fmt.Sprintf(`<clipPath id="mp%.0f"><rect x="%.1f" y="%.1f" width="%.1f" height="%.1f"/></clipPath>`, x0, cx-side/2, cy-side/2, side, side))
	seeds := Seeds(cs, r.Ro+0.4e-3, 0, 1, k)
	lines := r.Lines(cs, seeds, 0.25e-3)
	var sb strings.Builder
	for _, l := range lines {
		sb.WriteString("M")
		for i, p := range l {
			x, y := px(p)
			if i > 0 {
				sb.WriteString("L")
			}
			fmt.Fprintf(&sb, "%.1f %.1f", x, y)
		}
	}
	f.Raw(fmt.Sprintf(`<g clip-path="url(#mp%.0f)"><path d="%s" fill="none" stroke="%s" stroke-width="1" stroke-opacity="0.85"/>`, x0, sb.String(), hex(rgbLine)))
	// arrowheads at each line's midpoint
	for _, l := range lines {
		if len(l) < 8 {
			continue
		}
		m := len(l) / 2
		x1, y1 := px(l[m-1])
		x2, y2 := px(l[m+1])
		arrow(f, x2, y2, math.Atan2(y2-y1, x2-x1), 5, hex(rgbLine))
	}
	f.Raw("</g>")
	// sections
	for s := 0; s < r.N(); s++ {
		p0, p1 := r.Span(s)
		const n = 8
		for i := 0; i < n; i++ {
			a := p0 + (p1-p0)*float64(i)/n
			b := p0 + (p1-p0)*float64(i+1)/n
			c := r.poleRGB(s, 0.5*(a+b), false)
			pt := func(rho, phi float64) string {
				x, y := px([3]float64{rho * math.Cos(phi), rho * math.Sin(phi), 0})
				return fmt.Sprintf("%.1f,%.1f", x, y)
			}
			f.Raw(fmt.Sprintf(`<polygon points="%s %s %s %s" fill="%[5]s" stroke="%[5]s" stroke-width="0.5"/>`,
				pt(r.Ri, a), pt(r.Ro, a), pt(r.Ro, b), pt(r.Ri, b), hex(tint(c, 0.25))))
		}
		// joint line
		xa, ya := px([3]float64{r.Ri * math.Cos(p0), r.Ri * math.Sin(p0), 0})
		xb, yb := px([3]float64{r.Ro * math.Cos(p0), r.Ro * math.Sin(p0), 0})
		f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>`, xa, ya, xb, yb, colSurface))
		// M arrow
		M := r.Mag(s)
		mid := 0.5 * (p0 + p1)
		rc := 0.5 * (r.Ri + r.Ro)
		ux, uy := M[0]/r.M0, M[1]/r.M0
		l := 0.32 * (r.Ro - r.Ri)
		if !r.Radial {
			l = 0.5 * rc * (p1 - p0) * 0.6
		}
		ax, ay := px([3]float64{rc*math.Cos(mid) - l*ux, rc*math.Sin(mid) - l*uy, 0})
		bx, by := px([3]float64{rc*math.Cos(mid) + l*ux, rc*math.Sin(mid) + l*uy, 0})
		f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.6"/>`, ax, ay, bx, by, colInk))
		arrow(f, bx, by, math.Atan2(by-ay, bx-ax), 5, colInk)
	}
}

func arrow(f *Fig, x, y, ang, s float64, col string) {
	p := func(da, r float64) string {
		return fmt.Sprintf("%.1f,%.1f", x+r*math.Cos(ang+da), y+r*math.Sin(ang+da))
	}
	f.Raw(fmt.Sprintf(`<polygon points="%.1f,%.1f %s %s" fill="%s"/>`, x+s*math.Cos(ang), y+s*math.Sin(ang), p(2.5, s), p(-2.5, s), col))
}
