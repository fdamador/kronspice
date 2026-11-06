// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	_ "embed"
	"encoding/json"
	"fmt"
	"math"
	"math/rand/v2"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// The atom in 3-D, at its two scales. Left: the atom (Bohr radii) — the 1s
// cloud sampled from |ψ|², the 32 spin-current loops of the 21 cm circuit,
// the proton at the centre. Right: the particles (femtometres) — the
// electron's light on the (2, 1) double loop of photon_electron (twice
// round the ring of radius λ_C/4π while once round the tube of radius
// R/2) and the proton light ring of nucleon_circuits, at true relative
// size. The same geometry is written as JSON for the interactive view.

// view projects a point: azimuth az about z, then elevation el about the
// screen's x axis. Returns screen x, y (up) and depth (toward the viewer).
func view(p [3]float64, az, el float64) (x, y, d float64) {
	ca, sa := math.Cos(az), math.Sin(az)
	x1 := ca*p[0] - sa*p[1]
	y1 := sa*p[0] + ca*p[1]
	ce, se := math.Cos(el), math.Sin(el)
	return x1, ce*p[2] - se*y1, se*p[2] + ce*y1
}

// CloudPoints samples n points of the 1s density (radius in units of a):
// r ~ r²e^(−2r) (a sum of three exponentials over 2), direction uniform.
// The seed is fixed so the picture is the same on every run.
func CloudPoints(n int) [][3]float64 {
	rng := rand.New(rand.NewPCG(1, 32))
	pts := make([][3]float64, n)
	for i := range pts {
		r := -(math.Log(1-rng.Float64()) + math.Log(1-rng.Float64()) + math.Log(1-rng.Float64())) / 2
		c := 2*rng.Float64() - 1
		ph := 2 * math.Pi * rng.Float64()
		s := math.Sqrt(1 - c*c)
		pts[i] = [3]float64{r * s * math.Cos(ph), r * s * math.Sin(ph), r * c}
	}
	return pts
}

// DoubleLoop is the (2, 1) torus knot of radius R, tube R/2, n points.
func DoubleLoop(R float64, n int) [][3]float64 {
	pts := make([][3]float64, n+1)
	for i := range pts {
		u := 4 * math.Pi * float64(i) / float64(n)
		th := u / 2
		rr := R + R/2*math.Cos(th)
		pts[i] = [3]float64{rr * math.Cos(u), rr * math.Sin(u), R / 2 * math.Sin(th)}
	}
	return pts
}

// circle3 is a horizontal circle of radius rho at height z.
func circle3(rho, z float64, n int) [][3]float64 {
	pts := make([][3]float64, n+1)
	for i := range pts {
		t := 2 * math.Pi * float64(i) / float64(n)
		pts[i] = [3]float64{rho * math.Cos(t), rho * math.Sin(t), z}
	}
	return pts
}

func figAtom3D(dir string, _, _ Spectrum) error {
	const az, el = 0.55, 0.38
	f := NewFig(1000, 600, "Fig. 7 — Light hydrogen in 3-D, at the atom's scale and at the particles' scale",
		"left: 1s cloud (sampled from |ψ|²) with the 21 cm spin-current loops round the proton; right: the electron's double loop of light and the proton ring, true relative size")
	// Left panel: the atom.
	cx, cy, s := 260.0, 330.0, 105.0 // pixels per a0
	type dot struct{ x, y, d float64 }
	var dots []dot
	for _, p := range CloudPoints(2500) {
		x, y, d := view(p, az, el)
		if math.Abs(x) < 2.3 && math.Abs(y) < 2.3 {
			dots = append(dots, dot{x, y, d})
		}
	}
	sort.Slice(dots, func(i, j int) bool { return dots[i].d < dots[j].d })
	for _, q := range dots {
		op := 0.25 + 0.35*(q.d+2)/4
		f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="1.6" fill="%s" fill-opacity="%.2f"/>`, cx+s*q.x, cy-s*q.y, colBlue, math.Max(0.15, math.Min(0.6, op))))
	}
	a0 := A0 * Recoil
	loops := SpinLoops(Measured().Moment(), a0, 8)
	imax := 0.0
	for _, l := range loops {
		imax = math.Max(imax, l.I)
	}
	for _, l := range loops {
		pts := ""
		for _, p := range circle3(l.Rho/a0, l.Z/a0, 72) {
			x, y, _ := view(p, az, el)
			pts += fmt.Sprintf("%.1f,%.1f ", cx+s*x, cy-s*y)
		}
		f.Raw(fmt.Sprintf(`<polyline fill="none" stroke="%s" stroke-width="%.1f" stroke-opacity="0.85" points="%s"/>`, colOrange, 0.6+2.2*l.I/imax, pts))
	}
	// the axis and the proton
	x1, y1, _ := view([3]float64{0, 0, -2.1}, az, el)
	x2, y2, _ := view([3]float64{0, 0, 2.1}, az, el)
	f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" stroke-dasharray="3 3"/>`, cx+s*x1, cy-s*y1, cx+s*x2, cy-s*y2, colInk2))
	f.Text(cx+s*x2+6, cy-s*y2, "z (spin axis)", "start", colInk2, 12, 400)
	f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="4" fill="%s" stroke="%s" stroke-width="2"/>`, cx, cy, colInk, colSurface))
	f.Text(cx+10, cy+18, "proton", "start", colInk, 12, 600)
	f.Text(cx-s*2.2, 560, "scale bar:", "start", colInk2, 12, 400)
	f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="555" x2="%.1f" y2="555" stroke="%s" stroke-width="2"/>`, cx-s*2.2+70, cx-s*2.2+70+s, colInk))
	f.Text(cx-s*2.2+70+s+8, 560, "1 a0 = 52,918 fm", "start", colInk2, 12, 400)
	f.Legend(cx-s*2.2, 585, []string{"1s cloud", "spin-current loops (width ∝ current)"}, []string{colBlue, colOrange})

	// Right panel: the particles, in fm.
	px, py, ps := 760.0, 300.0, 0.62 // pixels per fm
	wv := WvdM()
	pts := ""
	for _, p := range DoubleLoop(wv.R*1e15, 400) {
		x, y, _ := view(p, az, 1.05)
		pts += fmt.Sprintf("%.1f,%.1f ", px+ps*x, py-ps*y)
	}
	f.Raw(fmt.Sprintf(`<polyline fill="none" stroke="%s" stroke-width="2.5" stroke-linejoin="round" points="%s"/>`, colOrange, pts))
	ring := ""
	for _, p := range circle3(ProtonRingRadius()*1e15, 0, 48) {
		x, y, _ := view(p, az, 1.05)
		ring += fmt.Sprintf("%.2f,%.2f ", px+ps*x, py-ps*y)
	}
	f.Raw(fmt.Sprintf(`<polyline fill="none" stroke="%s" stroke-width="2" points="%s"/>`, colBlue, ring))
	f.Raw(fmt.Sprintf(`<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>`, px, py, colBlue))
	f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>`, px, py+4, px, py+168, colInk2))
	f.Text(px, py+184, "proton ring, 0.21 fm (a dot at this scale)", "middle", colInk, 12, 400)
	f.Text(px, 112, "electron: light twice round a ring of 193 fm (H23)", "middle", colInk, 13, 600)
	f.Raw(fmt.Sprintf(`<line x1="%.1f" y1="555" x2="%.1f" y2="555" stroke="%s" stroke-width="2"/>`, px-60, px-60+100*ps, colInk))
	f.Text(px-60+100*ps+8, 560, "100 fm", "start", colInk2, 12, 400)
	f.Text(px, 590, "the whole right panel is 0.006 a0 wide: a point in the left panel", "middle", colInk2, 12, 400)
	f.Raw(fmt.Sprintf(`<line x1="520" y1="90" x2="520" y2="570" stroke="%s" stroke-width="1"/>`, colGrid))
	return f.Save(filepath.Join(dir, "fig7_atom3d.svg"))
}

// Scene3D is the geometry of the interactive view.
type Scene3D struct {
	Cloud      [][3]float64 `json:"cloud"`      // a0
	Cloud2p    [][3]float64 `json:"cloud2p"`    // a0, the 2p_z state
	Loops      [][3]float64 `json:"loops"`      // ρ, z [a0], relative current
	DoubleLoop [][3]float64 `json:"doubleLoop"` // fm
	ElectronR  float64      `json:"electronRingFm"`
	ProtonR    float64      `json:"protonRingFm"`
	BohrFm     float64      `json:"bohrFm"`
}

// atom3dTemplate is the interactive 3-D page (three.js from cdnjs); its
// __SCENE__ is replaced by the scene of the run.
//
//go:embed atom3d.html.tmpl
var atom3dTemplate string

// WriteScene writes the interactive view's data, and the interactive page
// atom3d.html beside it with the data inlined (a browser opening a local
// file cannot fetch the JSON).
func WriteScene(path string) error {
	a0 := A0 * Recoil
	sc := Scene3D{Cloud: CloudPoints(6000), Cloud2p: CloudPoints2p(6000), DoubleLoop: DoubleLoop(WvdM().R*1e15, 600),
		ElectronR: WvdM().R * 1e15, ProtonR: ProtonRingRadius() * 1e15, BohrFm: A0 * 1e15}
	loops := SpinLoops(Measured().Moment(), a0, 8)
	imax := 0.0
	for _, l := range loops {
		imax = math.Max(imax, l.I)
	}
	for _, l := range loops {
		sc.Loops = append(sc.Loops, [3]float64{l.Rho / a0, l.Z / a0, l.I / imax})
	}
	for _, cl := range [][][3]float64{sc.Cloud, sc.Cloud2p} {
		for i := range cl {
			for k := range cl[i] {
				cl[i][k] = math.Round(cl[i][k]*1e3) / 1e3
			}
		}
	}
	for i := range sc.DoubleLoop {
		for k := range sc.DoubleLoop[i] {
			sc.DoubleLoop[i][k] = math.Round(sc.DoubleLoop[i][k]*1e3) / 1e3
		}
	}
	b, err := json.Marshal(sc)
	if err != nil {
		return err
	}
	if err := os.WriteFile(path, b, 0o644); err != nil {
		return err
	}
	page := strings.Replace(atom3dTemplate, "__SCENE__", string(b), 1)
	return os.WriteFile(filepath.Join(filepath.Dir(path), "atom3d.html"), []byte(page), 0o644)
}
