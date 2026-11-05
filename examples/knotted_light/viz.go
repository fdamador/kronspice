// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// 3-D field lines for viz/render.py and the interactive viewer: magnetic
// field lines of the exact hopfion, traced from fixed seeds at several
// times (run with -viz viz). The solver's field agrees with the exact one
// through the energy and helicity (main.go); lines are traced on the exact
// field because a traced line is too sensitive to grid error to follow
// for long (README).

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
)

// Line is one traced field line.
type Line struct {
	Torus  int          `json:"torus"`  // which nest of seeds (0 = nearest the axis)
	Closed bool         `json:"closed"` // came back to its seed
	Points [][3]float64 `json:"points"` // in L0, downsampled
}

// Frame is the knot at one time.
type Frame struct {
	T       float64 `json:"t"`       // in L0/c
	RMS     float64 `json:"rms"`     // rms radius of the energy, in L0
	Linking float64 `json:"linking"` // Gauss linking of the first two lines
	Lines   []Line  `json:"lines"`
}

// VizRadii are the seed radii on the x–y plane (one nest of lines each)
// and VizAngles the number of seeds round each.
var VizRadii = []float64{0.3, 0.55, 0.8}

const VizAngles = 8

func writeViz(dir string, times []float64) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	h := Hopfion{}
	sign := h.Sign()
	var frames []Frame
	for _, t := range times {
		tt := t
		b := func(p [3]float64) [3]float64 {
			F := h.F(p[0], p[1], p[2], tt)
			return [3]float64{sign * imag(F[0]), sign * imag(F[1]), sign * imag(F[2])}
		}
		fr := Frame{T: t}
		var first [][][3]float64
		for ti, r := range VizRadii {
			for a := 0; a < VizAngles; a++ {
				phi := 2 * math.Pi * (float64(a) + 0.5) / VizAngles // half-step off the y-axis, whose line passes through infinity
				// seeds scale with the knot's size √(1 + t²), so the same family
				// of lines is followed as the knot grows
				sc := math.Sqrt(1 + t*t)
				seed := [3]float64{sc * r * math.Cos(phi), sc * r * math.Sin(phi), 0}
				pts, closed := trace(b, seed, 0.01*math.Sqrt(1+t*t), 80*math.Sqrt(1+t*t))
				if len(first) < 2 && a == 0 && closed {
					first = append(first, pts)
				}
				keep := len(pts)/240 + 1
				var ds [][3]float64
				for i := 0; i < len(pts); i += keep {
					p := pts[i]
					ds = append(ds, [3]float64{round(p[0]), round(p[1]), round(p[2])})
				}
				fr.Lines = append(fr.Lines, Line{Torus: ti, Closed: closed, Points: ds})
			}
		}
		fr.Linking = Linking(first[0], first[1])
		// rms radius of the exact energy on a coarse grid
		var u, ur float64
		for i := -60; i <= 60; i++ {
			for j := -60; j <= 60; j++ {
				for l := -60; l <= 60; l++ {
					x, y, z := float64(i)/10, float64(j)/10, float64(l)/10
					F := h.F(x, y, z, t)
					e := real(F[0])*real(F[0]) + imag(F[0])*imag(F[0]) + real(F[1])*real(F[1]) + imag(F[1])*imag(F[1]) + real(F[2])*real(F[2]) + imag(F[2])*imag(F[2])
					u += e
					ur += e * (x*x + y*y + z*z)
				}
			}
		}
		fr.RMS = math.Sqrt(ur / u)
		frames = append(frames, fr)
	}
	f, err := os.Create(filepath.Join(dir, "knot_lines.json"))
	if err != nil {
		return err
	}
	defer f.Close()
	return json.NewEncoder(f).Encode(frames)
}

func round(x float64) float64 { return math.Round(x*1e4) / 1e4 }
