// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"image"
	"image/color"
	"image/gif"
	"math"
	"os"
	"sort"
)

// A small rasterizer for the rotating 3D view (standard library only):
// polygons filled by scanline at 2× resolution, box-filtered down, and
// mapped to a palette of shaded pole colours.

type canvas struct {
	w, h int
	pix  []float64 // RGB
}

func newCanvas(w, h int) *canvas {
	c := &canvas{w: w, h: h, pix: make([]float64, 3*w*h)}
	for i := range c.pix {
		c.pix[i] = 252
	}
	return c
}

// fill paints a polygon (even–odd rule, pixel centres) with alpha.
func (c *canvas) fill(pts [][2]float64, rgb [3]float64, alpha float64) {
	ymin, ymax := math.Inf(1), math.Inf(-1)
	for _, p := range pts {
		ymin, ymax = math.Min(ymin, p[1]), math.Max(ymax, p[1])
	}
	y0, y1 := int(math.Max(0, math.Ceil(ymin-0.5))), int(math.Min(float64(c.h-1), math.Floor(ymax-0.5)))
	xs := make([]float64, 0, 8)
	for y := y0; y <= y1; y++ {
		yc := float64(y) + 0.5
		xs = xs[:0]
		for i := range pts {
			a, b := pts[i], pts[(i+1)%len(pts)]
			if (a[1] <= yc) != (b[1] <= yc) {
				xs = append(xs, a[0]+(yc-a[1])/(b[1]-a[1])*(b[0]-a[0]))
			}
		}
		sort.Float64s(xs)
		for k := 0; k+1 < len(xs); k += 2 {
			xa := int(math.Max(0, math.Ceil(xs[k]-0.5)))
			xb := int(math.Min(float64(c.w-1), math.Floor(xs[k+1]-0.5)))
			for x := xa; x <= xb; x++ {
				i := 3 * (y*c.w + x)
				for ch := 0; ch < 3; ch++ {
					c.pix[i+ch] += alpha * (rgb[ch] - c.pix[i+ch])
				}
			}
		}
	}
}

// line paints a segment of width w as a thin quad.
func (c *canvas) line(a, b [2]float64, w float64, rgb [3]float64, alpha float64) {
	dx, dy := b[0]-a[0], b[1]-a[1]
	l := math.Hypot(dx, dy)
	if l == 0 {
		return
	}
	nx, ny := -dy/l*w/2, dx/l*w/2
	ex, ey := dx/l*w/2, dy/l*w/2 // round-ish caps by extending
	c.fill([][2]float64{{a[0] - ex + nx, a[1] - ey + ny}, {b[0] + ex + nx, b[1] + ey + ny}, {b[0] + ex - nx, b[1] + ey - ny}, {a[0] - ex - nx, a[1] - ey - ny}}, rgb, alpha)
}

// ringPalette: shades of the N, S, tinted-face and line colours, and their
// blends toward the background.
func ringPalette() color.Palette {
	pal := color.Palette{color.RGBA{252, 252, 251, 255}}
	add := func(c [3]float64) {
		pal = append(pal, color.RGBA{uint8(math.Round(c[0])), uint8(math.Round(c[1])), uint8(math.Round(c[2])), 255})
	}
	for _, base := range [][3]float64{rgbN, rgbS, tint(rgbN, 0.45), tint(rgbS, 0.45), rgbLine} {
		for i := 0; i < 32; i++ { // shading 0.6 … 1.0
			f := 0.6 + 0.4*float64(i)/31
			add([3]float64{base[0] * f, base[1] * f, base[2] * f})
		}
		for i := 1; i < 18; i++ { // blends toward the background (edges)
			add(tint(base, float64(i)/18))
		}
	}
	return pal
}

// frame renders one view to a paletted image of w×h (rendered at 2×).
func (v RingView) frame(cam Camera, prims []Prim, w, h int, pal color.Palette) *image.Paletted {
	c := newCanvas(2*w, 2*h)
	scale := 0.3 * float64(2*w) / v.Ring.Ro
	for _, q := range cam.Paint(prims, scale, float64(w), float64(h)) {
		if q.Line {
			c.line(q.XY[0], q.XY[1], 2.4, q.RGB, 0.9)
		} else {
			c.fill(q.XY, q.RGB, 1)
		}
	}
	img := image.NewPaletted(image.Rect(0, 0, w, h), pal)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			var s [3]float64
			for dy := 0; dy < 2; dy++ {
				for dx := 0; dx < 2; dx++ {
					i := 3 * ((2*y+dy)*c.w + 2*x + dx)
					for ch := 0; ch < 3; ch++ {
						s[ch] += c.pix[i+ch] / 4
					}
				}
			}
			img.SetColorIndex(x, y, uint8(pal.Index(color.RGBA{uint8(s[0]), uint8(s[1]), uint8(s[2]), 255})))
		}
	}
	return img
}

// RingGIF writes a turntable animation of the ring and its field lines.
func RingGIF(path string, v RingView) error {
	const frames, w, h = 48, 480, 400
	pal := ringPalette()
	prims := v.Ring.Scene(v.Lines)
	g := &gif.GIF{}
	imgs := make([]*image.Paletted, frames)
	parallel(frames, func(i int) {
		cam := Camera{Az: -0.9 + 2*math.Pi*float64(i)/frames, El: 0.5}
		imgs[i] = v.frame(cam, prims, w, h, pal)
	})
	for _, im := range imgs {
		g.Image = append(g.Image, im)
		g.Delay = append(g.Delay, 8)
	}
	fo, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := gif.EncodeAll(fo, g); err != nil {
		fo.Close()
		return err
	}
	return fo.Close()
}
