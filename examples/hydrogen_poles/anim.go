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

// A small rasterizer for the animations (standard library only): polygons
// filled by scanline at 2× resolution, box-filtered down, mapped to a
// palette of the figures' colours in tints and shades.

var (
	rgbBg     = [3]float64{252, 252, 251}
	rgbN      = [3]float64{235, 104, 52}
	rgbS      = [3]float64{42, 120, 214}
	rgbInk    = [3]float64{11, 11, 11}
	rgbMuted  = [3]float64{154, 152, 143}
	tintSteps = 8
	shadeLo   = 0.55
	shadeN    = 14
)

type canvas struct {
	w, h int
	pix  []float64
}

func newCanvas(w, h int) *canvas {
	c := &canvas{w: w, h: h, pix: make([]float64, 3*w*h)}
	for i := 0; i < w*h; i++ {
		copy(c.pix[3*i:3*i+3], rgbBg[:])
	}
	return c
}

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

func (c *canvas) line(a, b [2]float64, w float64, rgb [3]float64) {
	dx, dy := b[0]-a[0], b[1]-a[1]
	l := math.Hypot(dx, dy)
	if l == 0 {
		return
	}
	nx, ny := -dy/l*w/2, dx/l*w/2
	ex, ey := dx/l*w/2, dy/l*w/2
	c.fill([][2]float64{{a[0] - ex + nx, a[1] - ey + ny}, {b[0] + ex + nx, b[1] + ey + ny}, {b[0] + ex - nx, b[1] + ey - ny}, {a[0] - ex - nx, a[1] - ey - ny}}, rgb, 1)
}

func (c *canvas) polyline(pts [][2]float64, w float64, rgb [3]float64) {
	for i := 1; i < len(pts); i++ {
		c.line(pts[i-1], pts[i], w, rgb)
	}
}

// arrow draws a vertical arrow from (x, y0) to (x, y1) (pixels).
func (c *canvas) arrow(x, y0, y1, w float64, rgb [3]float64) {
	d := math.Copysign(1, y1-y0)
	head := 3.2 * w
	c.line([2]float64{x, y0}, [2]float64{x, y1 - d*head}, w, rgb)
	c.fill([][2]float64{{x, y1}, {x - 1.6*head/1.8, y1 - d*head}, {x + 1.6*head/1.8, y1 - d*head}}, rgb, 1)
}

func (c *canvas) disc(x, y, r float64, rgb [3]float64) {
	var pts [][2]float64
	for k := 0; k < 24; k++ {
		a := 2 * math.Pi * float64(k) / 24
		pts = append(pts, [2]float64{x + r*math.Cos(a), y + r*math.Sin(a)})
	}
	c.fill(pts, rgb, 1)
}

func tint(c [3]float64, t float64) [3]float64 {
	return [3]float64{rgbBg[0] + (c[0]-rgbBg[0])*t, rgbBg[1] + (c[1]-rgbBg[1])*t, rgbBg[2] + (c[2]-rgbBg[2])*t}
}

func shade(c [3]float64, f float64) [3]float64 { return [3]float64{c[0] * f, c[1] * f, c[2] * f} }

// animPalette: N and S in tints × shades, grey and ink ramps, background.
func animPalette() color.Palette {
	pal := color.Palette{}
	add := func(c [3]float64) {
		pal = append(pal, color.RGBA{uint8(math.Round(c[0])), uint8(math.Round(c[1])), uint8(math.Round(c[2])), 255})
	}
	add(rgbBg)
	for _, base := range [][3]float64{rgbN, rgbS} {
		for t := 1; t <= tintSteps; t++ {
			for s := 0; s < shadeN; s++ {
				add(shade(tint(base, float64(t)/float64(tintSteps)), shadeLo+(1-shadeLo)*float64(s)/float64(shadeN-1)))
			}
		}
	}
	for i := 1; i <= 30; i++ {
		add(tint(rgbInk, float64(i)/30))
	}
	return pal
}

// quantTint and quantShade snap to the palette's levels.
func quantTint(t float64) float64 {
	return math.Max(1, math.Round(t*float64(tintSteps))) / float64(tintSteps)
}
func quantShade(f float64) float64 {
	k := math.Round((f - shadeLo) / (1 - shadeLo) * float64(shadeN-1))
	return shadeLo + (1-shadeLo)*math.Max(0, math.Min(float64(shadeN-1), k))/float64(shadeN-1)
}

func (c *canvas) paletted(pal color.Palette) *image.Paletted {
	w, h := c.w/2, c.h/2
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

func writeGIF(path string, imgs []*image.Paletted, delay int) error {
	g := &gif.GIF{}
	for _, im := range imgs {
		g.Image = append(g.Image, im)
		g.Delay = append(g.Delay, delay)
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

// AtomGIF: spheres round three levels (2p½, 2p3⁄2 m = 3⁄2, 2p3⁄2 m = ½), each
// coloured by B_r on it (orange N, blue S, deeper = stronger, scaled to the
// sphere's own peak), while the radius sweeps 1 → 32 a0 and back and the
// view turns. A bar at the foot shows the radius on a log scale, with a
// tick at 15.6 a0 where the m = ½ level's four zones end.
func AtomGIF(path string, fs []AtomField) error {
	const frames, W, H = 72, 720, 300
	pal := animPalette()
	imgs := make([]*image.Paletted, frames)
	logLo, logHi := math.Log(1), math.Log(32)
	parallel(frames, func(fr int) {
		ph := float64(fr) / frames
		s := 0.5 - 0.5*math.Cos(2*math.Pi*ph)
		R := math.Exp(logLo + (logHi-logLo)*s)
		az := 2 * math.Pi * ph
		el := 0.32
		c := newCanvas(2*W, 2*H)
		for k, f := range fs {
			cx, cy, rad := 2*(120+240*float64(k)), 2*128.0, 2*92.0
			// B_r(θ) on this sphere, scaled to its peak
			const nt = 72
			br := make([]float64, nt)
			mx := 0.0
			for i := range br {
				br[i] = f.Br(R, math.Pi*(float64(i)+0.5)/nt)
				mx = math.Max(mx, math.Abs(br[i]))
			}
			const np = 64
			for i := 0; i < nt; i++ {
				t0, t1 := math.Pi*float64(i)/nt, math.Pi*float64(i+1)/nt
				base := rgbN
				if br[i] < 0 {
					base = rgbS
				}
				tt := quantTint(0.25 + 0.75*math.Sqrt(math.Abs(br[i])/mx))
				for j := 0; j < np; j++ {
					p0, p1 := 2*math.Pi*float64(j)/np, 2*math.Pi*float64(j+1)/np
					pt := func(th, p float64) ([2]float64, [3]float64) {
						v := [3]float64{math.Sin(th) * math.Cos(p), math.Sin(th) * math.Sin(p), math.Cos(th)}
						ca, sa := math.Cos(az), math.Sin(az)
						x1 := v[0]*ca + v[1]*sa
						y1 := -v[0]*sa + v[1]*ca
						ce, se := math.Cos(el), math.Sin(el)
						return [2]float64{cx + rad*x1, cy - rad*(v[2]*ce+y1*se)}, [3]float64{x1, y1, v[2]}
					}
					a, n := pt(0.5*(t0+t1), 0.5*(p0+p1))
					_ = a
					// facing: n · (0, −cos el, sin el)
					if -n[1]*math.Cos(el)+n[2]*math.Sin(el) <= 0 {
						continue
					}
					light := math.Max(0, -0.4*n[0]-0.5*n[1]*math.Cos(el)+0.75*(n[2]*math.Cos(el)+n[1]*math.Sin(el)))
					col := shade(tint(base, tt), quantShade(0.62+0.38*light))
					q0, _ := pt(t0, p0)
					q1, _ := pt(t0, p1)
					q2, _ := pt(t1, p1)
					q3, _ := pt(t1, p0)
					c.fill([][2]float64{q0, q1, q2, q3}, col, 1)
				}
			}
			// spin axis stubs
			ce := math.Cos(el)
			c.line([2]float64{cx, cy - rad*ce}, [2]float64{cx, cy - rad*ce - 2*26}, 3, rgbInk)
			c.line([2]float64{cx, cy + rad*ce}, [2]float64{cx, cy + rad*ce + 2*26}, 3, rgbInk)
		}
		// radius bar (log scale 1 … 32 a0), tick at 15.6 a0
		x0, x1, y := 2*120.0, 2*600.0, 2*278.0
		c.line([2]float64{x0, y}, [2]float64{x1, y}, 4, tint(rgbInk, 0.35))
		for _, r := range []float64{1, 2, 5, 10, 20, 32} {
			x := x0 + (x1-x0)*math.Log(r)/logHi
			c.line([2]float64{x, y - 8}, [2]float64{x, y + 8}, 3, tint(rgbInk, 0.5))
		}
		xc := x0 + (x1-x0)*math.Log(15.6)/logHi
		c.line([2]float64{xc, y - 16}, [2]float64{xc, y + 16}, 3, rgbN)
		c.disc(x0+(x1-x0)*math.Log(R)/logHi, y, 11, rgbInk)
		imgs[fr] = c.paletted(pal)
	})
	return writeGIF(path, imgs, 9)
}

// BondGIF: two hydrogen atoms brought together from 7 a0 to the bond and
// back. Left, the triplet (spins ↑↑): its spin density in orange. Right,
// the singlet (↑↓): electron density in grey, spin density zero. Arrows are
// the electron spins (long) and proton spins (short, ortho / para). Below,
// the Heitler–London energies with the current separation marked.
func BondGIF(path string) error {
	const frames, W, H = 64, 640, 470
	pal := animPalette()
	imgs := make([]*image.Paletted, frames)
	parallel(frames, func(fr int) {
		ph := float64(fr) / frames
		s := 0.5 - 0.5*math.Cos(2*math.Pi*ph)
		R := 7 - (7-1.4011)*s
		c := newCanvas(2*W, 2*H)
		for p, trip := range []bool{true, false} {
			m := Molecule{R: R, Triplet: trip}
			x0, y0, side := 2*(20+320*float64(p)), 2*10.0, 2*280.0
			ext := 5.0
			sc := side / (2 * ext)
			const n = 110
			cell := side / n
			peak := 0.0
			vals := make([]float64, n*n)
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					x := -ext + (float64(i)+0.5)*2*ext/n
					z := -ext + (float64(j)+0.5)*2*ext/n
					v := m.Density(math.Abs(x), z)
					vals[i*n+j] = v
					peak = math.Max(peak, v)
				}
			}
			base := tint(rgbInk, 0.55)
			if trip {
				base = rgbN
			}
			for i := 0; i < n; i++ {
				for j := 0; j < n; j++ {
					t := math.Sqrt(vals[i*n+j] / peak)
					if t < 0.04 {
						continue
					}
					var col [3]float64
					if trip {
						col = tint(base, quantTint(math.Min(1, t)))
					} else {
						col = tint(rgbInk, math.Round(math.Min(1, t)*0.5*30)/30)
					}
					xa, ya := x0+float64(i)*cell, y0+float64(n-1-j)*cell
					c.fill([][2]float64{{xa, ya}, {xa + cell + 0.5, ya}, {xa + cell + 0.5, ya + cell + 0.5}, {xa, ya + cell + 0.5}}, col, 1)
				}
			}
			// frame
			c.polyline([][2]float64{{x0, y0}, {x0 + side, y0}, {x0 + side, y0 + side}, {x0, y0 + side}, {x0, y0}}, 2, tint(rgbInk, 0.25))
			cx, cy := x0+side/2, y0+side/2
			for k, zp := range []float64{-R / 2, R / 2} {
				y := cy - sc*zp
				// electron spin arrows beside each atom: ↑↑ triplet, ↑↓ singlet
				dir := 1.0
				if !trip && k == 1 {
					dir = -1
				}
				c.arrow(cx-2*38, y+dir*2*22, y-dir*2*22, 5, rgbInk)
				// proton spin arrow at the nucleus (short)
				c.disc(cx, y, 2*5, rgbInk)
				c.arrow(cx+2*22, y+dir*2*11, y-dir*2*11, 3.5, tint(rgbInk, 0.6))
			}
		}
		// energy panel
		ax0, ax1, ay0, ay1 := 2*70.0, 2*620.0, 2*305.0, 2*455.0
		rx := func(r float64) float64 { return ax0 + (ax1-ax0)*(r-1)/(7.2-1) }
		ey := func(e float64) float64 { return ay1 - (ay1-ay0)*(e+5)/(11) }
		c.line([2]float64{ax0, ey(0)}, [2]float64{ax1, ey(0)}, 2, tint(rgbInk, 0.35))
		c.line([2]float64{ax0, ay0}, [2]float64{ax0, ay1}, 2, tint(rgbInk, 0.35))
		var ps, pt [][2]float64
		for i := 0; i <= 200; i++ {
			r := 1 + 6.2*float64(i)/200
			es := (HLEnergy(r, 1) + 1) * Hartree
			et := math.Min(6, (HLEnergy(r, -1)+1)*Hartree)
			ps = append(ps, [2]float64{rx(r), ey(es)})
			if (HLEnergy(r, -1)+1)*Hartree <= 6 {
				pt = append(pt, [2]float64{rx(r), ey(et)})
			}
		}
		c.polyline(ps, 4, rgbS)
		c.polyline(pt, 4, rgbN)
		c.disc(rx(R), ey((HLEnergy(R, 1)+1)*Hartree), 10, rgbS)
		if e := (HLEnergy(R, -1) + 1) * Hartree; e <= 6 {
			c.disc(rx(R), ey(e), 10, rgbN)
		}
		imgs[fr] = c.paletted(pal)
	})
	return writeGIF(path, imgs, 10)
}
