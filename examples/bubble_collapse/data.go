// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"bytes"
	_ "embed"
	"encoding/json"
	"errors"
	"math"
	"math/cmplx"
	"os"

	"github.com/fdamador/kronspice/analysis"
)

// writeData runs every part and writes the series behind the report
// (go run ./examples/bubble_collapse -data report_data.json): the cycle and
// its ledger, the collapse close-up with the flash, the chain sweeps, the
// reversal paths, the flash spectrum, the polarization on the probe
// sphere, |E|² slices from the field solver and the knot's field lines
// from t = −3 to 3 L0/c for the animation.
func writeData(path string) error {
	out, err := reportData()
	if err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return json.NewEncoder(f).Encode(out)
}

//go:embed viz/report_template.html
var reportTemplate []byte

// writeReport writes the interactive report (viz/report.html): the
// template with the data series embedded. It opens in a browser on its own
// (three.js and the fonts load from cdnjs and Google Fonts).
func writeReport(path string) error {
	out, err := reportData()
	if err != nil {
		return err
	}
	js, err := json.Marshal(out)
	if err != nil {
		return err
	}
	if bytes.Contains(js, []byte("</")) || !bytes.Contains(reportTemplate, []byte("__DATA__")) {
		return errors.New("bubble_collapse: report template or data unusable")
	}
	return os.WriteFile(path, bytes.Replace(reportTemplate, []byte("__DATA__"), js, 1), 0o644)
}

// reportData runs every part and collects the report's series.
func reportData() (map[string]any, error) {
	d0, err := readDrive(0, 1.3*atm)
	if err != nil {
		return nil, err
	}
	b := DefaultBubble(d0.FPeak)
	cy, s, err := b.Run(6)
	if err != nil {
		return nil, err
	}
	T := 1 / b.F
	tmin := math.Mod(cy.TRmin, T)
	out := map[string]any{}
	r6 := func(x float64) float64 { return math.Round(x*1e6) / 1e6 }

	// cycle (µs, µm, K) and cumulative ledger (nJ)
	type pt = []float64
	var cyc, led []pt
	stride := len(s)/700 + 1
	for i, x := range s {
		near := math.Abs(x.T-tmin) < 300e-9
		if i%stride != 0 && !(near && i%3 == 0) && i != len(s)-1 {
			continue
		}
		cyc = append(cyc, pt{r6(x.T * 1e6), r6(x.R * 1e6), math.Round(x.Temp)})
		led = append(led, pt{r6(x.T * 1e6), r6(x.In * 1e9), r6(x.Visc * 1e9), r6(x.Rad * 1e9), r6(x.Stored * 1e9)})
	}
	var km []pt
	k := 0
	_, _, _, _, err = b.runKM(6, func(t, R float64) {
		if k++; k%3 == 0 || math.Abs(t-tmin) < 300e-9 {
			km = append(km, pt{r6(t * 1e6), r6(R * 1e6)})
		}
	})
	if err != nil {
		return nil, err
	}
	L := b.FlashEstimate(s)
	var zoom []pt // ns from the minimum, µm, K, radiated sound W, light µW
	for i, x := range s {
		dt := x.T - tmin
		if math.Abs(dt) < 4e-9 {
			zoom = append(zoom, pt{r6(dt * 1e9), r6(x.R * 1e6), math.Round(x.Temp), x.Prad, L.Power[i] * 1e6})
		}
	}
	ls, spec := flashSpectrum(b, s)
	tfit := PlanckFit(ls, spec)
	var sp []pt
	mx, bmx := 0.0, 0.0
	bb := make([]float64, len(ls))
	for i, l := range ls {
		mx = math.Max(mx, spec[i])
		bb[i] = 1 / (math.Pow(l, 5) * (math.Exp(hPlanck*cLight/(l*kB*tfit)) - 1))
		bmx = math.Max(bmx, bb[i])
	}
	for i, l := range ls {
		sp = append(sp, pt{l * 1e9, r6(spec[i] / mx), r6(bb[i] / bmx)})
	}
	out["bubble"] = map[string]any{
		"R0": b.R0 * 1e6, "Pa_atm": b.Pa / atm, "f_kHz": b.F / 1e3, "T_us": T * 1e6, "tmin_us": tmin * 1e6,
		"Rmax": cy.Rmax * 1e6, "Rmin": cy.Rmin * 1e6, "Tmax": cy.Tmax, "Pmax": cy.Pmax, "Mach": cy.MaxMach,
		"In_nJ": cy.In * 1e9, "Visc_nJ": cy.Viscous * 1e9, "Rad_nJ": cy.Radiated * 1e9, "Residual": cy.Residual / cy.In,
		"Burst": cy.BurstFrac, "Collapse_nJ": cy.CollapseWork * 1e9,
		"cycle": cyc, "km": km, "ledger": led, "zoom": zoom,
		"light": map[string]any{"E_fJ": L.Energy * 1e15, "FWHM_ps": L.FWHM * 1e12, "photons": L.Photons, "ion": L.PeakIon,
			"Tfit": tfit, "uv_ps": b.BandWidth(s, 250e-9, 300e-9) * 1e12, "red_ps": b.BandWidth(s, 600e-9, 700e-9) * 1e12},
		"spectrum": sp,
	}

	// chain sweeps, pressure per volt [kPa/V]
	d1, err := readDrive(cy.Y1, 1.3*atm)
	if err != nil {
		return nil, err
	}
	sweep := func(y complex128) ([]pt, error) {
		c, f0, err := chain(y)
		if err != nil {
			return nil, err
		}
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 401, FStart: 0.998 * f0, FStop: 1.002 * f0}, analysis.Options{})
		if err != nil {
			return nil, err
		}
		p, _ := r.NodeVoltage("p")
		var o []pt
		for i := range p {
			o = append(o, pt{r6(r.Freq[i] / 1e3), r6(cmplx.Abs(p[i]) / 1e3)})
		}
		return o, nil
	}
	s0, err := sweep(0)
	if err != nil {
		return nil, err
	}
	s1, err := sweep(cy.Y1)
	if err != nil {
		return nil, err
	}
	out["chain"] = map[string]any{"without": s0, "with": s1,
		"rows": [][]float64{{d0.FPeak, d0.PerVolt, d0.EMF, d0.Source, d0.Mode, 0}, {d1.FPeak, d1.PerVolt, d1.EMF, d1.Source, d1.Mode, d1.Bubble}}}

	// reversal paths (ns since R_max, µm)
	var tMax, rm float64
	for _, x := range s {
		if x.T < tmin && x.R > rm {
			rm, tMax = x.R, x.T
		}
	}
	tStart := 5*T + tMax
	half := cy.TRmin - tStart
	rev := map[string]any{"span_ns": 2 * half * 1e9}
	for name, l := range map[string]Losses{"lossless": {}, "full": {true, true}, "viscous": {true, false}, "radiation": {false, true}} {
		var fw, bw []pt
		n := 0
		dR, _, err := b.Reversal(l, tStart, half, func(pass int, t, R float64) {
			if n++; n%2 != 0 {
				return
			}
			if pass == 0 {
				fw = append(fw, pt{r6(t * 1e9), r6(R * 1e6)})
			} else {
				bw = append(bw, pt{r6(t * 1e9), r6(R * 1e6)})
			}
		})
		if err != nil {
			return nil, err
		}
		rev[name] = map[string]any{"forward": fw, "back": bw, "dR": dR}
	}
	out["reversal"] = rev

	// field solver
	hs := DefaultHelicity()
	h, err := hs.RunAll()
	if err != nil {
		return nil, err
	}
	ps := hs.probes()
	var ang []pt
	for _, p := range ps {
		r := math.Sqrt(p.x*p.x + p.y*p.y + p.z*p.z)
		ang = append(ang, pt{r6(math.Acos(p.z / r)), r6(math.Atan2(p.y, p.x))})
	}
	var runs []map[string]any
	for i, r := range h.Runs {
		loc := make([]float64, len(r.ProbeV))
		w := make([]float64, len(r.ProbeV))
		imax := 0.0
		for p := range r.ProbeV {
			imax = math.Max(imax, r.ProbeI[p])
		}
		for p := range r.ProbeV {
			loc[p] = math.Round(1000*r.ProbeV[p]/r.ProbeI[p]) / 1000
			w[p] = math.Round(1000*r.ProbeI[p]/imax) / 1000
		}
		m := map[string]any{"name": r.Name, "pc": r.Pc, "residual": r.Residual, "local": loc, "weight": w}
		if !math.IsNaN(r.PcExact) {
			m["pcExact"] = r.PcExact
		}
		if i <= 3 { // slices for the hopfion, the two dipoles and the first thermal seed
			var snaps [][]int
			for _, sn := range r.Snaps {
				q := make([]int, len(sn))
				big := float32(0)
				for _, v := range sn {
					if v > big {
						big = v
					}
				}
				for j, v := range sn {
					q[j] = int(math.Round(255 * math.Sqrt(float64(v/big))))
				}
				snaps = append(snaps, q)
			}
			var st []float64
			for _, t := range r.SnapT {
				st = append(st, r6(t*em3c()/hs.Lambda0))
			}
			m["snaps"], m["snapT"] = snaps, st
		}
		runs = append(runs, m)
	}
	n := 0
	if len(h.Runs) > 0 && len(h.Runs[0].Snaps) > 0 {
		n = int(math.Round(math.Sqrt(float64(len(h.Runs[0].Snaps[0])))))
	}
	out["helicity"] = map[string]any{"probes": ang, "runs": runs, "thermalMean": h.ThermalMean, "thermalRMS": h.ThermalRMS,
		"sliceN": n, "halfLambda": float64(n-1) / 2 / float64(hs.Cells), "core": hs.Core, "probeR": hs.Probe}

	// knot frames
	out["knot"] = knotFrames()

	return out, nil
}

func em3c() float64 { return cLight }

// knotFrames traces the exact hopfion's B lines from t = −3 to 3 L0/c,
// seeds on three rings scaled with the knot's size √(1 + t²) (as
// knotted_light's viz does), each loop resampled to 96 points from its seed.
func knotFrames() map[string]any {
	h := Hopfion{}
	sign := h.Sign()
	radii := []float64{0.3, 0.55, 0.8}
	const angles, npts = 4, 96
	var frames []map[string]any
	for it := -15; it <= 15; it++ {
		t := float64(it) / 5
		sc := math.Sqrt(1 + t*t)
		bf := func(p [3]float64) [3]float64 {
			F := h.F(p[0], p[1], p[2], t)
			return [3]float64{sign * imag(F[0]), sign * imag(F[1]), sign * imag(F[2])}
		}
		var lines [][][3]float64
		for ri, r := range radii {
			for a := 0; a < angles; a++ {
				phi := 2*math.Pi*(float64(a)+0.5)/angles + float64(ri)*0.4
				seed := [3]float64{sc * r * math.Cos(phi), sc * r * math.Sin(phi), 0}
				line := traceLoop(bf, seed, 0.01*sc, 60*sc)
				lines = append(lines, resample(line, npts))
			}
		}
		for _, l := range lines {
			for j := range l {
				for c := 0; c < 3; c++ {
					l[j][c] = math.Round(l[j][c]*1e4) / 1e4
				}
			}
		}
		frames = append(frames, map[string]any{"t": t, "lines": lines})
	}
	return map[string]any{"frames": frames, "radii": radii}
}

// traceLoop follows a field line with RK4 on the unit direction and stops
// at its closest return to the seed after leaving it (as knotted_light's
// traceGap).
func traceLoop(b func(p [3]float64) [3]float64, p0 [3]float64, step, maxLen float64) [][3]float64 {
	dir := func(p [3]float64) [3]float64 {
		v := b(p)
		n := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
		if n == 0 {
			return [3]float64{}
		}
		return [3]float64{v[0] / n, v[1] / n, v[2] / n}
	}
	add := func(p, d [3]float64, s float64) [3]float64 {
		return [3]float64{p[0] + s*d[0], p[1] + s*d[1], p[2] + s*d[2]}
	}
	dist := func(p [3]float64) float64 {
		return math.Sqrt((p[0]-p0[0])*(p[0]-p0[0]) + (p[1]-p0[1])*(p[1]-p0[1]) + (p[2]-p0[2])*(p[2]-p0[2]))
	}
	p := p0
	line := [][3]float64{p}
	best, at := math.Inf(1), -1
	scale := maxLen / 60
	for s := 0.0; s < maxLen; s += step {
		k1 := dir(p)
		k2 := dir(add(p, k1, step/2))
		k3 := dir(add(p, k2, step/2))
		k4 := dir(add(p, k3, step))
		for c := 0; c < 3; c++ {
			p[c] += step / 6 * (k1[c] + 2*k2[c] + 2*k3[c] + k4[c])
		}
		line = append(line, p)
		if s < scale {
			continue
		}
		d := dist(p)
		if d < best {
			best, at = d, len(line)-1
		}
		if best < 0.3*scale && d > best+0.3*scale {
			break
		}
	}
	if at < 0 {
		return line
	}
	return line[:at+1]
}

// resample spaces a polyline into n points by arc length.
func resample(line [][3]float64, n int) [][3]float64 {
	cum := make([]float64, len(line))
	for i := 1; i < len(line); i++ {
		a, b := line[i-1], line[i]
		cum[i] = cum[i-1] + math.Sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))
	}
	out := make([][3]float64, n)
	j := 0
	for k := 0; k < n; k++ {
		s := cum[len(cum)-1] * float64(k) / float64(n-1)
		for j < len(cum)-2 && cum[j+1] < s {
			j++
		}
		w := 0.0
		if cum[j+1] > cum[j] {
			w = (s - cum[j]) / (cum[j+1] - cum[j])
		}
		for c := 0; c < 3; c++ {
			out[k][c] = line[j][c] + w*(line[j+1][c]-line[j][c])
		}
	}
	return out
}
