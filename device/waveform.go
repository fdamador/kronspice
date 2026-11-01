// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"sort"

	"github.com/fdamador/kronspice/circuit"
)

// wave is a source waveform with SPICE3 defaults resolved for one analysis.
// Definitions follow the SPICE3 user's manual (Quarles et al., "SPICE3
// Version 3f3 User's Manual", §4.1 Independent Sources).
type wave struct {
	kind circuit.WaveKind
	a    [7]float64 // resolved arguments (PULSE, SIN, EXP)
	pwl  []float64  // t0 v0 t1 v1 ... (PWL)
}

// resolveWave fills SPICE3 defaults. step and stop are the transient print
// step and stop time; both are 0 outside a transient analysis, in which case
// only t = 0 is ever evaluated.
func resolveWave(w circuit.Waveform, step, stop float64) wave {
	r := wave{kind: w.Kind}
	if step <= 0 {
		step = 1e-9
	}
	if stop <= 0 {
		stop = 1e-6
	}
	arg := func(i int, def float64) float64 {
		if i < len(w.Args) {
			return w.Args[i]
		}
		return def
	}
	switch w.Kind {
	case circuit.WavePulse:
		// PULSE(v1 v2 td tr tf pw per); tr, tf default TSTEP; pw, per default TSTOP.
		r.a = [7]float64{arg(0, 0), arg(1, 0), arg(2, 0), arg(3, step), arg(4, step), arg(5, stop), arg(6, stop)}
		if r.a[3] == 0 {
			r.a[3] = step
		}
		if r.a[4] == 0 {
			r.a[4] = step
		}
	case circuit.WaveSin:
		// SIN(vo va freq td theta phase); freq defaults 1/TSTOP; phase in degrees.
		r.a = [7]float64{arg(0, 0), arg(1, 0), arg(2, 1/stop), arg(3, 0), arg(4, 0), arg(5, 0)}
	case circuit.WaveExp:
		// EXP(v1 v2 td1 tau1 td2 tau2); tau1, tau2 default TSTEP; td2 defaults td1+TSTEP.
		td1 := arg(2, 0)
		r.a = [7]float64{arg(0, 0), arg(1, 0), td1, arg(3, step), arg(4, td1+step), arg(5, step)}
		if r.a[3] == 0 {
			r.a[3] = step
		}
		if r.a[5] == 0 {
			r.a[5] = step
		}
	case circuit.WavePWL:
		r.pwl = append([]float64(nil), w.Args...)
	}
	return r
}

// value evaluates the waveform at time t.
func (w *wave) value(t float64) float64 {
	a := &w.a
	switch w.kind {
	case circuit.WavePulse:
		v1, v2, td, tr, tf, pw, per := a[0], a[1], a[2], a[3], a[4], a[5], a[6]
		if t < td {
			return v1
		}
		tt := t - td
		if per > 0 {
			tt = math.Mod(tt, per)
		}
		switch {
		case tt < tr:
			return v1 + (v2-v1)*tt/tr
		case tt < tr+pw:
			return v2
		case tt < tr+pw+tf:
			return v2 + (v1-v2)*(tt-tr-pw)/tf
		default:
			return v1
		}
	case circuit.WaveSin:
		vo, va, f, td, theta, ph := a[0], a[1], a[2], a[3], a[4], a[5]*math.Pi/180
		if t < td {
			return vo + va*math.Sin(ph)
		}
		tt := t - td
		return vo + va*math.Exp(-tt*theta)*math.Sin(2*math.Pi*f*tt+ph)
	case circuit.WaveExp:
		v1, v2, td1, tau1, td2, tau2 := a[0], a[1], a[2], a[3], a[4], a[5]
		switch {
		case t < td1:
			return v1
		case t < td2:
			return v1 + (v2-v1)*(1-math.Exp(-(t-td1)/tau1))
		default:
			return v1 + (v2-v1)*(1-math.Exp(-(t-td1)/tau1)) + (v1-v2)*(1-math.Exp(-(t-td2)/tau2))
		}
	case circuit.WavePWL:
		p := w.pwl
		n := len(p) / 2
		if t <= p[0] {
			return p[1]
		}
		if t >= p[2*(n-1)] {
			return p[2*n-1]
		}
		// First point with time > t.
		k := sort.Search(n, func(i int) bool { return p[2*i] > t })
		t0, v0, t1, v1 := p[2*(k-1)], p[2*(k-1)+1], p[2*k], p[2*k+1]
		if t1 == t0 {
			return v1
		}
		return v0 + (v1-v0)*(t-t0)/(t1-t0)
	}
	return 0
}

// slope returns the time derivative of the waveform at t. At a corner it is
// the slope of the segment that ends there (the left derivative), the one a
// time step landing on the corner has just integrated over.
func (w *wave) slope(t float64) float64 {
	a := &w.a
	switch w.kind {
	case circuit.WavePulse:
		v1, v2, td, tr, tf, pw, per := a[0], a[1], a[2], a[3], a[4], a[5], a[6]
		if t <= td {
			return 0
		}
		tt := t - td
		if per > 0 {
			tt = math.Mod(tt, per)
			if tt == 0 {
				tt = per // the end of the previous period
			}
		}
		switch {
		case tt <= tr:
			return (v2 - v1) / tr
		case tt <= tr+pw:
			return 0
		case tt <= tr+pw+tf:
			return (v1 - v2) / tf
		}
		return 0
	case circuit.WaveSin:
		va, f, td, theta, ph := a[1], a[2], a[3], a[4], a[5]*math.Pi/180
		if t <= td {
			return 0
		}
		tt, om := t-td, 2*math.Pi*f
		return va * math.Exp(-tt*theta) * (om*math.Cos(om*tt+ph) - theta*math.Sin(om*tt+ph))
	case circuit.WaveExp:
		v1, v2, td1, tau1, td2, tau2 := a[0], a[1], a[2], a[3], a[4], a[5]
		if t <= td1 {
			return 0
		}
		d := (v2 - v1) / tau1 * math.Exp(-(t-td1)/tau1)
		if t > td2 {
			d += (v1 - v2) / tau2 * math.Exp(-(t-td2)/tau2)
		}
		return d
	case circuit.WavePWL:
		p := w.pwl
		n := len(p) / 2
		if t <= p[0] || t > p[2*(n-1)] {
			return 0
		}
		k := sort.Search(n, func(i int) bool { return p[2*i] >= t }) // segment k−1 → k ends at or after t
		t0, v0, t1, v1 := p[2*(k-1)], p[2*(k-1)+1], p[2*k], p[2*k+1]
		if t1 == t0 {
			return 0
		}
		return (v1 - v0) / (t1 - t0)
	}
	return 0
}

// nextBreakpoint returns the first waveform corner strictly after t, or +Inf.
// The transient engine lands a time point exactly on each breakpoint.
func (w *wave) nextBreakpoint(t float64) float64 {
	const rel = 1e-12
	after := func(b float64) bool { return b > t+rel*math.Max(math.Abs(t), 1e-15) }
	a := &w.a
	switch w.kind {
	case circuit.WavePulse:
		td, tr, tf, pw, per := a[2], a[3], a[4], a[5], a[6]
		edges := [4]float64{0, tr, tr + pw, tr + pw + tf}
		if per <= 0 {
			for _, e := range edges {
				if b := td + e; after(b) {
					return b
				}
			}
			return math.Inf(1)
		}
		k := math.Floor((t - td) / per)
		if k < 0 {
			k = 0
		}
		for kk := k; kk <= k+1; kk++ {
			base := td + kk*per
			for _, e := range edges {
				if b := base + e; after(b) {
					return b
				}
			}
		}
		return td + (k+2)*per
	case circuit.WaveSin:
		if td := a[3]; after(td) {
			return td
		}
	case circuit.WaveExp:
		if after(a[2]) {
			return a[2]
		}
		if after(a[4]) {
			return a[4]
		}
	case circuit.WavePWL:
		p := w.pwl
		n := len(p) / 2
		k := sort.Search(n, func(i int) bool { return after(p[2*i]) })
		if k < n {
			return p[2*k]
		}
	}
	return math.Inf(1)
}
