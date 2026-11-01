// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func books(r *Result) float64 {
	in := r.ESource + r.EPump
	return math.Abs(in-r.EOh-(r.EStored1-r.EStored0)) / math.Abs(in)
}

// Paper Fig. 4: constant C at resonance, amplitude U0/R = 0.25 A, all from
// the source.
func TestClassical(t *testing.T) {
	f0 := W0 / (2 * math.Pi)
	r, err := Run(Case{Drive: drive(W0), DriveWave: sinWave(0, U0, W0, 0), Fd: f0, TEnd: 0.03})
	if err != nil {
		t.Fatal(err)
	}
	if ipk, _ := r.PeakLast(5 / f0); math.Abs(ipk-0.25) > 0.0025 {
		t.Errorf("amplitude %g A, want 0.25", ipk)
	}
	if b := books(r); b > 1e-4 {
		t.Errorf("books %g", b)
	}
}

// §3.1: tuned pumping grows without limit; the pump pays for it.
func TestTunedGrowth(t *testing.T) {
	f0 := W0 / (2 * math.Pi)
	s, ds := sinScale(0.2, 2*W0, 0)
	r, err := Run(Case{Drive: drive(W0), DriveWave: sinWave(0, U0, W0, 0), Fd: f0, Scale: s, DScale: ds, ScaleWave: sinWave(1, 0.2, 2*W0, 0), TEnd: 0.03})
	if err != nil {
		t.Fatal(err)
	}
	if ipk, _ := r.PeakLast(5 / f0); ipk < 50 {
		t.Errorf("peak %g A, expected growth", ipk)
	}
	if b := books(r); b > 2e-3 {
		t.Errorf("books %g", b)
	}
	if r.EPump < 0.99*(r.EOh+r.EStored1-r.EStored0) {
		t.Errorf("pump %g J does not pay for R %g J + stored %g J", r.EPump, r.EOh, r.EStored1-r.EStored0)
	}
}

// §3.1 at 0.97 f0: bounded (below the detuned threshold once the mean
// elastance shift is counted), as the paper reports.
func TestDetunedBounded(t *testing.T) {
	if m, mc := ElastanceIndex(0.2), Threshold((MeanShift(0.2)-0.97)/0.97); m >= mc {
		t.Fatalf("m %g ≥ m_c %g", m, mc)
	}
	f0 := W0 / (2 * math.Pi)
	w := 0.97 * W0
	s, ds := sinScale(0.2, 2*w, -math.Pi/2)
	r, err := Run(Case{Drive: drive(w), DriveWave: sinWave(0, U0, w, 0), Fd: 0.97 * f0, Scale: s, DScale: ds, ScaleWave: sinWave(1, 0.2, 2*w, -math.Pi/2), TEnd: 0.6})
	if err != nil {
		t.Fatal(err)
	}
	peak := func(t0, t1 float64) float64 {
		var p float64
		for k := range r.T {
			if r.T[k] > t0 && r.T[k] <= t1 {
				p = math.Max(p, math.Abs(r.I[k]))
			}
		}
		return p
	}
	p4, p6 := peak(0.39, 0.4), peak(0.59, 0.6)
	if p6 > 1.02*p4 || p6 < 10 {
		t.Errorf("peaks %g A at 0.4 s, %g A at 0.6 s: expected a bounded ~16 A", p4, p6)
	}
	if b := books(r); b > 1e-3 {
		t.Errorf("books %g", b)
	}
}

// §3.2: pumping at the drive frequency needs a deep pump (a0 ≈ 0.7);
// a0 = 0.2 does not grow, a0 = 0.8 does, paid by the pump.
func TestSecondTongue(t *testing.T) {
	f0 := W0 / (2 * math.Pi)
	wr := 0.95 * W0
	Tr := 2 * math.Pi / wr
	for _, a := range []float64{0.2, 0.8} {
		sq, dsq, pulse := Square(a, Tr)
		r, err := Run(Case{Fd: 0.95 * f0, Scale: sq, DScale: dsq, ScaleWave: pulse, I0: 1e-3, TEnd: 0.006, PerPeriod: 2000})
		if err != nil {
			t.Fatal(err)
		}
		ipk, _ := r.PeakLast(2 / (0.95 * f0))
		if grow := ipk > 1e-3; grow != (a > 0.5) {
			t.Errorf("a0 %g: peak %g A after a 1 mA start", a, ipk)
		}
		if a > 0.5 && r.EPump < 0.95*(r.EOh+r.EStored1-r.EStored0) {
			t.Errorf("pump %g J vs R %g J + stored %g J", r.EPump, r.EOh, r.EStored1-r.EStored0)
		}
	}
}
