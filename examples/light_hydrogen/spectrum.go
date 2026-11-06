// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The spectrum from the network's frequency response. Every s state of the
// atom is a resonance of Kron's network (dual form, ω² = 1/(E + V0)), so an
// AC sweep of the whole network shows the atom's levels as poles of its
// driving-point impedance. Each pole is pinned by bisection on the sign of
// the reactance; two networks (h = 0.1 and 0.05 a0) are extrapolated to
// h → 0 (the scheme's error is O(h²)); the levels' differences are the
// spectral lines. Non-relativistic Coulomb levels depend on n only, so the
// s levels fix the Lyman (np → 1s) and Balmer (n → 2) lines.

// NLevels is how many levels the spectrum resolves (n = 1…5).
const NLevels = 5

// Spectrum is one atom's frequency response and levels.
type Spectrum struct {
	Z        float64
	Nodes    int
	X, LogZ  []float64 // sweep: x = 1/√(−2E) (= n/Z at the poles), log10|Z_in| [Ω]
	Circuit  []float64 // levels from the h = 0.1 network's poles [Ha]
	Fine     []float64 // the h = 0.05 network's poles [Ha]
	Matrix   []float64 // the h = 0.1 network's matrix eigenvalues [Ha]
	Extrap   []float64 // circuit poles, Richardson h → 0 [Ha]
	SweepPts int
}

// acProbe is an AC solver for one network at single energies.
type acProbe struct {
	ckt   *circuit.Circuit
	drive string
	v0    float64
}

func (p acProbe) f(e float64) float64 { return 1 / math.Sqrt(e+p.v0) / (2 * math.Pi) }

// reactance is the sign-carrying part of the driving-point response at
// energy e: Im V/|V| (V = Z_in for a 1 A drive).
func (p acProbe) reactance(e float64) (float64, error) {
	f := p.f(e)
	res, err := analysis.AC(p.ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		return 0, err
	}
	v, ok := res.NodeVoltage(p.drive)
	if !ok {
		return 0, fmt.Errorf("no node %s", p.drive)
	}
	return imag(v[0]) / math.Max(math.Abs(imag(v[0]))+math.Abs(real(v[0])), 1e-300), nil
}

// pole pins the resonance in [lo, hi] (reactance − at lo, + at hi).
func (p acProbe) pole(lo, hi float64) (float64, error) {
	for i := 0; i < 100 && hi-lo > 1e-14; i++ {
		m := 0.5 * (lo + hi)
		s, err := p.reactance(m)
		if err != nil {
			return 0, err
		}
		if s < 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi), nil
}

func newProbe(z, h float64) (acProbe, Ladder, error) {
	l := NewLadder(z, h, 130/z)
	ckt, drive, v0, err := l.Network()
	return acProbe{ckt, drive, v0}, l, err
}

// RunSpectrum sweeps the network for coupling z and resolves its first
// NLevels levels.
func RunSpectrum(z float64, points int) (Spectrum, error) {
	sp := Spectrum{Z: z, SweepPts: points}
	p, l, err := newProbe(z, 0.1)
	if err != nil {
		return sp, err
	}
	sp.Nodes = len(l.D)
	// One sweep over E from −0.889 to −0.0149 Ha (x = 0.75 … 5.8).
	eLo, eHi := -1/(2*0.75*0.75), -1/(2*5.8*5.8)
	res, err := analysis.AC(p.ckt, analysis.ACSpec{Variation: "lin", Points: points, FStart: p.f(eHi), FStop: p.f(eLo)}, analysis.Options{})
	if err != nil {
		return sp, err
	}
	v, ok := res.NodeVoltage(p.drive)
	if !ok {
		return sp, fmt.Errorf("no node %s", p.drive)
	}
	// Walk upward in E (downward in frequency); a pole is a − → + change.
	var prevE, prevS float64
	first := true
	for i := len(res.Freq) - 1; i >= 0; i-- {
		w := 2 * math.Pi * res.Freq[i]
		e := 1/(w*w) - p.v0
		sp.X = append(sp.X, 1/math.Sqrt(-2*e))
		sp.LogZ = append(sp.LogZ, math.Log10(math.Max(math_abs(v[i]), 1e-300)))
		s := imag(v[i])
		if !first && prevS < 0 && s > 0 && len(sp.Circuit) < NLevels {
			e0, err := p.pole(prevE, e)
			if err != nil {
				return sp, err
			}
			sp.Circuit = append(sp.Circuit, e0)
		}
		prevE, prevS, first = e, s, false
	}
	if len(sp.Circuit) < NLevels {
		return sp, fmt.Errorf("found %d of %d poles", len(sp.Circuit), NLevels)
	}
	for k := 0; k < NLevels; k++ {
		sp.Matrix = append(sp.Matrix, l.Eigen(k))
	}
	// The finer network: each pole lies within 0.3 % of the coarse one.
	q, _, err := newProbe(z, 0.05)
	if err != nil {
		return sp, err
	}
	for k, e := range sp.Circuit {
		lo, hi := e*1.004, e*0.996
		if k > 0 {
			lo = math.Max(lo, 0.5*(e+sp.Circuit[k-1]))
		}
		if k+1 < NLevels {
			hi = math.Min(hi, 0.5*(e+sp.Circuit[k+1]))
		}
		sl, err := q.reactance(lo)
		if err != nil {
			return sp, err
		}
		sh, err := q.reactance(hi)
		if err != nil {
			return sp, err
		}
		if sl >= 0 || sh <= 0 {
			return sp, fmt.Errorf("level %d: no pole bracketed in the fine network", k+1)
		}
		e2, err := q.pole(lo, hi)
		if err != nil {
			return sp, err
		}
		sp.Fine = append(sp.Fine, e2)
		sp.Extrap = append(sp.Extrap, (4*e2-e)/3)
	}
	return sp, nil
}

func math_abs(c complex128) float64 { return math.Hypot(real(c), imag(c)) }

// Line is one spectral line n_up → n_lo.
type Line struct {
	Name       string
	Up, Lo     int
	MeasuredNm float64 // vacuum wavelength, NIST ASD, rounded to 0.01 nm
	Circuit    float64 // from the network's poles, h → 0 [nm]
	Bohr       float64 // exact −Z²/2n² levels [nm]
}

// MeasuredLines are the Lyman and Balmer lines the five levels reach
// (vacuum wavelengths, NIST Atomic Spectra Database, fine-structure
// components merged, to 0.01 nm).
var MeasuredLines = []Line{
	{Name: "Lyman α", Up: 2, Lo: 1, MeasuredNm: 121.57},
	{Name: "Lyman β", Up: 3, Lo: 1, MeasuredNm: 102.57},
	{Name: "Lyman γ", Up: 4, Lo: 1, MeasuredNm: 97.25},
	{Name: "Lyman δ", Up: 5, Lo: 1, MeasuredNm: 94.97},
	{Name: "Balmer α (Hα)", Up: 3, Lo: 2, MeasuredNm: 656.46},
	{Name: "Balmer β (Hβ)", Up: 4, Lo: 2, MeasuredNm: 486.27},
	{Name: "Balmer γ (Hγ)", Up: 5, Lo: 2, MeasuredNm: 434.17},
}

// HaJ is the hartree for the hydrogen reduced mass [J].
var HaJ = 2 * RyJ / Recoil

// nm converts an energy difference in hartree to a vacuum wavelength.
func nm(de float64) float64 { return H * C0 / (de * HaJ) * 1e9 }

// Lines turns a spectrum's levels into the Lyman and Balmer lines.
func (sp Spectrum) Lines() []Line {
	out := make([]Line, len(MeasuredLines))
	for i, l := range MeasuredLines {
		l.Circuit = nm(sp.Extrap[l.Up-1] - sp.Extrap[l.Lo-1])
		bohr := func(n int) float64 { return -sp.Z * sp.Z / (2 * float64(n*n)) }
		l.Bohr = nm(bohr(l.Up) - bohr(l.Lo))
		out[i] = l
	}
	return out
}

// OneSTwoS is the 1S–2S frequency of the network's levels [Hz].
func (sp Spectrum) OneSTwoS() float64 { return (sp.Extrap[1] - sp.Extrap[0]) * HaJ / H }
