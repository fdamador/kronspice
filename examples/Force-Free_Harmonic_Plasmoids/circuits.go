// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// kronspice analysis of the plasmoid.
//
// 1. Spectrum. A plasma sphere of fixed radius and density has a set of
//    l = 1 leaky resonances: TE overtones (interior argument near the
//    zeros of j1; the first is Nachamkin's αR = 4.4934) and TM overtones
//    (zeros of d(x·j1)/dx). Each is a parallel RLC tank (Foster form):
//    L·C = 1/ω², R = Q/(ωC). The tanks in series, driven by an AC current
//    source, give the plasmoid's harmonic spectrum; kronspice AC finds each
//    peak and its half-power width, which must return the mode's ω and Q.
// 2. Ring-down. One tank released with stored energy: kronspice transient
//    gives the energy decay rate ω/Q and the ledger (heat = energy lost).
// 3. His eq. 49 as a circuit. A one-dimensional plasma line (Kron's field
//    network: series L = μ0Δ, shunt C = ε0Δ) loaded at every node by the
//    electrons' kinetic inductance L_k/Δ, L_k = m/(n·e²), shorted at both
//    ends a half wavelength π/α apart, resonates at ω² = ωp² + α²c² (to
//    the lattice's (αΔ)² error): Nachamkin's charge-density condition is
//    the resonance of a plasma-loaded line.
// Frequencies are scaled by a reference ω (normalized circuit units) so
// that element values stay of order one; results are scaled back.

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// TMZeros returns the first count zeros of d(x·j1(x))/dx.
func TMZeros(count int) []float64 {
	f := func(x float64) float64 { return j1(x) + x*j1p(x) }
	var z []float64
	step := 0.01
	for x := 0.5; len(z) < count; x += step {
		if f(x)*f(x+step) <= 0 {
			lo, hi := x, x+step
			for i := 0; i < 100; i++ {
				m := 0.5 * (lo + hi)
				if f(lo)*f(m) <= 0 {
					hi = m
				} else {
					lo = m
				}
			}
			z = append(z, 0.5*(lo+hi))
		}
	}
	return z
}

// Spectrum returns the first n TE and n TM l = 1 resonances of a plasma
// sphere (radius R, plasma frequency wp), TE first.
func Spectrum(R, wp float64, n int) []Mode {
	s := Sphere{R: R, Wp: wp}
	var out []Mode
	for _, x := range J1Zeros(float64(n+1) * math.Pi)[:n] {
		out = append(out, s.FindModeFrom(true, x))
	}
	for _, x := range TMZeros(n) {
		out = append(out, s.FindModeFrom(false, x))
	}
	return out
}

// Peak is one resonance read back from the kronspice AC sweep.
type Peak struct {
	Mode           Mode
	FPeak, FMode   float64 // [Hz]
	QKron, QMode   float64
	FErr, QErr     float64 // relative
	ImpedanceAtTop float64 // |V/I| at the peak ÷ the tank's R
}

// SpectrumAC builds the Foster chain of the modes and sweeps each peak in
// kronspice (2001 points over ±6 half-widths).
func SpectrumAC(modes []Mode) ([]Peak, error) {
	ref := real(modes[0].Omega)
	b := circuit.NewBuilder("plasmoid harmonic spectrum (Foster chain)")
	b.Add(circuit.Element{Name: "IGEN", Kind: circuit.CurrentSource, ACMag: 1}, "0", "n0")
	for i, m := range modes {
		w := real(m.Omega) / ref
		a, c := fmt.Sprintf("n%d", i), fmt.Sprintf("n%d", i+1)
		if i == len(modes)-1 {
			c = "0"
		}
		b.Capacitor(fmt.Sprintf("C%d", i), a, c, 1)
		b.Inductor(fmt.Sprintf("L%d", i), a, c, 1/(w*w))
		b.Resistor(fmt.Sprintf("R%d", i), a, c, m.Q/w)
	}
	ckt, err := b.Build()
	if err != nil {
		return nil, err
	}
	var out []Peak
	for i, m := range modes {
		w := real(m.Omega) / ref
		f0 := w / (2 * math.Pi)
		hw := f0 / (2 * m.Q)
		r, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 2001, FStart: f0 - 6*hw, FStop: f0 + 6*hw}, analysis.Options{})
		if err != nil {
			return nil, err
		}
		a, c := fmt.Sprintf("n%d", i), fmt.Sprintf("n%d", i+1)
		va, _ := r.NodeVoltage(a)
		var vc []complex128
		if i == len(modes)-1 {
			vc = make([]complex128, len(va))
		} else {
			vc, _ = r.NodeVoltage(c)
		}
		mag := make([]float64, len(va))
		best := 0
		for k := range va {
			mag[k] = cmplx.Abs(va[k] - vc[k])
			if mag[k] > mag[best] {
				best = k
			}
		}
		fp := parabolaPeak(r.Freq, mag, best)
		half := mag[best] / math.Sqrt2
		lo := crossing(r.Freq, mag, best, -1, half)
		hi := crossing(r.Freq, mag, best, +1, half)
		p := Peak{Mode: m, FPeak: fp * ref, FMode: real(m.Omega) / (2 * math.Pi), QMode: m.Q}
		p.QKron = fp / (hi - lo)
		p.FErr = math.Abs(p.FPeak/p.FMode - 1)
		p.QErr = math.Abs(p.QKron/p.QMode - 1)
		p.ImpedanceAtTop = mag[best] / (m.Q / w)
		out = append(out, p)
	}
	return out, nil
}

func parabolaPeak(f, y []float64, k int) float64 {
	if k == 0 || k == len(y)-1 {
		return f[k]
	}
	d := y[k-1] - 2*y[k] + y[k+1]
	if d == 0 {
		return f[k]
	}
	return f[k] + 0.5*(y[k-1]-y[k+1])/d*(f[1]-f[0])
}

func crossing(f, y []float64, k, dir int, level float64) float64 {
	for i := k; i+dir >= 0 && i+dir < len(y); i += dir {
		if y[i+dir] < level {
			t := (y[i] - level) / (y[i] - y[i+dir])
			return f[i] + t*(f[i+dir]-f[i])
		}
	}
	return f[k]
}

// RingDown releases one mode's tank (normalized: ω = 1, C = 1, R = Q)
// from 1 V and runs periods cycles in kronspice. It returns Q from the
// fitted energy decay and the ledger mismatch (heat vs energy lost).
type RingDown struct {
	Q, QFit, Mismatch, Heat, Lost float64
}

func RunRingDown(m Mode, periods int) (RingDown, error) {
	rd := RingDown{Q: m.Q}
	ckt, err := circuit.NewBuilder("plasmoid mode ring-down").
		CapacitorIC("C1", "a", "0", 1, 1).
		Inductor("L1", "a", "0", 1).
		Resistor("R1", "a", "0", m.Q).BuildTransient()
	if err != nil {
		return rd, err
	}
	T := 2 * math.Pi
	step := T / 256
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: float64(periods) * T, MaxStep: step, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return rd, fmt.Errorf("kronspice: %w", err)
	}
	l := wf.Ledger
	W := l.StoredEnergy
	// least-squares slope of ln W over the run (oscillation averages out
	// over whole periods: ½Cv² + ½Li² is exactly the stored energy)
	var sx, sy, sxx, sxy, n float64
	for k, t := range wf.Time {
		y := math.Log(W[k])
		sx += t
		sy += y
		sxx += t * t
		sxy += t * y
		n++
	}
	slope := (n*sxy - sx*sy) / (n*sxx - sx*sx)
	rd.QFit = -1 / slope // energy rate = ω/Q with ω = 1
	rd.Heat, rd.Lost = l.EnergyDissipated, l.StoredStart-l.StoredEnd
	rd.Mismatch = math.Abs(rd.Heat-rd.Lost) / l.StoredStart
	return rd, nil
}

// Ladder is the plasma-line twin of eq. 49.
type Ladder struct {
	KOverAlpha              float64 // Nachamkin's k/α
	N                       int     // cells over the half wavelength π/α
	WKron, WPaper, WLattice float64 // resonance ÷ αc: kronspice, eq. 49, exact lattice
}

// RunLadder builds the line in units α = c = μ0 = ε0 = 1 (length π, N
// cells), sets the electrons' kinetic inductance from eq. 49 for the
// given k/α, and finds the lowest resonance with kronspice AC.
func RunLadder(kOverAlpha float64, N int) (Ladder, error) {
	ld := Ladder{KOverAlpha: kOverAlpha, N: N, WPaper: kOverAlpha}
	wp2 := kOverAlpha*kOverAlpha - 1 // ωp² from eq. 49 (α = c = 1)
	dx := math.Pi / float64(N)
	lk := 1 / wp2 // m/(n e²) in these units
	b := circuit.NewBuilder("plasma-loaded line (eq. 49 twin)")
	node := func(i int) string {
		if i == 0 || i == N {
			return "0"
		}
		return fmt.Sprintf("x%d", i)
	}
	// a 1 nΩ resistor in series with every inductor gives the DC operating
	// point a unique solution (inductor-only loops), and no visible damping
	for i := 0; i < N; i++ {
		m := fmt.Sprintf("s%d", i)
		b.Inductor(fmt.Sprintf("LS%d", i), node(i), m, dx)
		b.Resistor(fmt.Sprintf("RS%d", i), m, node(i+1), 1e-9)
	}
	for i := 1; i < N; i++ {
		m := fmt.Sprintf("p%d", i)
		b.Capacitor(fmt.Sprintf("CS%d", i), node(i), "0", dx)
		b.Inductor(fmt.Sprintf("LP%d", i), node(i), m, lk/dx)
		b.Resistor(fmt.Sprintf("RP%d", i), m, "0", 1e-9)
	}
	// light damping so the resonance has a finite peak
	b.Resistor("RD", node(N/2), "0", 1e6)
	b.Add(circuit.Element{Name: "IGEN", Kind: circuit.CurrentSource, ACMag: 1}, "0", node(N/2))
	ckt, err := b.Build()
	if err != nil {
		return ld, err
	}
	sd := 2 / dx * math.Sin(dx/2)
	ld.WLattice = math.Sqrt(wp2 + sd*sd)
	f0 := ld.WLattice / (2 * math.Pi)
	span := 1e-3 * f0
	for it := 0; it < 3; it++ {
		r, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 801, FStart: f0 - span, FStop: f0 + span}, analysis.Options{})
		if err != nil {
			return ld, err
		}
		v, _ := r.NodeVoltage(node(N / 2))
		best := 0
		mag := make([]float64, len(v))
		for k := range v {
			mag[k] = cmplx.Abs(v[k])
			if mag[k] > mag[best] {
				best = k
			}
		}
		f0 = parabolaPeak(r.Freq, mag, best)
		span /= 50
	}
	ld.WKron = 2 * math.Pi * f0
	return ld, nil
}
