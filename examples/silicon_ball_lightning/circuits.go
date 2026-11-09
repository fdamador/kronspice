// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// kronspice circuits.
//
// 1. The welding arc. Paiva et al. (PRL 98, 048501, 2007) made luminous
//    silicon balls with a 60 Hz arc at 20–25 V and 100–140 A: the
//    operating point of an AC stick welder. A welder is a transformer
//    whose leakage reactance sets the current: an open-circuit voltage
//    V_oc behind a series inductance L, feeding the arc. The arc is an
//    S-type element whose voltage hardly depends on its current
//    (v ≈ V_arc·sign(i)), entered as a current-controlled PWL element
//    (circuit.Builder.NonlinearCC). L is set so the arc current is
//    120 A rms; the ledger gives the arc power, the energy the arc spends
//    making vapor.
// 2. The power line. Cen et al. saw the ball's light oscillate at 99.4 Hz
//    with a 35 kV line 20 m away. The line couples to a ball of radius a at
//    distance d through a mutual capacitance C_m ≈ 4πε0·a²/d (a ≪ d), and
//    the ball has its own capacitance 4πε0·a to ground; any power the line
//    drives into the ball's resistance R goes as i², so it pulses at twice
//    the line frequency. kronspice gives the power for each R and its
//    largest value, against the ≈ 100 W the ball shines.

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

const eps0 = 8.8541878128e-12

// Weld is the welding-arc run.
type Weld struct {
	Voc, Varc, F   float64 // open-circuit voltage [V rms], arc voltage [V], frequency [Hz]
	L              float64 // leakage inductance [H]
	Irms, Vrms     float64 // arc current and voltage, rms
	Power          float64 // mean arc power [W]
	LedgerMismatch float64 // (source energy − arc energy − Δstored) ÷ source energy
}

func arcCurve(varc float64) []float64 {
	// v ≈ V_arc·sign(i), with a 0.5 A linear region through zero and a
	// slight positive slope (5 mΩ) above it
	return []float64{-1e4, -varc - 50, -0.5, -varc, 0.5, varc, 1e4, varc + 50}
}

// RunWeld simulates 20 cycles of the welder with leakage inductance L and
// measures the last 10.
func RunWeld(voc, varc, f, L float64) (Weld, error) {
	w := Weld{Voc: voc, Varc: varc, F: f, L: L}
	b := circuit.NewBuilder("AC stick welder and arc")
	b.VoltageWave("VW", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, voc * math.Sqrt2, f}})
	b.Inductor("LW", "s", "a", L)
	b.NonlinearCC("ARC", "a", "0", arcCurve(varc))
	ckt, err := b.BuildTransient()
	if err != nil {
		return w, err
	}
	T := 1 / f
	step := T / 400
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: 20 * T, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return w, fmt.Errorf("kronspice: %w", err)
	}
	ia, _ := wf.ElementCurrent("Vs_ARC")
	va, _ := wf.NodeVoltage("a")
	var s2i, s2v, sp, n float64
	for k, t := range wf.Time {
		if t < 10*T || k == 0 {
			continue
		}
		dt := t - wf.Time[k-1]
		s2i += ia[k] * ia[k] * dt
		s2v += va[k] * va[k] * dt
		sp += va[k] * ia[k] * dt
		n += dt
	}
	w.Irms, w.Vrms, w.Power = math.Sqrt(s2i/n), math.Sqrt(s2v/n), math.Abs(sp/n)
	l := wf.Ledger
	w.LedgerMismatch = l.IntegrationMismatch / math.Max(l.EnergyIn, 1e-30)
	return w, nil
}

// TuneWeld finds the leakage inductance for an arc current of irms (bisection
// in log L) and returns that run.
func TuneWeld(voc, varc, f, irms float64) (Weld, error) {
	lo, hi := 1e-5, 1e-1
	var w Weld
	var err error
	for it := 0; it < 40; it++ {
		m := math.Sqrt(lo * hi)
		w, err = RunWeld(voc, varc, f, m)
		if err != nil {
			return w, err
		}
		if w.Irms > irms {
			lo = m
		} else {
			hi = m
		}
		if hi/lo < 1+1e-4 {
			break
		}
	}
	return w, nil
}

// Line is the power-line coupling to a ball.
type Line struct {
	Vphase, F, A, D  float64 // phase voltage to ground [V rms], frequency, ball radius, distance [m]
	Cm, Cb           float64 // mutual and ball capacitance [F]
	BestR, BestPower float64 // the resistance that takes most power, and that power [W]
	RippleHz         float64 // frequency of the power pulsation from the transient [Hz]
}

// LinePower is the mean power into resistance R (kronspice AC, peak phasors).
func LinePower(l Line, R float64) (float64, error) {
	b := circuit.NewBuilder("line to ball")
	b.Add(circuit.Element{Name: "VL", Kind: circuit.VoltageSource, ACMag: l.Vphase * math.Sqrt2}, "line", "0")
	b.Capacitor("CM", "line", "ball", l.Cm)
	b.Capacitor("CB", "ball", "0", l.Cb)
	b.Resistor("RB", "ball", "0", R)
	ckt, err := b.Build()
	if err != nil {
		return 0, err
	}
	r, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: l.F, FStop: l.F}, analysis.Options{})
	if err != nil {
		return 0, err
	}
	v, _ := r.NodeVoltage("ball")
	m := real(v[0])*real(v[0]) + imag(v[0])*imag(v[0])
	return 0.5 * m / R, nil
}

// RunLine sweeps R over 1 Ω – 1 TΩ for the largest power, then runs a
// transient at that R and reads the power pulsation's frequency from its
// peaks.
func RunLine(vphase, f, a, d float64) (Line, error) {
	l := Line{Vphase: vphase, F: f, A: a, D: d}
	l.Cm = 4 * math.Pi * eps0 * a * a / d
	l.Cb = 4 * math.Pi * eps0 * a
	for e := 0.0; e <= 12; e += 0.05 {
		R := math.Pow(10, e)
		p, err := LinePower(l, R)
		if err != nil {
			return l, err
		}
		if p > l.BestPower {
			l.BestPower, l.BestR = p, R
		}
	}
	b := circuit.NewBuilder("line to ball, transient")
	b.VoltageWave("VL", "line", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vphase * math.Sqrt2, f}})
	b.Capacitor("CM", "line", "ball", l.Cm)
	b.Capacitor("CB", "ball", "0", l.Cb)
	b.Resistor("RB", "ball", "0", l.BestR)
	ckt, err := b.BuildTransient()
	if err != nil {
		return l, err
	}
	T := 1 / f
	step := T / 400
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: 12 * T, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return l, fmt.Errorf("kronspice: %w", err)
	}
	v, _ := wf.NodeVoltage("ball")
	var peaks []float64
	for k := 1; k+1 < len(v); k++ {
		if wf.Time[k] < 2*T {
			continue
		}
		p0, p1, p2 := v[k-1]*v[k-1], v[k]*v[k], v[k+1]*v[k+1]
		if p1 > p0 && p1 >= p2 {
			peaks = append(peaks, wf.Time[k])
		}
	}
	if len(peaks) > 2 {
		l.RippleHz = float64(len(peaks)-1) / (peaks[len(peaks)-1] - peaks[0])
	}
	return l, nil
}
