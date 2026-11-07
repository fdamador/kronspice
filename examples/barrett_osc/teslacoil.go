// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

// A spark-gap Tesla coil whose top-load is Barrett's floating terminal and
// receives a weak signal, with an ordinary and a Tesla bifilar secondary.
//
// Secondary: the 400-turn solenoid of examples/pulsed_bifilar (0.5 mm
// copper, 10 cm across, 126 m) in its lumped equivalent: L = 4.853 mH for
// both windings, self-capacitance 3.84 pF ordinary (f1 1.167 MHz) and
// 3.014 nF bifilar (f1 41.6 kHz), wire resistance 15 Ω (ordinary, skin
// effect near 525 kHz) and 10.8 Ω (bifilar near 41.6 kHz), a 15 pF top-load
// (declared). With the top-load the secondary rings at 525 kHz (ordinary)
// and 41.5 kHz (bifilar).
//
// Primary (declared): 10 kV DC through 10 kΩ charges Cp = 20 nF; a rotary
// spark gap (dry air, 1 cm: holds 35 kV, so only its trigger fires it,
// 1000 breaks/s; arc 30 V + 0.2 Ω, hot channel recovering over 100 µs) closes
// Cp onto the primary coil Lp, tuned to each secondary (Lp = 1/ω²Cp), k = 0.15.
//
// Signal: 1 V from 50 Ω in series with the top-load's capacitance to ground
// (the Thevenin form of a small antenna), at 0.9 of the secondary's
// frequency, placed between the spark's lines (multiples of 1 kHz) so that
// the sidebands fs ± k·1 kHz the gap makes are on bins of their own.
//
// Runs per winding: sparks with the signal at phase 0 and at phase π,
// sparks alone, and the signal alone. Half the difference of the two
// phases keeps only what is linear in the signal (a small shift of the
// spark timing by the signal, which is second order, cancels); subtracting
// the signal-alone response leaves what the switching does to the signal:
// the multiplication.
type TeslaCoil struct {
	Bifilar        bool
	Ls, Rsec       float64
	Cself, Ctop    float64
	Vdc, Rch, Cp   float64
	K              float64
	Fb             float64 // break rate [Hz]
	Vs, Rsig       float64
	Fs             float64 // signal frequency [Hz]
	Dt             float64 // step [s]
	Warm, Window   float64 // [s]
	GapD, Recovery float64
	Dwell          float64 // rotary gap: electrodes face each other this long per firing [s]
	PhiS           float64 // signal phase [rad]
	// Timed: the rotary gap as a timed switch (arc resistance during the
	// dwell, open otherwise, no arc voltage): the coil is then a linear
	// periodically switched circuit and its mixing is exactly linear in the
	// signal. Otherwise the arc model of gap.go (strikes at thresholds,
	// re-strikes while hot), which makes every switching decision depend on
	// tiny differences.
	Timed bool
}

// Secondary constants from the pulsed_bifilar distributed model.
const (
	tcLs        = 4.8532e-3
	tcCselfOrd  = 3.8353e-12
	tcCselfBif  = 3.0143e-9
	tcRsecOrd   = 15.0
	tcRsecBif   = 10.8
	tcCtop      = 15e-12
	tcBreakRate = 1000.0
)

// NewTeslaCoil returns the declared coil for one winding.
func NewTeslaCoil(bifilar bool) TeslaCoil {
	t := TeslaCoil{Bifilar: bifilar, Ls: tcLs, Ctop: tcCtop, Vdc: 10e3, Rch: 10e3, Cp: 20e-9, K: 0.15,
		Fb: tcBreakRate, Vs: 1, Rsig: 50, Warm: 3e-3, Window: 10e-3, GapD: 0.01, Recovery: 100e-6, Dwell: 100e-6}
	if bifilar {
		t.Cself, t.Rsec, t.Dt = tcCselfBif, tcRsecBif, 50e-9
		t.Warm = 10e-3 // its ring-down lasts about a millisecond (Q ≈ 117 at 41.5 kHz)
	} else {
		t.Cself, t.Rsec, t.Dt = tcCselfOrd, tcRsecOrd, 10e-9
	}
	// Signal at 0.9·f_sec, on the 100 Hz grid of the window, between the
	// spark's 1 kHz lines (a whole number of 100 Hz with an odd tenth).
	f := 0.9 * t.SecondaryF()
	t.Fs = (math.Floor(f/1000) + 0.5) * 1000
	return t
}

// SecondaryF is the secondary's resonance with its top-load [Hz].
func (t TeslaCoil) SecondaryF() float64 {
	return 1 / (2 * math.Pi * math.Sqrt(t.Ls*(t.Cself+t.Ctop)))
}

// Lp is the primary inductance tuned to the secondary.
func (t TeslaCoil) Lp() float64 {
	w := 2 * math.Pi * t.SecondaryF()
	return 1 / (w * w * t.Cp)
}

// TCRun is one run.
type TCRun struct {
	T, Top   []float64 // window: time, top-load voltage
	Peak     float64   // largest |top-load voltage| in the window [V]
	Residual float64   // ledger [J]
	Supplied float64   // [J] over the run
	GapHeat  float64   // [J]
	Strikes  int
}

// Run steps the coil with sparks and/or the signal.
func (t TeslaCoil) Run(sparks, signal bool) (TCRun, error) {
	var r TCRun
	n := newTC()
	n.V("VDC", "hv", "0", func(float64) float64 { return t.Vdc }, t.Rch)
	n.C("CP", "hv", "0", t.Cp, t.Vdc) // starts charged
	gap := AirGap(t.GapD)
	gap.Recovery, gap.Hold = t.Recovery, 20e-9
	if t.Timed {
		gap.Varc = 0
	}
	g := n.Gap("GAP", "hv", "pa", gap)
	g.forced = t.Timed
	lp := n.L("LP", "pa", "pb", t.Lp())
	n.R("RP", "pb", "0", 0.05)
	ls := n.L("LS", "0", "sm", t.Ls)
	n.R("RSEC", "sm", "top", t.Rsec)
	n.C("CSELF", "0", "top", t.Cself, 0)
	n.C("CTOP", "top", "ant", t.Ctop, 0)
	amp := 0.0
	if signal {
		amp = t.Vs
	}
	n.V("VSIG", "ant", "0", func(x float64) float64 { return amp * math.Sin(2*math.Pi*t.Fs*x+t.PhiS) }, t.Rsig)
	n.Couple(lp, ls, t.K)
	if t.Timed && sparks {
		g.on, g.sgn = true, 1 // the first firing is at t = 0
	}
	if err := n.Prepare(); err != nil {
		return r, err
	}
	steps := int(math.Round((t.Warm + t.Window) / t.Dt))
	w0 := int(math.Round(t.Warm / t.Dt))
	nextFire, lastFire, prevI := 0.0, -1.0, 0.0
	for k := 0; k < steps; k++ {
		tt := float64(k) * t.Dt
		if sparks && tt >= nextFire-1e-15 {
			g.armed = true
			lastFire = nextFire
			nextFire += 1 / t.Fb
			if t.Timed {
				g.on, g.sgn = true, 1
				g.strk++
			}
		}
		// The rotor carries the electrodes apart after the dwell: the arc
		// goes out at its next current zero and the gap cannot re-strike
		// until the next firing.
		if lastFire >= 0 && tt > lastFire+t.Dwell && g.on {
			if !t.Timed || g.i*prevI <= 0 {
				g.on, g.tOff = false, -1
			}
		}
		prevI = g.i
		if err := n.Step(t.Dt); err != nil {
			return r, err
		}
		if k >= w0 {
			v := n.Volt("top")
			r.T = append(r.T, tt+t.Dt-t.Warm)
			r.Top = append(r.Top, v)
			r.Peak = math.Max(r.Peak, math.Abs(v))
		}
	}
	r.Residual, r.Supplied, r.GapHeat, r.Strikes = n.Residual(), n.supplied, n.gapHeat, g.strk
	return r, nil
}

// TCResult is one winding's three runs and the mixing.
type TCResult struct {
	Coil                       TeslaCoil
	Both, Flip, Sparks, Signal TCRun
	Coherent                   []float64          // ½(both(0) − both(π)): the part linear in the signal
	Mix                        []float64          // coherent − signal alone
	MixLines                   map[int]complex128 // sidebands fs ± k·Fb, k = −5…5
	SigLine                    float64            // |signal line| in the signal-only run [V]
	MixRMS, SigRMS, SparkPk    float64
}

// RunTesla runs one winding.
func RunTesla(bifilar, timed bool) (TCResult, error) {
	t := NewTeslaCoil(bifilar)
	t.Timed = timed
	var r TCResult
	r.Coil = t
	var err error
	if r.Both, err = t.Run(true, true); err != nil {
		return r, err
	}
	tf := t
	tf.PhiS = math.Pi
	if r.Flip, err = tf.Run(true, true); err != nil {
		return r, err
	}
	if r.Sparks, err = t.Run(true, false); err != nil {
		return r, err
	}
	if r.Signal, err = t.Run(false, true); err != nil {
		return r, err
	}
	r.Mix = make([]float64, len(r.Both.Top))
	r.Coherent = make([]float64, len(r.Both.Top))
	var sm, ss float64
	for j := range r.Mix {
		r.Coherent[j] = (r.Both.Top[j] - r.Flip.Top[j]) / 2
		r.Mix[j] = r.Coherent[j] - r.Signal.Top[j]
		sm += r.Mix[j] * r.Mix[j]
		ss += r.Signal.Top[j] * r.Signal.Top[j]
	}
	r.MixRMS, r.SigRMS = math.Sqrt(sm/float64(len(r.Mix))), math.Sqrt(ss/float64(len(r.Mix)))
	r.SparkPk = r.Sparks.Peak
	f0 := 1 / t.Window
	ks := int(math.Round(t.Fs / f0))
	kb := int(math.Round(t.Fb / f0))
	r.MixLines = map[int]complex128{}
	for k := -5; k <= 5; k++ {
		r.MixLines[k] = dftReal(r.Mix, ks+k*kb)
	}
	r.SigLine = 2 * cmplx.Abs(dftReal(r.Signal.Top, ks))
	return r, nil
}

func printTesla() ([]TCResult, error) {
	fmt.Println("Spark-gap Tesla coil, top-load = Barrett's floating terminal receiving a weak signal")
	fmt.Println("  (secondary from examples/pulsed_bifilar: 400 turns, L 4.85 mH; 15 pF top-load; 10 kV, 20 nF, rotary gap 1000/s, k 0.15)")
	var out []TCResult
	fmt.Println("  rotary gap as a timed switch (0.2 Ω for 100 µs per firing): a linear, periodically switched coil")
	for _, bif := range []bool{false, true} {
		r, err := RunTesla(bif, true)
		if err != nil {
			return nil, err
		}
		t := r.Coil
		name := "ordinary solenoid"
		if bif {
			name = "Tesla bifilar solenoid"
		}
		fmt.Printf("  %s: self-C %.3g pF, rings at %.1f kHz, primary Lp %.1f µH; signal %.1f kHz\n",
			name, t.Cself*1e12, t.SecondaryF()/1e3, t.Lp()*1e6, t.Fs/1e3)
		fmt.Printf("    top-load peak with sparks %.3g kV; signal alone at the top %.3g V\n", r.SparkPk/1e3, r.SigLine)
		fmt.Printf("    sparks × signal (coherent part − signal alone): rms %.3g V = %.2f × the signal's own rms at the top\n", r.MixRMS, r.MixRMS/r.SigRMS)
		fmt.Print("    sidebands fs + k·1 kHz [V], k = −3…3:")
		for k := -3; k <= 3; k++ {
			fmt.Printf(" %.2e", 2*cmplx.Abs(r.MixLines[k]))
		}
		fmt.Printf("\n    ledger residual %.1e J of %.3g J supplied; gap heat %.3g J; firings %d in %.0f ms\n",
			r.Both.Residual, r.Both.Supplied, r.Both.GapHeat, r.Both.Strikes, (t.Warm+t.Window)*1e3)
		out = append(out, r)
	}
	fmt.Println("  the gap switches the coil 1000 times a second, so the weak signal at the top-load comes out")
	fmt.Println("  multiplied by the switching: sidebands at fs ± k·1 kHz. This is real mixing, made by the gap (a")
	fmt.Println("  time-varying element), of the kind Barrett's Tesla experiments would have seen.")
	ra, err := RunTesla(false, false)
	if err != nil {
		return nil, err
	}
	fmt.Println("  with the arc model (strikes at thresholds, re-strikes while the channel is hot), ordinary secondary:")
	fmt.Printf("    phase-flip difference rms %.3g V (%.0f× the signal at the top); %d strikes in %.0f ms: the 1 V signal\n",
		ra.MixRMS, ra.MixRMS/ra.SigRMS, ra.Both.Strikes, (ra.Coil.Warm+ra.Coil.Window)*1e3)
	fmt.Println("    flips individual re-strikes in a ring-down of hundreds of kV, so the runs diverge. That is the")
	fmt.Println("    coil's sensitivity to any small disturbance; a real gap's random jitter (ns–µs) swamps it, so it")
	fmt.Println("    is not usable mixing.")
	fmt.Println()
	return out, nil
}
