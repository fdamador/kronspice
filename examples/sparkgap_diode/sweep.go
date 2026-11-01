// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// Frequency sweep: every rectifying device in the same bridge, so the
// devices can be compared on one chart of amps per volt (mean load current
// ÷ source peak voltage) against frequency.
//
// Common test circuit (declared): full bridge, source resistance 0.1 Ω,
// load 1 Ω, no filter capacitor, 30 cycles of which the first 20 settle.
// A perfect diode bridge would give (2/π)/1.1 Ω = 0.579 A/V at every
// frequency; that line is the reference.

// SweepFreqs are the frequencies swept [Hz].
var SweepFreqs = []float64{50, 100, 300, 1e3, 3e3, 10e3, 30e3, 100e3}

// SweepRow is one point of the sweep.
type SweepRow struct {
	F, Vp      float64 // frequency [Hz], source peak [V]
	AperV      float64 // mean load current ÷ source peak [A/V]
	Efficiency float64
	Forward    int // spark-gap strikes forward / backward (gaps only)
	Backward   int
}

// SweepDevice is a device swept across frequency.
type SweepDevice struct {
	Name string
	Vp   float64
	Rows []SweepRow
}

// IdealAperV is the perfect bridge's amps per volt in the sweep circuit.
func IdealAperV() float64 { return 2 / math.Pi / 1.1 }

func sweepRect(vp, f float64) Rectifier {
	return Rectifier{Name: "sweep", Bridge: true, Vp: vp, F: f, Rs: 0.1, Rload: 1, Cycles: 30, Skip: 20}
}

// SweepCell sweeps a bridge of four cells (one per arm) at peak vp.
func SweepCell(c Cell, vp float64) (SweepDevice, error) {
	d := SweepDevice{Name: c.Name(), Vp: vp}
	for _, f := range SweepFreqs {
		r := sweepRect(vp, f)
		r.Arm = Arm{D: c, N: 1}
		res, err := r.Run()
		if err != nil {
			return d, fmt.Errorf("%s at %g Hz: %w", c.Name(), f, err)
		}
		d.Rows = append(d.Rows, SweepRow{F: f, Vp: vp, AperV: res.Idc / vp, Efficiency: res.Efficiency})
	}
	return d, nil
}

// SweepGap sweeps a bridge of four polarity gaps with the breakdown set for
// each frequency (AtFrequency, extrapolated relation). With vp > 0 the peak
// is fixed; with vp = 0 it is set at each frequency to the middle of the
// window where a bridge can work, between 2·Vpos (two gaps in series
// forward) and Vneg (each blocked gap holds the source): (2 + ratio)/2·Vpos.
// Amps per volt is then reported against that frequency's peak.
func SweepGap(g PolarGap, vp float64) SweepDevice {
	name := fmt.Sprintf("water polarity gap (ratio %.0f), %.0f kV fixed", g.NegRatio, vp/1e3)
	if vp == 0 {
		name = fmt.Sprintf("water polarity gap (ratio %.0f), voltage raised with frequency", g.NegRatio)
	}
	d := SweepDevice{Name: name, Vp: vp}
	for _, f := range SweepFreqs {
		gf := g.AtFrequency(f)
		v := vp
		if v == 0 {
			v = (2 + g.NegRatio) / 2 * gf.Vpos()
		}
		r := sweepRect(v, f)
		res := r.RunGaps(gf)
		d.Rows = append(d.Rows, SweepRow{F: f, Vp: v, AperV: res.Idc / v, Efficiency: res.Efficiency, Forward: res.Forward, Backward: res.Backward})
		continue
	}
	return d
}

// SweepGapDefault is the water gap used in the sweep: polarity ratio 3
// (the ratio that rectified in Part 3), hot-channel recovery 1 ms
// (declared).
func SweepGapDefault() PolarGap {
	g := DefaultGap()
	g.NegRatio = 3
	g.Recovery = 1e-3
	return g
}

// Sweep runs every device.
func Sweep() ([]SweepDevice, error) {
	var out []SweepDevice
	for _, a := range []Electrode{Pb, FeAgCoat, CuAgTips} {
		d, err := SweepCell(DefaultCell(a, Al), 100)
		if err != nil {
			return nil, err
		}
		out = append(out, d)
	}
	out = append(out, SweepGap(SweepGapDefault(), 6e3), SweepGap(SweepGapDefault(), 0))
	return out, nil
}
