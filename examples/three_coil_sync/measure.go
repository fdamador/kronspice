// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"strings"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/output"
)

// Result is one design's steady-state behaviour over the measuring window
// [tMeasure, tStop].
type Result struct {
	Design Design
	E      Electrical
	FR0    float64 // unloaded Royer tank frequency [Hz]

	POut          float64 // mean power into the load [W]
	PInHV, PInRoy float64 // mean power from the HV generator and the Royer supplies [W]
	Eta           float64 // POut / (PInHV + PInRoy)
	FRoyer, FOut  float64 // measured frequencies of the tank voltage and the load current [Hz]
	Locked        bool    // Royer within 0.2 % of the HV generator
	Beat          float64 // (max − min)/mean of the per-cycle output power over the window
	VdsPk, IdPk   float64 // MOSFET peak drain voltage and current
	VHVPk, VC3Pk  float64 // peak voltage at coil 1's HV terminal and across C3
	VOutPk        float64 // peak load voltage
	VgPk          float64 // peak HV generator (amplifier) output voltage
	I3rms, I2rms  float64 // rms current in coil 3 and in one coil-2 strand
	IChoke        float64 // mean supply current
	PDiss         map[string]float64
	Residual      float64 // ledger: max relative residual
	Mismatch      float64 // ledger: integration mismatch / energy in
	Steps         int
	Wave          *output.Waveforms
	Err           error
}

// Simulate runs the transient and measures it. keep retains the waveforms.
func Simulate(ctx context.Context, d Design, keep bool) Result {
	return SimulateStep(ctx, d, keep, maxStep)
}

// SimulateStep is Simulate with a chosen largest time step [s].
func SimulateStep(ctx context.Context, d Design, keep bool, step float64) Result {
	r := Result{Design: d, FR0: d.RoyerFreq0()}
	c, e, err := d.Build()
	r.E = e
	if err != nil {
		r.Err = err
		return r
	}
	if _, ok := ctx.Deadline(); !ok {
		var cancel context.CancelFunc
		ctx, cancel = context.WithTimeout(ctx, 120*time.Second)
		defer cancel()
	}
	opt := analysis.TranOptions{Options: analysis.Options{Context: ctx}}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: step, Stop: tStop, MaxStep: step}, opt)
	if err != nil {
		r.Err = err
		return r
	}
	if err := r.measure(wf); err != nil {
		r.Err = err
		return r
	}
	if keep {
		r.Wave = wf
	}
	return r
}

func elem(wf *output.Waveforms, name string) int {
	for i, n := range wf.ElementNames {
		if strings.EqualFold(n, name) {
			return i
		}
	}
	return -1
}

func node(wf *output.Waveforms, name string) []float64 {
	v, ok := wf.NodeVoltage(name)
	if !ok {
		return make([]float64, len(wf.Time))
	}
	return v
}

// window returns the first index at or after t.
func window(t []float64, from float64) int {
	for i, x := range t {
		if x >= from {
			return i
		}
	}
	return len(t)
}

// mean is the time average of y over [i0, end) by the trapezoidal rule.
func mean(t, y []float64, i0 int) float64 {
	s := 0.0
	for k := i0 + 1; k < len(t); k++ {
		s += 0.5 * (y[k] + y[k-1]) * (t[k] - t[k-1])
	}
	return s / (t[len(t)-1] - t[i0])
}

func rms(t, y []float64, i0 int) float64 {
	sq := make([]float64, len(y))
	for k, v := range y {
		sq[k] = v * v
	}
	return math.Sqrt(mean(t, sq, i0))
}

func peak(y []float64, i0 int) float64 {
	m := 0.0
	for _, v := range y[i0:] {
		m = math.Max(m, math.Abs(v))
	}
	return m
}

// crossings returns the times of rising zero crossings of y after i0
// (linear interpolation).
func crossings(t, y []float64, i0 int) []float64 {
	var c []float64
	for k := i0 + 1; k < len(t); k++ {
		if y[k-1] < 0 && y[k] >= 0 {
			c = append(c, t[k-1]+(t[k]-t[k-1])*(-y[k-1])/(y[k]-y[k-1]))
		}
	}
	return c
}

// freq is the mean frequency from rising zero crossings, 0 if fewer than 3.
func freq(c []float64) float64 {
	if len(c) < 3 {
		return 0
	}
	return float64(len(c)-1) / (c[len(c)-1] - c[0])
}

func (r *Result) measure(wf *output.Waveforms) error {
	t := wf.Time
	i0 := window(t, tMeasure)
	if i0 >= len(t)-10 {
		return fmt.Errorf("%w: measuring window empty", analysis.ErrNoMeasurement)
	}
	P := func(name string) []float64 {
		if k := elem(wf, name); k >= 0 {
			return wf.P[k]
		}
		return make([]float64, len(t))
	}
	I := func(name string) []float64 {
		if k := elem(wf, name); k >= 0 {
			return wf.I[k]
		}
		return make([]float64, len(t))
	}
	r.POut = mean(t, P("RL"), i0)
	r.PInHV = -mean(t, P("Vg"), i0)
	if r.Design.Sync {
		r.PInHV = -mean(t, P("Eg"), i0) - mean(t, P("Es"), i0)
	}
	r.PInRoy = -mean(t, P("Vcc"), i0) - mean(t, P("Vgate"), i0)
	if in := r.PInHV + r.PInRoy; in > 0 {
		r.Eta = r.POut / in
	}
	r.PDiss = map[string]float64{}
	for _, n := range []string{"Rg", "R1p", "R1s", "R3", "R2a", "R2b", "RC2", "Rs", "M1", "M2", "D1", "D2", "RG1", "RG2"} {
		r.PDiss[n] = mean(t, P(n), i0)
	}

	c1, c2 := node(wf, "c1"), node(wf, "c2")
	tank := make([]float64, len(t))
	for k := range t {
		tank[k] = c1[k] - c2[k]
	}
	iL := I("RL")
	r.FRoyer = freq(crossings(t, tank, i0))
	r.FOut = freq(crossings(t, iL, i0))
	if r.Design.Sync {
		r.Design.FHV = r.FRoyer
	}
	if r.Design.FHV > 0 && r.FRoyer > 0 && !r.Design.HVOff && !r.Design.RoyerOff {
		r.Locked = math.Abs(r.FRoyer-r.Design.FHV)/r.Design.FHV < 2e-3
	}
	// Per-cycle output power (cycles of the load current).
	if cr := crossings(t, iL, i0); len(cr) >= 3 {
		p := P("RL")
		lo, hi, sum := math.Inf(1), math.Inf(-1), 0.0
		n := 0
		for j := 1; j < len(cr); j++ {
			a, b := window(t, cr[j-1]), window(t, cr[j])
			if b-a < 4 {
				continue
			}
			s := 0.0
			for k := a + 1; k < b; k++ {
				s += 0.5 * (p[k] + p[k-1]) * (t[k] - t[k-1])
			}
			s /= t[b-1] - t[a]
			lo, hi, sum, n = math.Min(lo, s), math.Max(hi, s), sum+s, n+1
		}
		if n > 0 && sum > 0 {
			r.Beat = (hi - lo) / (sum / float64(n))
		}
	}
	r.VdsPk = math.Max(peak(c1, i0), peak(c2, i0))
	r.IdPk = math.Max(peak(I("M1"), i0), peak(I("M2"), i0))
	r.VHVPk = peak(node(wf, "hv"), i0)
	vb3, vout := node(wf, "b3"), node(wf, "out")
	vc3 := make([]float64, len(t))
	for k := range t {
		vc3[k] = vb3[k] - vout[k]
	}
	r.VC3Pk = peak(vc3, i0)
	r.VOutPk = peak(vout, i0)
	r.VgPk = peak(node(wf, "hvg"), i0)
	r.I3rms = rms(t, I("L3"), i0)
	r.I2rms = rms(t, I("L2a"), i0)
	r.IChoke = mean(t, I("Lchoke"), i0)
	r.Residual = wf.Ledger.MaxRelativeResidual
	if wf.Ledger.EnergyIn > 0 {
		r.Mismatch = wf.Ledger.IntegrationMismatch / wf.Ledger.EnergyIn
	}
	r.Steps = wf.Stats.Accepted
	return nil
}
