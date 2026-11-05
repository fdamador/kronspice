// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render.py: the resonance curves, ring-downs and the
// zero-point run, all as kronspice computed them (run with -viz viz).

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Sweep is an AC sweep of the E-ring's loop current over f0·(1 ± half).
func (e ERing) Sweep(half float64, pts int) (freq, imag []float64, err error) {
	ckt, err := e.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, 0)
	if err != nil {
		return nil, nil, err
	}
	f0 := e.Omega0() / (2 * math.Pi)
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: pts, FStart: f0 * (1 - half), FStop: f0 * (1 + half)}, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	i, _ := res.ElementCurrent("LS0")
	for k := range i {
		imag = append(imag, cmplx.Abs(i[k]))
	}
	return res.Freq, imag, nil
}

// Spectrum is the dual ring's differential impedance |Z| at (a0, b0).
func (d DualRing) Spectrum(fStart, fStop float64, pts int) (freq, z []float64, err error) {
	ckt, err := d.Build()
	if err != nil {
		return nil, nil, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: pts, FStart: fStart, FStop: fStop}, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	for k := range res.Freq {
		z = append(z, cmplx.Abs(va[k]-vb[k]))
	}
	return res.Freq, z, nil
}

// Trace is one zero-point realization: scaled time, loop current and
// stored energy (in ħ), every `keep`-th point.
func (z ZPFRing) Trace(keep int) (t, i, e []float64, err error) {
	rng := rand.New(rand.NewSource(z.Seed))
	v := z.noise(rng, 1/z.Q)
	ramp := z.Ramp * z.Q
	for k := range v {
		if tt := float64(k) * z.Dt; tt < ramp {
			v[k] *= 0.5 * (1 - math.Cos(math.Pi*tt/ramp))
		}
	}
	ckt, err := z.Build(v)
	if err != nil {
		return nil, nil, nil, err
	}
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: z.Step, Stop: float64(z.Pts-1) * z.Dt, MaxStep: z.Step}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return nil, nil, nil, err
	}
	il, _ := wf.ElementCurrent("LS0")
	for k := 0; k < len(wf.Time); k += keep {
		t, i, e = append(t, wf.Time[k]), append(i, il[k]), append(e, wf.Ledger.StoredEnergy[k])
	}
	return t, i, e, nil
}

func writeCSV(dir, name, header string, cols ...[]float64) error {
	f, err := os.Create(filepath.Join(dir, name))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, header)
	for k := range cols[0] {
		for j, c := range cols {
			if j > 0 {
				fmt.Fprint(f, ",")
			}
			fmt.Fprintf(f, "%.9g", c[k])
		}
		fmt.Fprintln(f)
	}
	return nil
}

// writeViz writes every CSV render.py reads.
func writeViz(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	// 1. E-ring resonance curves at four permittivities, and Q against εr.
	var ers, qf, qac []float64
	for _, er := range []float64{30, 100, 300, 1000} {
		e := DefaultERing()
		e.EpsR = er
		f, i, err := e.Sweep(0.2, 8001)
		if err != nil {
			return err
		}
		f0 := e.Omega0() / (2 * math.Pi)
		rel := make([]float64, len(f))
		peak := 0.0
		for k := range i {
			peak = math.Max(peak, i[k])
		}
		norm := make([]float64, len(i))
		for k := range f {
			rel[k] = f[k]/f0 - 1
			norm[k] = i[k] / peak
		}
		if err := writeCSV(dir, fmt.Sprintf("ering_er%g.csv", er), "f_over_f0_minus_1,current_over_peak", rel, norm); err != nil {
			return err
		}
		ac, err := e.Resonance()
		if err != nil {
			return err
		}
		ers, qf, qac = append(ers, er), append(qf, e.Q()), append(qac, ac.Q)
	}
	if err := writeCSV(dir, "ering_q.csv", "eps_r,q_formula,q_kronspice_ac", ers, qf, qac); err != nil {
		return err
	}
	// 2. The inverted ring lines' spectra.
	for _, crossed := range []bool{false, true} {
		d := DefaultDualRing(crossed)
		f, z, err := d.Spectrum(10e6, 200e6, 19001)
		if err != nil {
			return err
		}
		peaks, _, err := d.Modes(4, 10e6, 200e6)
		if err != nil {
			return err
		}
		name := map[bool]string{false: "plain", true: "crossed"}[crossed]
		if err := writeCSV(dir, "dual_"+name+".csv", "freq_hz,z_ohm", f, z); err != nil {
			return err
		}
		if err := writeCSV(dir, "dual_"+name+"_modes.csv", "freq_hz", peaks); err != nil {
			return err
		}
	}
	// 3. Ring-downs at εr 100 (Q 40) and 1000 (Q 1140), time in periods.
	for _, c := range []struct {
		er     float64
		cycles int
	}{{100, 60}, {1000, 1500}} {
		e := DefaultERing()
		e.EpsR = c.er
		rd, err := e.RingDown(c.cycles)
		if err != nil {
			return err
		}
		T := 2 * math.Pi / e.Omega0()
		var tp, in []float64
		peak := 0.0
		for _, x := range rd.I {
			peak = math.Max(peak, math.Abs(x))
		}
		for k := 0; k < len(rd.Time); k += 2 {
			tp, in = append(tp, rd.Time[k]/T), append(in, rd.I[k]/peak)
		}
		if err := writeCSV(dir, fmt.Sprintf("ringdown_er%g.csv", c.er), "time_periods,current_over_peak", tp, in); err != nil {
			return err
		}
	}
	// 4. One zero-point realization of the bench-Q ring (scaled: ω0 = 1).
	t, i, en, err := DefaultZPFRing().Trace(3)
	if err != nil {
		return err
	}
	return writeCSV(dir, "zpf_trace.csv", "time_s_scaled,current,stored_energy_hbar", t, i, en)
}
