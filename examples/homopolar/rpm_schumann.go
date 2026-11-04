// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// The shafts' rpm oscillating at the first Schumann frequency:
//
//	rpm(t) = 1768·(1 ± h·sin(2π·7.8·t)),
//
// 1768 rpm being a 4-pole induction motor at full load. "+" on both shafts
// is in phase; "+" on one and "−" on the other is antiphase. The speed drives
// (VFDs) impose the profile; the field's equations are then linear.
//
// For one self-excited loop at a held speed, L_s·İ = (−R + M_s·ω(t))·I, so
//
//	I(t) = I(0)·exp(∫(−R + M_s·ω)/L_s dt),
//
// and with ω = ω0·(1 + h·sin Ωt) the sine integrates to a bounded term:
//
//	I(t) ∝ exp((M_s·ω0 − R)/L_s · t) · exp(−β·cos Ωt),   β = M_s·ω0·h/(L_s·Ω).
//
// The field grows or decays at the MEAN speed's rate, whatever h; the
// oscillation only ripples it. exp(−β·cos θ) = I0(β) + 2·Σ (−1)ⁿ·Iₙ(β)·cos nθ,
// so the n-th harmonic relative to the DC is 2·Iₙ(β)/I0(β) (modified Bessel
// functions): an exact twin for the engine at the self-excitation speed.

package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

const (
	rpmBase  = 1768.0
	fSchu    = 7.8
	bSatIron = 1.6 // declared: a usable flux density for electrical steel [T]
)

// besselI is the modified Bessel function Iₙ(x) by its series.
func besselI(n int, x float64) float64 {
	var s, term float64
	term = math.Pow(x/2, float64(n))
	for k := 1; k <= n; k++ {
		term /= float64(k)
	}
	for k := 0; k < 60; k++ {
		s += term
		term *= (x / 2) * (x / 2) / (float64(k+1) * float64(n+k+1))
	}
	return s
}

// rpmMachine is the stacked machine with a load; gain is the iron factor.
func rpmMachine(gain float64) Dual {
	m := DefaultDual()
	m.Stacked, m.RLoad, m.Windage, m.IronGain = true, 50e-6, 5e-3, gain
	return m
}

// rpmRun holds both shafts on rpm(t) = 1768·(1 + h·sin) (shaft 2 with the
// given phase) and returns the RK4 trace.
func rpmRun(m Dual, h, phase2, tEnd float64) DualTrace {
	w0 := rpmBase * 2 * math.Pi / 60
	m.Speed = &[2]float64{w0, w0}
	m.SpeedAmp, m.SpeedFreq, m.SpeedPhase2 = h, fSchu, phase2
	m.I0 = [2]float64{1000, 1000}
	step := 2e-4
	if m.CParallel > 0 {
		step = math.Min(step, m.RLoad*m.CParallel/5)
	}
	return m.RunRK4(tEnd, step)
}

// ripple returns the mean growth rate of I1, and its 7.8 Hz, 15.6 Hz and
// 23.4 Hz components relative to its DC part, the growth taken out first.
func ripple(tr DualTrace) (g float64, h1, h2, h3 float64) {
	g = rmsGrowth(tr.T, tr.I1, 1)
	y := make([]float64, len(tr.I1))
	for k := range y {
		y[k] = tr.I1[k] * math.Exp(-g*tr.T[k])
	}
	dc := mean(tr.T, y, 1)
	h1 = cmplx.Abs(amplitude(tr.T, y, 1, fSchu)) / math.Abs(dc)
	h2 = cmplx.Abs(amplitude(tr.T, y, 1, 2*fSchu)) / math.Abs(dc)
	h3 = cmplx.Abs(amplitude(tr.T, y, 1, 3*fSchu)) / math.Abs(dc)
	return
}

func runRPMSchumann() error {
	w0 := rpmBase * 2 * math.Pi / 60
	fmt.Printf("  F. rpm(t) = %.0f·(1 ± h·sin(2π·%.1f·t)): the shafts follow the Schumann sine (speed drives)\n", rpmBase, fSchu)

	// The copper machine at 1768 rpm.
	a := rpmMachine(1)
	l11, l12, ks, r := a.loop()
	ls, ms := l11+l12, (1+ks)*a.Coupling()
	fmt.Printf("     copper machine (air coils): self-excites at %.0f rpm; at %.0f rpm the field decays at (−R + M_s·ω0)/L_s = %+.2f 1/s\n",
		r/ms*60/(2*math.Pi), rpmBase, (-r+ms*w0)/ls)
	fmt.Print("       growth with the sine, in phase:")
	for _, h := range []float64{0, 0.5, 1} {
		g, _, _, _ := ripple(rpmRun(a, h, 0, 3))
		fmt.Printf("  h = %.1f: %+.2f 1/s", h, g)
	}
	fmt.Println()
	fmt.Println("       the sine averages out: the mean speed decides, and 1768 rpm is a fifth of the onset. Nothing builds.")

	// The same disks with iron: the coupling raised until 1768 rpm is the onset.
	k := r / (ms * w0)
	b := rpmMachine(k)
	l11, l12, ks, r = b.loop()
	ls, ms = l11+l12, (1+ks)*b.Coupling()
	iSat := bSatIron * math.Pi * (b.Disk.Ro*b.Disk.Ro - b.Disk.Ri*b.Disk.Ri) / (2 * math.Pi * ms)
	fmt.Printf("     with an iron magnetic circuit of gain k = %.2f (declared) the onset is 1768 rpm: L_s = %.0f µH, R/L_s = %.2f 1/s;\n",
		k, ls*1e6, r/ls)
	fmt.Printf("     the mean field on the disk reaches %.1f T (iron's usable limit, declared) at %.1f kA per loop\n", bSatIron, iSat/1e3)
	fmt.Println("     in phase, at the onset the mean growth is zero and the 7.8 Hz sine ripples the field:")
	fmt.Println("       h      β = (R/L_s)·h/Ω   growth [1/s]   7.8 Hz / DC (exact 2I1/I0)   15.6 Hz (2I2/I0)   23.4 Hz (2I3/I0)")
	om := 2 * math.Pi * fSchu
	for _, h := range []float64{0.01, 0.1, 0.5, 1} {
		g, h1, h2, h3 := ripple(rpmRun(b, h, 0, 6))
		beta := (r / ls) * h / om
		i0 := besselI(0, beta)
		fmt.Printf("       %-5.2f  %.4f            %+.4f        %.5f (%.5f)              %.6f (%.6f)  %.7f (%.7f)\n",
			h, beta, g, h1, 2*besselI(1, beta)/i0, h2, 2*besselI(2, beta)/i0, h3, 2*besselI(3, beta)/i0)
	}
	{
		q := b
		q.Speed = &[2]float64{w0, w0}
		q.SpeedAmp, q.SpeedFreq = 1, fSchu
		q.I0 = [2]float64{1000, 1000}
		ct, err := q.Run(6, 2e-4)
		if err != nil {
			return err
		}
		g := rmsGrowth(ct.T, ct.I1, 1)
		y := make([]float64, len(ct.I1))
		for j := range y {
			y[j] = ct.I1[j] * math.Exp(-g*ct.T[j])
		}
		fmt.Printf("       circuit at h = 1: 7.8 Hz / DC = %.5f, ledger mismatch %.0e of the drive work\n",
			cmplx.Abs(amplitude(ct.T, y, 1, fSchu))/math.Abs(mean(ct.T, y, 1)), ct.Mismatch/ct.In)
	}
	fmt.Printf("       the loop's L/R time (%.2f s) is long against 1/(2π·7.8) s: the coils smooth the sine (β ≪ 1).\n", ls/r)
	fmt.Println("       (harmonics above the fundamental sit near the analysis floor of ~1e-4.)")

	// What comes out: DC amps and AC amps, at a DC level just under the iron's limit.
	const idc = 10e3
	fmt.Printf("     out of each loop at %.0f kA DC (the level is set by how far the field has built; the speed decides growth):\n", idc/1e3)
	fmt.Println("       h      DC [kA]   AC at 7.8 Hz [A peak / A rms]   at 15.6 Hz [A peak]   load: DC [kW]   AC [W]   AC share")
	for _, h := range []float64{0.1, 0.5, 1} {
		beta := (r / ls) * h / om
		i0 := besselI(0, beta)
		a1, a2 := idc*2*besselI(1, beta)/i0, idc*2*besselI(2, beta)/i0
		pdc, pac := b.RLoad*idc*idc, 0.5*b.RLoad*(a1*a1+a2*a2)
		fmt.Printf("       %-5.1f  %6.1f    %6.0f / %6.0f                  %6.1f                %6.2f          %6.1f   %.2f %%\n",
			h, idc/1e3, a1, a1/math.Sqrt2, a2, pdc/1e3, pac, 100*pac/(pdc+pac))
	}
	fmt.Println("       DC and AC leave on the same conductors: the load gets both. A capacitor across the load shunts some")
	fmt.Println("       AC but cannot separate it (below); an AC tap needs a branch that blocks DC, e.g. a series capacitor to")
	fmt.Printf("       an AC load, whose reactance at 7.8 Hz must be well under that load: C ≫ 1/(2π·7.8·R_ac) (%.0f F for a\n",
		10/(2*math.Pi*fSchu*b.RLoad))
	fmt.Println("       50 µΩ AC load), or a transformer, whose core the DC would saturate unless its flux is cancelled.")
	fmt.Println("       A larger AC share needs a larger β: more R/L_s (it also raises the onset speed) or a slower sine.")

	fmt.Println("     antiphase (one shaft 1 + h·sin, the other 1 − h·sin), at the iron machine's onset:")
	for _, h := range []float64{0.1, 0.5, 1} {
		g, _, _, _ := ripple(rpmRun(b, h, 180, 6))
		fmt.Printf("       h = %.1f: growth %+.4f 1/s\n", h, g)
	}
	fmt.Println("       the two disks' speeds then sum to a constant: the symmetric mode sees no oscillation at all,")
	fmt.Println("       and the stray self-coupling (KS) of each disk ripples the field only at second order.")

	fmt.Println("     with a capacitor across each load, in phase, h = 1, at the onset:")
	for _, f0 := range []float64{fSchu, fSchu / 2} {
		c := 1 / (math.Pow(2*math.Pi*f0, 2) * ls)
		for _, tau := range []float64{0, 0.87} {
			q := b
			q.CParallel, q.ESR = c, tau/c
			g, h1, _, _ := ripple(rpmRun(q, 1, 0, 6))
			lab := "ideal"
			if tau > 0 {
				lab = fmt.Sprintf("ESR %.1f mΩ", 1e3*tau/c)
			}
			fmt.Printf("       tank at %.1f Hz: C = %.1f F (%s), Q = %.3f: growth %+.3f 1/s, 7.8 Hz / DC %.4f\n",
				f0, c, lab, q.RLoad/math.Sqrt(ls/c), g, h1)
		}
	}
	fmt.Println("       the tank's Q = R_L/√(L/C) is far below 1 at these sizes: it neither locks in nor pumps.")

	// The drives.
	j := b.Inertia()
	for _, h := range []float64{0.1, 1} {
		tq := j * w0 * h * om
		fmt.Printf("     drives at h = %.1f: inertial torque ±%.0f N·m per shaft, reactive power up to %.0f kW, swapped with the rotor\n",
			h, tq, tq*w0*(1+h)/1e3)
	}
	fmt.Println("     every 1/7.8 s: four-quadrant (regenerative) drives. Rim at the peak speed: well under the copper limit.")
	fmt.Println("  So rpm = 1768·(1 ± sin 7.8 Hz) puts a 7.8 Hz ripple on the field, with Bessel-function harmonics; it")
	fmt.Println("  cannot make the field build. Building is set by the mean speed against the onset alone.")
	fmt.Println()
	return nil
}
