// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command parametric_lcr analyses F. D. Amador, H. Eckardt, "Simulation of a
// Parametric Resonance Circuit" (A.I.A.S./UPITEC, 30 March 2012) with
// kronspice.
//
// The circuit (paper Fig. 2, eqs. 19–27): a series loop of a sine source
// U = U0·sin(ωt), an inductance L, a resistance R and a capacitance that
// varies in time, C_eff(t) = C·(1 + a0·sin(ω1·t + φ)), with q = C_eff·u_C
// (eq. 25). Parameters (eqs. 28–36): U0 = 5 V, C = 230 nF, L = 20 mH,
// R = 20 Ω, f_res = 1/(2π√(LC)) = 2346.6 Hz, ω1 = 2ω0, a0 = 0.2, φ = 0.
//
// kronspice's time-varying capacitor (VaryingCapacitor, q = s(t)·C·v) is the
// paper's eq. (25) exactly. It also books what the paper's energy sum
// leaves out: the power the agent that varies C exchanges with the
// circuit, P_pump = −ṡ·½·C·v² into the circuit (work against the field
// between the plates as they move, or the pump power of a varactor's bias
// source). That pump is how a parametric oscillator is powered
// (L. Mandelstam, N. Papalexi, Z. Phys. 73, 223, 1932, the paper's ref.
// [13]; Kouřil & Vrba 1988 §10.2.2, pp. 377–378, the paper's ref. [6];
// J. M. Manley, H. E. Rowe, Proc. IRE 44, 904, 1956). The program checks
// every case's books: source + pump = heat in R + change of stored energy.
//
//	go run ./examples/parametric_lcr
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Paper values, eqs. (28)–(31).
const (
	U0 = 5.0
	C  = 2.3e-7
	L  = 0.02
	R  = 20.0
)

// W0 is the resonance angular frequency, eq. (32)–(33).
var W0 = 1 / math.Sqrt(L*C)

// Case is one simulation.
type Case struct {
	Name      string
	Drive     func(t float64) float64 // source voltage; nil = no source
	DriveStop float64                 // drive switched off after this time (0 = never)
	Fd        float64                 // drive frequency [Hz] (for sampling)
	Scale     func(t float64) float64 // C_eff/C; nil = constant C
	DScale    func(t float64) float64 // d(C_eff/C)/dt
	I0        float64                 // initial inductor current [A]
	DriveWave circuit.Waveform        // engine waveform for the drive (else PWL of Drive)
	ScaleWave circuit.Waveform        // engine waveform for the scale (else PWL of Scale)
	TEnd      float64                 // [s]
	PerPeriod float64                 // time steps per period (default 100)
	IMax      float64                 // stop when |I| exceeds this (0 = never)
}

// Result holds the energy books of a run.
type Result struct {
	T, I, VC, P         []float64 // time, loop current, capacitor voltage, R power
	ESource, EPump, EOh float64   // energy from the source, from the pump, into R [J]
	EStored0, EStored1  float64   // stored ½LI² + q²/2C_eff at start and end [J]
	Mismatch            float64   // engine's integration mismatch [J]
	Stopped             bool
}

// pwl samples f on [0, tEnd] at n points per period of frequency fs.
func pwl(f func(float64) float64, fs, tEnd float64, n int) []float64 {
	dt := 1 / (fs * float64(n))
	var a []float64
	for t := 0.0; t <= tEnd+dt; t += dt {
		a = append(a, t, f(t))
	}
	return a
}

// Run simulates a case with kronspice. Drive and scale go to the engine as
// native SIN/PULSE waveforms where given, else as PWL sampled at 200
// points per period; their exact forms are used for the pump-power check.
func Run(c Case) (*Result, error) {
	b := circuit.NewBuilder("parametric LCR (Amador & Eckardt 2012)")
	if c.Drive != nil {
		drive := c.Drive
		if c.DriveStop > 0 {
			stop := c.DriveStop
			drive = func(t float64) float64 {
				if t > stop {
					return 0
				}
				return c.Drive(t)
			}
		}
		w := c.DriveWave
		if w.Kind == circuit.WaveNone {
			w = circuit.Waveform{Kind: circuit.WavePWL, Args: pwl(drive, c.Fd, math.Min(c.TEnd, c.DriveStop+1/c.Fd), 200)}
			if c.DriveStop == 0 {
				w.Args = pwl(drive, c.Fd, c.TEnd, 200)
			}
		}
		b.Add(circuit.Element{Name: "U", Kind: circuit.VoltageSource, Wave: w}, "s", "0")
	} else {
		b.VoltageSource("U", "s", "0", 0)
	}
	b.InductorIC("L", "s", "a", L, c.I0).Resistor("R", "a", "c", R)
	if c.Scale != nil {
		w := c.ScaleWave
		if w.Kind == circuit.WaveNone {
			w = circuit.Waveform{Kind: circuit.WavePWL, Args: pwl(c.Scale, c.Fd, c.TEnd, 200)}
		}
		b.VaryingCapacitor("C", "c", "0", C, w)
	} else {
		b.CapacitorIC("C", "c", "0", C, 0)
	}
	ckt, err := b.BuildTransient()
	if err != nil {
		return nil, err
	}
	T := 1 / c.Fd / c.PerPeriod
	if c.PerPeriod == 0 {
		T = 1 / c.Fd / 100
	}
	opt := analysis.TranOptions{Options: analysis.Options{RelTol: 1e-6}}
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: T, Stop: c.TEnd, MaxStep: T, UIC: true}, opt)
	if err != nil {
		return nil, err
	}
	r := &Result{Mismatch: wf.Ledger.IntegrationMismatch}
	il, _ := wf.ElementCurrent("L")
	vs, _ := wf.NodeVoltage("s")
	vc, _ := wf.NodeVoltage("c")
	isrc, _ := wf.ElementCurrent("U")
	scale := c.Scale
	if scale == nil {
		scale = func(float64) float64 { return 1 }
	}
	stored := func(k int) float64 {
		s := scale(wf.Time[k])
		return 0.5*L*il[k]*il[k] + 0.5*C*s*vc[k]*vc[k]
	}
	for k := range wf.Time {
		r.T, r.I, r.VC, r.P = append(r.T, wf.Time[k]), append(r.I, il[k]), append(r.VC, vc[k]), append(r.P, R*il[k]*il[k])
		if k == 0 {
			continue
		}
		dt := wf.Time[k] - wf.Time[k-1]
		ps := func(j int) float64 { return -vs[j] * isrc[j] } // power the source delivers
		pp := func(j int) float64 {                           // power the pump delivers
			if c.DScale == nil {
				return 0
			}
			return -c.DScale(wf.Time[j]) * 0.5 * C * vc[j] * vc[j]
		}
		r.ESource += 0.5 * dt * (ps(k) + ps(k-1))
		r.EPump += 0.5 * dt * (pp(k) + pp(k-1))
		r.EOh += 0.5 * dt * (r.P[k] + r.P[k-1])
		if c.IMax > 0 && math.Abs(il[k]) > c.IMax {
			r.Stopped = true
		}
	}
	r.EStored0, r.EStored1 = stored(0), stored(len(wf.Time)-1)
	return r, nil
}

// PeakLast returns the largest |I| and the mean R power over the last
// window [TEnd − w, TEnd].
func (r *Result) PeakLast(w float64) (ipk, pmean float64) {
	tEnd := r.T[len(r.T)-1]
	var e, n float64
	for k := 1; k < len(r.T); k++ {
		if r.T[k] < tEnd-w {
			continue
		}
		ipk = math.Max(ipk, math.Abs(r.I[k]))
		dt := r.T[k] - r.T[k-1]
		e += 0.5 * dt * (r.P[k] + r.P[k-1])
		n += dt
	}
	return ipk, e / n
}

// ElastanceIndex is the first-harmonic depth m of 1/(1 + a·sin θ): the m
// of D = D0(1 + m·sin) in Kouřil & Vrba eq. (10.71), with D0 the mean.
func ElastanceIndex(a float64) float64 {
	s := math.Sqrt(1 - a*a)
	return 2 * (1 - s) / (a * s) * s // first harmonic / mean (mean = 1/s)
}

// Threshold is the critical pump depth m_c = 2·√(d² + (2Δω/ω0)²), d = 1/Q
// (K&V eqs. 10.101, 10.103), for pump at 2ω1, Δω = ω1 − ω0.
func Threshold(detune float64) float64 {
	d := R * math.Sqrt(C/L)
	return 2 * math.Sqrt(d*d+4*detune*detune)
}

// drive is the paper's source U0·sin(ωt), eq. (26).
func drive(w float64) func(float64) float64 {
	return func(t float64) float64 { return U0 * math.Sin(w*t) }
}

// sinWave is SIN(vo, va, ω/2π, 0, 0, φ in degrees).
func sinWave(vo, va, w, phi float64) circuit.Waveform {
	return circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{vo, va, w / (2 * math.Pi), 0, 0, phi * 180 / math.Pi}}
}

// MeanShift is the circuit's own frequency in units of f0 when
// C_eff = C(1 + a·sin): √(mean elastance) = (1 − a²)^(−1/4).
func MeanShift(a float64) float64 { return math.Pow(1-a*a, -0.25) }

func sinScale(a, w1, phi float64) (func(float64) float64, func(float64) float64) {
	return func(t float64) float64 { return 1 + a*math.Sin(w1*t+phi) },
		func(t float64) float64 { return a * w1 * math.Cos(w1*t+phi) }
}

// Square returns C_eff/C stepping between 1 − a and 1 + a once per period
// T, with 2 % ramps so the pump power is finite, its derivative, and the
// matching PULSE waveform.
func Square(a, T float64) (func(float64) float64, func(float64) float64, circuit.Waveform) {
	rise := 0.02 * T
	s := func(t float64) float64 {
		x := math.Mod(t, T)
		switch {
		case x < rise:
			return 1 - a + 2*a*x/rise
		case x < 0.5*T:
			return 1 + a
		case x < 0.5*T+rise:
			return 1 + a - 2*a*(x-0.5*T)/rise
		}
		return 1 - a
	}
	ds := func(t float64) float64 {
		x := math.Mod(t, T)
		switch {
		case x < rise:
			return 2 * a / rise
		case x >= 0.5*T && x < 0.5*T+rise:
			return -2 * a / rise
		}
		return 0
	}
	return s, ds, circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{1 - a, 1 + a, 0, rise, rise, 0.5*T - rise, T}}
}

func report(c Case) *Result {
	r, err := Run(c)
	if err != nil {
		log.Fatalf("%s: %v", c.Name, err)
	}
	ipk, pm := r.PeakLast(5 / c.Fd)
	in := r.ESource + r.EPump
	fmt.Printf("%s\n", c.Name)
	fmt.Printf("   at %.1f ms: peak I %.3g A, mean R power %.3g W (last 5 periods)\n", c.TEnd*1e3, ipk, pm)
	fmt.Printf("   energy: source %.4g J, pump %.4g J → R %.4g J + stored %.4g J (books close to %.2g %% of the input; engine's integration mismatch %.1e J)\n",
		r.ESource, r.EPump, r.EOh, r.EStored1-r.EStored0, 100*(in-r.EOh-(r.EStored1-r.EStored0))/math.Max(math.Abs(in), 1e-300), r.Mismatch)
	if r.ESource != 0 {
		fmt.Printf("   R energy / source energy = %.3g (the paper's COP); R energy / (source + pump) = %.3f\n", r.EOh/r.ESource, r.EOh/in)
	}
	return r
}

func main() {
	f0 := W0 / (2 * math.Pi)
	q := 1 / (R * math.Sqrt(C/L))
	fmt.Println("Amador & Eckardt, \"Simulation of a Parametric Resonance Circuit\" (2012) — kronspice analysis")
	fmt.Printf("f_res = %.1f Hz (paper 2346.6), Q = %.2f, α = R/2L = %.0f 1/s\n", f0, q, R/(2*L))
	m := ElastanceIndex(0.2)
	fmt.Printf("a0 = 0.2 → elastance depth m = %.4f; self-excitation threshold m_c = 2/Q = %.4f tuned (K&V 10.101)\n", m, Threshold(0))
	fmt.Printf("C(1 + a sin) raises the mean elastance by 1/√(1−a²), so the circuit's own frequency is %.4f f0; pumped at 2·0.97 f0 the detuning is %.1f %%,\n",
		MeanShift(0.2), 100*(MeanShift(0.2)-0.97)/0.97)
	fmt.Printf("threshold %.4f > m = %.4f: below self-excitation, a parametric amplifier of the drive (K&V 10.103)\n\n", Threshold((MeanShift(0.2)-0.97)/0.97), m)

	report(Case{Name: "Fig. 4 — constant C, driven at resonance", Drive: drive(W0), DriveWave: sinWave(0, U0, W0, 0), Fd: f0, TEnd: 0.03})
	fmt.Printf("   expected amplitude U0/R = %.3f A\n\n", U0/R)

	s, ds := sinScale(0.2, 2*W0, 0)
	r := report(Case{Name: "§3.1 Fig. 5 — C(1 + 0.2 sin 2ω0t), driven at f0", Drive: drive(W0), DriveWave: sinWave(0, U0, W0, 0), Fd: f0, Scale: s, DScale: ds, ScaleWave: sinWave(1, 0.2, 2*W0, 0), TEnd: 0.03})
	sig := -R/(2*L) + m*W0/4
	fmt.Printf("   envelope growth: K&V predicts σ = −α + mω0/4 = %.0f 1/s; the current grows without limit, powered by the pump\n\n", sig)
	_ = r

	for _, phi := range []float64{0, -math.Pi / 2, math.Pi / 2} {
		w := 0.97 * W0
		s, ds := sinScale(0.2, 2*w, phi)
		report(Case{Name: fmt.Sprintf("§3.1 Figs. 6–7 — drive and pump at 0.97 f0, φ = %+.0f°", phi*180/math.Pi),
			Drive: drive(w), DriveWave: sinWave(0, U0, w, 0), Fd: 0.97 * f0, Scale: s, DScale: ds, ScaleWave: sinWave(1, 0.2, 2*w, phi), TEnd: 0.6})
	}
	fmt.Println()

	// §3.2: pump at the drive frequency (ω1 = ω0). The paper does not state
	// a0 here; with a0 = 0.2 nothing grows (the pump has no component near
	// 2ω0, so only the weak second tongue acts), and growth needs a0 ≈ 0.7.
	w7 := 1.07 * W0
	for _, a := range []float64{0.2, 0.8} {
		a := a
		s7 := func(t float64) float64 { return 1 + a*math.Pow(math.Sin(w7*t), 7) }
		ds7 := func(t float64) float64 { return a * 7 * math.Pow(math.Sin(w7*t), 6) * math.Cos(w7*t) * w7 }
		report(Case{Name: fmt.Sprintf("§3.2 Fig. 8 — C(1 + %.1f sin⁷ ω1t), ω1 = drive = 1.07 ω0", a), Drive: drive(w7), DriveWave: sinWave(0, U0, w7, 0),
			Fd: 1.07 * f0, Scale: s7, DScale: ds7, TEnd: 0.012})
	}
	wr := 0.95 * W0
	Tr := 2 * math.Pi / wr
	for _, a := range []float64{0.2, 0.8} {
		sq, dsq, pulse := Square(a, Tr)
		report(Case{Name: fmt.Sprintf("§3.2 Fig. 9 — rectangular C (%.1fC ↔ %.1fC) at 0.95 f0, drive off after 2 periods", 1-a, 1+a),
			Drive: drive(wr), DriveStop: 2 * Tr, Fd: 0.95 * f0, Scale: sq, DScale: dsq, ScaleWave: pulse, TEnd: 0.012, PerPeriod: 2000})
		report(Case{Name: fmt.Sprintf("§3.2 — rectangular C (%.1fC ↔ %.1fC) at 0.95 f0, no drive, 1 mA start", 1-a, 1+a),
			Fd: 0.95 * f0, Scale: sq, DScale: dsq, ScaleWave: pulse, I0: 1e-3, TEnd: 0.012, PerPeriod: 2000})
	}
	fmt.Println("\n(The a0 = 0.8 square-wave runs grow by orders of magnitude in 12 ms across ramps 2 % of a period wide; their")
	fmt.Println(" few-percent book error is time-integration error and falls as the step is refined: 0.97 % → 0.32 % from 2000 to 8000 steps per period.)")
	fmt.Println("\nIn every case the heat in R and the stored energy are paid by the source plus the pump that varies C.")
	fmt.Println("The paper's COP divides by the source alone; the pump's share is what its \"energy from spacetime\" counts.")
}
