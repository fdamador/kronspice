// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Studies on the layer-1 machine (CLAUDE.md in this folder, next tasks 2–4):
//
//   - convergence: the reversal-interval distribution, not the count, is
//     what an integrator of a chaotic system can be asked to reproduce;
//   - the (μ, A) plane against the geomagnetic polarity record;
//   - the expansions of paper107 (Evans, Amador & Eckardt, "Spin connection
//     resonance in the Faraday disk generator," reference/paper107.pdf):
//     a time-varying drive, independently speed-controlled disks, a feedback
//     coil around each disk, Tesla's opposed fields, and the rotating-frame
//     terms of its §14.4, each with the standard prediction.

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"sort"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// RevStat is a long run summarized on the fly: reversal times, the peaks
// of |I1| (local maxima), and the RMS of I1, all after tSkip.
type RevStat struct {
	Rev       []float64
	PeakT     []float64
	Peak      []float64
	RMS       float64
	Diverged  bool
	Intervals []float64
}

// StatRun integrates the machine by RK4 with step h without storing the
// trace (runs of 10⁷ steps).
func (p Params) StatRun(tEnd, h, tSkip, floor float64) RevStat {
	var r RevStat
	s := [4]float64{p.I1, p.I2, p.W1, p.W2}
	n := int(math.Round(tEnd / h))
	last := 0.0
	var a0, a1 float64 // |I1| at the two previous steps
	var sum2 float64
	var cnt int
	for k := 0; k < n; k++ {
		t := float64(k) * h
		s = rk4(p.rhs, t, s, h)
		t += h
		x := s[0]
		if math.IsNaN(x) || math.Abs(x) > 1e6 {
			r.Diverged = true
			break
		}
		ax := math.Abs(x)
		if t >= tSkip {
			sum2 += x * x
			cnt++
			if a1 > a0 && a1 >= ax && a1 > floor {
				r.PeakT, r.Peak = append(r.PeakT, t-h), append(r.Peak, a1)
			}
			if ax >= floor {
				sg := math.Copysign(1, x)
				if last != 0 && sg != last {
					r.Rev = append(r.Rev, t)
				}
				last = sg
			}
		}
		a0, a1 = a1, ax
	}
	if cnt > 0 {
		r.RMS = math.Sqrt(sum2 / float64(cnt))
	}
	for k := 1; k < len(r.Rev); k++ {
		r.Intervals = append(r.Intervals, r.Rev[k]-r.Rev[k-1])
	}
	return r
}

// meanCV returns the mean, the coefficient of variation and the largest
// value of a sample.
func meanCV(x []float64) (mean, cv, mx float64) {
	if len(x) < 2 {
		return math.NaN(), math.NaN(), math.NaN()
	}
	var s, s2 float64
	for _, v := range x {
		s += v
		s2 += v * v
		mx = math.Max(mx, v)
	}
	n := float64(len(x))
	mean = s / n
	cv = math.Sqrt(math.Max(0, s2/n-mean*mean)) / mean
	return
}

// ks2 is the two-sample Kolmogorov–Smirnov statistic and its asymptotic
// p-value (Kolmogorov distribution with the effective size nm/(n+m)).
func ks2(a, b []float64) (d, p float64) {
	x, y := append([]float64(nil), a...), append([]float64(nil), b...)
	sort.Float64s(x)
	sort.Float64s(y)
	i, j := 0, 0
	for i < len(x) && j < len(y) {
		if x[i] <= y[j] {
			i++
		} else {
			j++
		}
		d = math.Max(d, math.Abs(float64(i)/float64(len(x))-float64(j)/float64(len(y))))
	}
	ne := float64(len(x)*len(y)) / float64(len(x)+len(y))
	return d, kolmogorovQ((math.Sqrt(ne) + 0.12 + 0.11/math.Sqrt(ne)) * d)
}

// kolmogorovQ is P(K > λ) for the Kolmogorov distribution.
func kolmogorovQ(l float64) float64 {
	if l < 1e-3 {
		return 1
	}
	var s float64
	for k := 1; k <= 100; k++ {
		s += 2 * math.Pow(-1, float64(k-1)) * math.Exp(-2*float64(k*k)*l*l)
	}
	return math.Min(1, math.Max(0, s))
}

// dExp is the one-sample KS distance of x from an exponential with the
// sample's own mean.
func dExp(x []float64) float64 {
	y := append([]float64(nil), x...)
	sort.Float64s(y)
	m, _, _ := meanCV(y)
	n := float64(len(y))
	var d float64
	for i, v := range y {
		f := 1 - math.Exp(-v/m)
		d = math.Max(d, math.Max(math.Abs(float64(i+1)/n-f), math.Abs(f-float64(i)/n)))
	}
	return d
}

// ksExp tests x against an exponential distribution whose mean is fitted
// from x: the statistic, and a parametric-bootstrap p-value (the mean is
// estimated, so the standard KS table does not apply; Lilliefors 1969).
// The seed makes it reproducible.
func ksExp(x []float64, boot int, seed uint64) (d, p float64) {
	d = dExp(x)
	rng := rand.New(rand.NewPCG(seed, 0x9e3779b97f4a7c15))
	var hit int
	sim := make([]float64, len(x))
	for b := 0; b < boot; b++ {
		for i := range sim {
			sim[i] = rng.ExpFloat64()
		}
		if dExp(sim) >= d {
			hit++
		}
	}
	return d, float64(hit+1) / float64(boot+1)
}

// precursor is the median, over reversals, of the last |I1| peak before
// each reversal divided by the median of all peaks: below 1 means the
// field weakens before it reverses, above 1 that it swings wider.
func precursor(r RevStat) float64 {
	if len(r.Rev) == 0 || len(r.Peak) == 0 {
		return math.NaN()
	}
	var ratios []float64
	all := append([]float64(nil), r.Peak...)
	sort.Float64s(all)
	med := all[len(all)/2]
	j := 0
	for _, tr := range r.Rev {
		for j+1 < len(r.PeakT) && r.PeakT[j+1] < tr {
			j++
		}
		if r.PeakT[j] < tr {
			ratios = append(ratios, r.Peak[j]/med)
		}
	}
	if len(ratios) == 0 {
		return math.NaN()
	}
	sort.Float64s(ratios)
	return ratios[len(ratios)/2]
}

// withMuA returns Rikitake's machine with μ = R and A = ω1 − ω2 (equal
// torques keep A constant), started at (x, y, z) = (1, 0.5, 1).
func withMuA(mu, a float64) Params {
	p := Rikitake()
	p.R = mu
	p.W2 = p.W1 - a
	return p
}

func runConvergence() error {
	fmt.Println("Task 2 — convergence of the reversal statistics (μ = 1, A = 1)")
	const tEnd, skip, floor = 20000.0, 50.0, 1e-2
	hs := []float64{0.004, 0.002, 0.001, 0.0005}
	var ref RevStat
	runs := make([]RevStat, len(hs))
	for i := len(hs) - 1; i >= 0; i-- {
		runs[i] = Rikitake().StatRun(tEnd, hs[i], skip, floor)
	}
	ref = runs[len(hs)-1]
	fmt.Printf("  RK4, t = %.0f:\n", tEnd)
	for i, h := range hs {
		r := runs[i]
		m, cv, mx := meanCV(r.Intervals)
		d, pv := ks2(r.Intervals, ref.Intervals)
		fmt.Printf("    h = %-6g %4d reversals, mean %.2f, CV %.2f, max %.0f; KS vs h = %g: D = %.3f, p = %.2f\n",
			h, len(r.Rev), m, cv, mx, hs[len(hs)-1], d, pv)
	}
	fmt.Println("  kronspice circuit (trapezoidal, LTE control, max step dt), t = 4000:")
	for _, dt := range []float64{0.02, 0.01, 0.005} {
		tr, _, err := Rikitake().RunCircuit(4000, dt)
		if err != nil {
			return err
		}
		rev := Reversals(tr, skip, floor)
		var iv []float64
		for k := 1; k < len(rev); k++ {
			iv = append(iv, rev[k]-rev[k-1])
		}
		m, cv, _ := meanCV(iv)
		d, pv := ks2(iv, ref.Intervals)
		fmt.Printf("    dt = %-5g %4d reversals, mean %.2f, CV %.2f; KS vs RK4 h = 0.0005: D = %.3f, p = %.2f\n", dt, len(rev), m, cv, d, pv)
	}
	fmt.Println("  counts at a fixed time differ run to run (chaos); the interval distributions agree.")
	fmt.Println()
	return nil
}

func runMuA() error {
	fmt.Println("Task 3 — the (μ, A) plane: reversal statistics against an exponential (Poisson) and the polarity record")
	const tEnd, h, skip, floor = 20000.0, 0.002, 200.0, 1e-2
	fmt.Println("    μ     A    reversals  mean    CV    max/mean   KS D vs exp  p(boot)  last peak/median peak")
	seed := uint64(1)
	for _, mu := range []float64{0.5, 1, 2} {
		for _, a := range []float64{0.5, 1, 2, 4} {
			r := withMuA(mu, a).StatRun(tEnd, h, skip, floor)
			seed++
			if r.Diverged || len(r.Intervals) < 10 {
				fmt.Printf("    %-4g  %-4g %5d       (no sustained reversals; RMS I1 %.3g)\n", mu, a, len(r.Rev), r.RMS)
				continue
			}
			m, cv, mx := meanCV(r.Intervals)
			d, pv := ksExp(r.Intervals, 400, seed)
			fmt.Printf("    %-4g  %-4g %5d     %6.2f  %5.2f   %6.1f      %.3f       %.4f    %.2f\n",
				mu, a, len(r.Rev), m, cv, mx/m, d, pv, precursor(r))
		}
	}
	fmt.Println("  polarity record (declared): Cenozoic mean polarity interval of order 0.2–0.3 Myr; the Cretaceous Normal")
	fmt.Println("  Superchron ≈ 38 Myr (≈ 121–84 Ma; Ogg, Geologic Time Scale 2020, ch. 5), so longest/mean ≳ 100;")
	fmt.Println("  the dipole weakens before reversing (Valet & Meynadier, Nature 366, 234, 1993; Valet et al., Nature 435, 802, 2005).")
	fmt.Println("  Rikitake: intervals have a floor at the oscillation period (no short chrons), longest/mean ≲ 10, and the")
	fmt.Println("  field swings wider, not weaker, before it reverses. Self-consistent 3-D dynamos (Glatzmaier & Roberts,")
	fmt.Println("  Nature 377, 203, 1995) are the benchmark the two-disk model does not reach.")
	fmt.Println()
	return nil
}

// runDrive is paper107's time-varying rotation (task 4a): both torques
// modulated, Gᵢ = 1 + GMod·sin(Wd·t), swept across the natural frequency.
func runDrive() error {
	_, per := Equilibrium(1)
	wn := 2 * math.Pi / per
	fmt.Printf("Task 4a — modulated drive torque G = 1 + 0.3·sin(ωd·t), swept across the natural frequency %.3f\n", wn)
	const tEnd, h, skip, floor = 10000.0, 0.002, 200.0, 1e-2
	base := Rikitake().StatRun(tEnd, h, skip, floor)
	fmt.Printf("    unmodulated: RMS I1 %.3f, %.4f reversals per unit time\n", base.RMS, float64(len(base.Rev))/(tEnd-skip))
	fmt.Println("    ωd/ωn   RMS I1   reversals/time   mean interval")
	for _, x := range []float64{0.25, 0.5, 0.75, 0.9, 1, 1.1, 1.25, 1.5, 2, 3} {
		p := Rikitake()
		p.GMod, p.Wd = 0.3, x*wn
		r := p.StatRun(tEnd, h, skip, floor)
		m, _, _ := meanCV(r.Intervals)
		fmt.Printf("    %-5.2f   %.3f    %.4f           %.2f\n", x, r.RMS, float64(len(r.Rev))/(tEnd-skip), m)
	}
	// Twin: the circuit with the modulated torque sources against RK4.
	p := Rikitake()
	p.GMod, p.Wd = 0.3, wn
	tr, mis, err := p.RunCircuit(10, 0.01)
	if err != nil {
		return err
	}
	ref := p.RunRK4(10, 0.01)
	var dmax float64
	for k := 1; k < len(tr.T); k++ { // k = 0 is the imposed UIC point
		j := int(math.Round(tr.T[k] / 0.01))
		dmax = math.Max(dmax, math.Abs(tr.I1[k]-ref.I1[j]))
	}
	fmt.Printf("  circuit (SIN torque sources) vs RK4 at ωd = ωn, t ≤ 10: |ΔI1| ≤ %.2g; engine ledger mismatch %.1g J\n", dmax, mis)
	fmt.Println("  standard prediction: a driven chaotic oscillator; any resonance is the system's own (no new term).")
	fmt.Println()
	return nil
}

// speedRun drives both shafts at prescribed speeds (task 4b) and returns
// the circuit's late growth rate of the field and the twin linear RK4 one.
func speedRun(p Params, w1, w2 circuit.Waveform, tEnd float64) (circ, ref float64, err error) {
	p.Speed = &[2]circuit.Waveform{w1, w2}
	tr, _, err := p.RunCircuit(tEnd, 0.005)
	if err != nil {
		return
	}
	growth := func(t, x, y []float64) float64 { // slope of ln|I| over the last half
		n := len(t)
		a, b := n/2, n-1
		na := math.Hypot(x[a], y[a])
		nb := math.Hypot(x[b], y[b])
		return (math.Log(nb) - math.Log(na)) / (t[b] - t[a])
	}
	circ = growth(tr.T, tr.I1, tr.I2)
	// Linear ODE with ω(t) prescribed.
	wave := func(w circuit.Waveform, t float64) float64 {
		if w.Kind == circuit.WaveSin {
			return w.Args[0] + w.Args[1]*math.Sin(2*math.Pi*w.Args[2]*t+w.Args[5]*math.Pi/180)
		}
		return w.Args[0]
	}
	f := func(t float64, s [4]float64) [4]float64 {
		o1, o2 := wave(w1, t), wave(w2, t)
		return [4]float64{
			(-p.R*s[0] + p.Wind2*p.M*o1*s[1] + p.K*p.M*o1*s[0]) / p.L,
			(-p.R*s[1] + p.M*o2*s[0] + p.K*p.M*o2*s[1]) / p.L,
		}
	}
	s := [4]float64{p.I1, p.I2}
	const hh = 1e-3
	n := int(math.Round(tEnd / hh))
	var ts, xs, ys []float64
	for k := 0; k <= n; k++ {
		ts, xs, ys = append(ts, float64(k)*hh), append(xs, s[0]), append(ys, s[1])
		s = rk4(f, float64(k)*hh, s, hh)
	}
	ref = growth(ts, xs, ys)
	return
}

func constW(w float64) circuit.Waveform {
	return circuit.Waveform{Kind: circuit.WaveNone, Args: []float64{w}}
}

func runSpeedControl() error {
	fmt.Println("Task 4b — each disk on its own speed drive (paper107 §14.5.2.2): the field equations become linear")
	fmt.Println("  standard prediction: L·İ = (−R·1 + M·[[0, ω1],[ω2, 0]])·I grows at (−R + M·√(ω1ω2))/L when ω1ω2 > 0,")
	fmt.Println("  decays at −R/L (oscillating at M·√|ω1ω2|/L) when ω1ω2 < 0: self-excitation needs ω1ω2 > (R/M)².")
	fmt.Println("    ω1     ω2     exact     circuit   RK4")
	for _, w := range [][2]float64{{2, 2}, {2, 0.8}, {1, 1}, {0.5, 1.5}, {3, 0.2}, {2, -2}} {
		exact := -1.0
		if w[0]*w[1] > 0 {
			exact = -1 + math.Sqrt(w[0]*w[1])
		}
		tEnd := math.Min(20, 6/math.Max(math.Abs(exact), 0.05)) // stop near e⁶ of growth
		c, r, err := speedRun(Rikitake(), constW(w[0]), constW(w[1]), tEnd)
		if err != nil {
			return err
		}
		fmt.Printf("    %-5g  %-5g  %+.4f   %+.4f   %+.4f\n", w[0], w[1], exact, c, r)
	}
	fmt.Println("  prescribed profiles ω1 = 1 + 0.5·sin(ωm·t), ω2 = 1 + 0.5·sin(ωm·t + φ) (mean speeds at threshold):")
	fmt.Println("    ωm    φ      circuit   RK4     (in phase the sum mode's mean growth is ⟨ω⟩ − 1 = 0: bounded; a 30-unit window reads it within ±0.03)")
	for _, c := range []struct{ wm, ph float64 }{{1, 0}, {1, 180}, {0.2, 180}, {2, 180}, {4, 180}} {
		w1 := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, c.wm / (2 * math.Pi), 0, 0, 0}}
		w2 := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, c.wm / (2 * math.Pi), 0, 0, c.ph}}
		g, r, err := speedRun(Rikitake(), w1, w2, 60)
		if err != nil {
			return err
		}
		fmt.Printf("    %-4g  %-5g  %+.4f   %+.4f\n", c.wm, c.ph, g, r)
	}
	fmt.Println("  opposed modulation lowers ⟨√(ω1ω2)⟩ below the mean speed: the field decays at the old threshold.")
	fmt.Println()
	return nil
}

// runFeedback is paper107's feedback inductor (§14.5.2.3, task 4c): each
// disk also in the field of a coil carrying its own current (K·M).
func runFeedback() error {
	fmt.Println("Task 4c — a feedback coil around each disk (K·M; K > 0 aids the cross field, K < 0 opposes it)")
	const tEnd, h, skip, floor = 20000.0, 0.002, 200.0, 1e-2
	fmt.Println("    K      reversals  mean interval  CV    RMS I1   mean ω1   (Bullard alone: self-excites when K·M·ω > R)")
	for _, k := range []float64{-0.5, -0.25, 0, 0.25, 0.5, 0.75} {
		p := Rikitake()
		p.K = k
		r := p.StatRun(tEnd, h, skip, floor)
		m, cv, _ := meanCV(r.Intervals)
		// mean ω1 from a short stored run
		tr := p.RunRK4(2000, 0.01)
		var wsum float64
		for _, w := range tr.W1[len(tr.W1)/2:] {
			wsum += w
		}
		fmt.Printf("    %+.2f   %5d      %8.2f       %.2f   %.3f    %.3f\n", k, len(r.Rev), m, cv, r.RMS, wsum/float64(len(tr.W1)-len(tr.W1)/2))
	}
	p := Rikitake()
	p.K = 0.5
	tr, mis, err := p.RunCircuit(10, 0.01)
	if err != nil {
		return err
	}
	ref := p.RunRK4(10, 0.01)
	var dmax float64
	for k := 1; k < len(tr.T); k++ {
		dmax = math.Max(dmax, math.Abs(tr.I1[k]-ref.I1[int(math.Round(tr.T[k]/0.01))]))
	}
	fmt.Printf("  circuit (machine element with G[0][0] = −K·M) vs RK4 at K = 0.5, t ≤ 10: |ΔI1| ≤ %.2g; ledger mismatch %.1g J\n", dmax, mis)
	fmt.Println("  K > 0: the field locks to one polarity (no reversals): each loop now feeds itself as well as the other.")
	fmt.Println("  K < 0: the coil's I² torque drives the shaft forward, the field weakens and reversals become rare.")
	fmt.Println("  the feedback coil changes how the drive's energy splits between field and heat; it adds no input:")
	fmt.Println("  its EMF K·M·ω·I is paid by the back torque K·M·I² on the shaft (the ledger closes in every case).")
	fmt.Println()
	return nil
}

// runRotatingFrame is task 4e: the machine element under dω/dt ≠ 0 and
// the size of the rotating-frame terms of paper107 §14.4 for a laboratory
// disk.
func runRotatingFrame() error {
	fmt.Println("Task 4e — rotating-frame terms (paper107 §14.4) and the machine element under dω/dt ≠ 0")
	// One disk, field coil fed by a 1 A source, shaft driven at
	// ω(t) = α0·sin(ω0·t) (paper's α0 = 5, ω0 = 10 in its units).
	const a0, w0 = 5.0, 10.0
	sp := circuit.NewMachineSpec(2)
	sp.L0[0][0], sp.L0[1][1] = 0.1, 1
	sp.G[0][1] = 1
	b := circuit.NewBuilder("Faraday disk, prescribed speed").Domain(exchange.Mechanical, "w")
	b.MachineModel("disk", sp)
	b.Machine("N1", "disk", [][2]string{{"a", "b"}, {"f", "0"}}, "w", "0")
	b.Resistor("RL", "b", "a", 1).Resistor("RA", "a", "0", 1e6)
	b.Add(circuit.Element{Name: "IF", Kind: circuit.CurrentSource, Value: 1, HasDC: true}, "f", "0")
	b.Add(circuit.Element{Name: "VW", Kind: circuit.VoltageSource, Value: 0, HasDC: true, Channel: exchange.Mechanical,
		Wave: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, a0, w0 / (2 * math.Pi), 0, 0, 0}}}, "w", "0")
	c, err := b.Build()
	if err != nil {
		return err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: 1e-3, Stop: 3, MaxStep: 1e-3}, analysis.TranOptions{})
	if err != nil {
		return err
	}
	il, _ := wf.ElementCurrent("RL")
	// Exact: 0.1·i' = −i + ω(t)·1, i(0) = 0, ω = a0·sin(w0 t):
	// i = A·(sin − τw0·cos)/(1 + (τw0)²) + A·τw0/(1 + (τw0)²)·e^(−t/τ), τ = 0.1.
	tau := 0.1
	x := tau * w0
	var dmax, imax float64
	for k, t := range wf.Time {
		ex := a0 * ((math.Sin(w0*t) - x*math.Cos(w0*t)) + x*math.Exp(-t/tau)) / (1 + x*x)
		// RL's current is booked b → a; the loop current runs a → b inside the disk, b → a in RL.
		dmax = math.Max(dmax, math.Abs(il[k]-ex))
		imax = math.Max(imax, math.Abs(ex))
	}
	fmt.Printf("  disk at ω(t) = 5·sin(10t), field 1 A: loop current vs exact solution within %.1g of its peak; engine ledger mismatch %.1g J\n",
		dmax/imax, wf.Ledger.IntegrationMismatch)
	fmt.Println("  the speed EMF M·ω(t)·I_f follows the instantaneous speed: nothing in the lumped machine needs dω/dt.")
	// Sizes of the inertial terms for a laboratory copper disk.
	const (
		me = 9.1093837015e-31
		qe = 1.602176634e-19
		B  = 1.0
		r  = 0.1
	)
	om := 2 * math.Pi * 50
	al := om * 10 // ω(t) = ω0·cos(10 t): peak dω/dt
	wc := qe * B / me
	fmt.Printf("  laboratory disk: B = %.0f T, r = %.0f cm, Ω = 2π·50 rad/s, peak dΩ/dt = %.0f rad/s², electron cyclotron ω_c = %.3g 1/s\n", B, 100*r, al, wc)
	emf := 0.5 * B * om * r * r
	st := me / qe * al * r * r / 2
	cf := me / qe * om * om * r * r / 2
	fmt.Printf("    motional EMF ½BΩr² = %.4g V\n", emf)
	fmt.Printf("    angular (Stewart–Tolman) EMF (m/e)·(dΩ/dt)·r²/2 = %.3g V (ratio %.1g)\n", st, st/emf)
	fmt.Printf("    centrifugal EMF (m/e)·Ω²r²/2 = %.3g V (ratio Ω/ω_c = %.1g)\n", cf, om/wc)
	fmt.Printf("    Coriolis force on a drift electron 2mΩv vs evB: ratio 2Ω/ω_c = %.1g\n", 2*om/wc)
	fmt.Println("  these are the standard inertial terms of the rotating frame; the lumped machine omits them at ~1e-9.")
	fmt.Println("  paper107 states no further, measurable term; its ECE reading stays in the register (H15, H16),")
	fmt.Println("  where the Maxwell–Lorentz prediction above is the standard counterpart.")
	fmt.Println()
	return nil
}

// runTesla is task 4d: Tesla's two-disk machine (US 406,968) has opposed
// field magnets and co-rotating disks, so the two EMFs add around the belt;
// in the cross-excited Rikitake machine the same opposition is case 4.
func runTesla() error {
	fmt.Println("Task 4d — Tesla's opposed fields (US 406,968; paper107 §14.5.2.1) in the cross-excited machine")
	p := Rikitake()
	p.Wind2 = -1
	// Linear growth about I = 0 at the initial speeds: eigenvalues −R ± i·M·√(ω1ω2).
	lam := complex(-p.R, 0) + cmplx.Sqrt(complex(-p.M*p.M*p.W1*1, 0))
	fmt.Printf("  co-rotating, coil 2 rewound: about I = 0, λ = −R ± i·M·√(ω1ω2): %.3f ± %.3fi at ω1 = ω2 = 1\n", real(lam), math.Abs(imag(lam)))
	tr := p.RunRK4(60, 0.01)
	fmt.Printf("  run: |I1| after t = 60 is %.2g; the shafts spin up without load (ω1 = %.1f)\n", math.Abs(tr.I1[len(tr.I1)-1]), tr.W1[len(tr.W1)-1])
	fmt.Println("  so with each disk excited by the other, opposed fields cannot self-excite; Tesla's magnets excite")
	fmt.Println("  his disks separately, and his opposed arrangement only makes the EMFs add along the belt (examples/homopolar).")
	fmt.Println("  Reversing a disk and its coil together (case 3) restores the loop gain: the sign of the loop matters, not the spin.")
	fmt.Println()
	return nil
}

func runStudies() error {
	for _, f := range []func() error{runConvergence, runMuA, runDrive, runSpeedControl, runFeedback, runTesla, runRotatingFrame} {
		if err := f(); err != nil {
			return err
		}
	}
	return nil
}
