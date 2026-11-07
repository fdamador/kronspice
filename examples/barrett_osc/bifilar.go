// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The Fig. 2G bench with Tesla bifilar pancakes, and with a Metglas core.
//
// Coils. The pump coil α and the shuttle coil β are the flat Archimedean
// spirals of examples/pulsed_bifilar (214 turns of 0.5 mm copper, 127 m),
// whose 20-section distributed model gives, for the same L = 6.84 mH:
//
//	ordinary winding  f1 = 757 kHz → self-capacitance 6.5 pF
//	Tesla bifilar     f1 = 34.9 kHz → self-capacitance 3.03 nF (470×)
//
// Here each coil is that lumped equivalent: L with its wire resistance
// (10.9 Ω) and the self-capacitance across its terminals, valid up to and
// around the first resonance. A bifilar pancake is a capacitor above 35 kHz,
// so the bench runs below it: pump 10 kHz, signal 6.3 kHz, idler 3.7 kHz.
//
// Core. The two pancakes lie on one axis (Fig. 2G), coupled k = 0.6 in air
// (declared). They are modelled as a transformer: each winding's leakage
// (1 − k)·L in series, and a shared magnetizing inductance k·L across an
// ideal 1:1 transformer (a VCVS and a CCCS, lossless), which keeps the
// windings isolated. A Metglas 2605SA1 ribbon backing (declared) multiplies
// the magnetizing inductance by MuEff and saturates it:
//
//	λ(i) = Lair·i + (Lcore − Lair)·Isat·tanh(i/Isat)
//
// (Lair = k·L, Lcore = MuEff·k·L), with a core-loss resistance across it
// (Q = 50 at the pump). The ribbon's magnetostrictive resonance (Metglas
// k_me ≈ 0.3, mechanical Q ≈ 300, declared) is the motional branch of a
// magnetoelastic resonator: a parallel tank (k_me²·Lcore, tuned to the
// ribbon's first longitudinal mode, Q) in series with the magnetizing
// inductance, the standard equivalent of a magnetostrictive resonator. A
// ribbon of length ℓ rings at v/2ℓ, v ≈ 4.0 km/s in annealed 2605SA1
// (declared), so 3.7 kHz needs ℓ ≈ 0.54 m.
//
// Signal and output. The shuttle β runs from G1 through the signal source
// (6.3 kHz, 50 Ω) to its store G2, a 10 nF capacitor with a 1 GΩ leak
// (at kHz a pF top-load would pass almost nothing). The load (50 Ω) is in
// the pump loop with α, as drawn.
type Bench struct {
	Bifilar        bool    // bifilar pancakes (else ordinary windings)
	Metglas        bool    // saturable Metglas core
	Ribbon         bool    // magnetostrictive ribbon resonance
	Bias           float64 // DC magnetizing current through the pump loop [A]
	RibbonF        float64 // ribbon resonance [Hz]
	L, Rw          float64 // coil inductance [H], wire resistance [Ω]
	Cself          float64 // set from Bifilar
	K              float64 // air coupling
	MuEff, Isat    float64 // core gain and saturation current of the magnetizing branch [A]
	Kme, Qmech     float64 // magnetoelastic coupling, mechanical Q
	F0             float64 // grid [Hz]
	MP, MS         int     // pump and signal multiples of F0
	Vp, Vs         float64 // amplitudes [V]
	PhiS           float64 // signal phase [deg]
	Rs, RL, Cstore float64 // source and load [Ω], shuttle store [F]
	Points         int     // time points per window (80 000 for the printed results)
}

// Pancake self-capacitance from the pulsed_bifilar distributed model.
const (
	pancakeL       = 6.842e-3  // H, both windings
	pancakeCOrd    = 6.455e-12 // F, ordinary flat spiral (f1 = 757 kHz)
	pancakeCBif    = 3.033e-9  // F, Tesla bifilar flat spiral (f1 = 34.9 kHz)
	pancakeRw      = 10.9      // Ω, 127 m of 0.5 mm copper
	metglasSoundV  = 4000.0    // m/s, longitudinal sound speed of annealed 2605SA1 (declared)
	metglasMuEff   = 4.0       // open-path pancake backing (declared)
	metglasIsat    = 0.02      // A, magnetizing current at saturation onset (declared)
	metglasKme     = 0.3       // magnetoelastic coupling (declared)
	metglasQmech   = 300.0     // ribbon mechanical Q (declared)
	benchCoreQ     = 50.0      // core loss Q at the pump (declared)
	benchStoreLeak = 1e9       // Ω
	benchBias      = 0.013     // A: DC bias at the core curve's steepest bend (≈ 0.66·Isat), declared
)

// DefaultBench is the bifilar air-core bench.
func DefaultBench() Bench {
	return Bench{Bifilar: true, L: pancakeL, Rw: pancakeRw, K: 0.6, MuEff: metglasMuEff, Isat: metglasIsat,
		Kme: metglasKme, Qmech: metglasQmech, F0: 100, MP: 100, MS: 63, Vp: 30, Vs: 1,
		Rs: 50, RL: 50, Cstore: 10e-9, RibbonF: 3700, Points: 80000}
}

func (b Bench) cself() float64 {
	if b.Bifilar {
		return pancakeCBif
	}
	return pancakeCOrd
}

// lm returns the unsaturated and air magnetizing inductances.
func (b Bench) lm() (lcore, lair float64) {
	lair = b.K * b.L
	lcore = lair
	if b.Metglas {
		lcore = b.MuEff * lair
	}
	return
}

// Lambda is the magnetizing flux linkage λ(i).
func (b Bench) Lambda(i float64) float64 {
	lc, la := b.lm()
	if !b.Metglas {
		return la * i
	}
	return la*i + (lc-la)*b.Isat*math.Tanh(i/b.Isat)
}

// Build makes the circuit.
func (b Bench) Build(pumpOn, sigOn bool) (*circuit.Circuit, error) {
	cb := circuit.NewBuilder("Fig. 2G bench")
	lc, _ := b.lm()
	leak := (1 - b.K) * b.L
	amp := func(on bool, a float64) float64 {
		if on {
			return a
		}
		return 0
	}
	src := func(name, p, n string, vo, a, f, ph float64) {
		cb.Add(circuit.Element{Name: name, Kind: circuit.VoltageSource, Value: vo, HasDC: true,
			Wave: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{vo, a, f, 0, 0, ph}}}, p, n)
	}
	fp, fs := float64(b.MP)*b.F0, float64(b.MS)*b.F0
	// Pump loop: G1 – VP – Rs – a – [α] – c – RL – G1.
	// A DC bias through the pump loop (Rs + Rw + RL) sets the core's
	// operating point; without it the symmetric core makes only odd orders.
	src("VP", "vp", "0", b.Bias*(b.Rs+b.Rw+b.RL), amp(pumpOn, b.Vp), fp, 0)
	cb.Resistor("RS", "vp", "a", b.Rs)
	// α winding a→c: Rw, leakage, then the transformer primary p1 (to c).
	cb.Resistor("RWA", "a", "a1", b.Rw)
	cb.Inductor("LLA", "a1", "p1", leak)
	cb.Capacitor("CSA", "a", "c", b.cself())
	cb.Resistor("RL", "c", "0", b.RL)
	// Magnetizing branch across p1–c: (motional tank) in series with λ(i).
	top := "p1"
	if b.Ribbon && b.Metglas {
		lmot := b.Kme * b.Kme * lc
		w := 2 * math.Pi * b.RibbonF
		cb.Inductor("LMOT", "p1", "mt", lmot)
		cb.Capacitor("CMOT", "p1", "mt", 1/(w*w*lmot))
		cb.Resistor("RMOT", "p1", "mt", b.Qmech*w*lmot)
		top = "mt"
	}
	if b.Metglas {
		var cv []float64
		for j := -600; j <= 600; j++ {
			i := 8 * b.Isat * float64(j) / 600
			cv = append(cv, i, b.Lambda(i))
		}
		cb.NonlinearInductor("LM", top, "c", cv, circuit.Waveform{})
	} else {
		cb.Inductor("LM", top, "c", lc)
	}
	cb.Resistor("RCORE", "p1", "c", benchCoreQ*2*math.Pi*fp*lc)
	// Ideal 1:1 transformer: secondary voltage = primary voltage (VCVS),
	// primary current = secondary current (CCCS sensing the secondary).
	cb.Add(circuit.Element{Name: "ET", Kind: circuit.VCVS, Value: 1, CtrlNodes: [2]int{cb.Node("p1"), cb.Node("c")}}, "s1", "s0")
	cb.Add(circuit.Element{Name: "VSENSE", Kind: circuit.VoltageSource}, "s1", "s2")
	cb.Add(circuit.Element{Name: "FT", Kind: circuit.CCCS, Value: 1, CtrlSource: "VSENSE"}, "p1", "c")
	// β winding: G1 – VS – Rs – g – [β] – w (store G2).
	src("VS", "vs", "0", 0, amp(sigOn, b.Vs), fs, b.PhiS)
	cb.Resistor("RSS", "vs", "g", b.Rs)
	cb.Resistor("RWB", "g", "s0", b.Rw)
	cb.Inductor("LLB", "s2", "w", leak)
	cb.Capacitor("CSB", "g", "w", b.cself())
	cb.Capacitor("CST", "w", "0", b.Cstore)
	cb.Resistor("RLEAKB", "w", "0", benchStoreLeak)
	return cb.Build()
}

// BenchRun is one run's window.
type BenchRun struct {
	T, VLoad, VLoadLinear []float64
	Load                  map[int]complex128 // load voltage spectrum
	IM                    []float64          // magnetizing current over the window
	Ledger                float64
}

// Run simulates warm windows and one analysis window 1/F0 of the load
// voltage.
func (b Bench) Run(pumpOn, sigOn bool, kmax int) (BenchRun, error) {
	var r BenchRun
	c, err := b.Build(pumpOn, sigOn)
	if err != nil {
		return r, err
	}
	const warm = 3
	per := b.Points
	T := 1 / b.F0
	dt := T / float64(per)
	stop := float64(warm+1) * T
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: stop, MaxStep: dt}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return r, err
	}
	r.Ledger = math.Abs(wf.Ledger.IntegrationMismatch) / math.Max(wf.Ledger.EnergyIn, 1e-300)
	v, _ := wf.NodeVoltage("c")
	t0 := stop - T
	r.VLoad = resample(wf.Time, v, t0, dt, per)
	r.T = make([]float64, per)
	for j := range r.T {
		r.T[j] = float64(j) * dt
	}
	if im, ok := wf.ElementCurrent("LM"); ok {
		r.IM = resample(wf.Time, im, t0, dt, per)
	}
	r.Load = map[int]complex128{}
	for k := 1; k <= kmax; k++ {
		r.Load[k] = dftReal(r.VLoad, k)
	}
	return r, nil
}

// ManleyRowe for the magnetizing inductance, from its flux curve: at each
// line, V = j2πf·Λ(f) with Λ the transform of λ(i(t)).
func (b Bench) ManleyRowe(r BenchRun, order int) (mr1, mr2 float64) {
	lam := make([]float64, len(r.IM))
	for j, i := range r.IM {
		lam[j] = b.Lambda(i)
	}
	pw := func(m, n int) float64 {
		k := m*b.MP + n*b.MS
		if k < 0 {
			k = -k
		}
		I := dftReal(r.IM, k)
		V := complex(0, 2*math.Pi*float64(k)*b.F0) * dftReal(lam, k)
		return 2 * real(V*cmplx.Conj(I))
	}
	var s1, a1, s2, a2 float64
	for m := -order; m <= order; m++ {
		for n := -order; n <= order; n++ {
			f := float64(m*b.MP + n*b.MS)
			if f == 0 {
				continue
			}
			P := pw(m, n)
			if m > 0 {
				s1 += float64(m) * P / f
				a1 += math.Abs(float64(m) * P / f)
			}
			if n > 0 {
				s2 += float64(n) * P / f
				a2 += math.Abs(float64(n) * P / f)
			}
		}
	}
	return s1 / a1, s2 / a2
}

// BenchCase is one configuration's summary.
type BenchCase struct {
	Name                  string
	B                     Bench
	Pump, Sig, Idler, IM3 float64 // load-voltage line amplitudes [V]: fp, fs, fp−fs, 2fp−fs
	MR1, MR2              float64
	Run                   BenchRun
}

// BenchCases runs the five configurations.
func BenchCases() ([]BenchCase, error) { return benchCases(DefaultBench(), nil) }

// benchCases runs the configurations whose index is in only (all if nil).
func benchCases(base Bench, only []int) ([]BenchCase, error) {
	cfg := []struct {
		name                 string
		bif, metglas, ribbon bool
		bias                 float64
	}{
		{"ordinary pancakes, air", false, false, false, 0},
		{"bifilar pancakes, air", true, false, false, 0},
		{"bifilar, Metglas core", true, true, false, 0},
		{"bifilar, Metglas, DC bias", true, true, false, benchBias},
		{"bifilar, Metglas, bias, ribbon 3.7 kHz", true, true, true, benchBias},
	}
	var out []BenchCase
	for idx, c := range cfg {
		if only != nil && !contains(only, idx) {
			continue
		}
		b := base
		b.Bifilar, b.Metglas, b.Ribbon, b.Bias = c.bif, c.metglas, c.ribbon, c.bias
		r, err := b.Run(true, true, 4*b.MP)
		if err != nil {
			return nil, fmt.Errorf("%s: %w", c.name, err)
		}
		a := func(k int) float64 { return 2 * cmplx.Abs(r.Load[k]) }
		bc := BenchCase{Name: c.name, B: b, Pump: a(b.MP), Sig: a(b.MS), Idler: a(b.MP - b.MS), IM3: a(2*b.MP - b.MS), Run: r}
		if b.Metglas {
			bc.MR1, bc.MR2 = b.ManleyRowe(r, 4)
		}
		out = append(out, bc)
	}
	return out, nil
}

func printBench() ([]BenchCase, error) {
	b := DefaultBench()
	fmt.Println("Fig. 2G bench with Tesla bifilar pancakes and a Metglas core (pancakes from examples/pulsed_bifilar)")
	fmt.Printf("  each pancake L %.2f mH, %.1f Ω; self-capacitance ordinary %.1f pF (f1 757 kHz), bifilar %.2f nF (f1 34.9 kHz)\n",
		pancakeL*1e3, pancakeRw, pancakeCOrd*1e12, pancakeCBif*1e9)
	fmt.Printf("  pump %.0f kHz %.0f V, signal %.1f kHz %.0f V (50 Ω each), idler %.1f kHz; coupling k %.1f; store %.0f nF; load %.0f Ω\n",
		float64(b.MP)*b.F0/1e3, b.Vp, float64(b.MS)*b.F0/1e3, b.Vs, float64(b.MP-b.MS)*b.F0/1e3, b.K, b.Cstore*1e9, b.RL)
	fmt.Printf("  Metglas (declared): magnetizing ×%.0f, saturating from %.0f mA, core Q %.0f; ribbon k_me %.1f, Q %.0f, %.2f m long for 3.7 kHz\n",
		b.MuEff, b.Isat*1e3, benchCoreQ, b.Kme, b.Qmech, metglasSoundV/(2*b.RibbonF))
	cs, err := BenchCases()
	if err != nil {
		return nil, err
	}
	fmt.Println("                                         load at fp     fs         idler fp−fs   2fp−fs     Manley–Rowe   ledger")
	for _, c := range cs {
		mr := "    —"
		if c.B.Metglas {
			mr = fmt.Sprintf("%.0e %.0e", c.MR1, c.MR2)
		}
		fmt.Printf("  %-38s %8.3f V   %8.2e   %8.2e V    %8.2e   %s   %.0e\n", c.Name, c.Pump, c.Sig, c.Idler, c.IM3, mr, c.Run.Ledger)
	}
	fmt.Println("  air-core pancakes, ordinary or bifilar, are linear: no idler (the bifilar winding changes the")
	fmt.Println("  self-capacitance 470×, so the lines differ, but adds no new frequency). An unbiased Metglas core")
	fmt.Println("  is symmetric and makes only odd orders (the four-wave line 2fp − fs); a DC bias makes the even")
	fmt.Println("  orders, the three-wave idler fp − fs. Manley–Rowe holds for the core's flux curve. A ribbon tuned")
	fmt.Println("  to the idler cuts the idler at the load (here 38×): a magnetostrictive resonance is an impedance")
	fmt.Println("  peak of the core branch, a trap at that frequency; to help the idler, tune it elsewhere.")
	fmt.Println()
	return cs, nil
}

func contains(s []int, x int) bool {
	for _, v := range s {
		if v == x {
			return true
		}
	}
	return false
}
