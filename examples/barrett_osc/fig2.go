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

// Fig. 2 of Barrett 1991: the simplified oscillator-shuttle-circuits A–G.
// Each has a pump source (~), the pump/circuit coil α, the shuttle coil β,
// and a load box "L" with its coil γ; G1 is the common earth, G2 and G3
// are floating terminals (charge stores). Read from the drawing:
//
//	A  pump p–G1; α p→q; load q–G1; β p→s; store Cs s→t; terminal G2 at t
//	B  pump p–G1; α and β in parallel p→q, coupled (the arrow); load q–G1
//	C  loop: load, pump, α; G1 at the pump/α junction; β from that
//	   junction up to the terminal G2, coupled to α (arrow i2)
//	D  as C, but β hangs from the far end of α down to G2
//	E  loop: load, pump, α, bottom wire; β (G2 to G3 = G1) lies along the
//	   bottom wire and is coupled to it
//	F  as E, with a capacitor to G1 on the bottom wire (the drawn G1)
//	G  loop: load, pump, pancake α to G1; pancake β on α (coupled) to G2
//
// A floating terminal is a top-load: capacitance Ct to its surroundings
// (G1). The shuttle receives the signal there, as the paper and patent say
// ("oscillator shuttles … are also receive … antennas"): the signal source
// sits in series with Ct, the Thevenin form of a small antenna. The output
// (the "idler") is the current in the load. Every coil has a series
// resistance (declared). Values are declared in Fig2Params: the paper gives
// none.
type Fig2Params struct {
	F0               float64 // base frequency: both tones are multiples of it [Hz]
	MP, MS           int     // pump and signal at MP·F0 and MS·F0
	Vp, Vs           float64 // pump and signal amplitudes [V]
	PhiS             float64 // signal phase [deg]
	Alpha, Beta, Gam float64 // α, β, γ [H]
	Wire             float64 // the bottom wire as an inductor in E, F [H]
	ESR              float64 // series resistance of every coil [Ω]
	RL               float64 // load [Ω]
	Rsrc             float64 // pump source resistance [Ω]
	Cs, Ct, Cg       float64 // A's store, terminal top-load, F's capacitor [F]
	K                float64 // coupling of coupled coils
	// Varactor, when > 0, makes the terminal store nonlinear (A only):
	// q(v) = Ct·(v + v²/(2·Varactor)), a χ⁽²⁾ element (Varactor in volts).
	Varactor float64
}

// DefaultFig2: pump 1 MHz, signal 630 kHz (both on a 10 kHz grid, so every
// mixing product m·fp + n·fs with |m|, |n| < 63 falls on its own bin).
// leak is the insulation leakage of a floating terminal to its
// surroundings [Ω] (declared; it gives the floating node its DC level).
const leak = 1e9

var DefaultFig2 = Fig2Params{F0: 10e3, MP: 100, MS: 63, Vp: 10, Vs: 1,
	Alpha: 10e-6, Beta: 30e-6, Gam: 5e-6, Wire: 1e-6, ESR: 2, RL: 50, Rsrc: 10,
	Cs: 1e-9, Ct: 10e-12, Cg: 1e-9, K: 0.3}

type fig2Build struct {
	b   *circuit.Builder
	p   Fig2Params
	n   int
	sig bool
}

func (f *fig2Build) node() string { f.n++; return fmt.Sprintf("x%d", f.n) }

// coil adds inductor name (with its ESR) from a to c.
func (f *fig2Build) coil(name, a, c string, h float64) {
	m := f.node()
	f.b.Inductor(name, a, m, h)
	f.b.Resistor("R"+name, m, c, f.p.ESR)
}

func (f *fig2Build) src(name, a, c string, amp, freq, phase float64) {
	f.b.Add(circuit.Element{Name: name, Kind: circuit.VoltageSource,
		Wave:  circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, amp, freq, 0, 0, phase}},
		ACMag: amp, ACPhase: phase}, a, c)
}

// pump adds the pump source from a (+) to c, with its source resistance
// at a. (An ideal pump would pin its node and cut the shuttle off from the
// load in A and B.)
func (f *fig2Build) pump(a, c string, on bool) {
	amp := 0.0
	if on {
		amp = f.p.Vp
	}
	m := f.node()
	f.src("VP", m, c, amp, float64(f.p.MP)*f.p.F0, 0)
	f.b.Resistor("RP", m, a, f.p.Rsrc)
}

// load adds the load box (RL with coil γ across it) from a to c.
func (f *fig2Build) load(a, c string) {
	f.b.Resistor("RL", a, c, f.p.RL)
	f.coil("LG", a, c, f.p.Gam)
}

// terminal adds a floating terminal at node t: the received signal in
// series with the top-load to G1.
func (f *fig2Build) terminal(t string) {
	amp := 0.0
	if f.sig {
		amp = f.p.Vs
	}
	m := "tc" // the store's top plate
	f.src("VS", t, m, amp, float64(f.p.MS)*f.p.F0, f.p.PhiS)
	f.b.Resistor("RLEAK", m, "0", leak) // insulation leakage to G1 (declared)
	if f.p.Varactor > 0 {
		var q []float64
		for v := -0.8 * f.p.Varactor; v <= 0.8*f.p.Varactor+1e-9; v += f.p.Varactor / 400 {
			q = append(q, v, f.p.Ct*(v+v*v/(2*f.p.Varactor)))
		}
		f.b.NonlinearCapacitor("CT", m, "0", q, circuit.Waveform{})
		return
	}
	f.b.Capacitor("CT", m, "0", f.p.Ct)
}

// Fig2 builds variant v ('A'–'G') with the pump and/or signal on.
func Fig2(v byte, p Fig2Params, pumpOn, sigOn bool) (*circuit.Circuit, error) {
	f := &fig2Build{b: circuit.NewBuilder(fmt.Sprintf("Barrett 1991 Fig. 2%c", v)), p: p, sig: sigOn}
	switch v {
	case 'A':
		f.pump("p", "0", pumpOn)
		f.coil("LA", "p", "q", p.Alpha)
		f.load("q", "0")
		f.coil("LB", "p", "s", p.Beta)
		f.b.Capacitor("CS", "s", "t", p.Cs)
		f.terminal("t")
	case 'B':
		f.pump("p", "0", pumpOn)
		f.coil("LA", "p", "q", p.Alpha)
		// β carries the received signal in series (B has no terminal).
		amp := 0.0
		if sigOn {
			amp = p.Vs
		}
		f.src("VS", "p", "pb", amp, float64(p.MS)*p.F0, p.PhiS)
		f.coil("LB", "pb", "q", p.Beta)
		f.b.Mutual("K1", "LA", "LB", p.K)
		f.load("q", "0")
	case 'C', 'D', 'G':
		f.pump("a", "0", pumpOn)
		f.coil("LA", "0", "u", p.Alpha)
		f.load("u", "a")
		if v == 'G' {
			p.K = math.Min(2*p.K, 0.9) // pancakes wound on each other: tighter (declared)
		}
		from := "0"
		if v == 'D' {
			from = "u"
		}
		f.coil("LB", from, "w", p.Beta)
		f.b.Mutual("K1", "LA", "LB", p.K)
		f.terminal("w")
	case 'E', 'F':
		f.pump("a", "0", pumpOn)
		f.coil("LA", "0", "u", p.Alpha)
		f.coil("LW", "u", "m", p.Wire)
		f.load("m", "a")
		if v == 'F' {
			f.b.Capacitor("CG", "m", "0", p.Cg)
		}
		f.coil("LB", "w", "0", p.Beta)
		f.b.Mutual("K1", "LW", "LB", p.K)
		f.terminal("w")
	default:
		return nil, fmt.Errorf("no Fig. 2%c", v)
	}
	return f.b.Build()
}

// Spectrum is the load current at the tones and at mixing products.
type Spectrum struct {
	Line      map[[2]int]complex128 // coefficient of e^{+j2π(m·fp+n·fs)t} by (m, n)
	Fund      float64               // largest of |I(fp)|, |I(fs)| [A]
	MixMax    float64               // largest mixing product / Fund
	MixAt     [2]int
	Ledger    float64   // relative ledger mismatch
	TerminalV []float64 // varactor voltage over the window
}

// Run simulates warm windows then one analysis window of 1/F0.
func (p Fig2Params) Run(v byte, pumpOn, sigOn bool, order int) (Spectrum, error) {
	return p.runOut(v, pumpOn, sigOn, order, "RL")
}

func (p Fig2Params) runOut(v byte, pumpOn, sigOn bool, order int, out string) (Spectrum, error) {
	var s Spectrum
	c, err := Fig2(v, p, pumpOn, sigOn)
	if err != nil {
		return s, err
	}
	const warm, per = 3, 20000 // windows; points per window
	T := 1 / p.F0
	stop := float64(warm+1) * T
	dt := T / per
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: stop, MaxStep: dt}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return s, err
	}
	l := wf.Ledger
	s.Ledger = math.Abs(l.IntegrationMismatch) / math.Max(l.EnergyIn, 1e-30)
	i, ok := wf.ElementCurrent(out)
	if !ok {
		return s, fmt.Errorf("no current %s", out)
	}
	t0 := stop - T
	x := resample(wf.Time, i, t0, dt, per)
	s.Line = map[[2]int]complex128{}
	fp, fs := p.MP, p.MS
	for m := -order; m <= order; m++ {
		for n := -order; n <= order; n++ {
			k := m*fp + n*fs
			if k <= 0 {
				continue
			}
			s.Line[[2]int{m, n}] = dftReal(x, k)
		}
	}
	s.Fund = math.Max(cmplx.Abs(s.Line[[2]int{1, 0}]), cmplx.Abs(s.Line[[2]int{0, 1}]))
	for mn, c := range s.Line {
		if mn == [2]int{1, 0} || mn == [2]int{0, 1} {
			continue
		}
		if r := cmplx.Abs(c) / s.Fund; r > s.MixMax {
			s.MixMax, s.MixAt = r, mn
		}
	}
	// Terminal store voltage and current over the window (for Manley–Rowe).
	if vt, ok := wf.NodeVoltage("tc"); ok {
		s.TerminalV = resample(wf.Time, vt, t0, dt, per)
	}
	return s, nil
}

// ACLine returns the load-current phasor of a source alone (AC analysis,
// pump or signal), as the e^{+jωt} coefficient of the transient line.
func (p Fig2Params) ACLine(v byte, pump bool) (complex128, error) {
	c, err := Fig2(v, p, pump, !pump)
	if err != nil {
		return 0, err
	}
	f := float64(p.MS) * p.F0
	if pump {
		f = float64(p.MP) * p.F0
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		return 0, err
	}
	i, ok := r.ElementCurrent("RL")
	if !ok {
		return 0, fmt.Errorf("no RL current")
	}
	// A sin(ωt+φ) source has AC phasor A∠φ; its e^{+jωt} coefficient is
	// (A/2)∠(φ−90°). The circuit is linear, so the output scales the same.
	return i[0] / 2 * cmplx.Exp(complex(0, -math.Pi/2)), nil
}

// resample interpolates (t, y) at t0 + j·dt, j = 0..m−1.
func resample(t, y []float64, t0, dt float64, m int) []float64 {
	out := make([]float64, m)
	k := 0
	for j := range out {
		tj := t0 + float64(j)*dt
		for k < len(t)-2 && t[k+1] < tj {
			k++
		}
		u := (tj - t[k]) / (t[k+1] - t[k])
		out[j] = y[k] + u*(y[k+1]-y[k])
	}
	return out
}

// dftReal is the coefficient of e^{+j2πkt/T} of x sampled over one window
// T, under a Hann window normalised to unit gain. A tone on the grid (bin k,
// neighbours k ± 1 empty) comes out exact; what is not periodic in the
// window (the slow charge-up of a floating terminal, the tail of a ringing
// that has not died out) leaks into other bins only as 1/Δk³.
func dftReal(x []float64, k int) complex128 {
	var acc complex128
	m := len(x)
	for j, v := range x {
		w := 1 - math.Cos(2*math.Pi*float64(j)/float64(m))
		acc += complex(v*w, 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(k*j)/float64(m)))
	}
	return acc / complex(float64(m), 0)
}

// ManleyRowe returns the two Manley–Rowe sums for the powers flowing into
// the varactor, each relative to the sum of its terms' magnitudes:
// Σ_{m≥0,n} m·P_mn/(m·fp+n·fs) and Σ_{m,n≥0} n·P_mn/(m·fp+n·fs), with P_mn
// the mean power into it at m·fp + n·fs (P_mn = P_−m,−n). The current at
// each line is taken from the element's own charge, I = j2πf·Q(f) with
// Q the transform of q(v(t)), which is exact for a lossless element; the
// sampled current carries the integrator's small phase error, which is
// large next to the real power of a nearly reactive line.
func (p Fig2Params) ManleyRowe(s Spectrum, order int) (mr1, mr2 float64) {
	q := make([]float64, len(s.TerminalV))
	for j, v := range s.TerminalV {
		q[j] = p.Ct * (v + v*v/(2*p.Varactor))
	}
	tv := s.TerminalV
	pw := func(m, n int) float64 {
		k := m*p.MP + n*p.MS
		if k < 0 {
			k = -k
		}
		V := dftReal(tv, k)
		I := complex(0, 2*math.Pi*float64(k)*p.F0) * dftReal(q, k)
		return 2 * real(V*cmplx.Conj(I))
	}
	var s1, a1, s2, a2 float64
	for m := -order; m <= order; m++ {
		for n := -order; n <= order; n++ {
			f := float64(m*p.MP + n*p.MS)
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

func printFig2() error {
	p := DefaultFig2
	fmt.Println("Fig. 2 — the oscillator-shuttle-circuits A–G as drawn (coils, top-loads, load; no nonlinear part)")
	fmt.Printf("  pump %.0f kHz %.0f V, signal %.0f kHz %.0f V at the floating terminal; α %.0f µH, β %.0f µH, γ %.0f µH,\n",
		float64(p.MP)*p.F0/1e3, p.Vp, float64(p.MS)*p.F0/1e3, p.Vs, p.Alpha*1e6, p.Beta*1e6, p.Gam*1e6)
	fmt.Printf("  coil ESR %.0f Ω, load %.0f Ω, A's store %.0f nF, top-load %.0f pF, F's capacitor %.0f nF, coupling k %.1f (declared)\n",
		p.ESR, p.RL, p.Cs*1e9, p.Ct*1e12, p.Cg*1e9, p.K)
	fmt.Println("        load |I(fp)|   load |I(fs)|   tran vs AC (fp, fs)   largest mixing product   ledger")
	for _, v := range []byte("ABCDEFG") {
		s, err := p.Run(v, true, true, 3)
		if err != nil {
			return err
		}
		ap, err := p.ACLine(v, true)
		if err != nil {
			return err
		}
		as, err := p.ACLine(v, false)
		if err != nil {
			return err
		}
		ip, is := s.Line[[2]int{1, 0}], s.Line[[2]int{0, 1}]
		ep, es := cmplx.Abs(ip-ap)/cmplx.Abs(ap), cmplx.Abs(is-as)/cmplx.Abs(as)
		fmt.Printf("  %c   %9.3e A   %9.3e A   %8.1e  %8.1e     %8.1e at (%+d,%+d)    %.0e\n",
			v, 2*cmplx.Abs(ip), 2*cmplx.Abs(is), ep, es, s.MixMax, s.MixAt[0], s.MixAt[1], s.Ledger)
	}
	fmt.Println("  every circuit gives the pump and signal lines that AC analysis gives for each source alone")
	fmt.Println("  (superposition) and no mixing product above round-off: coils and capacitors, however their")
	fmt.Println("  grounds float, make a linear network, and a linear network makes no idler. Bench test: drive")
	fmt.Println("  the pump, feed a signal to G2, and look for fp − fs in the load; the standard prediction is only")
	fmt.Println("  what the parts' own nonlinearity gives (ferrite cores, ceramic dielectrics; air-core coils and")
	fmt.Println("  film or vacuum capacitors give least).")
	fmt.Println()

	q := VaractorFig2()
	fmt.Printf("Fig. 2A with a varactor as the terminal store (χ⁽²⁾ element: q = Ct·(v + v²/2·%.0f V)), pump %.0f V, signal %.1f V\n",
		q.Varactor, q.Vp, q.Vs)
	base, err := q.Run('A', true, true, 4)
	if err != nil {
		return err
	}
	k := q.MP - q.MS
	vi, vpump := dftReal(base.TerminalV, k), dftReal(base.TerminalV, q.MP)
	fmt.Printf("  at the varactor terminal (the shuttle's transmitting end): idler fp − fs = %.0f kHz at %.3e of the pump line\n",
		float64(k)*q.F0/1e3, cmplx.Abs(vi)/cmplx.Abs(vpump))
	mr1, mr2 := q.ManleyRowe(base, 4)
	fmt.Printf("  Manley–Rowe sums into the varactor (relative, harmonics to 4th): %.1e, %.1e\n", mr1, mr2)
	ex, err := q.idlerExponents()
	if err != nil {
		return err
	}
	fmt.Printf("  idler ∝ Vs^%.3f · Vp^%.3f (standard χ⁽²⁾: Es*·Ep, exponents 1 and 1; the paper's Eq. 24 / patent\n", ex[0], ex[1])
	fmt.Println("  Eq. 15: |Es|²·Ep, exponent 2 in the signal)")
	dphi, err := q.idlerPhaseShift(60)
	if err != nil {
		return err
	}
	fmt.Printf("  signal phase +60° → idler phase %+.2f° (conjugate: −60°): the phase conjugation the patent claims\n", dphi)
	fmt.Println("  comes from the nonlinear element, as in any parametric mixer, and obeys Manley–Rowe.")
	fmt.Println()
	return nil
}

// VaractorFig2 is Fig. 2A with the terminal store a varactor, at a pump
// small enough for second-order mixing to dominate (declared).
func VaractorFig2() Fig2Params {
	q := DefaultFig2
	q.Varactor, q.Vp, q.Vs = 20, 2, 0.2
	return q
}

// terminalIdler is the idler line (fp − fs) of the varactor terminal's
// voltage: the shuttle's terminal is the transmitting antenna (patent: "the
// line including C2 and L3 … transmit antennas").
func (p Fig2Params) terminalIdler() (complex128, error) {
	s, err := p.Run('A', true, true, 1)
	if err != nil {
		return 0, err
	}
	return dftReal(s.TerminalV, p.MP-p.MS), nil
}

// idlerExponents returns d log|idler| / d log Vs and d log Vp (doubling each).
func (p Fig2Params) idlerExponents() ([2]float64, error) {
	var ex [2]float64
	b, err := p.terminalIdler()
	if err != nil {
		return ex, err
	}
	q := p
	q.Vs *= 2
	s, err := q.terminalIdler()
	if err != nil {
		return ex, err
	}
	ex[0] = math.Log2(cmplx.Abs(s) / cmplx.Abs(b))
	q = p
	q.Vp *= 2
	s, err = q.terminalIdler()
	if err != nil {
		return ex, err
	}
	ex[1] = math.Log2(cmplx.Abs(s) / cmplx.Abs(b))
	return ex, nil
}

// idlerPhaseShift returns the change of the idler's phase [deg] when the
// signal's phase moves by d degrees.
func (p Fig2Params) idlerPhaseShift(d float64) (float64, error) {
	b, err := p.terminalIdler()
	if err != nil {
		return 0, err
	}
	q := p
	q.PhiS += d
	s, err := q.terminalIdler()
	if err != nil {
		return 0, err
	}
	return cmplx.Phase(s/b) * 180 / math.Pi, nil
}
