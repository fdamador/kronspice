// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// PaschenAir is the breakdown voltage of air at pressure p [Torr] and gap
// d [cm] (Raizer, Gas Discharge Physics, 1991, Table 4.1 and eq. 7.3,
// γ = 0.01 declared), as in examples/water_lamp.
func PaschenAir(p, d float64) float64 {
	pd := p * d
	return 365 * pd / math.Log(15*pd/math.Log(1+1/0.01))
}

// AirGap is a symmetric dry-air spark gap of d [m] at 1 atm: arc 30 V +
// 0.2 Ω, out below 0.1 A (declared).
func AirGap(d float64) PolarGap {
	return PolarGap{Fixed: PaschenAir(760, d*100), NegRatio: 1, Tip: 1e-3, Sigma: 1e-12, Varc: 30, Rarc: 0.2, Iext: 0.1, Floor: 330}
}

// bridge adds a four-arm bridge between AC nodes ac1/ac2 and output p/n.
// arm(name, anode, cathode) adds one arm.
func bridge(arm func(name, a, k string), ac1, ac2, p, n string) {
	arm("DA", ac1, p)
	arm("DB", ac2, p)
	arm("DC", n, ac1)
	arm("DD", n, ac2)
}

// cellArm returns an arm of m cells in series for the solver.
func cellArm(net *Net, c Cell, m int, group string) func(name, a, k string) {
	return func(name, a, k string) {
		prev := a
		for i := 0; i < m; i++ {
			next := k
			if i < m-1 {
				next = fmt.Sprintf("%s_%d", name, i)
			}
			net.Cell(fmt.Sprintf("%s%d", name, i), prev, next, c, group)
			prev = next
		}
	}
}

func run(n *Net, stop, dt float64, each func(t float64)) error {
	n.Prepare()
	steps := int(math.Round(stop / dt))
	for s := 0; s < steps; s++ {
		t := float64(s) * dt
		if err := n.Step(t, dt); err != nil {
			return fmt.Errorf("t = %g s: %w", t+dt, err)
		}
		if each != nil {
			each(t + dt)
		}
	}
	return nil
}

func sine(vp, f float64) func(float64) float64 {
	return func(t float64) float64 { return vp * math.Sin(2*math.Pi*f*t) }
}

// ---------------------------------------------------------------------------
// A. Mains DC supply: four lead–aluminium cells (the 1890s Pollak bridge).

// Supply is the reference bridge: 100 V peak, 60 Hz, source 0.5 Ω, 2,200 µF,
// 50 Ω load.
type Supply struct{ Vdc, Ripple, Efficiency, KronVdc float64 }

// RunSupply runs the supply in the solver and in kronspice.
func RunSupply() (Supply, Ledger, error) {
	const f, vp, cycles, skip = 60.0, 100.0, 20, 10
	net := NewNet()
	net.V("VAC", "s", "0", sine(vp, f), 0.5)
	c := DefaultCell(Pb, Al)
	bridge(cellArm(net, c, 1, "bridge"), "s", "0", "p", "n")
	net.R("RREF", "n", "0", 1e9)
	net.C("CL", "p", "n", 2200e-6, 0)
	load := net.R("RL", "p", "n", 50)
	var s Supply
	per := 1 / f
	var T, sv, lo, hi, eLoad, eIn float64
	lo, hi = math.Inf(1), math.Inf(-1)
	dt := per / 2000
	src := net.el[0]
	err := run(net, cycles*per, dt, func(t float64) {
		if t <= skip*per {
			return
		}
		vo := load.v
		T += dt
		sv += vo * dt
		lo, hi = math.Min(lo, vo), math.Max(hi, vo)
		eLoad += vo * vo / 50 * dt
		eIn += -src.src(t) * src.i * dt
	})
	if err != nil {
		return s, Ledger{}, err
	}
	s.Vdc, s.Ripple, s.Efficiency = sv/T, hi-lo, eLoad/eIn
	// kronspice twin
	b := circuit.NewBuilder("supply")
	b.VoltageWave("VAC", "s0", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vp, f}})
	b.Resistor("RS", "s0", "s", 0.5)
	bridge(func(name, a, k string) { c.Add(b, name, a, k) }, "s", "0", "p", "n")
	b.Resistor("RREF", "n", "0", 1e9)
	b.Capacitor("CL", "p", "n", 2200e-6)
	b.Resistor("RL", "p", "n", 50)
	ckt, err := b.BuildTransient()
	if err != nil {
		return s, Ledger{}, err
	}
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: per / 400, Stop: cycles * per, MaxStep: per / 400}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return s, Ledger{}, err
	}
	vpk, _ := wv.NodeVoltage("p")
	vnk, _ := wv.NodeVoltage("n")
	var kt, ks float64
	for i := 1; i < len(wv.Time); i++ {
		if wv.Time[i] <= skip*per {
			continue
		}
		d := wv.Time[i] - wv.Time[i-1]
		kt += d
		ks += (vpk[i] - vnk[i]) * d
	}
	s.KronVdc = ks / kt
	return s, net.Ledger(), nil
}

// ---------------------------------------------------------------------------
// B. Surge-protected supply: a lightning-type surge on the mains line.

// SurgeCase is one run of the surge test.
type SurgeCase struct {
	Name        string
	PeakArm     float64 // largest voltage across a bridge cell (either way) [V]
	BridgeSpark float64 // energy the bridge cells took above their normal conduction [J]
	Arrester    float64 // energy absorbed by the air-gap arrester [J]
	Clamp       float64 // energy absorbed by the aluminium–aluminium clamp [J]
	PeakOut     float64 // highest DC output [V]
	Strikes     int
	Ledger      Ledger
	Time, VArm  []float64 // largest |voltage| across any bridge cell
}

// Surge: 1.2/50 µs open-circuit voltage wave of peak vs from a 2 Ω source
// (the IEC 61000-4-5 combination-wave impedance; the generator's 8/20 µs
// current shape is not modelled — declared), added to the mains at its
// peak. Double exponential fit: τ1 = 0.4074 µs, τ2 = 68.22 µs, k = 1.037.
func surgeWave(vs, t0 float64) func(float64) float64 {
	return func(t float64) float64 {
		v := 100.0 // mains at its peak (constant over the 300 µs run)
		if t > t0 {
			x := t - t0
			v += vs * 1.037 * (math.Exp(-x/68.22e-6) - math.Exp(-x/0.4074e-6))
		}
		return v
	}
}

// RunSurge runs the supply hit by a 4 kV surge, unprotected and protected
// by a 0.2 mm air gap (arrester) and an aluminium–aluminium cell (clamp).
func RunSurge(protected bool, vs float64) (SurgeCase, error) {
	sc := SurgeCase{Name: "unprotected"}
	if protected {
		sc.Name = "protected: air-gap arrester + aluminium–aluminium clamp"
	}
	net := NewNet()
	net.V("VLINE", "in", "0", surgeWave(vs, 5e-6), 2)
	var arr *part
	if protected {
		arr = net.Gap("ARR", "in", "0", AirGap(0.2e-3), "arrester")
	}
	net.L("LS", "in", "x", 10e-6)
	net.R("RLINE", "x", "ac1", 0.5)
	if protected {
		net.Cell("CLAMP", "ac1", "0", DefaultCell(Al, Al), "clamp")
	}
	c := DefaultCell(Pb, Al)
	bridge(cellArm(net, c, 1, "bridge"), "ac1", "0", "p", "n")
	net.R("RREF", "n", "0", 1e9)
	cl := net.C("CL", "p", "n", 2200e-6, 88)
	net.R("RL", "p", "n", 50)
	var arms []*part
	for _, p := range net.el {
		if p.k == kCell && p.group == "bridge" {
			arms = append(arms, p)
		}
	}
	// the run before the surge sets the bridge's normal heat
	normal := 0.0
	stop, dt := 300e-6, 2e-9
	err := run(net, stop, dt, func(t float64) {
		for _, a := range arms {
			sc.PeakArm = math.Max(sc.PeakArm, math.Abs(a.v))
		}
		sc.PeakOut = math.Max(sc.PeakOut, cl.v)
		if t <= 5e-6 {
			normal = 0
			for _, a := range arms {
				normal += a.heat
			}
		}
		if int(math.Round(t/dt))%50 == 0 {
			sc.Time = append(sc.Time, t)
			worst := 0.0
			for _, a := range arms {
				worst = math.Max(worst, math.Abs(a.v))
			}
			sc.VArm = append(sc.VArm, worst)
		}
	})
	if err != nil {
		return sc, err
	}
	l := net.Ledger()
	sc.Ledger = l
	// bridge heat beyond what it would have taken over the run at the pre-surge rate
	sc.BridgeSpark = l.Heat["bridge"] - normal*stop/5e-6
	sc.Arrester, sc.Clamp = l.Heat["arrester"], l.Heat["clamp"]
	if arr != nil {
		sc.Strikes = arr.strikesF + arr.strikesB
	}
	return sc, nil
}

// ---------------------------------------------------------------------------
// C. Spark pulser: a water spark-gap bridge charges a capacitor from a 6 kV,
// 60 Hz transformer; an output air gap dumps it into a load.

// Pulser is the pulser's result.
type Pulser struct {
	Name               string
	Pulses             int
	Rate               float64 // pulses per second
	PeakAmps           float64
	PerPulse           float64 // mean energy into the load per pulse [J]
	LoadPower, InPower float64 // [W]
	BridgeGaps, OutGap float64 // heat [J]
	Rs                 float64 // charging (source) resistance [Ω]
	Vfire              float64 // output gap self-breakdown [V]
	Rload              float64 // [Ω]
	OnTime             float64 // total time the output gap conducted [s]
	BackStrikes        int     // bridge gaps breaking down backwards
	BackEnergy         float64 // energy those backward arcs dissipated [J]
	Ledger             Ledger
	Time, VStore, IOut []float64
}

// PulserDesign chooses the output switch and load.
type PulserDesign struct {
	Name      string
	Triggered bool    // output gap fired at each mains zero crossing (else self-breaking)
	Rload     float64 // [Ω]
}

// PulserDesigns: two self-breaking 1 mm gaps (ringing and overdamped
// load), and a 2 mm gap triggered at each zero crossing of the mains — the
// timing of the synchronous rotary spark gap of Tesla-coil practice — into
// the overdamped load. 1 µH with 0.1 µF is critically damped at 6.3 Ω.
var PulserDesigns = []PulserDesign{
	{"self-breaking 1 mm air gap, 1 Ω load (rings)", false, 1},
	{"self-breaking 1 mm air gap, 8 Ω load (overdamped)", false, 8},
	{"2 mm air gap fired at each mains zero crossing, 8 Ω load", true, 8},
}

// RunPulser runs a design for the given number of mains cycles.
func RunPulser(d PulserDesign, cycles int) (Pulser, error) {
	const f, vp, rs = 60.0, 6000.0, 1000.0
	r := Pulser{Name: d.Name, Rload: d.Rload, Rs: rs}
	net := NewNet()
	src := net.V("VT", "s", "0", sine(vp, f), rs)
	g := DefaultGap()
	g.NegRatio = 3
	g.Recovery = 1e-3
	bridge(func(name, a, k string) { net.Gap(name, a, k, g, "bridge gaps") }, "s", "0", "p", "n")
	net.R("RREF", "n", "0", 1e8)
	store := net.C("CS", "p", "n", 0.1e-6, 0)
	out := AirGap(1e-3)
	if d.Triggered {
		out = AirGap(2e-3)
	}
	r.Vfire = out.Vpos()
	og := net.Gap("GOUT", "p", "o1", out, "output gap")
	lo := net.L("LO", "o1", "o2", 1e-6)
	net.R("RLOAD", "o2", "n", d.Rload)
	var bridgeGaps []*part
	for _, p := range net.el {
		if p.k == kGap && p.group == "bridge gaps" {
			bridgeGaps = append(bridgeGaps, p)
		}
	}
	net.Prepare()
	per := 1 / f
	slow, fast := per/20000, 5e-9
	t, stop := 0.0, float64(cycles)*per
	lastSample, nextTrig := -1.0, per/2
	for t < stop {
		dt := slow
		if og.on || og.v > 0.97*r.Vfire {
			dt = fast
		}
		if d.Triggered && t+dt >= nextTrig {
			og.armed = true
			nextTrig += per / 2
			dt = fast
		}
		if err := net.Step(t, dt); err != nil {
			return r, fmt.Errorf("t = %g s: %w", t+dt, err)
		}
		t += dt
		if og.on {
			r.OnTime += dt
		}
		r.PeakAmps = math.Max(r.PeakAmps, math.Abs(lo.i))
		if t-lastSample >= 10*slow { // an 8 µs trace for the picture
			lastSample = t
			if t < 3*per {
				r.Time = append(r.Time, t)
				r.VStore = append(r.VStore, store.v)
				r.IOut = append(r.IOut, lo.i)
			}
		}
	}
	l := net.Ledger()
	r.Ledger = l
	r.Pulses = og.strikesF + og.strikesB
	r.Rate = float64(r.Pulses) / stop
	if r.Pulses > 0 {
		r.PerPulse = l.Heat["RLOAD"] / float64(r.Pulses)
	}
	r.LoadPower, r.InPower = l.Heat["RLOAD"]/stop, src.supplied/stop
	r.BridgeGaps, r.OutGap = l.Heat["bridge gaps"], l.Heat["output gap"]
	for _, b := range bridgeGaps {
		r.BackStrikes += b.strikesB
		r.BackEnergy += b.heatBack
	}
	return r, nil
}

// ---------------------------------------------------------------------------
// D. Voltage doubler (Greinacher) of water cells, in kronspice.

// Doubler is a doubler run.
type Doubler struct {
	PerPosition        int
	Vdc, Efficiency    float64
	ReversePeak        float64
	EnergyIn, Mismatch float64
}

// RunDoubler: 100 V peak, 60 Hz, source 0.5 Ω, pump and output capacitors
// 2,200 µF, 200 Ω load, m cells in series at each position.
func RunDoubler(m int) (Doubler, error) {
	const f, vp, cycles, skip = 60.0, 100.0, 30, 20
	d := Doubler{PerPosition: m}
	c := DefaultCell(Pb, Al)
	b := circuit.NewBuilder("doubler")
	b.VoltageWave("VAC", "s0", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vp, f}})
	b.Resistor("RS", "s0", "s", 0.5)
	b.Capacitor("C1", "s", "x", 2200e-6)
	arm := func(name, a, k string) {
		prev := a
		for i := 0; i < m; i++ {
			next := k
			if i < m-1 {
				next = fmt.Sprintf("%s_%d", name, i)
			}
			c.Add(b, fmt.Sprintf("%s%d", name, i), prev, next)
			prev = next
		}
	}
	arm("D1", "0", "x")
	arm("D2", "x", "out")
	b.Capacitor("C2", "out", "0", 2200e-6)
	b.Resistor("RL", "out", "0", 200)
	ckt, err := b.BuildTransient()
	if err != nil {
		return d, err
	}
	per := 1 / f
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: per / 400, Stop: cycles * per, MaxStep: per / 400}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return d, err
	}
	vo, _ := wv.NodeVoltage("out")
	vx, _ := wv.NodeVoltage("x")
	vs, _ := wv.NodeVoltage("s")
	v0, _ := wv.NodeVoltage("s0")
	var T, sv, pl, ps float64
	for i := 1; i < len(wv.Time); i++ {
		d.ReversePeak = math.Max(d.ReversePeak, vx[i]) // D1 blocked: x above ground
		if wv.Time[i] <= skip*per {
			continue
		}
		dt := wv.Time[i] - wv.Time[i-1]
		T += dt
		sv += vo[i] * dt
		pl += vo[i] * vo[i] / 200 * dt
		ps += v0[i] * (v0[i] - vs[i]) / 0.5 * dt
	}
	d.Vdc, d.Efficiency = sv/T, pl/ps
	d.EnergyIn, d.Mismatch = wv.Ledger.EnergyIn, wv.Ledger.IntegrationMismatch
	return d, nil
}

// ---------------------------------------------------------------------------
// E. A varactor with memory: the aluminium oxide film as a tuning capacitor.

// Varactor is one tuning point.
type Varactor struct {
	Vform, Bias float64
	FilmC       float64 // [F]
	F0, Q       float64 // kronspice AC resonance and quality factor
	F0Formula   float64 // 1/(2π√(L·C))
}

// varactorL is the tank inductor [H] (declared).
const varactorL = 1e-3

// RunVaractor forms the aluminium electrode's oxide at vform (the film is
// then about 1.3 nm per volt thick, so its capacitance goes as 1/vform),
// biases the cell at bias (positive = blocking), and finds the resonance of
// a 1 mH tank across it with a kronspice AC sweep.
func RunVaractor(vform, bias float64) (Varactor, error) {
	c := DefaultCell(Al, Pb) // anode aluminium: blocking for positive bias
	c.Vspark = vform         // the film sparks above its formation voltage
	v := Varactor{Vform: vform, Bias: bias, FilmC: c.FilmC()}
	v.F0Formula = 1 / (2 * math.Pi * math.Sqrt(varactorL*v.FilmC))
	b := circuit.NewBuilder("varactor")
	b.Add(circuit.Element{Name: "VAC", Kind: circuit.VoltageSource, ACMag: 1}, "src", "0")
	b.Resistor("RSRC", "src", "t", 100e3)
	c.Add(b, "X", "t", "0")
	b.Inductor("LT", "t", "m", varactorL)
	b.Capacitor("CB", "m", "0", 1) // AC short; lets the bias reach the cell through LT
	b.VoltageSource("VB", "bias", "0", bias)
	b.Resistor("RB", "bias", "m", 1e3)
	ckt, err := b.Build()
	if err != nil {
		return v, err
	}
	lo, hi := v.F0Formula/4, v.F0Formula*4
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "dec", Points: 2000, FStart: lo, FStop: hi}, analysis.Options{})
	if err != nil {
		return v, err
	}
	vt, _ := res.NodeVoltage("t")
	k := 0
	for i := range vt {
		if cmplxAbs(vt[i]) > cmplxAbs(vt[k]) {
			k = i
		}
	}
	v.F0 = res.Freq[k]
	// Q from the −3 dB width
	top := cmplxAbs(vt[k]) / math.Sqrt2
	i1, i2 := k, k
	for i1 > 0 && cmplxAbs(vt[i1]) > top {
		i1--
	}
	for i2 < len(vt)-1 && cmplxAbs(vt[i2]) > top {
		i2++
	}
	if i1 > 0 && i2 < len(vt)-1 && i2 > i1 && k > 0 && k < len(vt)-1 {
		v.Q = v.F0 / (res.Freq[i2] - res.Freq[i1])
	} else {
		v.F0 = 0 // no resonance in the sweep: the conducting cell shorts the tank
	}
	return v, nil
}

func cmplxAbs(z complex128) float64 { return math.Hypot(real(z), imag(z)) }

// ---------------------------------------------------------------------------
// F. Frequency doubler: a bridge's output repeats twice per mains cycle.

// Harmonics is the output's spectrum relative to the source peak.
type Harmonics struct {
	Name       string
	F          float64
	DC, H2, H4 float64 // amplitudes ÷ source peak
	H2Fraction float64 // share of the output's AC power at 2f
}

func spectrum(name string, f float64, t, v []float64, vp float64) Harmonics {
	h := Harmonics{Name: name, F: f}
	var T, a0, total float64
	var c2, s2, c4, s4 float64
	for i := 1; i < len(t); i++ {
		dt := t[i] - t[i-1]
		w := 2 * math.Pi * f * t[i]
		T += dt
		a0 += v[i] * dt
		c2 += v[i] * math.Cos(2*w) * dt
		s2 += v[i] * math.Sin(2*w) * dt
		c4 += v[i] * math.Cos(4*w) * dt
		s4 += v[i] * math.Sin(4*w) * dt
		total += v[i] * v[i] * dt
	}
	h.DC = a0 / T / vp
	h.H2 = 2 * math.Hypot(c2, s2) / T / vp
	h.H4 = 2 * math.Hypot(c4, s4) / T / vp
	ac := total/T - (a0/T)*(a0/T)
	h.H2Fraction = 0.5 * (h.H2 * vp) * (h.H2 * vp) / ac
	return h
}

// RunDoubling: the lead–aluminium bridge (100 V, 60 Hz) and the water
// spark-gap bridge (6 kV, 60 Hz) into resistive loads, last 10 cycles.
func RunDoubling() ([]Harmonics, error) {
	const f = 60.0
	per := 1 / f
	var out []Harmonics
	// cells
	{
		net := NewNet()
		net.V("VAC", "s", "0", sine(100, f), 0.1)
		bridge(cellArm(net, DefaultCell(Pb, Al), 1, "bridge"), "s", "0", "p", "n")
		net.R("RREF", "n", "0", 1e9)
		load := net.R("RL", "p", "n", 10)
		var t, v []float64
		dt := per / 4000
		err := run(net, 12*per, dt, func(tt float64) {
			if tt > 2*per {
				t, v = append(t, tt), append(v, load.v)
			}
		})
		if err != nil {
			return nil, err
		}
		out = append(out, spectrum("lead–aluminium cell bridge, 100 V", f, t, v, 100))
	}
	// water spark gaps
	{
		g := DefaultGap()
		g.NegRatio = 3
		g.Recovery = 1e-3
		net := NewNet()
		net.V("VAC", "s", "0", sine(6000, f), 10)
		bridge(func(name, a, k string) { net.Gap(name, a, k, g, "gaps") }, "s", "0", "p", "n")
		net.R("RREF", "n", "0", 1e9)
		load := net.R("RL", "p", "n", 1000)
		var t, v []float64
		dt := per / 20000
		err := run(net, 12*per, dt, func(tt float64) {
			if tt > 2*per {
				t, v = append(t, tt), append(v, load.v)
			}
		})
		if err != nil {
			return nil, err
		}
		out = append(out, spectrum("water spark-gap bridge, 6 kV", f, t, v, 6000))
	}
	return out, nil
}

// IdealDoubling is a perfect full-wave bridge: DC 2/π, 2f 4/(3π), 4f 4/(15π).
func IdealDoubling() Harmonics {
	h := Harmonics{Name: "perfect diode bridge", DC: 2 / math.Pi, H2: 4 / (3 * math.Pi), H4: 4 / (15 * math.Pi)}
	h.H2Fraction = 0.5 * h.H2 * h.H2 / (0.5 - h.DC*h.DC) // |sin| has mean square ½
	return h
}
