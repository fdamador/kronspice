// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 5a: implosion and explosion of the E-ring (H28 × H30), as a circuit.
//
// The ring's radiation resistance R_rad is where its energy leaves: power
// R_rad·i² goes out as magnetic-dipole radiation. Here it is treated as a
// port to the outside, a line of impedance R_rad. A wave arriving on that
// line is a Thévenin source V_s = 2·v_inc behind R_rad (the standard
// equivalent of a matched line carrying an incident wave), so
//
//	v_inc = V_s/2,  v_ref = v_inc − R_rad·i,  P_inc = v_inc²/R_rad,  P_ref = v_ref²/R_rad,
//
// with i the loop current leaving the source. P_inc − P_ref = V_s·i −
// R_rad·i² is the power that enters the ring. v_ref is the outgoing wave:
// the reflection while the ring is driven, the explosion after.
//
// The implosion is the explosion run backwards: the ring-down's outgoing
// wave reversed in time is a wave growing as e^{+αt} at the ring's own
// frequency. For a one-port resonator a drive e^{st} is reflected with
// amplitude 1 − 2R_rad/Z(s), Z(s) = sL + 1/(sC) + R_rad + R_d; it vanishes
// at sL + 1/(sC) = R_rad − R_d, i.e. at the growth rate
// α_cpa = (R_rad − R_d)/(2L): the time-reversed ring-down of the ring with
// its internal loss reversed in sign, so it is absorbed with no reflection
// at all ("coherent perfect absorption", the time-reversed laser:
// Y. D. Chong, L. Ge, H. Cao, A. D. Stone, Phys. Rev. Lett. 105, 053901,
// 2010). The literal reversal of the lossy ring-down, α = (R_rad +
// R_d)/(2L), leaves the amplitude reflection R_d/(R_rad + R_d).
//
// Standard prediction: every joule held at the peak came in on the
// declared incoming wave; the stored energy then leaves as the ordinary
// ring-down (Q = 40.5); a drive that is not the time-reversed ring-down is
// partly reflected (a resonator at resonance with R_rad ≫ R_d reflects a
// steady wave almost completely). The run books all of it.

import (
	"math"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// BuildPort makes the E-ring loop with R_rad as a port: section 0 holds
// the incident-wave sources (in series, + toward q0), the port resistor
// RP = R_rad from q0 to s0 and the dielectric loss R_d from s0 to n1. The drive's SIN
// segments become one source each (a single 0 V source if there are none),
// so the only breakpoints are the segments' starts.
func (e ERing) BuildPort(d Drive) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("E-ring with its radiation port (fdamador_neutrino Part 5a)")
	s := strconv.Itoa
	n := func(k int) string { return "n" + s(k%e.N) }
	cs, ls := e.C()*float64(e.N), e.L()/float64(e.N)
	for k := 0; k < e.N; k++ {
		b.Capacitor("CS"+s(k), n(k), "m"+s(k), cs)
		if k == 0 {
			b.Inductor("LS0", "m0", "p0", ls)
		} else {
			b.Inductor("LS"+s(k), "m"+s(k), n(k+1), ls)
		}
		b.Resistor("RB"+s(k), "m"+s(k), "0", 1e9)
		b.Resistor("RN"+s(k), n(k), "0", 1e9)
	}
	// The dielectric loss is lumped in section 0 (series elements commute
	// in a loop this small against the wavelength): 24 resistors of
	// R_d/24 = 1.2 mΩ beside the 1 GΩ bleeders made the DC matrix too
	// ill-conditioned for the dense solver's pivot test.
	b.Resistor("RS0", "s0", "n1", e.RDiel())
	b.ResistorTo("RP", "q0", "s0", e.RRad(), exchange.Magnetic)
	addSources(b, d.Segs, "q0", "p0")
	return b.BuildTransient()
}

// addSources puts one SIN source per segment (value 2·v_inc) in series
// from plus to minus.
func addSources(b *circuit.Builder, segs []sinSeg, plus, minus string) {
	if len(segs) == 0 {
		b.Add(circuit.Element{Name: "VP" + plus, Kind: circuit.VoltageSource, Channel: exchange.Magnetic,
			ChannelNote: "passive radiation port"}, plus, minus)
		return
	}
	from := plus
	for k, sg := range segs {
		to := minus
		if k < len(segs)-1 {
			to = "x" + plus + "_" + strconv.Itoa(k)
		}
		b.Add(circuit.Element{Name: "VP" + plus + "_" + strconv.Itoa(k), Kind: circuit.VoltageSource, Channel: exchange.Magnetic,
			ChannelNote: "incident wave on the radiation port (declared)",
			Wave:        circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 2 * sg.VA, sg.F, sg.TD, sg.Theta, 0}}}, from, to)
		from = to
	}
}

// AlphaTotal is the ring-down's amplitude decay rate (R_rad + R_d)/(2L).
func (e ERing) AlphaTotal() float64 { return (e.RRad() + e.RDiel()) / (2 * e.L()) }

// AlphaCPA is the growth rate absorbed with no reflection, (R_rad − R_d)/(2L).
func (e ERing) AlphaCPA() float64 { return (e.RRad() - e.RDiel()) / (2 * e.L()) }

// sinSeg is v_inc = VA·e^{−θ(t−TD)}·sin(2πF(t−TD)) for t ≥ TD (a SPICE
// SIN source; θ < 0 grows).
type sinSeg struct{ TD, VA, F, Theta float64 }

func (g sinSeg) v(t float64) float64 {
	if t < g.TD {
		return 0
	}
	return g.VA * math.Exp(-g.Theta*(t-g.TD)) * math.Sin(2*math.Pi*g.F*(t-g.TD))
}

// piece is a segment from td0 to tp, ended at tp by a second segment that
// cancels it exactly; F·(tp − td0) must be whole, so the sine is at a zero
// there and neither segment jumps.
func piece(td0, tp, va, f, theta float64) []sinSeg {
	return []sinSeg{{td0, va, f, theta}, {tp, -va * math.Exp(-theta*(tp-td0)), f, theta}}
}

// Drive is an incident wave v_inc(t) on the port, zero outside [Start, Stop].
type Drive struct {
	Name        string
	Start, Stop float64
	Segs        []sinSeg
}

// V is v_inc(t) [V].
func (d Drive) V(t float64) float64 {
	var s float64
	for _, g := range d.Segs {
		s += g.v(t)
	}
	return s
}

// PortRun is one driven run of the ring with its port.
type PortRun struct {
	Name                string
	EInc                float64 // incident energy [J]
	ERefDrive           float64 // outgoing energy while the drive is on (the reflection) [J]
	EOutAfter           float64 // outgoing energy after the drive (the explosion) [J]
	HeatDiel            float64 // dielectric heat (R_d) [J]
	PeakStored          float64 // largest stored energy [J]
	StoredAtStop        float64 // stored energy when the drive ends [J]
	StoredEnd           float64
	Absorbed            float64 // 1 − ERefDrive/EInc
	SourceWork, RadHeat float64 // kronspice ledger: work of VP, heat of RP [J]
	LedgerRel           float64 // |work − heat(all) − stored_end| / work
	WaveRel             float64 // |E_inc − E_out(all) − heat_d − stored_end| / E_inc
	Mismatch            float64 // kronspice time-stepping mismatch / work
	Time, Stored, PInc  []float64
	PRef                []float64
}

// RunPort drives the port with d and follows the ring to stop.
func (e ERing) RunPort(d Drive, stop float64) (PortRun, error) {
	T := 2 * math.Pi / e.Omega0()
	ckt, err := e.BuildPort(d)
	if err != nil {
		return PortRun{}, err
	}
	step := T / 200
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return PortRun{}, err
	}
	vq, _ := wf.NodeVoltage("q0")
	vs0, _ := wf.NodeVoltage("s0")
	vp, _ := wf.NodeVoltage("p0")
	R := e.RRad()
	r := PortRun{Name: d.Name}
	l := wf.Ledger
	var prevInc, prevRef float64
	for k := range wf.Time {
		i := (vq[k] - vs0[k]) / R // loop current out of the source, through RP
		vsrc := vq[k] - vp[k]
		vin := vsrc / 2
		vref := vin - R*i
		pin, pref := vin*vin/R, vref*vref/R
		if k > 0 {
			h := wf.Time[k] - wf.Time[k-1]
			r.EInc += 0.5 * (pin + prevInc) * h
			eo := 0.5 * (pref + prevRef) * h
			if wf.Time[k] <= d.Stop+1e-15 {
				r.ERefDrive += eo
			} else {
				r.EOutAfter += eo
			}
		}
		prevInc, prevRef = pin, pref
		w := l.StoredEnergy[k]
		r.PeakStored = math.Max(r.PeakStored, w)
		if wf.Time[k] <= d.Stop+1e-15 {
			r.StoredAtStop = w
		}
		if k%4 == 0 {
			r.Time, r.Stored, r.PInc, r.PRef = append(r.Time, wf.Time[k]), append(r.Stored, w), append(r.PInc, pin), append(r.PRef, pref)
		}
	}
	r.SourceWork = sourceWork(wf.Time, wf.ElementNames, wf.P)
	r.RadHeat = elementEnergy(wf.Time, wf.ElementNames, wf.P, "RP")
	r.HeatDiel = elementEnergy(wf.Time, wf.ElementNames, wf.P, "RS0")
	for k := 0; k < e.N; k++ {
		r.HeatDiel += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RB"+strconv.Itoa(k))
		r.HeatDiel += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RN"+strconv.Itoa(k))
	}
	r.StoredEnd = l.StoredEnd
	r.Absorbed = 1 - r.ERefDrive/r.EInc
	r.LedgerRel = math.Abs(r.SourceWork-r.RadHeat-r.HeatDiel-(l.StoredEnd-l.StoredStart)) / r.SourceWork
	r.WaveRel = math.Abs(r.EInc-r.ERefDrive-r.EOutAfter-r.HeatDiel-r.StoredEnd) / r.EInc
	r.Mismatch = math.Abs(l.IntegrationMismatch) / r.SourceWork
	return r, nil
}

// Growth is the drive A·e^{a(t−tp)}·sin(ω_d(t−tp)), ω_d = √(ω0² − a²),
// over a whole number of periods (about `window`) ending at tp, where it
// is at a zero crossing and the ring holds its peak.
func (e ERing) Growth(name string, a, window, tp, amp float64) Drive {
	f := math.Sqrt(e.Omega0()*e.Omega0()-a*a) / (2 * math.Pi)
	K := math.Round(window * f)
	td0 := tp - K/f
	return Drive{Name: name, Start: td0, Stop: tp, Segs: piece(td0, tp, amp*math.Exp(-a*K/f), f, -a)}
}

// Burst is a sine at f0 of constant amplitude (theta = 0) or decaying at
// rate theta, over a whole number of periods ending at tp.
func (e ERing) Burst(name string, length, tp, theta float64) Drive {
	f := e.Omega0() / (2 * math.Pi)
	K := math.Round(length * f)
	td0 := tp - K/f
	return Drive{Name: name, Start: td0, Stop: tp, Segs: piece(td0, tp, 1, f, theta)}
}

// driveEnergy is ∫ v_inc²/R_rad dt of a drive, sampled finely.
func (e ERing) driveEnergy(d Drive) float64 {
	T := 2 * math.Pi / e.Omega0()
	h := T / 400
	var s float64
	for t := d.Start; t < d.Stop; t += h {
		v := d.V(t + h/2)
		s += v * v * h
	}
	return s / e.RRad()
}

// scaled returns d with its amplitude set so its incident energy is target.
func (e ERing) scaled(d Drive, target float64) Drive {
	k := math.Sqrt(target / e.driveEnergy(d))
	segs := make([]sinSeg, len(d.Segs))
	for i, g := range d.Segs {
		g.VA *= k
		segs[i] = g
	}
	d.Segs = segs
	return d
}

// ImplosionSet is Part 5a: the time-reversed drives and the comparisons,
// all with the same incident energy, then the explosion.
type ImplosionSet struct {
	TauE       float64 // energy decay time Q/ω0 [s]
	AlphaTot   float64
	AlphaCPA   float64
	RefCPA     float64 // predicted amplitude reflection for exact reversal, R_d/(R_rad+R_d)
	Runs       []PortRun
	WindowTauE float64
}

// Implosions runs every drive. Window: 8 energy decay times (the growing
// drive starts at e⁻⁸ of its final power); the explosion is followed for
// as long again.
func (e ERing) Implosions() (ImplosionSet, error) {
	T := 2 * math.Pi / e.Omega0()
	tauE := e.Q() / e.Omega0()
	win := 8 * tauE
	tp := math.Ceil((win+2*T)/T) * T
	stop := tp + win
	set := ImplosionSet{TauE: tauE, AlphaTot: e.AlphaTotal(), AlphaCPA: e.AlphaCPA(), WindowTauE: 8,
		RefCPA: e.RDiel() / (e.RRad() + e.RDiel())}
	ref := e.Growth("time-reversed ring-down, CPA rate (R_rad−R_d)/2L", e.AlphaCPA(), win, tp, 1)
	ein := e.driveEnergy(ref)
	drives := []Drive{
		ref,
		e.scaled(e.Growth("time-reversed ring-down, exact rate (R_rad+R_d)/2L", e.AlphaTotal(), win, tp, 1), ein),
		e.scaled(e.Growth("growing at twice the rate", 2*e.AlphaCPA(), win, tp, 1), ein),
		e.scaled(e.Growth("growing at half the rate", 0.5*e.AlphaCPA(), win, tp, 1), ein),
		e.scaled(e.Burst("the ring-down itself (decaying, not reversed)", win, tp, e.AlphaTotal()), ein),
		e.scaled(e.Burst("constant burst at f0, 2 τ_E long", 2*tauE, tp, 0), ein),
		e.scaled(e.Burst("constant burst at f0, 8 τ_E long", win, tp, 0), ein),
	}
	for _, d := range drives {
		r, err := e.RunPort(d, stop)
		if err != nil {
			return set, err
		}
		set.Runs = append(set.Runs, r)
	}
	return set, nil
}

// TrainCycle is one cycle of the pulsing ring.
type TrainCycle struct {
	EInc, EOut, Heat, StoredStart, StoredEnd float64
	Peak                                     float64
}

// TrainRun is the pulsing ring: an implosion every Period, each drive the
// CPA growth over Window ending at the cycle's peak.
type TrainRun struct {
	Name           string
	PeriodT        float64 // cycle length in ring periods
	Cycles         []TrainCycle
	LedgerRel      float64
	SourceWork     float64
	NetPerCycle    float64 // E_inc − E_out of the last cycle [J]
	HeatLastCycle  float64
	PeakLastCycle  float64
	SingleShotPeak float64
}

// Train drives `cycles` implosions, each 3 τ_E of CPA growth, one every
// `periods` ring periods (whole periods of the drive frequency). With
// sign = −1 every drive after the first is reversed in sign, so it meets
// the free ring left from the cycle before in antiphase.
func (e ERing) Train(name string, periods float64, cycles int, sign float64) (TrainRun, error) {
	T := 2 * math.Pi / e.Omega0()
	tauE := e.Q() / e.Omega0()
	a := e.AlphaCPA()
	f := math.Sqrt(e.Omega0()*e.Omega0()-a*a) / (2 * math.Pi)
	Td := 1 / f
	K := math.Round(3 * tauE * f)
	win := K * Td
	P := periods * Td
	first := math.Ceil((win+T)/Td) * Td
	stop := first + float64(cycles-1)*P + 0.5*P
	var peaks []float64
	drive := Drive{Name: name, Start: 0, Stop: stop}
	for c := 0; c < cycles; c++ {
		tp := first + float64(c)*P
		peaks = append(peaks, tp)
		sg := 1.0
		if c > 0 {
			sg = sign
		}
		drive.Segs = append(drive.Segs, piece(tp-win, tp, sg*math.Exp(-a*win), f, -a)...)
	}
	ckt, err := e.BuildPort(drive)
	if err != nil {
		return TrainRun{}, err
	}
	step := T / 200
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return TrainRun{}, err
	}
	vq, _ := wf.NodeVoltage("q0")
	vs0, _ := wf.NodeVoltage("s0")
	vp, _ := wf.NodeVoltage("p0")
	R := e.RRad()
	tr := TrainRun{Name: name, PeriodT: periods}
	bounds := []float64{0}
	for c := 1; c < cycles; c++ {
		bounds = append(bounds, peaks[c]-win-0.01*Td)
	}
	bounds = append(bounds, stop+1)
	cyc := make([]TrainCycle, cycles)
	var pPrevIn, pPrevOut float64
	pd := make([]float64, len(wf.Time)) // dielectric heat power
	for k := 0; k < e.N; k++ {
		for _, nm := range []string{"RS", "RB", "RN"} {
			for j, n := range wf.ElementNames {
				if n == nm+strconv.Itoa(k) {
					for t := range pd {
						pd[t] += wf.P[j][t]
					}
				}
			}
		}
	}
	ci := 0
	for k := range wf.Time {
		for ci+1 < cycles && wf.Time[k] >= bounds[ci+1] {
			cyc[ci].StoredEnd = wf.Ledger.StoredEnergy[k]
			ci++
			cyc[ci].StoredStart = wf.Ledger.StoredEnergy[k]
		}
		i := (vq[k] - vs0[k]) / R
		vin := (vq[k] - vp[k]) / 2
		vref := vin - R*i
		pin, pout := vin*vin/R, vref*vref/R
		if k > 0 {
			h := wf.Time[k] - wf.Time[k-1]
			cyc[ci].EInc += 0.5 * (pin + pPrevIn) * h
			cyc[ci].EOut += 0.5 * (pout + pPrevOut) * h
			cyc[ci].Heat += 0.5 * (pd[k] + pd[k-1]) * h
		}
		pPrevIn, pPrevOut = pin, pout
		cyc[ci].Peak = math.Max(cyc[ci].Peak, wf.Ledger.StoredEnergy[k])
	}
	cyc[cycles-1].StoredEnd = wf.Ledger.StoredEnd
	tr.Cycles = cyc
	l := wf.Ledger
	tr.SourceWork = sourceWork(wf.Time, wf.ElementNames, wf.P)
	var heat float64
	for j := range wf.ElementNames {
		if wf.ElementNames[j][0] == 'R' {
			heat += elementEnergy(wf.Time, wf.ElementNames, wf.P, wf.ElementNames[j])
		}
	}
	tr.LedgerRel = math.Abs(tr.SourceWork-heat-(l.StoredEnd-l.StoredStart)) / tr.SourceWork
	last := cyc[cycles-2] // the last full cycle
	tr.NetPerCycle = last.EInc - last.EOut
	tr.HeatLastCycle = last.Heat
	tr.PeakLastCycle = last.Peak
	tr.SingleShotPeak = cyc[0].Peak
	return tr, nil
}

// sourceWork is the work done by every incident-wave source (VP…) [J].
func sourceWork(time []float64, names []string, p [][]float64) float64 {
	var w float64
	for _, n := range names {
		if len(n) > 2 && n[:2] == "VP" {
			w -= elementEnergy(time, names, p, n)
		}
	}
	return w
}
