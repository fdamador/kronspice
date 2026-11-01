// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// A low-voltage, high-current pickup on the pulse-driven ordinary solenoid.
//
// The ordinary coil rings as a free half-wave line: its ends carry the
// voltage maxima and its middle the current maximum. A pickup of a few
// turns of thick wire wound round the middle, tuned with a series
// capacitor to the coil's resonance, takes energy out magnetically at low
// voltage and high current. The tuning capacitor is the one that makes the
// pickup's own inductance resonate at f1, C = 1/((2πf1)²·L_pickup) (the
// series LC of the pickup loop: its reactance cancels at f1, so only the
// load and the wire resistance remain). The pickup's inductance and its
// coupling to every coil section and to the primary come from the same
// coaxial-circle formulas as the coil.
//
// Compared: the tuned pickup at the middle; the same pickup untuned (no
// capacitor); the tuned pickup moved to the high-voltage top end; each
// over a range of load resistances. The ledger books the load's power
// against the pulse source's.

import (
	"fmt"
	"math"
	"math/cmplx"
	"runtime/debug"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/material"
)

// Declared pickup geometry.
const (
	PickupTurns  = 4
	PickupRadius = 0.065 // m: wound round the 5 cm solenoid, inside its 10 cm can
	PickupWireR  = 1e-3  // m: 2 mm copper
	PickupPitch  = 3e-3  // m
	PickupGround = "u0"  // the pickup loop's low end is the circuit ground reference
)

// Pickup is one pickup arrangement.
type Pickup struct {
	Name    string
	Centre  float64 // axial position of the pickup's middle [m]
	Tuned   bool
	CFactor float64 // tuning capacitor as a multiple of the resonant value (1 = correct)
	RLoad   float64 // [Ω]
}

func (pk Pickup) z(i int) float64 {
	return pk.Centre + (float64(i)-float64(PickupTurns-1)/2)*PickupPitch
}

// PickupL is the pickup's self inductance [H].
func PickupL() float64 {
	pk := Pickup{}
	l := 0.0
	for i := 0; i < PickupTurns; i++ {
		for j := 0; j < PickupTurns; j++ {
			if i == j {
				l += em.Mu0 * PickupRadius * (math.Log(8*PickupRadius/PickupWireR) - 2)
			} else {
				l += em.MutualCircles(PickupRadius, PickupRadius, pk.z(i)-pk.z(j))
			}
		}
	}
	return l
}

// PickupR is the pickup wire's resistance at f [Ω] (skin effect as
// WireResistance, for the thicker wire).
func PickupR(f float64) float64 {
	length := PickupTurns * 2 * math.Pi * PickupRadius
	area := math.Pi * PickupWireR * PickupWireR
	rdc := material.Copper.Resistivity * length / area
	delta := math.Sqrt(material.Copper.Resistivity / (math.Pi * f * em.Mu0))
	return rdc * math.Max(1, 2*PickupWireR/(4*delta))
}

// TuningC is the capacitor that resonates the pickup at f [F].
func TuningC(f float64) float64 {
	w := 2 * math.Pi * f
	return 1 / (w * w * PickupL())
}

// BuildWithPickup is Build plus the pickup loop: LPU u0 → u1, wire
// resistance u1 → u2, tuning capacitor u2 → u3 (or a short when untuned),
// load u3 → u0, u0 grounded; K elements couple the pickup to every coil
// section and to the primary.
func (c *Coil) BuildWithPickup(src circuit.Element, f1 float64, pk Pickup) (*circuit.Circuit, error) {
	return c.build(src, c.SectionResistance(f1), func(b *circuit.Builder) {
		g := c.Geom
		lpu := PickupL()
		b.Inductor("LPU", PickupGround, "u1", lpu)
		b.Resistor("RPU", "u1", "u2", PickupR(f1))
		if pk.Tuned {
			b.Capacitor("CT", "u2", "u3", pk.CFactor*TuningC(f1))
		} else {
			b.Resistor("RSHORT", "u2", "u3", 1e-6)
		}
		b.Resistor("RLOAD", "u3", PickupGround, pk.RLoad)
		b.Resistor("RGND", PickupGround, "0", 1e-6)
		// Coupling to each section.
		ms := make([]float64, NumSections())
		for t := 0; t < g.N(); t++ {
			p := c.Pos[t]
			for i := 0; i < PickupTurns; i++ {
				ms[sectionOf(t, g.N())] += em.MutualCircles(PickupRadius, g.R[p], g.Z[p]-pk.z(i))
			}
		}
		for s, m := range ms {
			b.Mutual(fmt.Sprintf("KU%d", s), "LPU", "L"+strconv.Itoa(s), m/math.Sqrt(lpu*c.L[s][s]))
		}
		mp := 0.0
		for i := 0; i < PickupTurns; i++ {
			for j := 0; j < PrimTurns; j++ {
				mp += em.MutualCircles(PickupRadius, g.PrimRadius, pk.z(i)-(g.PrimTop-float64(j)*PrimPitch))
			}
		}
		b.Mutual("KUP", "LPU", "LP", mp/math.Sqrt(lpu*c.Lp))
	})
}

// PickupAC returns the load current per volt of drive at f [A/V] and the
// frequency (searched from 0.7 to 1.4 × f1, the coil's main mode: a tightly coupled tuned pickup splits it in two) where it is largest.
func (c *Coil) PickupAC(f1 float64, pk Pickup) (fBest, iBest float64, err error) {
	ckt, err := c.BuildWithPickup(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, f1, pk)
	if err != nil {
		return 0, 0, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 2801, FStart: 0.7 * f1, FStop: 1.4 * f1}, analysis.Options{})
	if err != nil {
		return 0, 0, err
	}
	i, ok := res.ElementCurrent("RLOAD")
	if !ok {
		return 0, 0, fmt.Errorf("no load current")
	}
	for k, f := range res.Freq {
		if a := cmplx.Abs(i[k]); a > iBest {
			fBest, iBest = f, a
		}
	}
	return fBest, iBest, nil
}

// PickupRun is one pulse-driven run with a pickup.
type PickupRun struct {
	Pickup
	Rate                float64 // pulse rate [Hz]
	ILoad, VLoad, PLoad float64 // peak load current [A], peak load voltage [V], mean load power [W] (last 10 %)
	PIn                 float64 // mean pulse-source power (last 10 %) [W]
	PSource             float64 // of which heats the source's own 50 Ω (last 10 %) [W]
	VCoil               float64 // peak coil-end voltage (last 10 %) [V]
	EIn, EDiss, EStore  float64 // ledger over the run [J]
	ELoad, Mismatch     float64
}

// PulsePickup drives the coil with one-sign pulses at `rate` (width 10 %
// of 1/f1) for `cycles` periods and measures the pickup's load.
func (c *Coil) PulsePickup(f1, rate float64, cycles int, pk Pickup) (PickupRun, error) {
	src := circuit.Element{Kind: circuit.VoltageSource, Wave: circuit.Waveform{Kind: circuit.WavePulse,
		Args: []float64{0, PulseVolts, 0, PulseEdge / f1, PulseEdge / f1, 0.1 / f1, 1 / rate}}}
	ckt, err := c.BuildWithPickup(src, f1, pk)
	if err != nil {
		return PickupRun{}, err
	}
	stop := float64(cycles) / f1
	step := 1 / f1 / 50
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, TranOpts)
	if err != nil {
		return PickupRun{}, err
	}
	r := PickupRun{Pickup: pk, Rate: rate}
	il, _ := w.ElementCurrent("RLOAD")
	top, _ := w.NodeVoltage("n" + strconv.Itoa(NumSections()))
	base, _ := w.NodeVoltage("n0")
	var pl, ps, prs []float64
	for i, n := range w.ElementNames {
		switch n {
		case "RLOAD":
			pl = w.P[i]
		case "VP":
			ps = w.P[i]
		case "RS":
			prs = w.P[i]
		}
	}
	var span, eload float64
	for k := 1; k < len(w.Time); k++ {
		dt := w.Time[k] - w.Time[k-1]
		eload += 0.5 * (pl[k] + pl[k-1]) * dt
		if w.Time[k] < 0.9*stop {
			continue
		}
		span += dt
		r.PLoad += 0.5 * (pl[k] + pl[k-1]) * dt
		r.PIn += -0.5 * (ps[k] + ps[k-1]) * dt // a source's absorbed power is negative when it delivers
		r.PSource += 0.5 * (prs[k] + prs[k-1]) * dt
		r.ILoad = math.Max(r.ILoad, math.Abs(il[k]))
		r.VCoil = math.Max(r.VCoil, math.Abs(top[k]-base[k]))
	}
	r.PLoad /= span
	r.PIn /= span
	r.PSource /= span
	r.VLoad = r.ILoad * pk.RLoad
	r.ELoad = eload
	r.EIn, r.EDiss, r.EStore, r.Mismatch = w.Ledger.EnergyIn, w.Ledger.EnergyDissipated, w.Ledger.StoredEnd-w.Ledger.StoredStart, w.Ledger.IntegrationMismatch
	return r, nil
}

// PickupPower is the mean load power per square volt of drive amplitude at
// the frequency where the load current is largest [W/V²], with that
// frequency.
func (c *Coil) PickupPower(f1 float64, pk Pickup) (f, p float64, err error) {
	f, i, err := c.PickupAC(f1, pk)
	return f, 0.5 * i * i * pk.RLoad, err
}

// LoadGrid and CGrid are the load resistances [Ω] and tuning-capacitor
// factors searched.
var (
	LoadGrid = []float64{0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 25}
	CGrid    = []float64{0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.5, 4}
)

// BestPickup searches the load (and, when cFactors has more than one
// entry, the tuning capacitor) for the largest load power per V² of
// drive.
func (c *Coil) BestPickup(f1 float64, base Pickup, cFactors []float64) (Pickup, float64, float64, error) {
	best, bestP, bestF := base, -1.0, 0.0
	for _, cf := range cFactors {
		for _, rl := range LoadGrid {
			pk := base
			pk.CFactor, pk.RLoad = cf, rl
			f, p, err := c.PickupPower(f1, pk)
			if err != nil {
				return pk, 0, 0, err
			}
			if p > bestP {
				best, bestP, bestF = pk, p, f
			}
		}
	}
	return best, bestF, bestP, nil
}

// RunPickups is the secondary-coil study on the ordinary solenoid.
func RunPickups() error {
	c := NewCoil(Solenoid(), false)
	f1, _, q, err := c.Resonance(1e4, 1e8)
	if err != nil {
		return err
	}
	mid, top := 0.5*c.Geom.Z[c.Geom.N()-1], c.Geom.Z[c.Geom.N()-1]
	fmt.Printf("\n== Low-voltage, high-current pickup on the ordinary solenoid (f1 = %.4g Hz) ==\n", f1)
	fmt.Printf("   pickup: %d turns of 2 mm copper, radius %.1f cm, L = %.3f µH, wire %.3f Ω at f1\n",
		PickupTurns, PickupRadius*100, PickupL()*1e6, PickupR(f1))
	fmt.Printf("   resonant tuning capacitor C = 1/((2π f1)² L) = %.4g nF\n", TuningC(f1)*1e9)
	cases := []struct {
		name string
		base Pickup
		cf   []float64
	}{
		{"middle, tuned C = 1/(ω1²L)", Pickup{Name: "middle-tuned", Centre: mid, Tuned: true}, []float64{1}},
		{"middle, capacitor optimized", Pickup{Name: "middle-optimized", Centre: mid, Tuned: true}, CGrid},
		{"middle, untuned (no capacitor)", Pickup{Name: "middle-untuned", Centre: mid}, []float64{1}},
		{"top (high-voltage end), tuned C = 1/(ω1²L)", Pickup{Name: "top-tuned", Centre: top, Tuned: true}, []float64{1}},
	}
	cycles := int(math.Min(600, math.Max(60, 5*q/math.Pi+20)))
	fmt.Printf("   one-sign pulses: %.0f V, %.0f Ω, width 10 %% of 1/f1, at the coupled system's best frequency, %d cycles\n", PulseVolts, SourceOhms, cycles)
	fmt.Printf("   %-44s %8s %9s %9s %10s %10s %10s %10s %12s %8s %10s\n", "pickup", "C [nF]", "load [Ω]", "rate/f1", "I load [A]", "V load [V]", "P load [W]", "P in [W]", "past 50 Ω [W]", "eff [%]", "coil [V]")
	for _, cs := range cases {
		pk, f, _, err := c.BestPickup(f1, cs.base, cs.cf)
		if err != nil {
			return err
		}
		r, err := c.PulsePickup(f1, f, cycles, pk)
		if err != nil {
			return err
		}
		cnf := 0.0
		if pk.Tuned {
			cnf = pk.CFactor * TuningC(f1) * 1e9
		}
		fmt.Printf("   %-44s %8.3g %9.3g %9.4f %10.3f %10.3f %10.4f %10.4f %12.4f %8.1f %10.0f\n", cs.name, cnf, pk.RLoad, f/f1,
			r.ILoad, r.VLoad, r.PLoad, r.PIn, r.PIn-r.PSource, 100*r.PLoad/(r.PIn-r.PSource), r.VCoil)
		fmt.Printf("   %-44s ledger: in %.4g J = heat %.4g J (load %.4g J) + stored %.3g J (+ time-stepping error %.2g J)\n", "", r.EIn, r.EDiss, r.ELoad, r.EStore, r.Mismatch)
		debug.FreeOSMemory()
	}
	fmt.Println("   eff = load power ÷ the power that reaches the primary (past the source's own 50 Ω).")
	fmt.Println("   The load's power is always less than the pulse source's: the pickup converts voltage to current, it does not add energy.")
	return nil
}
