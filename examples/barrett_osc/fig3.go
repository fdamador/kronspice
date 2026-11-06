// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// Fig. 3 of Barrett 1991: three capacitors charged −, +, − joined by switches
// in a chain to a line L that runs to G; closing the switches in turn (the
// arrows point toward L) sends out a − + − pulse (the drawn waveform). The
// switches are spark gaps: the gap model and solver are the ones in
// examples/sparkgap_diodebridge (copied in gap.go, sim.go).
//
// Netlist (declared values in FrozenWave):
//
//	G1  bus–n1   C1 n1–0 (−V)
//	G2  n1–n2    C2 n2–0 (+V)
//	G3  n2–n3    C3 n3–0 (−V)
//	L   bus–g    R  g–0  (the load G)
//
// The gaps are dry-air gaps that hold the static voltage across them (2V
// between neighbours) and are fired by triggers at Trig[k]. A gap that has
// just conducted re-strikes at a reduced voltage while its channel is hot
// (Recovery), so the gaps nearer L carry the later charges through.
type FrozenWave struct {
	V        float64    // charge voltage [V]
	C        float64    // each capacitor [F]
	L, R     float64    // line inductance [H], load G [Ω]
	GapD     float64    // air gap [m]
	Recovery float64    // hot-channel recovery time [s]
	Trig     [3]float64 // trigger times of G1, G2, G3 [s]
	Step     float64    // solver step [s]
	Stop     float64    // [s]
}

// DefaultFrozenWave: 10 kV on 1 nF each, 1 µH into 50 Ω, 1 cm air gaps
// (Paschen ≈ 31 kV, so they hold 20 kV), triggers 100 ns apart (about half
// a period of 1 µH with 1 nF), hot channel recovering over 1 µs.
var DefaultFrozenWave = FrozenWave{V: 10e3, C: 1e-9, L: 1e-6, R: 50, GapD: 0.01, Recovery: 1e-6,
	Trig: [3]float64{0, 100e-9, 200e-9}, Step: 0.01e-9, Stop: 600e-9}

// FrozenResult is one run.
type FrozenResult struct {
	T, Vload       []float64
	Lobes          []float64 // signed peak of each lobe of the load voltage, in order [V]
	Stored0        float64   // Σ ½CV² at the start [J]
	Load, GapHeat  float64   // energy into R and dissipated in the gaps [J]
	Left, NumLoss  float64   // stored at the end, backward-Euler damping [J]
	Residual       float64   // ledger residual [J]
	Strikes        [3]int
	BreakdownVolts float64
}

// Run steps the circuit.
func (f FrozenWave) Run() (FrozenResult, error) {
	var r FrozenResult
	gap := AirGap(f.GapD)
	gap.Recovery = f.Recovery
	gap.Hold = 5e-9
	r.BreakdownVolts = gap.Vpos()
	n := NewNet()
	g := [3]*part{
		n.Gap("G1", "bus", "n1", gap, "gaps"),
		n.Gap("G2", "n1", "n2", gap, "gaps"),
		n.Gap("G3", "n2", "n3", gap, "gaps"),
	}
	n.C("C1", "n1", "0", f.C, -f.V)
	n.C("C2", "n2", "0", f.C, +f.V)
	n.C("C3", "n3", "0", f.C, -f.V)
	n.L("L", "bus", "g", f.L)
	n.R("R", "g", "0", f.R).group = "load"
	n.Prepare()
	steps := int(math.Round(f.Stop / f.Step))
	for k := 0; k < steps; k++ {
		t := float64(k) * f.Step
		for j := range g {
			if t <= f.Trig[j] && f.Trig[j] < t+f.Step {
				g[j].armed = true
			}
		}
		if err := n.Step(t, f.Step); err != nil {
			return r, err
		}
		r.T = append(r.T, t+f.Step)
		r.Vload = append(r.Vload, n.Volt("g"))
	}
	l := n.Ledger()
	r.Stored0, r.Left, r.NumLoss, r.Residual = l.Stored0, l.Stored, l.NumLoss, l.Residual
	r.Load, r.GapHeat = l.Heat["load"], l.Heat["gaps"]
	for j := range g {
		r.Strikes[j] = g[j].strikesF + g[j].strikesB
	}
	r.Lobes = lobes(r.Vload, 0.05)
	return r, nil
}

// lobes returns the signed peak of each excursion of x that exceeds frac of
// max|x|, in time order.
func lobes(x []float64, frac float64) []float64 {
	m := 0.0
	for _, v := range x {
		m = math.Max(m, math.Abs(v))
	}
	var out []float64
	cur := 0.0
	for _, v := range x {
		if math.Abs(v) < frac*m {
			continue
		}
		if cur != 0 && math.Signbit(v) != math.Signbit(cur) {
			out = append(out, cur)
			cur = 0
		}
		if math.Abs(v) > math.Abs(cur) {
			cur = v
		}
	}
	if cur != 0 {
		out = append(out, cur)
	}
	return out
}

func printFig3() error {
	f := DefaultFrozenWave
	r, err := f.Run()
	if err != nil {
		return err
	}
	fmt.Println("Fig. 3 — three charged capacitors (−, +, −) switched in turn by spark gaps through L into G")
	fmt.Printf("  %.0f kV on %.0f nF each, L %.0f µH, G %.0f Ω, %.0f mm air gaps (breakdown %.1f kV), triggers at %.0f/%.0f/%.0f ns\n",
		f.V/1e3, f.C*1e9, f.L*1e6, f.R, f.GapD*1e3, r.BreakdownVolts/1e3, f.Trig[0]*1e9, f.Trig[1]*1e9, f.Trig[2]*1e9)
	fmt.Print("  load-voltage lobes in order [kV]:")
	for _, v := range r.Lobes {
		fmt.Printf(" %+.2f", v/1e3)
	}
	fmt.Printf("\n  strikes per gap: %v\n", r.Strikes)
	fmt.Printf("  ledger: stored %.2f mJ = load %.2f mJ + gaps %.2f mJ + left %.3f mJ + solver damping %.3f mJ (residual %.1e J)\n",
		r.Stored0*1e3, r.Load*1e3, r.GapHeat*1e3, r.Left*1e3, r.NumLoss*1e3, r.Residual)
	fmt.Println("  standard prediction: the drawn − + − pulse is what sequentially switched charged capacitors give")
	fmt.Println("  (a lumped frozen-wave generator, the impulse source of later ultrawideband radar); the load gets")
	fmt.Println("  less than the ½CV² stored: when a gap fires, the next capacitor shares charge with the one before")
	fmt.Println("  it with no inductance between them, and that loss goes to the arcs (why the + and second − lobes")
	fmt.Println("  are smaller; transmission-line segments in place of capacitors avoid it).")
	fmt.Println()
	return nil
}

// AirGap is a symmetric dry-air spark gap of d [m] at 1 atm: arc 30 V +
// 0.2 Ω, out below 0.1 A (declared). Copied from
// examples/sparkgap_diodebridge/apps.go.
func AirGap(d float64) PolarGap {
	return PolarGap{Fixed: PaschenAir(760, d*100), NegRatio: 1, Tip: 1e-3, Sigma: 1e-12, Varc: 30, Rarc: 0.2, Iext: 0.1, Floor: 330}
}

// PaschenAir is the breakdown voltage of air at pressure p [Torr] and gap
// d [cm] (Raizer, Gas Discharge Physics, 1991, Table 4.1 and eq. 7.3,
// γ = 0.01 declared). Copied from examples/sparkgap_diodebridge/apps.go.
func PaschenAir(p, d float64) float64 {
	pd := p * d
	return 365 * pd / math.Log(15*pd/math.Log(1+1/0.01))
}
