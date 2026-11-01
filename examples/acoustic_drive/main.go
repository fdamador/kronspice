// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command acoustic_drive is the electrical → mechanical → acoustic chain of
// a single-bubble sonoluminescence flask (CLAUDE.md §3h, acoustically
// driven mode): a 50 Ω generator, an optional series matching inductor, a
// piezoelectric ring transducer (clamped capacitance and force factor, its
// own mechanical resonance and loss), the coupling into the flask's
// breathing mode, and the mode itself (a 3 cm water or heavy-water sphere).
//
// It reports the center pressure per volt, the drive and electrical power
// needed for a 1.3 atm acoustic amplitude (the SBSL range is about
// 1.2–1.5 atm at 20–40 kHz, Brenner, Hilgenfeldt & Lohse, Rev. Mod. Phys.
// 74, 425 (2002)), and where every watt goes: transducer heat, heat in the
// liquid from the mode's losses, and source power, from the energy ledger.
//
// The bubble itself (Keller–Miksis dynamics, collapse, light) is milestone
// M9p; this chain ends at the pressure the bubble would be driven by.
//
//	go run ./examples/acoustic_drive
package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/fluid"
	"github.com/fdamador/kronspice/output"
)

// Assumed transducer and flask (declared, not measured: replace with the
// bench values; the coupling area and Q in particular need a hydrophone
// calibration of the real flask).
const (
	radius = 0.03   // flask radius [m]
	q      = 1500.0 // breathing-mode quality factor (assumed)
	sEff   = 2e-4   // effective coupling area of the transducer into the mode [m²] (assumed)
	c0     = 4e-9   // transducer clamped capacitance [F]
	phi    = 0.6    // force factor [N/V]
	mass   = 4e-3   // transducer moving mass [kg]
	bMech  = 3.0    // transducer mechanical loss [N·s/m]
	target = 1.3 * 101325
)

// chain builds the circuit; withL adds a series inductor tuning out C0 at f0.
func chain(m fluid.Medium, withL bool) (*circuit.Circuit, float64, error) {
	f0, Ca, Ma, Rpar, err := fluid.SphericalResonator(radius, m, q)
	if err != nil {
		return nil, 0, err
	}
	w0 := 2 * math.Pi * f0
	k := w0 * w0 * mass // transducer resonance placed at the mode frequency
	b := circuit.NewBuilder("acoustic drive "+m.Name).
		Domain(exchange.Mechanical, "v").Domain(exchange.Acoustic, "p").
		Add(circuit.Element{Name: "VGEN", Kind: circuit.VoltageSource, ACMag: 1}, "g", "0").
		Resistor("RGEN", "g", "a", 50)
	if withL {
		b.Inductor("LMATCH", "a", "e", 1/(w0*w0*c0)).Resistor("RLMATCH", "e", "x", 0.5)
	} else {
		b.Resistor("RWIRE", "a", "x", 0.5)
	}
	b.Capacitor("C0", "x", "0", c0).Gyrator("YPZT", "x", "0", "v", "0", phi).
		Capacitor("CMASS", "v", "0", mass).Inductor("LSPRING", "v", "0", 1/k).Resistor("RLOSS", "v", "0", 1/bMech).
		Gyrator("YCOUPLE", "v", "0", "p", "0", sEff).
		Capacitor("CMODE", "p", "0", Ca).Inductor("LMODE", "p", "0", Ma).Resistor("RMODE", "p", "0", Rpar)
	c, err := b.Build()
	return c, f0, err
}

func main() {
	fmt.Printf("Single-bubble sonoluminescence drive chain, %.0f mm flask, Q = %.0f (assumed), S_eff = %.1f cm²\n\n",
		2e3*radius, q, sEff*1e4)
	for _, m := range []fluid.Medium{fluid.Water, fluid.HeavyWater} {
		for _, withL := range []bool{false, true} {
			if err := report(m, withL); err != nil {
				fmt.Fprintln(os.Stderr, "acoustic_drive:", err)
				os.Exit(1)
			}
		}
	}
}

func report(m fluid.Medium, withL bool) error {
	c, f0, err := chain(m, withL)
	if err != nil {
		return err
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 2001, FStart: 0.97 * f0, FStop: 1.03 * f0}, analysis.Options{})
	if err != nil {
		return err
	}
	p, _ := r.NodeVoltage("p")
	best := 0
	for i := range p {
		if cmplx.Abs(p[i]) > cmplx.Abs(p[best]) {
			best = i
		}
	}
	perVolt := cmplx.Abs(p[best])
	drive := target / perVolt // generator amplitude [V] for the target pressure
	scale := drive * drive    // powers scale with the square of the drive
	l := r.Ledger
	heat := map[string]float64{}
	for e, name := range r.ElementNames {
		ie := r.I[e][best]
		var v complex128
		el := c.Elements[c.Index(name)]
		v = node(r, c, el.Nodes[0], best) - node(r, c, el.Nodes[1], best)
		if el.Kind == circuit.Resistor && name != "RGEN" {
			heat[name] = 0.5 * real(v*cmplx.Conj(ie)) * scale
		}
	}
	match := "no matching inductor"
	if withL {
		match = "series inductor tuning out C0"
	}
	fmt.Printf("%s, %s\n", m.Name, match)
	fmt.Printf("  mode f0 = %s (peak at %s), center pressure %s per volt of generator EMF\n",
		output.Eng(f0, "Hz"), output.Eng(r.Freq[best], "Hz"), output.Eng(perVolt, "Pa"))
	fmt.Printf("  for 1.3 atm: generator EMF %s amplitude, source power %s\n",
		output.Eng(drive, "V"), output.Eng(l.Inflow[best]*scale, "W"))
	fmt.Printf("  into the liquid (mode losses) %s, transducer loss %s, wiring %s, generator's own 50 Ω %s\n",
		output.Eng(heat["RMODE"], "W"), output.Eng(heat["RLOSS"], "W"),
		output.Eng(heat["RWIRE"]+heat["RLMATCH"], "W"),
		output.Eng(l.Dissipation[best]*scale-heat["RMODE"]-heat["RLOSS"]-heat["RWIRE"]-heat["RLMATCH"], "W"))
	fmt.Printf("  ledger: in %s, out as heat %s, residual %.2g of the scale\n\n",
		output.Eng(l.Inflow[best]*scale, "W"), output.Eng(l.Dissipation[best]*scale, "W"), l.MaxRelativeResidual)
	return nil
}

func node(r *output.ACResult, c *circuit.Circuit, n, k int) complex128 {
	if n == 0 {
		return 0
	}
	return r.V[n][k]
}
