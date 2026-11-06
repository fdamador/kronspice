// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/fluid"
	"github.com/fdamador/kronspice/output"
)

// The drive chain of examples/acoustic_drive (matched, water), with the
// same declared transducer and flask values, and the bubble added at the
// centre node p as its fundamental admittance Y₁ = G + jB: a resistor 1/G
// and, for B > 0, a capacitor B/ω (for B < 0, an inductor −1/(ωB)). The
// flask mode is a narrow resonance, so the centre pressure is the
// fundamental; the bubble's harmonics (its sound bursts) see the mode's
// off-resonance impedance and are booked in the bubble's own ledger as
// radiated sound. This is a harmonic balance at the fundamental: exact for
// the power the bubble draws (only the fundamental of U does work against
// a sinusoidal pressure), first-order for the frequency shift.
const (
	radius = 0.03   // flask radius [m]
	qMode  = 1500.0 // breathing-mode Q (assumed, as in acoustic_drive)
	sEff   = 2e-4   // coupling area [m²] (assumed)
	c0     = 4e-9   // transducer clamped capacitance [F]
	phi    = 0.6    // force factor [N/V]
	mass   = 4e-3   // transducer moving mass [kg]
	bMech  = 3.0    // transducer mechanical loss [N·s/m]
)

// chain builds the matched water chain, with the bubble's admittance at
// the centre when y ≠ 0.
func chain(y complex128) (*circuit.Circuit, float64, error) {
	m := fluid.Water
	f0, Ca, Ma, Rpar, err := fluid.SphericalResonator(radius, m, qMode)
	if err != nil {
		return nil, 0, err
	}
	w0 := 2 * math.Pi * f0
	k := w0 * w0 * mass
	name := "acoustic drive with bubble"
	if y == 0 {
		name = "acoustic drive"
	}
	b := circuit.NewBuilder(name).
		Domain(exchange.Mechanical, "v").Domain(exchange.Acoustic, "p").
		Add(circuit.Element{Name: "VGEN", Kind: circuit.VoltageSource, ACMag: 1}, "g", "0").
		Resistor("RGEN", "g", "a", 50).
		Inductor("LMATCH", "a", "e", 1/(w0*w0*c0)).Resistor("RLMATCH", "e", "x", 0.5).
		Capacitor("C0", "x", "0", c0).Gyrator("YPZT", "x", "0", "v", "0", phi).
		Capacitor("CMASS", "v", "0", mass).Inductor("LSPRING", "v", "0", 1/k).Resistor("RLOSS", "v", "0", 1/bMech).
		Gyrator("YCOUPLE", "v", "0", "p", "0", sEff).
		Capacitor("CMODE", "p", "0", Ca).Inductor("LMODE", "p", "0", Ma).Resistor("RMODE", "p", "0", Rpar)
	if y != 0 {
		b.Resistor("RBUBBLE", "p", "0", 1/real(y))
		if imag(y) > 0 {
			b.Capacitor("CBUBBLE", "p", "0", imag(y)/w0)
		} else if imag(y) < 0 {
			b.Inductor("LBUBBLE", "p", "0", -1/(w0*imag(y)))
		}
	}
	c, err := b.Build()
	return c, f0, err
}

// Drive is the chain read at its peak for a target centre pressure.
type Drive struct {
	F0, FPeak float64 // mode and peak frequency [Hz]
	PerVolt   float64 // centre pressure per volt [Pa/V]
	EMF       float64 // generator amplitude for the target [V]
	Source    float64 // source power at that drive [W]
	Mode      float64 // heat in the liquid from the mode's losses [W]
	Bubble    float64 // power into the bubble element [W]
	Residual  float64 // ledger relative residual
}

func readDrive(y complex128, target float64) (Drive, error) {
	c, f0, err := chain(y)
	if err != nil {
		return Drive{}, err
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 4001, FStart: 0.995 * f0, FStop: 1.005 * f0}, analysis.Options{})
	if err != nil {
		return Drive{}, err
	}
	p, _ := r.NodeVoltage("p")
	best := 0
	for i := range p {
		if cmplx.Abs(p[i]) > cmplx.Abs(p[best]) {
			best = i
		}
	}
	d := Drive{F0: f0, FPeak: r.Freq[best], PerVolt: cmplx.Abs(p[best])}
	d.EMF = target / d.PerVolt
	scale := d.EMF * d.EMF
	d.Source = r.Ledger.Inflow[best] * scale
	d.Residual = r.Ledger.MaxRelativeResidual
	heat := func(name string) float64 {
		i := c.Index(name)
		if i < 0 {
			return 0
		}
		for e, n := range r.ElementNames {
			if n == name {
				el := c.Elements[i]
				v := nodeV(r, el.Nodes[0], best) - nodeV(r, el.Nodes[1], best)
				return 0.5 * real(v*cmplx.Conj(r.I[e][best])) * scale
			}
		}
		return 0
	}
	d.Mode = heat("RMODE")
	d.Bubble = heat("RBUBBLE")
	return d, nil
}

func nodeV(r *output.ACResult, n, k int) complex128 {
	if n == 0 {
		return 0
	}
	return r.V[n][k]
}
