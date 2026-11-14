// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// Part 8b: the spinning core as a battery. In the membrane paradigm a
// spinning horizon threaded by field lines is a generator (Blandford &
// Znajek, MNRAS 179, 433, 1977; Thorne, Price & Macdonald 1986, ch. IV):
// the rotation drives an EMF ∝ Ω_H·Ψ round a circuit whose internal
// resistance is the horizon (its ~377 Ω per square, summed over the
// geometry into R_H) and whose external load is whatever the field lines
// feed — jets, or here the photon ring. The load's voltage is ∝ Ω_F·Ψ
// (Ω_F the field lines' angular velocity), so Ω_F/Ω_H = R_L/(R_H + R_L), and
// the power to the load is greatest at the match R_L = R_H, Ω_F = Ω_H/2,
// with half the power heating the horizon.
//
// In kronspice: VSPIN (the spin's EMF, booked as the spin column), RH, RL.
// Two views of the same core on the Smith chart:
//   - sources off (analysis.Smith, as a network analyser sees it): Z = +R_H,
//     inside the unit circle — the horizon is a resistor that absorbs;
//   - at the operating point, seen from the load: Z = V/I_in = −R_L, outside
//     the circle — the core supplies power, paid by the spin. Every point
//     meets Γ_core·Γ_load = 1 (Kurokawa's condition of a steady state).

// SpinPoint is the circuit at one load.
type SpinPoint struct {
	RLoverRH        float64
	OmegaF          float64 // Ω_F/Ω_H
	PLoad, PHorizon float64 // [W] for EMF 1 V, R_H = 1 Ω
	PSpin           float64 // power the EMF delivers (the spin column)
	Residual        float64 // kronspice ledger, relative
	GammaLoad       complex128
	GammaCore       complex128 // of Z = −R_L against Z0 = R_H
}

// SpinSweep runs the circuit in kronspice for each R_L/R_H.
func SpinSweep(ratios []float64) ([]SpinPoint, error) {
	out := make([]SpinPoint, 0, len(ratios))
	for i, x := range ratios {
		b := circuit.NewBuilder("spinning core " + strconv.Itoa(i))
		b.VoltageSource("VSPIN", "e", "0", 1) // EMF ∝ Ω_H·Ψ, normalised to 1 V
		b.Resistor("RH", "e", "m", 1)         // the horizon, normalised to 1 Ω
		b.Resistor("RL", "m", "0", x)         // the load: jets, or the photon ring
		c, err := b.Build()
		if err != nil {
			return nil, err
		}
		r, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			return nil, err
		}
		p := SpinPoint{RLoverRH: x}
		for _, e := range r.Elements {
			switch e.Name {
			case "vspin", "VSPIN":
				p.PSpin = -e.PowerAbsorbed
			case "rh", "RH":
				p.PHorizon = e.PowerAbsorbed
			case "rl", "RL":
				p.PLoad = e.PowerAbsorbed
			}
		}
		vm, _ := r.Voltage("m")
		p.OmegaF = vm // load voltage / EMF
		p.Residual = r.Ledger.RelativeResidual()
		p.GammaLoad = output.SmithGamma(complex(x, 0), 1)
		p.GammaCore = output.SmithGamma(complex(-x, 0), 1)
		out = append(out, p)
	}
	return out, nil
}

// CoreSmith measures the core alone (VSPIN + RH, port m–0) with kronspice's
// Smith sweep: sources off, so it sees +R_H.
func CoreSmith() (*output.SmithResult, error) {
	b := circuit.NewBuilder("spinning core, port m-0")
	b.VoltageSource("VSPIN", "e", "0", 1)
	b.Resistor("RH", "e", "m", 1)
	b.Resistor("RBLEED", "m", "0", 1e9)
	c, err := b.Build()
	if err != nil {
		return nil, err
	}
	return analysis.Smith(c, "m", "0", 1, analysis.ACSpec{Variation: "lin", Points: 3, FStart: 1, FStop: 3}, analysis.Options{})
}

// Superradiance for the photon ring. Its wave comes back to itself only
// after two trips (photon_electron's 720°), so its azimuthal number is
// m = ½; a spinning absorber amplifies it when ω < m·Ω (Zel'dovich, JETP
// Lett. 14, 180, 1971): Ω > 2ω_C.
const RingM = 0.5

// SuperradianceOmega is the least core spin that amplifies the ring mode [rad/s].
func SuperradianceOmega() float64 { return OmegaC() / RingM }

// RimSpeed is the core's rim speed at that spin for a core of radius rc [m/s];
// it equals c·rc/R, since ω_C·R = c/2.
func RimSpeed(rc float64) float64 { return SuperradianceOmega() * rc }

// RingQ is the photon ring's own Q from photon_electron's field-solver
// run (README §3): 5.8 % of the energy left within 2R after a quarter
// period, W ∝ exp(−ωt/Q) ⇒ Q = (π/2)/ln(1/0.058).
func RingQ() float64 { return (math.Pi / 2) / math.Log(1/0.058) }

// SpinLifetime is how long a core of spin ħ/2 turning at Ω can keep the
// ring topped up through a matched load (half the drawn power heats the
// horizon): E_rot = ħΩ/4, drawn at 2·ω_C·m_e c²/Q [s].
func SpinLifetime(omega float64) float64 {
	erot := Hbar * omega / 4
	pRing := OmegaC() * Me * C0 * C0 / RingQ()
	return erot / (2 * pRing)
}
