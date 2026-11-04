// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Layer 2 of the Earth machine: the global atmospheric electric circuit
// (C. T. R. Wilson 1920; M. J. Rycroft, S. Israelsson, C. Price, "The global
// atmospheric electric circuit, solar activity and climate change," J. Atmos.
// Sol.-Terr. Phys. 62, 1563, 2000).
//
// Thunderstorms and electrified clouds act as current generators that carry
// positive charge up to the ionosphere, an equipotential at about 250 kV; the
// return current flows down through the weakly conducting fair-weather air
// everywhere else. The air is a leaky dielectric: each spherical shell of it
// is a resistor in parallel with a capacitor (Maxwell's total current
// ℭ = 𝔎 + ∂𝔇/∂t, CLAUDE.md §3c). Stacking the shells gives a Kron ladder
// from the ground (node 0) to the ionosphere (node "ion"), and the storms
// are a current source in parallel with the whole ladder.
//
// Conductivity rises with height, σ(z) = σ0·exp(z/H). The model is
// calibrated on three measured anchors, declared below: the ionospheric
// potential V_I, the fair-weather air–earth current density J and the
// fair-weather surface field E0. σ0 = J/E0 and H follow; everything else
// (global current and resistance, the Earth's charge, where the balancing
// positive charge sits, the relaxation time, the response to the daily
// thunderstorm cycle) is computed.

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

const (
	eps0   = 8.8541878128e-12 // F/m
	mu0    = 1.25663706212e-6 // H/m
	cLight = 299792458.0      // m/s
	rEarth = 6.371e6          // m, mean radius
)

// Atmosphere is the declared input of layer 2.
type Atmosphere struct {
	VI     float64 // ionospheric potential [V]; ≈ 250 kV (Markson 2007, Rycroft et al. 2000)
	J      float64 // fair-weather air–earth current density [A/m²]; ≈ 2 pA/m² (Rycroft et al. 2000)
	E0     float64 // fair-weather surface field magnitude [V/m]; ≈ 120 V/m (Chalmers 1967, typical)
	Top    float64 // height of the ionospheric equipotential [m]
	Layers int     // shells in the ladder
}

// StandardAtmosphere returns the declared anchors.
func StandardAtmosphere() Atmosphere {
	return Atmosphere{VI: 250e3, J: 2e-12, E0: 120, Top: 60e3, Layers: 60}
}

// Calibrate returns σ0 = J/E0 and the scale height H for which a flat
// column of σ0·exp(z/H) up to Top has V_I = J·∫dz/σ.
func (a Atmosphere) Calibrate() (sigma0, H float64) {
	sigma0 = a.J / a.E0
	H = a.VI / a.E0
	for k := 0; k < 50; k++ { // H = V_I / (E0·(1 − e^(−Top/H)))
		H = a.VI / (a.E0 * (1 - math.Exp(-a.Top/H)))
	}
	return
}

// Shell is one layer of the ladder.
type Shell struct {
	Z0, Z1 float64 // bottom and top heights [m]
	R, C   float64 // global resistance [Ω] and capacitance [F]
}

// Shells discretizes the atmosphere, thin near the ground where the
// resistance is.
func (a Atmosphere) Shells() []Shell {
	s0, H := a.Calibrate()
	const stretch = 8.0
	z := func(k int) float64 {
		return a.Top * (math.Exp(stretch*float64(k)/float64(a.Layers)) - 1) / (math.Exp(stretch) - 1)
	}
	out := make([]Shell, a.Layers)
	for k := range out {
		z0, z1 := z(k), z(k+1)
		r0, r1 := rEarth+z0, rEarth+z1
		rm := 0.5 * (r0 + r1)
		area := 4 * math.Pi * rm * rm
		out[k] = Shell{
			Z0: z0, Z1: z1,
			R: H / (s0 * area) * (math.Exp(-z0/H) - math.Exp(-z1/H)), // ∫dz/σ over the shell area
			C: 4 * math.Pi * eps0 * r0 * r1 / (r1 - r0),              // spherical shell capacitor
		}
	}
	return out
}

// GlobalCurrent is the storm current that balances the fair-weather
// return current: J times the Earth's surface.
func (a Atmosphere) GlobalCurrent() float64 { return a.J * 4 * math.Pi * rEarth * rEarth }

// node names the junction at the bottom of shell k (k = 0 is the ground).
func (a Atmosphere) node(k int) string {
	switch {
	case k == 0:
		return "0"
	case k == a.Layers:
		return "ion"
	}
	return fmt.Sprintf("z%d", k)
}

// Circuit builds the ladder with the storm generator. storm is its
// waveform (nil: DC at the global current); acMag sets its small-signal
// amplitude.
func (a Atmosphere) Circuit(storm *circuit.Waveform, acMag float64) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("Global atmospheric electric circuit")
	for k, s := range a.Shells() {
		lo, hi := a.node(k), a.node(k+1)
		b.Resistor(fmt.Sprintf("R%d", k), hi, lo, s.R)
		b.Capacitor(fmt.Sprintf("C%d", k), hi, lo, s.C)
	}
	// Thunderstorms: a convective heat engine (latent heat) that lifts
	// positive charge from the ground to the ionosphere. SPICE convention:
	// the current is drawn from n+ (ground) and delivered into n− (ion).
	e := circuit.Element{Name: "ISTORM", Kind: circuit.CurrentSource, Value: a.GlobalCurrent(), HasDC: true,
		ACMag: acMag, Channel: exchange.Thermal}
	if storm != nil {
		e.Wave = *storm
	}
	b.Add(e, "0", "ion")
	return b.Build()
}

// AtmosphereReport is the computed state of layer 2.
type AtmosphereReport struct {
	Sigma0, H         float64 // calibrated profile
	I, R              float64 // global current [A] and resistance [Ω]
	VI, VIExact       float64 // circuit and exact ladder (I·ΣR) ionospheric potential
	QEarth, QEarthGeo float64 // Earth's charge from the circuit and −ε0·E0·4πR²
	QSpace, QIon      float64 // positive space charge in the air, charge at the top
	ZHalf             float64 // height below which half the space charge lies
	CGeo, CEff        float64 // geometric capacitance to Top, Q/V_I
	Tau               float64 // ε0/σ0
	Power             float64 // dissipated power [W]
	LedgerRes         float64 // relative residual of the engine's OP ledger
}

// Solve runs the operating point and books the charges.
func (a Atmosphere) Solve() (AtmosphereReport, error) {
	var r AtmosphereReport
	r.Sigma0, r.H = a.Calibrate()
	r.I = a.GlobalCurrent()
	c, err := a.Circuit(nil, 0)
	if err != nil {
		return r, err
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return r, err
	}
	shells := a.Shells()
	v := make([]float64, len(shells)) // voltage across each shell
	for k := range shells {
		hi, _ := op.Voltage(a.node(k + 1))
		lo, _ := op.Voltage(a.node(k))
		v[k] = hi - lo
		r.VIExact += r.I * shells[k].R
	}
	r.VI, _ = op.Voltage("ion")
	r.R = r.VIExact / r.I
	n := len(shells) - 1
	r.QEarth = -shells[0].C * v[0]
	r.QIon = shells[n].C * v[n]
	r.QSpace = -r.QEarth - r.QIon
	r.QEarthGeo = -eps0 * a.E0 * 4 * math.Pi * rEarth * rEarth
	// Junction k (bottom of shell k) holds C_{k−1}v_{k−1} − C_k v_k.
	var cum float64
	for k := 1; k <= n; k++ {
		cum += shells[k-1].C*v[k-1] - shells[k].C*v[k]
		if cum >= r.QSpace/2 {
			r.ZHalf = shells[k].Z0
			break
		}
	}
	rt := rEarth + a.Top
	r.CGeo = 4 * math.Pi * eps0 * rEarth * rt / a.Top
	r.CEff = -r.QEarth / r.VI
	r.Tau = eps0 / r.Sigma0
	r.Power = r.VI * r.I
	if op.Ledger != nil && op.Ledger.Scale > 0 {
		r.LedgerRes = op.Ledger.Residual / op.Ledger.Scale
	}
	return r, nil
}

// SwitchOn charges an uncharged Earth with the storms switched on at t = 0
// and returns V_I(t) from the circuit and from the exact ladder solution:
// a current source into series RC shells charges each shell independently,
// v_k(t) = I·R_k·(1 − exp(−t/R_kC_k)).
func (a Atmosphere) SwitchOn(tEnd, dt float64) (t, vi, exact []float64, err error) {
	c, err := a.Circuit(nil, 0)
	if err != nil {
		return
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: tEnd, MaxStep: dt, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return
	}
	t = wf.Time
	vi, _ = wf.NodeVoltage("ion")
	I := a.GlobalCurrent()
	shells := a.Shells()
	for _, tk := range t {
		var s float64
		for _, sh := range shells {
			s += I * sh.R * -math.Expm1(-tk/(sh.R*sh.C))
		}
		exact = append(exact, s)
	}
	return
}

// Carnegie returns the amplitude ratio and the lag [s] of V_I behind a
// daily modulation of the storm current (the Carnegie curve; R. G. Harrison,
// "The Carnegie curve," Surv. Geophys. 34, 209, 2013), from the engine's AC
// analysis at one cycle per day: |ΔV_I| / (R·|ΔI|) and −phase/ω.
func (a Atmosphere) Carnegie(rel float64) (ratio, lag float64, err error) {
	I := a.GlobalCurrent()
	c, err := a.Circuit(nil, rel*I)
	if err != nil {
		return
	}
	f := 1.0 / 86400
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		return
	}
	v, _ := ac.NodeVoltage("ion")
	var rsum float64
	for _, s := range a.Shells() {
		rsum += s.R
	}
	ratio = cmplx.Abs(v[0]) / (rsum * rel * I)
	lag = -cmplx.Phase(v[0]) / (2 * math.Pi * f)
	return
}

func runAtmosphere() error {
	a := StandardAtmosphere()
	fmt.Println("Layer 2 — the global atmospheric electric circuit as a Kron ladder of leaky shells")
	fmt.Printf("  anchors (declared, measured): V_I = %.0f kV, J = %.1f pA/m², E0 = %.0f V/m, ionosphere at %.0f km\n",
		a.VI/1e3, a.J*1e12, a.E0, a.Top/1e3)
	r, err := a.Solve()
	if err != nil {
		return err
	}
	fmt.Printf("  calibrated profile σ(z) = σ0·exp(z/H): σ0 = %.3g S/m, H = %.2f km (%d shells)\n", r.Sigma0, r.H/1e3, a.Layers)
	fmt.Printf("  global storm current I = J·4πR² = %.0f A; global resistance %.1f Ω\n", r.I, r.R)
	fmt.Printf("  V_I: circuit %.2f kV, exact ladder %.2f kV (rel. diff %.1g); flat column %.0f kV (sphericity −%.2f %%)\n",
		r.VI/1e3, r.VIExact/1e3, math.Abs(r.VI-r.VIExact)/r.VIExact, a.VI/1e3, 100*(1-r.VIExact/a.VI))
	fmt.Printf("  Earth's charge: circuit %.3g C, −ε0·E0·4πR² = %.3g C\n", r.QEarth, r.QEarthGeo)
	fmt.Printf("  balancing positive charge: %.3g C in the air (half of it below %.2f km), %.2g C at the ionosphere\n",
		r.QSpace, r.ZHalf/1e3, r.QIon)
	fmt.Printf("  capacitance: geometric shell to %.0f km %.3f F; effective Q/V_I = %.2f F\n", a.Top/1e3, r.CGeo, r.CEff)
	fmt.Printf("  relaxation: R·C_eff = %.0f s = ε0/σ0 = %.0f s (%.1f min) — without storms the charge would leak away in minutes\n",
		r.R*r.CEff, r.Tau, r.Tau/60)
	fmt.Printf("  power: V_I·I = %.0f MW from the storms, all dissipated as heat; engine ledger residual %.1g\n", r.Power/1e6, r.LedgerRes)

	t, vi, ex, err := a.SwitchOn(4000, 1)
	if err != nil {
		return err
	}
	var dmax float64
	t63 := math.NaN()
	for k := range t {
		dmax = math.Max(dmax, math.Abs(vi[k]-ex[k])/r.VIExact)
		if math.IsNaN(t63) && vi[k] >= (1-math.Exp(-1))*r.VIExact {
			t63 = t[k]
		}
	}
	fmt.Printf("  switch-on transient: V_I reaches 63 %% in %.0f s; circuit vs exact ladder within %.1g of V_I\n", t63, dmax)
	ratio, lag, err := a.Carnegie(0.15)
	if err != nil {
		return err
	}
	fmt.Printf("  Carnegie curve (±15 %% daily storm current): V_I follows with amplitude ratio %.6f and a lag of %.1f min\n\n", ratio, lag/60)
	return nil
}
