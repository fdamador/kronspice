// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The zero-frequency limit of Nachamkin's field is the spheromak: a ring of
// toroidal current (the "electron ring") inside a sphere of force-free
// magnetic field, B = B0·F with F the l = 1 Beltrami field. Its current
// J = αB/μ0 runs along B, so there is no Lorentz force: the ring does not
// pinch or burst. What ends it is resistance. With resistivity η the
// force-free field decays as one mode (∂B/∂t = (η/μ0)∇²B = −(ηα²/μ0)B):
//
//	τ_B = μ0/(η·α²),   energy lifetime τ_W = τ_B/2
//
// η is Spitzer's resistivity of a hydrogen plasma, η ≈ 5.2×10⁻⁵·lnΛ/T^{3/2}
// Ω·m with T in eV (NRL Plasma Formulary, parallel resistivity; lnΛ = 10
// declared). As a circuit, the ring is an inductor L = 2W/I² carrying its
// toroidal current I, in series with R = P/I², so L/R = τ_B.

// Spitzer returns η [Ω·m] for electron temperature t [eV].
func Spitzer(t float64) float64 { return 5.2e-5 * 10 / math.Pow(t, 1.5) }

// Spheromak is a force-free sphere of radius R and energy W.
type Spheromak struct {
	R, W float64
}

// B0 returns the field scale for energy W: W = (B0²/2μ0)·V·⟨|F|²⟩.
func (s Spheromak) B0() float64 {
	alpha := J1Zero / s.R
	// ⟨|F|²⟩ in the units of Beltrami (max-normalized meanSquare × max |F|²)
	const n = 200
	var sum, w float64
	for i := 0; i < n; i++ {
		r := (float64(i) + 0.5) / n * s.R
		for j := 0; j < n; j++ {
			th := (float64(j) + 0.5) / n * math.Pi
			fr, ft, fp := Beltrami(alpha, r, th)
			dv := r * r * math.Sin(th)
			sum += (fr*fr + ft*ft + fp*fp) * dv
			w += dv
		}
	}
	vol := 4.0 / 3 * math.Pi * s.R * s.R * s.R
	return math.Sqrt(2 * Mu0 * s.W / (vol * sum / w))
}

// RingCurrent is the toroidal current through a meridional half-plane:
// I = (α/μ0)·∫ B_φ dA.
func (s Spheromak) RingCurrent() float64 {
	alpha := J1Zero / s.R
	b0 := s.B0()
	const n = 400
	var sum float64
	for i := 0; i < n; i++ {
		r := (float64(i) + 0.5) / n * s.R
		for j := 0; j < n; j++ {
			th := (float64(j) + 0.5) / n * math.Pi
			_, _, fp := Beltrami(alpha, r, th)
			sum += fp * r * (s.R / n) * (math.Pi / n)
		}
	}
	return alpha / Mu0 * b0 * sum
}

// TauB is the field decay time at electron temperature t [eV].
func (s Spheromak) TauB(t float64) float64 {
	alpha := J1Zero / s.R
	return Mu0 / (Spitzer(t) * alpha * alpha)
}

// TempFor is the temperature [eV] that gives energy lifetime tw [s].
func (s Spheromak) TempFor(tw float64) float64 {
	alpha := J1Zero / s.R
	eta := Mu0 / (2 * tw * alpha * alpha)
	return math.Pow(5.2e-4/eta, 2.0/3)
}

// Decay is a kronspice run of the ring circuit.
type Decay struct {
	T, L, Rc, I0, Tau       float64
	Time, Energy            []float64
	Heat, Stored0, Mismatch float64
	Fitted                  float64 // energy lifetime fitted from the run
}

// RunDecay runs the ring as L (with its current) in series with R in
// kronspice for 4 field lifetimes.
func (s Spheromak) RunDecay(t float64) (Decay, error) {
	d := Decay{T: t, I0: s.RingCurrent(), Tau: s.TauB(t)}
	d.L = 2 * s.W / (d.I0 * d.I0)
	d.Rc = d.L / d.Tau
	ckt, err := circuit.NewBuilder("spheromak ring").
		InductorIC("LR", "a", "0", d.L, d.I0).
		Resistor("RR", "a", "0", d.Rc).BuildTransient()
	if err != nil {
		return d, err
	}
	step := d.Tau / 200
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: 4 * d.Tau, MaxStep: step, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return d, fmt.Errorf("kronspice: %w", err)
	}
	il, _ := wf.ElementCurrent("LR")
	for k := range wf.Time {
		d.Time = append(d.Time, wf.Time[k])
		d.Energy = append(d.Energy, 0.5*d.L*il[k]*il[k])
	}
	// energy lifetime from the run: time for W to fall by e
	for k := range d.Energy {
		if d.Energy[k] < d.Energy[0]/math.E {
			f := (d.Energy[k-1] - d.Energy[0]/math.E) / (d.Energy[k-1] - d.Energy[k])
			d.Fitted = d.Time[k-1] + f*(d.Time[k]-d.Time[k-1])
			break
		}
	}
	l := wf.Ledger
	d.Heat, d.Stored0 = l.EnergyDissipated, l.StoredStart
	d.Mismatch = l.IntegrationMismatch / l.StoredStart
	return d, nil
}
