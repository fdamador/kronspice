// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Part 12: the donut void — a hollow vortex ring in the ether with the
// photon running round inside its empty core.
//
// Ring radius R, hollow core radius a, aspect t = a/R. Energies:
//
//	light   E_L = ħc·√((x/a)² + (1/2R)²)   lowest TE₁₁ mode of the hollow
//	        core (x = 1.8412), wavenumber 1/2R along the ring: the wave
//	        turns over once per circuit (the crossed ring), so it carries
//	        ħ/2 about the axis for every R and a;
//	vortex  E_V = K·R·(ln(8R/a) − 2)       hollow vortex ring (Lamb §163),
//	        K = ½ρκ²: line energy times length, ∝ R at fixed shape;
//	ether   E_P = p∞·2π²R·a²               work to empty the core;
//	charge  U   = α·ħc·⟨1/d⟩/2             e spread over the core's wall
//	        (TorusInvD), its static field.
//
// Lengths in ℓ = √(ħc/K), energies in E0 = √(ħcK); the one parameter is
// Π = p∞·ħc/K². At fixed shape t every term is a power of the size λ = R/ℓ:
//
//	e(λ, t) = A(t)/λ + B(t)·λ + C(t)·λ³,
//	A = √(x²/t² + ¼) + α⟨1/d⟩(t)/2,  B = ln(8/t) − 2,  C = 2π²Π·t².
//
// The vortex term is the only one linear in λ: a constant-force string.

// XGuide is the first zero of J₁′: the TE₁₁ cutoff k_c·a of a round
// channel (Part 9's XTE11 is the sphere's root, a different number).
const XGuide = 1.8411837813

// Donut is the balance at one shape t and pressure Π.
type Donut struct {
	T, Pi          float64 // a/R, p∞ħc/K²
	Lambda         float64 // R/ℓ at the size balance
	EL, EV, EP, EU float64 // energies [E0]
	M              float64 // total [E0]
	G              float64 // g = M/E_L
	LongFrac       float64 // share of the light's energy that is motion round the ring
	D2             float64 // d²e/dλ² at the balance (> 0: stable in size)
}

// DonutAt solves the size balance −A/λ² + B + 3Cλ² = 0 at fixed shape.
// invD is ⟨1/d⟩ for the core wall at aspect t (pass 0 for no charge).
func DonutAt(t, pi, invD float64) Donut {
	AL := math.Sqrt(XGuide*XGuide/(t*t) + 0.25)
	AU := Alpha * invD / 2
	A := AL + AU
	B := math.Log(8/t) - 2
	C := 2 * math.Pi * math.Pi * pi * t * t
	var l float64
	if C == 0 {
		l = math.Sqrt(A / B)
	} else {
		u := (-B + math.Sqrt(B*B+12*C*A)) / (6 * C)
		l = math.Sqrt(u)
	}
	d := Donut{T: t, Pi: pi, Lambda: l, EL: AL / l, EU: AU / l, EV: B * l, EP: C * l * l * l}
	d.M = d.EL + d.EU + d.EV + d.EP
	d.G = d.M / d.EL
	d.LongFrac = 0.25 / (XGuide*XGuide/(t*t) + 0.25)
	d.D2 = 2*A/(l*l*l) + 6*C*l
	return d
}

// RadiusOverLbar is R/ƛ_C for the donut scaled to the electron's mass:
// ℓ = ħc/E0 and E0 = m c²/M, so R/ƛ = λ·M.
func (d Donut) RadiusOverLbar() float64 { return d.Lambda * d.M }

// ShapeScan returns the balance along t at fixed Π and whether the total
// energy falls all the way to the fattest shape scanned (no thin donut is
// held).
func ShapeScan(ts []float64, pi float64, invD func(float64) float64) (ds []Donut, fattens bool) {
	fattens = true
	for i, t := range ts {
		ds = append(ds, DonutAt(t, pi, invD(t)))
		if i > 0 && ds[i].M >= ds[i-1].M {
			fattens = false
		}
	}
	return ds, fattens
}

// ---------------------------------------------------------------- kronspice twin

// GuideRing is the crossed two-wire ring line with a shunt inductor Ld
// across each section: a waveguide in circuit form, with cutoff
// ω_c = 1/√(Ld·Cd) and the crossed ring's half-integer wave on top,
// ω² = ω_c² + 4sin²(θ/2)/(2L·Cd), θ = π/N for the lowest wave.
type GuideRing struct {
	N         int
	L, Cd, Ld float64
}

// Analytic is the ladder's lowest note [Hz].
func (g GuideRing) Analytic() float64 {
	s2 := math.Pow(math.Sin(math.Pi/(2*float64(g.N))), 2)
	w2 := 1/(g.Ld*g.Cd) + 4*s2/(2*g.L*g.Cd)
	return math.Sqrt(w2) / (2 * math.Pi)
}

// Cutoff is the guide's cutoff frequency [Hz].
func (g GuideRing) Cutoff() float64 { return 1 / (2 * math.Pi * math.Sqrt(g.Ld*g.Cd)) }

// Peak runs kronspice AC round the analytic note and returns the peak of
// the differential response.
func (g GuideRing) Peak() (float64, error) {
	b := circuit.NewBuilder("crossed guide ring")
	a := func(k int) string { return "a" + strconv.Itoa(k%g.N) }
	bb := func(k int) string { return "b" + strconv.Itoa(k%g.N) }
	for k := 0; k < g.N; k++ {
		ks := strconv.Itoa(k)
		na, nb := a(k+1), bb(k+1)
		if k == g.N-1 {
			na, nb = bb(0), a(0)
		}
		b.Inductor("LA"+ks, a(k), "ma"+ks, g.L)
		b.Resistor("RA"+ks, "ma"+ks, na, 1e-3)
		b.Inductor("LB"+ks, bb(k), "mb"+ks, g.L)
		b.Resistor("RB"+ks, "mb"+ks, nb, 1e-3)
		b.Capacitor("CD"+ks, a(k), bb(k), g.Cd)
		b.Inductor("LD"+ks, a(k), bb(k), g.Ld)
	}
	b.Resistor("RBA", a(0), "0", 1e6)
	b.Resistor("RBB", bb(0), "0", 1e6)
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "b0", "a0")
	c, err := b.Build()
	if err != nil {
		return 0, err
	}
	f0 := g.Analytic()
	res, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 4001, FStart: 0.99 * f0, FStop: 1.0107 * f0}, analysis.Options{})
	if err != nil {
		return 0, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	k, best := 0, 0.0
	for i := range res.Freq {
		d := va[i] - vb[i]
		if m := real(d)*real(d) + imag(d)*imag(d); m > best {
			best, k = m, i
		}
	}
	return res.Freq[k], nil
}
