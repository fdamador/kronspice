// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

const mu0 = 1.25663706212e-6

// stepsPerCycle of the fastest tone sets the time step (with the mechanics).
var stepsPerCycle = 60.0

// Material is a magnetostrictive core material in the anhysteretic,
// quadratic-magnetostriction model (declared values; see README):
//
//	M = Ms·tanh(Ht/a),  λ(Ht) = λm·tanh²(Ht/a),  Ht = H + Hbias,
//	σ = E·(ε − λ),      B = μ0·(Ht + M) + σ·λ′(Ht).
//
// Both relations come from one potential Φ(Ht, ε) = ½E(ε−λ)² −
// μ0·(Ht²/2 + Ms·a·ln cosh(Ht/a)) (σ = ∂Φ/∂ε, B = −∂Φ/∂Ht), so the
// magnetic and mechanical sides exchange energy exactly; the stored energy
// density is Φ + Ht·B.
type Material struct {
	Name           string
	Ms, A, LambdaM float64 // [A/m], [A/m], [–]
	E, Rho         float64 // [Pa], [kg/m³]
	Q              float64 // mechanical Q (Kelvin–Voigt damping set at the drive tones)
}

// Declared materials (sources in the README).
var (
	Terfenol = Material{Name: "Terfenol-D", Ms: 8.0e5, A: 8.7e4, LambdaM: 1500e-6, E: 30e9, Rho: 9250, Q: 10}
	Galfenol = Material{Name: "Galfenol (Fe81Ga19)", Ms: 1.19e6, A: 9.2e3, LambdaM: 250e-6, E: 60e9, Rho: 7800, Q: 100}
	// Terfenol-D/epoxy 1-3 composite, 51 % Terfenol: Ms and ρ by volume.
	TerfenolComposite = Material{Name: "Terfenol-D composite", Ms: 0.51 * 8.0e5, A: 1.0e5, LambdaM: 700e-6, E: 15e9, Rho: 5200, Q: 15}
)

func logcosh(x float64) float64 {
	x = math.Abs(x)
	return x + math.Log1p(math.Exp(-2*x)) - math.Ln2
}

// at returns M, M′, λ, λ′, λ″ at Ht.
func (m Material) at(ht float64) (M, dM, l, dl, ddl float64) {
	x := ht / m.A
	t := math.Tanh(x)
	s2 := 1 - t*t
	M, dM = m.Ms*t, m.Ms*s2/m.A
	l = m.LambdaM * t * t
	dl = m.LambdaM * 2 * t * s2 / m.A
	ddl = m.LambdaM * 2 * s2 * (1 - 3*t*t) / (m.A * m.A)
	return
}

// Bias is the operating point: Hbias = 0.66·a, where λ′ is largest.
func (m Material) Bias() float64 { return 0.66 * m.A }

// Small-signal properties at the bias: d = λ′, μr = 1 + M′, k² = E·d²/(μ0·μr).
func (m Material) Small() (d, mur, k float64) {
	_, dM, _, dl, _ := m.at(m.Bias())
	d, mur = dl, 1+dM
	k = math.Sqrt(m.E * d * d / (mu0 * mur))
	return
}

// Shape is the core's form.
type Shape int

const (
	Rod   Shape = iota // straight rod, flux closed by a return yoke (declared)
	Ring               // closed ring (bonded segments, or a moulded composite)
	Helix              // helical spring of wire, coils wound along the wire
)

func (s Shape) String() string { return [...]string{"rod", "ring", "helix"}[s] }

// Core is one variant: material, shape and dimensions.
type Core struct {
	Mat   Material
	Shape Shape
	Len   float64 // magnetic path = mechanical line length [m] (rod length, ring circumference, wire length)
	Area  float64 // cross-section [m²]
	R     float64 // ring radius [m]
	Segs  int     // mechanical segments (rod, helix)
	// Helix only: coil diameter, turns, pitch (for the drawing and the spring rate).
	CoilD float64
	Turns int
	Pitch float64
}

// Wave speed at constant field.
func (c Core) Speed() float64 { return math.Sqrt(c.Mat.E / c.Mat.Rho) }

// ModeF returns the n-th mode that the coils can feel: rod and helix wire
// (free–free) f = n·c/2ℓ for odd n (even modes have no net strain), ring
// breathing c/2πR (its circumferential modes have no net strain either).
func (c Core) ModeF(n int) float64 {
	if c.Shape == Ring {
		return c.Speed() / (2 * math.Pi * c.R)
	}
	return float64(n) * c.Speed() / (2 * c.Len)
}

// Cores are the five variants compared.
func Cores() []Core {
	ringR := 0.02
	wireL := 20 * math.Pi * 0.02
	return []Core{
		{Mat: Terfenol, Shape: Rod, Len: 0.05, Area: math.Pi * 0.005 * 0.005, Segs: 40},
		{Mat: Terfenol, Shape: Ring, Len: 2 * math.Pi * ringR, Area: 25e-6, R: ringR},
		{Mat: Galfenol, Shape: Ring, Len: 2 * math.Pi * ringR, Area: 25e-6, R: ringR},
		{Mat: TerfenolComposite, Shape: Helix, Len: wireL, Area: math.Pi * 1e-3 * 1e-3, Segs: 320, CoilD: 0.02, Turns: 20, Pitch: 0.004},
		{Mat: Galfenol, Shape: Helix, Len: wireL, Area: math.Pi * 1e-3 * 1e-3, Segs: 320, CoilD: 0.02, Turns: 20, Pitch: 0.004},
	}
}

// Rig is the electrical and acoustic set-up shared by all variants.
type Rig struct {
	N1, N2   float64 // primary and secondary turns, wound uniformly along the core
	Ra, Rb   float64 // port terminations [Ω]: primary side, secondary side (the load)
	Leak     float64 // leakage inductance of each winding as a fraction of its own magnetizing inductance
	Vdrive   float64 // drive source amplitude [V]
	F0       float64 // frequency grid [Hz]
	Fe       float64 // electrical drive [Hz]
	Force    float64 // each transducer's force amplitude [N]
	Warm     float64 // [s]
	Reverse  bool    // drive the secondary (in series with Rb), watch the primary
	Acoustic bool    // transducers on
	Electric bool    // drive on
}

// DefaultRig: 200 and 400 turns, 2 Ω and 20 Ω, 5 % leakage, 20 V at 5 kHz,
// 5 N per transducer.
func DefaultRig() Rig {
	return Rig{N1: 200, N2: 400, Ra: 2, Rb: 20, Leak: 0.05, Vdrive: 20, F0: 100, Fe: 5000, Force: 5,
		Warm: 10e-3, Acoustic: true, Electric: true}
}

// Tones returns the two transducer frequencies on the grid: rod modes 1 and
// 3; ring breathing and 1.37× it (the ring's one coupled mode, driven on and
// off resonance); helix wire odd modes nearest 20 and 40 kHz.
func (c Core) Tones(f0 float64) (fa, fb float64) {
	grid := func(f float64) float64 { return math.Round(f/f0) * f0 }
	switch c.Shape {
	case Rod:
		return grid(c.ModeF(1)), grid(c.ModeF(3))
	case Ring:
		return grid(c.ModeF(0)), grid(1.37 * c.ModeF(0))
	default:
		odd := func(f float64) float64 {
			n := int(math.Round(f / c.ModeF(1)))
			if n%2 == 0 {
				n++
			}
			return grid(c.ModeF(n))
		}
		return odd(20e3), odd(40e3)
	}
}

// Result is one run.
type Result struct {
	Core   Core
	Rig    Rig
	Fa, Fb float64
	// Window: output current (secondary forward, primary reverse) and its lines.
	T, Iout []float64
	Lines   map[string]complex128 // by name: fe, fa, fb, fa+fe, …
	// Ledger [J]
	In, Lost, Stored0, Stored1, Residual float64
	Hac                                  float64 // peak |H − 0| reached [A/m]
	Steps                                int
}

// Run integrates the coupled system for Warm + one window 1/F0.
func (c Core) Run(r Rig) (Result, error) {
	res := Result{Core: c, Rig: r}
	m := c.Mat
	fa, fb := c.Tones(r.F0)
	res.Fa, res.Fb = fa, fb
	hb := m.Bias()
	_, dM0, l0, dl0, _ := m.at(hb)
	g0 := mu0*(1+dM0) - m.E*dl0*dl0 // ∂B/∂H at the bias, unstressed
	lm1 := r.N1 * r.N1 * c.Area * g0 / c.Len
	lm2 := r.N2 * r.N2 * c.Area * g0 / c.Len
	ll1, ll2 := r.Leak*lm1, r.Leak*lm2
	wRef := 2 * math.Pi * fa
	eta := m.E / (wRef * m.Q) // Kelvin–Voigt viscosity giving Q at the first tone

	// Mechanical layout.
	ring := c.Shape == Ring
	ns := c.Segs
	if ring {
		ns = 1
	}
	dx := c.Len / float64(ns)
	nn := ns + 1 // nodes (rod/helix)
	if ring {
		nn = 1 // the radial coordinate w
	}
	mass := make([]float64, nn)
	if ring {
		mass[0] = m.Rho * c.Area * c.Len
	} else {
		for k := range mass {
			mass[k] = m.Rho * c.Area * dx
			if k == 0 || k == nn-1 {
				mass[k] /= 2
			}
		}
	}
	// State: currents, node positions and velocities (ring: w and ẇ).
	var i1, i2 float64
	pos := make([]float64, nn)
	vel := make([]float64, nn)
	// Start at the bias equilibrium: strain = λ(Hb), so σ = 0.
	if ring {
		pos[0] = l0 * c.R
	} else {
		for k := range pos {
			pos[k] = l0 * dx * float64(k)
		}
	}
	src := func(t float64) (v1, v2, Fa, Fb float64) {
		if r.Electric {
			v := r.Vdrive * math.Sin(2*math.Pi*r.Fe*t)
			if r.Reverse {
				v2 = v
			} else {
				v1 = v
			}
		}
		if r.Acoustic {
			Fa = r.Force * math.Sin(2*math.Pi*fa*t)
			Fb = r.Force * math.Sin(2*math.Pi*fb*t)
		}
		return
	}
	strain := func(j int) (e, ed float64) {
		if ring {
			return pos[0] / c.R, vel[0] / c.R
		}
		return (pos[j+1] - pos[j]) / dx, (vel[j+1] - vel[j]) / dx
	}
	vol := c.Area * dx
	if ring {
		vol = c.Area * c.Len
	}
	stored := func() float64 {
		ht := (r.N1*i1+r.N2*i2)/c.Len + hb
		M, _, l, dl, _ := m.at(ht)
		e := 0.5*ll1*i1*i1 + 0.5*ll2*i2*i2
		for j := 0; j < ns; j++ {
			eps, _ := strain(j)
			sig := m.E * (eps - l)
			B := mu0*(ht+M) + sig*dl
			phi := 0.5*m.E*(eps-l)*(eps-l) - mu0*(ht*ht/2+m.Ms*m.A*logcosh(ht/m.A))
			e += vol * (phi + ht*B)
		}
		for k := range vel {
			e += 0.5 * mass[k] * vel[k] * vel[k]
		}
		return e
	}
	// Step: the coil currents are stiff (a thin core's leakage time can be
	// nanoseconds), so they step by backward Euler; the mechanics by
	// symplectic Euler. The step is set by the mechanics and the fastest
	// tone. The ledger is first order in the step; its residual is reported.
	fmax := math.Max(fb, r.Fe)
	dt := 1 / (stepsPerCycle * fmax)
	if !ring {
		dt = math.Min(dt, 0.3*dx/c.Speed())
	}
	T := 1 / r.F0
	per := int(math.Ceil(T / dt))
	dt = T / float64(per)
	warm := int(math.Round(r.Warm / dt))
	res.Stored0 = stored()
	force := make([]float64, nn)
	sig := make([]float64, ns)
	t := 0.0
	for s := 0; s < warm+per; s++ {
		ht := (r.N1*i1+r.N2*i2)/c.Len + hb
		_, dM, l, dl, ddl := m.at(ht)
		_, _, Fa, Fb := src(t)
		var sAvg float64
		for j := 0; j < ns; j++ {
			e, ed := strain(j)
			sAvg += m.E * (e - l)
			sig[j] = m.E*(e-l) + eta*ed
		}
		sAvg /= float64(ns)
		// Mechanics (symplectic Euler).
		var pMech, pDamp float64
		if ring {
			F := Fa + Fb
			vel[0] += dt * (-2*math.Pi*c.Area*sig[0] + F) / mass[0]
			pos[0] += dt * vel[0]
			pMech = F * vel[0]
		} else {
			for k := range force {
				f := 0.0
				if k < ns {
					f += c.Area * sig[k]
				}
				if k > 0 {
					f -= c.Area * sig[k-1]
				}
				if k == 0 {
					f += Fa
				}
				if k == nn-1 {
					f -= Fb
				}
				force[k] = f
			}
			for k := range vel {
				vel[k] += dt * force[k] / mass[k]
				pos[k] += dt * vel[k]
			}
			pMech = Fa*vel[0] - Fb*vel[nn-1]
		}
		var edAvg float64
		for j := 0; j < ns; j++ {
			_, ed := strain(j)
			edAvg += ed
			pDamp += eta * ed * ed * vol
		}
		edAvg /= float64(ns)
		// Coils (backward Euler, G and λ′ frozen over the step).
		t += dt
		v1, v2, _, _ := src(t)
		G := mu0*(1+dM) + sAvg*ddl - m.E*dl*dl
		cc := m.E * dl * edAvg
		beta := c.Area * G / (c.Len * dt)
		a1, a2 := ll1/dt+r.Ra, ll2/dt+r.Rb
		m11, m12, m22 := a1+r.N1*r.N1*beta, r.N1*r.N2*beta, a2+r.N2*r.N2*beta
		b1 := v1 - r.Ra*i1 - r.N1*c.Area*cc
		b2 := v2 - r.Rb*i2 - r.N2*c.Area*cc
		det := m11*m22 - m12*m12
		d1 := (b1*m22 - m12*b2) / det
		d2 := (m11*b2 - m12*b1) / det
		bd := G*(r.N1*d1+r.N2*d2)/(c.Len*dt) + cc
		i1 += d1
		i2 += d2
		pBias := hb * c.Area * c.Len * bd
		res.In += dt * (v1*i1 + v2*i2 + pMech + pBias)
		res.Lost += dt * (r.Ra*i1*i1 + r.Rb*i2*i2 + pDamp)
		if math.IsNaN(i1) || math.IsInf(i1, 0) || math.IsNaN(pos[0]) {
			return res, fmt.Errorf("%s %s: integration blew up at t = %g", c.Mat.Name, c.Shape, t)
		}
		res.Hac = math.Max(res.Hac, math.Abs(r.N1*i1+r.N2*i2)/c.Len)
		if s >= warm {
			res.T = append(res.T, t)
			if r.Reverse {
				res.Iout = append(res.Iout, i1)
			} else {
				res.Iout = append(res.Iout, i2)
			}
		}
	}
	res.Steps = warm + per
	res.Stored1 = stored()
	res.Residual = res.In - res.Lost - (res.Stored1 - res.Stored0)
	res.Lines = map[string]complex128{}
	k := func(f float64) int { return int(math.Round(f / r.F0)) }
	for name, f := range map[string]float64{
		"fe": r.Fe, "2fe": 2 * r.Fe, "fa": fa, "fb": fb, "fa+fe": fa + r.Fe, "fa-fe": fa - r.Fe,
		"fb+fe": fb + r.Fe, "fb-fe": fb - r.Fe, "fa+fb": fa + fb, "fb-fa": fb - fa, "2fa": 2 * fa,
	} {
		res.Lines[name] = hann(res.Iout, k(f))
	}
	return res, nil
}

// hann is the Hann-windowed (unit gain) coefficient of e^{+j2πkt/T}.
func hann(x []float64, k int) complex128 {
	var acc complex128
	m := len(x)
	for j, v := range x {
		w := 1 - math.Cos(2*math.Pi*float64(j)/float64(m))
		acc += complex(v*w, 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(k*j)/float64(m)))
	}
	return acc / complex(float64(m), 0)
}

// Amp is a line's amplitude [A].
func (r Result) Amp(name string) float64 { return 2 * cmplx.Abs(r.Lines[name]) }
