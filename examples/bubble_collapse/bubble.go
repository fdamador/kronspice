// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
)

// The bubble as an equivalent circuit on the acoustic node of the flask
// (CLAUDE.md §1 non-goals: the bubble enters as a cited equivalent-circuit
// model; doc/domains.md for the acoustic domain: pressure across, volume
// velocity through).
//
// Across the bubble wall the circuit is, from the far liquid inwards:
//
//	p∞(t) ── [ M(V) ∥ R_rad(V) ] ── [ R_visc(V) ] ── wall ── gas C(V), surface 2σ/R
//
// with V = 4πR³/3 the bubble volume (the "charge" of the acoustic node),
// U = dV/dt the volume velocity through the branch, and
//
//	M(V)     = ρ/(4πR)          the liquid's inertance (kinetic energy 2πρR³Ṙ² = ½MU²)
//	R_rad(V) = ρc/(4πR²)        the radiation resistance of a pulsating sphere
//	R_visc   = μ/(πR³)          Newtonian viscosity at the wall, 4μṘ/R = R_visc·U
//	p_g(V)   = p_g0·((V0−Vh)/(V−Vh))^γ   argon, adiabatic, van der Waals hard core Vh
//
// The inertance in parallel with the radiation resistance is the exact
// radiation impedance of a pulsating sphere of fixed radius,
// Z = (ρc/4πR²)·jkR/(1 + jkR) (Kinsler, Frey, Coppens & Sanders,
// Fundamentals of Acoustics, 4th ed., §7.2, the pulsating sphere;
// the same as 1/Z = 1/(jωM) + 1/R_rad). Here R moves, so M and R_rad are
// configuration-dependent elements, and the inertance's equation is Kron's
// for an inductance that depends on a coordinate (Kron 1951, Epilogue
// eqs. 3–4; doc/M9c1.md): with momentum Π = M·U_M,
//
//	dΠ/dt = Δp + ½·(dM/dV)·U·U_M
//
// the ½M′ term being the Christoffel (speed-voltage) term. With R_rad → ∞
// and R_visc = 0 (incompressible, inviscid) this is exactly the
// Rayleigh–Plesset equation ρ(RR̈ + 3Ṙ²/2) = p_wall − p∞ (Lord Rayleigh,
// Phil. Mag. 34, 94, 1917; M. S. Plesset, J. Appl. Mech. 16, 277, 1949);
// R_rad brings in the liquid's compressibility to first order in Ṙ/c, the
// role the Keller–Miksis equation plays (J. B. Keller, M. Miksis, J. Acoust.
// Soc. Am. 68, 628, 1980), which is run beside it as the twin.
//
// The ledger closes by construction (derivation in README.md):
//
//	∫ −p_ac·U dt = Δ[½MU_M² + U_gas + 4πσR² + P0·V] + ∫R_visc·U² dt + ∫Δp²/R_rad dt
//
// in = Δ(stored) + viscous heat + sound radiated out, with
// p_ac = p∞ − P0 the flask mode's pressure at the centre and
// U_gas = p_g·(V − Vh)/(γ − 1).

// Bubble declares one bubble and its liquid.
type Bubble struct {
	R0    float64 // ambient radius [m]
	Pa    float64 // acoustic amplitude at the bubble [Pa]
	F     float64 // drive frequency [Hz]
	P0    float64 // static pressure [Pa]
	T0    float64 // liquid and gas temperature at rest [K]
	Rho   float64 // liquid density [kg/m³]
	C     float64 // liquid sound speed [m/s]
	Mu    float64 // liquid viscosity [Pa·s]
	Sigma float64 // surface tension [N/m]
	Gamma float64 // gas adiabatic exponent
	Core  float64 // R0/h, the van der Waals hard-core ratio
}

// Water's viscosity and surface tension at 300 K.
//
//	μ: IAPWS R12-08 (Huber et al., J. Phys. Chem. Ref. Data 38, 101, 2009),
//	   853.8 µPa·s at 300 K, 0.1 MPa.
//	σ: IAPWS R1-76(2014), σ = 235.8·τ^1.256·(1 − 0.625τ) mN/m,
//	   τ = 1 − T/647.096 K.
const muWater300 = 853.8e-6

func sigmaWater(T float64) float64 {
	tau := 1 - T/647.096
	return 235.8e-3 * math.Pow(tau, 1.256) * (1 - 0.625*tau)
}

// DefaultBubble is the argon bubble of Brenner, Hilgenfeldt & Lohse, Rev.
// Mod. Phys. 74, 425 (2002), §III (R0 = 4.5 µm, Pa = 1.3 atm), in water at
// 300 K, driven at the acoustic_drive flask's breathing mode, 25.0 kHz.
// R0/h = 8.86 is argon's van der Waals excluded volume (same §).
func DefaultBubble(f float64) Bubble {
	return Bubble{R0: 4.5e-6, Pa: 1.3 * 101325, F: f, P0: 101325, T0: 300,
		Rho: 996.556, C: 1501.51914, Mu: muWater300, Sigma: sigmaWater(300), Gamma: 5.0 / 3, Core: 8.86}
}

func (b Bubble) h() float64   { return b.R0 / b.Core }
func (b Bubble) pg0() float64 { return b.P0 + 2*b.Sigma/b.R0 }

// gas returns the gas pressure and internal energy at radius R.
func (b Bubble) gas(R float64) (p, u float64) {
	h := b.h()
	v := 4 * math.Pi / 3 * (R*R*R - h*h*h)
	v0 := 4 * math.Pi / 3 * (b.R0*b.R0*b.R0 - h*h*h)
	p = b.pg0() * math.Pow(v0/v, b.Gamma)
	return p, p * v / (b.Gamma - 1)
}

// Temperature is the adiabatic gas temperature at radius R.
func (b Bubble) Temperature(R float64) float64 {
	h := b.h()
	return b.T0 * math.Pow((b.R0*b.R0*b.R0-h*h*h)/(R*R*R-h*h*h), b.Gamma-1)
}

// Atoms is the number of argon atoms (p(V − Vh) = NkT at rest).
func (b Bubble) Atoms() float64 {
	h := b.h()
	return b.pg0() * 4 * math.Pi / 3 * (b.R0*b.R0*b.R0 - h*h*h) / (kB * b.T0)
}

const kB = 1.380649e-23

// Losses switches the two dissipative elements (for the reversibility run).
type Losses struct{ Viscous, Radiation bool }

// state: R, Π, then the integrals ∫in, ∫visc, ∫rad, ∫U cos ωt, ∫U sin ωt.
type state [7]float64

// circuitDeriv is the equivalent circuit's state equation; drive(t) is the
// acoustic pressure p_ac(t) at the centre.
func (b Bubble) circuitDeriv(t float64, y state, drive func(float64) float64, l Losses) (dy state, u float64) {
	R, Pi := y[0], y[1]
	M := b.Rho / (4 * math.Pi * R)
	dM := -b.Rho / (16 * math.Pi * math.Pi * R * R * R * R) // dM/dV
	um := Pi / M
	pg, _ := b.gas(R)
	pac := drive(t)
	a0 := pg - 2*b.Sigma/R - (b.P0 + pac)
	kv := 0.0
	if l.Viscous {
		kv = b.Mu / (math.Pi * R * R * R)
	}
	g := 0.0 // 1/R_rad
	if l.Radiation {
		g = 4 * math.Pi * R * R / (b.Rho * b.C)
	}
	u = (um + a0*g) / (1 + kv*g)
	dp := a0 - kv*u
	w := 2 * math.Pi * b.F
	dy[0] = u / (4 * math.Pi * R * R)
	dy[1] = dp + 0.5*dM*u*um
	dy[2] = -pac * u
	dy[3] = kv * u * u
	dy[4] = dp * dp * g
	dy[5] = u * math.Cos(w*t)
	dy[6] = u * math.Sin(w*t)
	return dy, u
}

// Stored is ½MU_M² + U_gas + 4πσR² + P0·V for a state.
func (b Bubble) Stored(y state) float64 {
	R, Pi := y[0], y[1]
	M := b.Rho / (4 * math.Pi * R)
	_, ug := b.gas(R)
	return 0.5*Pi*Pi/M + ug + 4*math.Pi*b.Sigma*R*R + b.P0*4*math.Pi/3*R*R*R
}

// kmDeriv is the Keller–Miksis twin in R, Ṙ, in the form of Brenner,
// Hilgenfeldt & Lohse 2002 (with the viscous part of dp_B/dt kept, its R̈
// moved to the left).
func (b Bubble) kmDeriv(t float64, y state, drive func(float64) float64) (dy state) {
	R, Rd := y[0], y[1]
	pg, _ := b.gas(R)
	h := b.h()
	dpg := -b.Gamma * pg * 3 * R * R * Rd / (R*R*R - h*h*h)
	pB := pg - 2*b.Sigma/R - 4*b.Mu*Rd/R
	dpB := dpg + 2*b.Sigma*Rd/(R*R) + 4*b.Mu*Rd*Rd/(R*R) // − 4μR̈/R moved left
	pinf := b.P0 + drive(t)
	c := b.C
	rhs := (1+Rd/c)*(pB-pinf)/b.Rho + R/(b.Rho*c)*dpB - 1.5*(1-Rd/(3*c))*Rd*Rd
	dy[0] = Rd
	dy[1] = rhs / ((1-Rd/c)*R + 4*b.Mu/(b.Rho*c))
	return dy
}

// ErrStep reports a collapse the integrator could not resolve.
var ErrStep = errors.New("bubble_collapse: step size underflow")

// dopri integrates y' = f(t, y) from t0 to t1 with the Dormand–Prince 5(4)
// pair (Hairer, Nørsett & Wanner, Solving ODEs I, §II.5), error measured on
// the first two components only (the rest are integrals carried along),
// calling obs after every accepted step.
func dopri(f func(float64, state) state, y state, t0, t1, rtol float64, scale [2]float64, obs func(t float64, y state)) (state, error) {
	const (
		a21                     = 1.0 / 5
		a31, a32                = 3.0 / 40, 9.0 / 40
		a41, a42, a43           = 44.0 / 45, -56.0 / 15, 32.0 / 9
		a51, a52, a53, a54      = 19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729
		a61, a62, a63, a64, a65 = 9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656
		b1, b3, b4, b5, b6      = 35.0 / 384, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84
		e1, e3, e4, e5, e6, e7  = 71.0 / 57600, -71.0 / 16695, 71.0 / 1920, -17253.0 / 339200, 22.0 / 525, -1.0 / 40
		c2, c3, c4, c5          = 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9
	)
	dir := 1.0
	if t1 < t0 {
		dir = -1
	}
	add := func(y state, h float64, ks ...struct {
		c float64
		k state
	}) state {
		for _, kc := range ks {
			for i := range y {
				y[i] += h * kc.c * kc.k[i]
			}
		}
		return y
	}
	type kc = struct {
		c float64
		k state
	}
	t := t0
	h := dir * 1e-12
	k1 := f(t, y)
	for dir*(t1-t) > 0 {
		if dir*(t+h-t1) > 0 {
			h = t1 - t
		}
		k2 := f(t+c2*h, add(y, h, kc{a21, k1}))
		k3 := f(t+c3*h, add(y, h, kc{a31, k1}, kc{a32, k2}))
		k4 := f(t+c4*h, add(y, h, kc{a41, k1}, kc{a42, k2}, kc{a43, k3}))
		k5 := f(t+c5*h, add(y, h, kc{a51, k1}, kc{a52, k2}, kc{a53, k3}, kc{a54, k4}))
		k6 := f(t+h, add(y, h, kc{a61, k1}, kc{a62, k2}, kc{a63, k3}, kc{a64, k4}, kc{a65, k5}))
		yn := add(y, h, kc{b1, k1}, kc{b3, k3}, kc{b4, k4}, kc{b5, k5}, kc{b6, k6})
		k7 := f(t+h, yn)
		err := 0.0
		for i := 0; i < 2; i++ {
			e := h * (e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k7[i])
			sc := rtol * (math.Abs(yn[i]) + scale[i])
			err = math.Max(err, math.Abs(e)/sc)
		}
		if err <= 1 {
			t += h
			y = yn
			k1 = k7
			if obs != nil {
				obs(t, y)
			}
		}
		fac := 0.9 * math.Pow(math.Max(err, 1e-10), -0.2)
		fac = math.Min(5, math.Max(0.2, fac))
		if math.IsNaN(err) {
			fac = 0.1
		}
		h *= fac
		if math.Abs(h) < 1e-20 {
			return y, ErrStep
		}
	}
	return y, nil
}

// Cycle is the bubble over one acoustic period.
type Cycle struct {
	Rmax, Rmin   float64    // [m]
	TRmin        float64    // time of the minimum within the run [s]
	MaxMach      float64    // max |Ṙ|/c
	Tmax, Pmax   float64    // peak adiabatic gas temperature [K] and pressure [Pa]
	In           float64    // work done by the sound field on the bubble [J]
	Viscous      float64    // viscous heat [J]
	Radiated     float64    // sound radiated by the bubble [J]
	DStored      float64    // change of stored energy [J]
	Residual     float64    // In − DStored − Viscous − Radiated [J]
	BurstFrac    float64    // fraction of Radiated emitted within ±BurstWin of the minimum
	Y1           complex128 // fundamental admittance −U₁/P₁ the bubble presents at the node [m³/(Pa·s)]
	CollapseWork float64    // P0·Vmax, the energy the collapse releases [J]
}

// BurstWin is the window around the minimum for BurstFrac.
const BurstWin = 10e-9

// drive returns the centre's acoustic pressure p_ac(t) = −Pa·sin ωt.
func (b Bubble) drive() func(float64) float64 {
	w := 2 * math.Pi * b.F
	return func(t float64) float64 { return -b.Pa * math.Sin(w*t) }
}

// Run integrates the equivalent circuit for n cycles from rest and returns
// the last one.
func (b Bubble) Run(n int) (Cycle, []Sample, error) {
	T := 1 / b.F
	dr := b.drive()
	l := Losses{true, true}
	f := func(t float64, y state) state { d, _ := b.circuitDeriv(t, y, dr, l); return d }
	y := state{b.R0, 0}
	scale := [2]float64{b.R0 * 1e-3, b.Rho / (4 * math.Pi * b.R0) * 4 * math.Pi * b.R0 * b.R0 * 1e-2}
	var err error
	y, err = dopri(f, y, 0, float64(n-1)*T, 1e-10, scale, nil)
	if err != nil {
		return Cycle{}, nil, err
	}
	start := y
	t0 := float64(n-1) * T
	cy := Cycle{Rmin: math.Inf(1)}
	var samples []Sample
	type radAt struct{ t, w float64 }
	var rad []radAt
	y, err = dopri(f, y, t0, t0+T, 1e-10, scale, func(t float64, y state) {
		d, u := b.circuitDeriv(t, y, dr, l)
		R := y[0]
		if R > cy.Rmax {
			cy.Rmax = R
		}
		if R < cy.Rmin {
			cy.Rmin, cy.TRmin = R, t
		}
		cy.MaxMach = math.Max(cy.MaxMach, math.Abs(d[0])/b.C)
		p, _ := b.gas(R)
		cy.Pmax = math.Max(cy.Pmax, p)
		samples = append(samples, Sample{T: t - t0, R: R, U: u, Temp: b.Temperature(R),
			In: y[2] - start[2], Visc: y[3] - start[3], Rad: y[4] - start[4], Stored: b.Stored(y) - b.Stored(start), Prad: d[4]})
		rad = append(rad, radAt{t, y[4]})
	})
	if err != nil {
		return Cycle{}, nil, err
	}
	cy.Tmax = b.Temperature(cy.Rmin)
	cy.In = y[2] - start[2]
	cy.Viscous = y[3] - start[3]
	cy.Radiated = y[4] - start[4]
	cy.DStored = b.Stored(y) - b.Stored(start)
	cy.Residual = cy.In - cy.DStored - cy.Viscous - cy.Radiated
	// radiated energy within the burst window
	var wa, wb float64
	wa, wb = start[4], y[4]
	for _, r := range rad {
		if r.t <= cy.TRmin-BurstWin {
			wa = r.w
		}
		if r.t <= cy.TRmin+BurstWin {
			wb = r.w
		}
	}
	cy.BurstFrac = (wb - wa) / cy.Radiated
	// Fundamental: p_ac = Re(P1 e^{jωt}) with P1 = j·Pa; the element's
	// current into it is −U, I1 = (2/T)∫ −U e^{−jωt} dt.
	cu := y[5] - start[5]
	su := y[6] - start[6]
	I1 := -complex(2/T*cu, -2/T*su)
	cy.Y1 = I1 / complex(0, b.Pa)
	cy.CollapseWork = b.P0 * 4 * math.Pi / 3 * cy.Rmax * cy.Rmax * cy.Rmax
	return cy, samples, nil
}

// Sample is one accepted step of the last cycle.
type Sample struct {
	T, R, U, Temp float64
	// cumulative ledger since the cycle began [J] and the radiated power [W]
	In, Visc, Rad, Stored, Prad float64
}

// RunKM integrates the Keller–Miksis twin and returns Rmax, Rmin and the
// time of the minimum within the last of n cycles.
func (b Bubble) RunKM(n int) (rmax, rmin, tmin, mach float64, err error) {
	return b.runKM(n, nil)
}

// runKM is RunKM with an observer of (time in the last cycle, R).
func (b Bubble) runKM(n int, obs func(t, R float64)) (rmax, rmin, tmin, mach float64, err error) {
	T := 1 / b.F
	dr := b.drive()
	f := func(t float64, y state) state { return b.kmDeriv(t, y, dr) }
	y := state{b.R0, 0}
	scale := [2]float64{b.R0 * 1e-3, 1e-2}
	y, err = dopri(f, y, 0, float64(n-1)*T, 1e-10, scale, nil)
	if err != nil {
		return
	}
	rmin = math.Inf(1)
	t0 := float64(n-1) * T
	_, err = dopri(f, y, t0, t0+T, 1e-10, scale, func(t float64, y state) {
		rmax = math.Max(rmax, y[0])
		if y[0] < rmin {
			rmin, tmin = y[0], t
		}
		mach = math.Max(mach, math.Abs(y[1])/b.C)
		if obs != nil {
			obs(t-t0, y[0])
		}
	})
	return
}

// Reversal runs the collapse forward from the bubble's maximum through its
// minimum to the next maximum-ish point (a span of 2·half), flips the
// momentum, runs forward again under the time-reversed drive, and returns
// how far from the starting radius and momentum it comes back, relative to
// R and to the peak momentum.
//
// trace, if not nil, receives every accepted step: pass 0 forward, pass 1
// after the flip, s the time since tStart in each pass.
func (b Bubble) Reversal(l Losses, tStart, half float64, trace func(pass int, s, R float64)) (dR, dPi float64, err error) {
	if trace == nil {
		trace = func(int, float64, float64) {}
	}
	dr := b.drive()
	f := func(t float64, y state) state { d, _ := b.circuitDeriv(t, y, dr, l); return d }
	// bring the bubble to tStart with the full model, then switch
	full := func(t float64, y state) state { d, _ := b.circuitDeriv(t, y, dr, Losses{true, true}); return d }
	scale := [2]float64{b.R0 * 1e-3, b.Rho / (4 * math.Pi * b.R0) * 4 * math.Pi * b.R0 * b.R0 * 1e-2}
	y0, err := dopri(full, state{b.R0, 0}, 0, tStart, 1e-11, scale, nil)
	if err != nil {
		return
	}
	y0[2], y0[3], y0[4], y0[5], y0[6] = 0, 0, 0, 0, 0
	piMax := 0.0
	t1 := tStart + 2*half
	y1, err := dopri(f, y0, tStart, t1, 1e-12, scale, func(t float64, y state) { piMax = math.Max(piMax, math.Abs(y[1])); trace(0, t-tStart, y[0]) })
	if err != nil {
		return
	}
	y1[1] = -y1[1]
	rev := func(s float64) float64 { return dr(t1 + tStart - s) }
	g := func(s float64, y state) state { d, _ := b.circuitDeriv(s, y, rev, l); return d }
	y2, err := dopri(g, y1, tStart, t1, 1e-12, scale, func(t float64, y state) { trace(1, t-tStart, y[0]) })
	if err != nil {
		return
	}
	return math.Abs(y2[0]-y0[0]) / y0[0], math.Abs(-y2[1]-y0[1]) / piMax, nil
}
