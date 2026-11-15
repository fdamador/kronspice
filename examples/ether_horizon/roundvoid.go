// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Part 9: a round void. In Parts 1 and 8 the void came out as a line,
// because swirl empties only the axis. A round void needs something to hold
// it open against the ether, and something to stop it closing. Here the
// light of the photon electron sits inside the void and pushes out
// (radiation pressure), and the nonlinear ether pushes in, with a pressure
// p∞ and a surface tension σ at the void's wall. This is the structure of
// the MIT bag model of hadrons, where quark fields inside a bubble push
// against a vacuum pressure B (Chodos, Jaffe, Johnson, Thorn & Weisskopf,
// Phys. Rev. D 9, 3471, 1974), and of Poincaré's stresses that hold a
// charged electron together (Poincaré 1906).
//
// 9a. The energy of the void of radius R holding one photon in a cavity
// mode of wave number x/R (x the mode's root):
//
//	E(R) = ħc·x/R + (4π/3)·p∞·R³ + 4π·σ·R².
//
// The light's energy falls as 1/R as the void grows (the photon number is
// conserved and its frequency scales as 1/R: Part 9b checks this in
// kronspice). dE/dR = 0 gives the radius; d²E/dR² > 0 says it is stable.

// Cavity mode roots for a spherical void with a perfectly reflecting wall:
// TM_{1,1}: tan x = x/(1 − x²)… root of [x·j1(x)]' = 0, x = 2.7437; TE_{1,1}:
// j1(x) = 0, x = 4.4934 (Jackson, Classical Electrodynamics, §8.9).
const (
	XTM11 = 2.743707269992269
	XTE11 = J1Zero
)

// VoidEnergy is E(R) [J].
func VoidEnergy(R, x, pInf, sigma float64) float64 {
	return Hbar*C0*x/R + 4*math.Pi/3*pInf*R*R*R + 4*math.Pi*sigma*R*R
}

// VoidForce is −dE/dR [N]: positive pushes the wall out.
func VoidForce(R, x, pInf, sigma float64) float64 {
	return Hbar*C0*x/(R*R) - 4*math.Pi*pInf*R*R - 8*math.Pi*sigma*R
}

// VoidStiffness is d²E/dR² [N/m]: positive is stable.
func VoidStiffness(R, x, pInf, sigma float64) float64 {
	return 2*Hbar*C0*x/(R*R*R) + 8*math.Pi*pInf*R + 8*math.Pi*sigma
}

// VoidRadius solves VoidForce = 0 by bisection [m].
func VoidRadius(x, pInf, sigma float64) float64 {
	lo, hi := 1e-30, 1.0
	for i := 0; i < 400; i++ {
		m := math.Sqrt(lo * hi)
		if VoidForce(m, x, pInf, sigma) > 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return math.Sqrt(lo * hi)
}

// PressureForElectron is the ether pressure (σ = 0) that makes the void's
// equilibrium energy m_e c²: p∞ = ħcx/(4πR⁴) at R = (4/3)·x·ƛ_C, where
// E = (4/3)·ħcx/R (the bag model's virial).
func PressureForElectron(x float64) (R, p float64) {
	lc := Hbar / (Me * C0)
	R = 4.0 / 3 * x * lc
	return R, Hbar * C0 * x / (4 * math.Pi * math.Pow(R, 4))
}

// TensionForElectron is the surface tension (p∞ = 0) that does the same:
// σ = ħcx/(8πR³) at R = (3/2)·x·ƛ_C, where E = (3/2)·ħcx/R.
func TensionForElectron(x float64) (R, sigma float64) {
	lc := Hbar / (Me * C0)
	R = 1.5 * x * lc
	return R, Hbar * C0 * x / (8 * math.Pi * R * R * R)
}

// ---------------------------------------------------------------- 9b

// BreathingTank is the void's light as a lossless LC tank whose L and C
// both scale with the void's radius, L = L0·s(t), C = C0·s(t), so its
// frequency is ω0/s — the cavity mode of a void of radius s·R0. kronspice
// runs it with a scaled inductor and a VaryingCapacitor while s grows from
// 1 to SEnd. Scaling L and C together is a pure change of time scale, so
// the stored energy goes exactly as 1/s, the photon's ħω with ω ∝ 1/R, and
// what it loses is the work done on the moving wall: radiation pressure.
type BreathingTank struct {
	L0, C0, I0     float64
	SEnd           float64 // final radius / initial radius
	Periods        float64 // ramp length in initial periods
	StepsPerPeriod float64 // time steps per initial period
}

// BreathingResult is what kronspice gives.
type BreathingResult struct {
	W0, W1         float64 // stored energy at the start and the end [J]
	Invariant      float64 // (W1·SEnd)/W0: 1 if W ∝ 1/R
	WallWork       float64 // energy given to the moving wall, ∫½ṡ(Cv² + Li²)/… dt [J]
	BooksResidual  float64 // (W0 − W1 − WallWork)/W0
	EngineMismatch float64 // kronspice's own integration mismatch [J]
	Freq0, Freq1   float64 // zero-crossing frequency at the start and the end [Hz]
}

// Run builds and solves the tank in kronspice transient. The radius grows
// smoothly, s(t) = 1 + (SEnd − 1)·(1 − e^{−(t − t1)/τ}) after t1 (SPICE's
// EXP waveform: two breakpoints only, so no backward-Euler restarts damp
// the tank), with τ = Periods initial periods.
func (bt BreathingTank) Run() (BreathingResult, error) {
	f0 := 1 / (2 * math.Pi * math.Sqrt(bt.L0*bt.C0))
	T0 := 1 / f0
	tau := bt.Periods * T0
	t1 := 20 * T0
	tStop := t1 + 14*tau + 20*T0*bt.SEnd
	s := func(t float64) float64 {
		if t <= t1 {
			return 1
		}
		return 1 + (bt.SEnd-1)*(1-math.Exp(-(t-t1)/tau))
	}
	ds := func(t float64) float64 {
		if t <= t1 {
			return 0
		}
		return (bt.SEnd - 1) / tau * math.Exp(-(t-t1)/tau)
	}
	w := circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{1, bt.SEnd, t1, tau, 2 * tStop, tau}}
	b := circuit.NewBuilder("breathing void")
	b.Add(circuit.Element{Name: "LV", Kind: circuit.NonlinearInductor, Value: bt.L0, Scale: w, IC: bt.I0, HasIC: true}, "a", "0")
	b.VaryingCapacitor("CV", "a", "0", bt.C0, w)
	c, err := b.BuildTransient()
	if err != nil {
		return BreathingResult{}, err
	}
	step := T0 / bt.StepsPerPeriod
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: step, Stop: tStop, MaxStep: step, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-7}})
	if err != nil {
		return BreathingResult{}, err
	}
	v, _ := wf.NodeVoltage("a")
	il, _ := wf.ElementCurrent("LV")
	stored := func(k int) float64 {
		sk := s(wf.Time[k])
		return 0.5*sk*bt.L0*il[k]*il[k] + 0.5*sk*bt.C0*v[k]*v[k]
	}
	wall := func(k int) float64 {
		return 0.5 * ds(wf.Time[k]) * (bt.L0*il[k]*il[k] + bt.C0*v[k]*v[k])
	}
	res := BreathingResult{EngineMismatch: wf.Ledger.IntegrationMismatch}
	avg := func(ta, tb float64) float64 {
		sum, cnt := 0.0, 0
		for k := range wf.Time {
			if wf.Time[k] >= ta && wf.Time[k] <= tb {
				sum += stored(k)
				cnt++
			}
		}
		return sum / float64(cnt)
	}
	te := wf.Time[len(wf.Time)-1]
	res.W0 = avg(2*T0, 18*T0)
	res.W1 = avg(te-18*T0*bt.SEnd, te-2*T0*bt.SEnd)
	for k := 1; k < len(wf.Time); k++ {
		res.WallWork += 0.5 * (wf.Time[k] - wf.Time[k-1]) * (wall(k) + wall(k-1))
	}
	res.Invariant = res.W1 * s(te) / res.W0
	res.BooksResidual = (res.W0 - res.W1 - res.WallWork) / res.W0
	freq := func(ta, tb float64) float64 {
		var first, last float64
		cnt := -1
		for k := 1; k < len(wf.Time); k++ {
			if wf.Time[k] < ta || wf.Time[k] > tb {
				continue
			}
			if v[k-1] < 0 && v[k] >= 0 {
				tc := wf.Time[k-1] + (wf.Time[k]-wf.Time[k-1])*(-v[k-1])/(v[k]-v[k-1])
				if cnt < 0 {
					first = tc
				}
				last = tc
				cnt++
			}
		}
		return float64(cnt) / (last - first)
	}
	res.Freq0 = freq(1*T0, 19*T0)
	res.Freq1 = freq(te-19*T0*bt.SEnd, te-T0*bt.SEnd)
	return res, nil
}

// ---------------------------------------------------------------- 9c

// The fountain round the void. Near the centre of the two-pole flow (Part
// 8a, l = 2) the ether comes in round the equator and leaves along the
// axis: an axisymmetric straining flow, potential φ = E·r²·P2(cos θ). With
// a round void of radius b at the centre the flow must go round it:
//
//	φ = E·(r² + (2/3)·b⁵/r³)·P2(cos θ)      (u_r = 0 at r = b),
//
// so the ether slides over the wall with speed u_θ = −(5/2)·E·b·sin 2θ:
// zero at the poles and the equator (stagnation), fastest at 45°.
// Bernoulli then presses hardest at the poles and the equator and least at
// 45°, and the wall gives way to that pattern. Linearized Young–Laplace for
// r = b·(1 + Σ ε_l·P_l(cos θ)): where the ether slides fastest its pressure
// is lowest and the wall bulges out, held by the tension:
// ε_l = q_l·b/(σ·(l − 1)(l + 2)) for l ≥ 2, q_l the P_l part of ½ρu²
// (Lamb, Hydrodynamics, §275).

// WallSpeed is |u_θ| on the void's wall, in units of E·b.
func WallSpeed(theta float64) float64 { return 2.5 * math.Abs(math.Sin(2*theta)) }

// legendre returns P_l(x).
func legendre(l int, x float64) float64 {
	p0, p1 := 1.0, x
	if l == 0 {
		return p0
	}
	for k := 2; k <= l; k++ {
		p0, p1 = p1, ((2*float64(k)-1)*x*p1-(float64(k)-1)*p0)/float64(k)
	}
	return p1
}

// ShapeModes returns ε_l (l = 2, 4) of the void's wall for Weber number
// We = ρ·(E·b)²·b/σ: the ½ρu² pattern projected on P_l, divided by the
// tension's stiffness. ε_0 (the size) is set by 9a and is left out; odd l
// vanish by the up–down symmetry.
func ShapeModes(we float64) map[int]float64 {
	out := map[int]float64{}
	for _, l := range []int{2, 4, 6} {
		n := 4000
		s := 0.0
		for i := 0; i < n; i++ {
			th := (float64(i) + 0.5) * math.Pi / float64(n)
			q := 0.5 * WallSpeed(th) * WallSpeed(th) // ½ρu²/(ρE²b²)
			s += q * legendre(l, math.Cos(th)) * math.Sin(th) * math.Pi / float64(n)
		}
		pl := (2*float64(l) + 1) / 2 * s // Legendre coefficient of the dynamic pressure
		// Bernoulli: wall pressure = H − ½ρu², so its P_l part is −pl.
		out[l] = -(-pl) * we / ((float64(l) - 1) * (float64(l) + 2))
	}
	return out
}

// WallRadius is the deformed wall r(θ)/b for Weber number we.
func WallRadius(theta, we float64) float64 {
	r := 1.0
	for l, e := range ShapeModes(we) {
		r += e * legendre(l, math.Cos(theta))
	}
	return r
}

// HollowTail is the radius of the hollow core the axis vortex opens in an
// ether of pressure pInf and density rho: the swirl's Bernoulli drop
// ½ρ(C0/s)² reaches p∞ at s = C0/√(2p∞/ρ) (Rankine's hollow vortex; Lamb §160).
func HollowTail(c0, pInf, rho float64) float64 { return c0 / math.Sqrt(2*pInf/rho) }

// Elongation is r(pole)/r(equator) of the wall at Weber number we.
func Elongation(we float64) float64 { return WallRadius(0, we) / WallRadius(math.Pi/2, we) }

// WeberFor finds the Weber number at which the elongation reaches e.
func WeberFor(e float64) float64 {
	lo, hi := 0.0, 1.0 // the linear theory holds only for small We
	for i := 0; i < 80; i++ {
		m := 0.5 * (lo + hi)
		if Elongation(m) < e {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi)
}
