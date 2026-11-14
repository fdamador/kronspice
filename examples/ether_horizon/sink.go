// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"math"
)

// Part 5: the sink. A horizon in the ether needs an inflow, and an inflow
// needs somewhere for the ether to go. Here the sink is a localized loss
// in the Gross–Pitaevskii ether, −iγ(r)ψ with γ = γ0 inside r < Rs — the
// standard model of a BEC with a localized dissipative defect (an electron
// beam that removes atoms: Barontini et al., PRL 110, 035302, 2013). The
// steady state is solved in Madelung form, units ħ = m = g = 1 and bulk
// density 1, so μ = 1, the bulk wave speed is 1 and ξ = 1/√2:
//
//	a'' + a'/r = 2a(½v² + n²/2r² + a² − μ)      (Bernoulli, with quantum pressure)
//	F' = −2γ r a²,   v = F/(r a²)              (continuity, F = r·ρ·v_r)
//
// a = √ρ, v the radial velocity (negative: inflow), n the winding. The
// solution regular on the axis (a ≈ c·rⁿ, F ≈ 0) is shot on c so that
// a → 1 far out. The local wave speed is a, so a horizon is where the
// inflow Mach number |v|/a passes 1.

// ErrNoSink is returned when the shooting cannot bracket the sink solution.
var ErrNoSink = errors.New("ether_horizon: sink solution not bracketed")

// Sink is a steady inflow into a localized loss.
type Sink struct {
	N        int
	Gamma0   float64 // loss rate inside the sink
	Rs       float64 // sink radius (units of 1/√(gρm)/ħ; ξ = 1/√2)
	R        []float64
	A, F, V  []float64 // amplitude √ρ, flux r·ρ·v, radial velocity
	FluxFar  float64   // 2π·F at the outer edge: ether particles per second per unit length flowing in
	LossRate float64   // ∫2γρ dA: ether removed per second per unit length
}

func (s *Sink) gamma(r float64) float64 {
	if r < s.Rs {
		return s.Gamma0
	}
	return 0
}

// integrate runs RK4 from r0 to rmax for a trial axis coefficient c and
// reports +1 if a overshoots 1 (c too large), −1 if a turns back towards 0
// (c too small), 0 if it reached rmax inside the band.
func (s *Sink) integrate(c, rmax, h float64, keep bool) int {
	n := float64(s.N)
	r := h
	step := 1
	a := c * math.Pow(r, n)
	ap := n * c * math.Pow(r, n-1)
	if s.N == 0 {
		ap = 0
	}
	F := -2 * s.gamma(r) * c * c * math.Pow(r, 2*n+2) / (2*n + 2)
	gs := 0.0 // the loss rate of the current step, taken at its midpoint
	deriv := func(r, a, ap, F float64) (float64, float64, float64) {
		v := 0.0
		if a > 0 {
			v = F / (r * a * a)
		}
		app := -ap/r + 2*a*(0.5*v*v+n*n/(2*r*r)+a*a-1)
		return ap, app, -2 * gs * r * a * a
	}
	if keep {
		s.R, s.A, s.F, s.V = s.R[:0], s.A[:0], s.F[:0], s.V[:0]
	}
	for r < rmax {
		if keep {
			s.R = append(s.R, r)
			s.A = append(s.A, a)
			s.F = append(s.F, F)
			s.V = append(s.V, F/(r*a*a))
		}
		gs = s.gamma(r + h/2) // whole steps inside or outside: the edge is a grid point
		k1a, k1p, k1f := deriv(r, a, ap, F)
		k2a, k2p, k2f := deriv(r+h/2, a+h/2*k1a, ap+h/2*k1p, F+h/2*k1f)
		k3a, k3p, k3f := deriv(r+h/2, a+h/2*k2a, ap+h/2*k2p, F+h/2*k2f)
		k4a, k4p, k4f := deriv(r+h, a+h*k3a, ap+h*k3p, F+h*k3f)
		a += h / 6 * (k1a + 2*k2a + 2*k3a + k4a)
		ap += h / 6 * (k1p + 2*k2p + 2*k3p + k4p)
		F += h / 6 * (k1f + 2*k2f + 2*k3f + k4f)
		step++
		r = float64(step) * h // no drift: the sink edge falls on the grid
		if a > 1.2 || math.IsNaN(a) {
			return +1
		}
		if (r > s.Rs+2 && ap < 0 && a < 0.95) || a <= 0 {
			return -1
		}
	}
	return 0
}

// SolveSink shoots the steady inflow out to rmax with step h.
func SolveSink(n int, gamma0, rs, rmax, h float64) (*Sink, error) {
	s := &Sink{N: n, Gamma0: gamma0, Rs: rs}
	h = rs / math.Round(rs/h) // the sink edge on a grid point
	lo, hi := 1e-40, 10.0
	if s.integrate(lo, rmax, h, false) != -1 || s.integrate(hi, rmax, h, false) != +1 {
		return nil, ErrNoSink
	}
	// Bisect to the separatrix; the last run that went furthest is kept.
	for i := 0; i < 200 && hi-lo > 1e-15*hi; i++ {
		m := 0.5 * (lo + hi)
		switch s.integrate(m, rmax, h, false) {
		case +1:
			hi = m
		case -1:
			lo = m
		default:
			lo, hi = m, m
		}
	}
	s.integrate(0.5*(lo+hi), rmax, h, true)
	// Keep the part before the shot peels off the separatrix: up to the
	// first point past the sink where a is within 2e-3 of its far value,
	// then a little further.
	end := len(s.R) - 1
	for i := range s.R {
		if s.R[i] > s.Rs+3 && math.Abs(s.A[i]*s.A[i]+0.5*s.V[i]*s.V[i]+float64(n*n)/(2*s.R[i]*s.R[i])-1) < 2e-3 {
			end = i
			break
		}
	}
	s.R, s.A, s.F, s.V = s.R[:end+1], s.A[:end+1], s.F[:end+1], s.V[:end+1]
	s.FluxFar = 2 * math.Pi * s.F[end]
	// Loss by the trapezoid rule on r·a², interval by interval inside the
	// sink (Rs is a grid point), plus the piece from 0 to the first point.
	for i := 0; i+1 < len(s.R); i++ {
		rm := 0.5 * (s.R[i] + s.R[i+1])
		if s.gamma(rm) == 0 {
			continue
		}
		g0 := s.R[i] * s.A[i] * s.A[i]
		g1 := s.R[i+1] * s.A[i+1] * s.A[i+1]
		s.LossRate += 2 * s.gamma(rm) * 2 * math.Pi * 0.5 * (g0 + g1) * (s.R[i+1] - s.R[i])
	}
	s.LossRate += -2 * math.Pi * s.F[0] // F(h) is the analytic loss inside r < h
	return s, nil
}

// MaxMach is the largest inflow Mach number |v|/a and its radius.
func (s *Sink) MaxMach() (m, r float64) {
	for i := range s.R {
		if s.A[i] <= 0 {
			continue
		}
		if x := math.Abs(s.V[i]) / s.A[i]; x > m {
			m, r = x, s.R[i]
		}
	}
	return
}

// Horizons returns the radii where the inflow Mach number crosses 1.
func (s *Sink) Horizons() []float64 {
	var out []float64
	for i := 1; i < len(s.R); i++ {
		m0, m1 := math.Abs(s.V[i-1])/s.A[i-1], math.Abs(s.V[i])/s.A[i]
		if (m0-1)*(m1-1) < 0 {
			out = append(out, s.R[i-1]+(s.R[i]-s.R[i-1])*(m0-1)/(m0-m1))
		}
	}
	return out
}

// CentreDensity is ρ at the inner end of the grid.
func (s *Sink) CentreDensity() float64 { return s.A[0] * s.A[0] }

// EnergyColumn is the energy the sink takes per second per unit length:
// each ether particle arriving from far away carries the Bernoulli
// constant μ = 1, so the column is μ·(loss rate).
func (s *Sink) EnergyColumn() float64 { return 1 * s.LossRate }

// MachAt is the inflow Mach number |v|/a at radius r (nearest grid point).
func (s *Sink) MachAt(r float64) float64 {
	i := int(math.Round((r - s.R[0]) / (s.R[1] - s.R[0])))
	if i < 0 {
		i = 0
	}
	if i >= len(s.R) {
		i = len(s.R) - 1
	}
	return math.Abs(s.V[i]) / s.A[i]
}

// OuterHorizon is the outermost Mach-1 crossing, 0 if there is none.
func (s *Sink) OuterHorizon() float64 {
	h := s.Horizons()
	if len(h) == 0 {
		return 0
	}
	return h[len(h)-1]
}

// CriticalGamma is the loss rate at which the inflow first reaches the
// local wave speed at the sink's edge (bisection in log γ); above it the
// horizon stands outside the sink, in loss-free ether.
func CriticalGamma(n int, rs, rmax, h float64) (float64, *Sink, error) {
	lo, hi := math.Log(1e-3), math.Log(100.0)
	var last *Sink
	for i := 0; i < 40; i++ {
		m := 0.5 * (lo + hi)
		s, err := SolveSink(n, math.Exp(m), rs, rmax, h)
		if err != nil {
			return 0, nil, err
		}
		if s.MachAt(rs) >= 1 {
			hi = m
			last = s
		} else {
			lo = m
		}
	}
	return math.Exp(hi), last, nil
}
