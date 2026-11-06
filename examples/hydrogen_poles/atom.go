// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
)

// Units: lengths in Bohr radii a0, magnetic moments in Bohr magnetons μB,
// magnetization in μB/a0³, current density in μB/a0⁴ (= A/m² × a0⁴/μB),
// fields in B0 = μ0·μB/(4π·a0³).

// Physical constants (CODATA 2018).
const (
	MuB  = 9.2740100783e-24 // J/T
	A0   = 5.29177210903e-11
	Mu0  = 1.25663706212e-6
	GS   = 2.00231930436256 // electron spin g-factor (magnitude)
	MuPr = 1.52103220230e-3 // proton magnetic moment in μB (2.79284734463 μN)
)

// B0 is the field unit μ0μB/(4πa0³) in tesla.
var B0 = Mu0 * MuB / (4 * math.Pi * A0 * A0 * A0)

// State is a hydrogen level |n l j m_j⟩ in the nonrelativistic Pauli
// form: ψ = c↑ R_nl Y_{l,m−½} χ↑ + c↓ R_nl Y_{l,m+½} χ↓.
type State struct {
	Name   string
	N, L   int
	Jq, M  float64 // j and m_j
	cu, cd float64 // Clebsch–Gordan coefficients
	mu, md int     // m_l of the two components
}

// NewState builds |n l j m⟩ with the Condon–Shortley coefficients for
// coupling l and ½ (e.g. Sakurai, Modern Quantum Mechanics, eq. 3.8.64).
func NewState(name string, n, l int, j, m float64) State {
	s := State{Name: name, N: n, L: l, Jq: j, M: m}
	s.mu, s.md = int(math.Round(m-0.5)), int(math.Round(m+0.5))
	L := float64(l)
	if j > L {
		s.cu = math.Sqrt((L + m + 0.5) / (2*L + 1))
		s.cd = math.Sqrt((L - m + 0.5) / (2*L + 1))
	} else {
		s.cu = -math.Sqrt((L - m + 0.5) / (2*L + 1))
		s.cd = math.Sqrt((L + m + 0.5) / (2*L + 1))
	}
	if abs(s.mu) > l {
		s.cu = 0
	}
	if abs(s.md) > l {
		s.cd = 0
	}
	return s
}

func abs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}

// States are the levels studied: the ground state and the n = 2 p levels.
func States() []State {
	return []State{
		NewState("1s½, m = +½", 1, 0, 0.5, 0.5),
		NewState("2p½, m = +½", 2, 1, 0.5, 0.5),
		NewState("2p3⁄2, m = +3⁄2", 2, 1, 1.5, 1.5),
		NewState("2p3⁄2, m = +½", 2, 1, 1.5, 0.5),
	}
}

// Radial is R_nl(r) (a0^−3/2).
func Radial(n, l int, r float64) float64 {
	switch {
	case n == 1 && l == 0:
		return 2 * math.Exp(-r)
	case n == 2 && l == 0:
		return (1 - r/2) * math.Exp(-r/2) / math.Sqrt2
	case n == 2 && l == 1:
		return r * math.Exp(-r/2) / (2 * math.Sqrt(6))
	}
	return 0
}

// Ylm is the spherical harmonic at φ = 0 (real there), l ≤ 1.
func Ylm(l, m int, th float64) float64 {
	switch {
	case l == 0:
		return 1 / math.Sqrt(4*math.Pi)
	case l == 1 && m == 0:
		return math.Sqrt(3/(4*math.Pi)) * math.Cos(th)
	case l == 1 && m == 1:
		return -math.Sqrt(3/(8*math.Pi)) * math.Sin(th)
	case l == 1 && m == -1:
		return math.Sqrt(3/(8*math.Pi)) * math.Sin(th)
	}
	return 0
}

// Spinor returns the two components at (r, θ, φ = 0).
func (s State) Spinor(r, th float64) (up, dn float64) {
	R := Radial(s.N, s.L, r)
	return s.cu * R * Ylm(s.L, s.mu, th), s.cd * R * Ylm(s.L, s.md, th)
}

// Density is |ψ|².
func (s State) Density(r, th float64) float64 {
	u, d := s.Spinor(r, th)
	return u*u + d*d
}

// Mag is the spin magnetization M = −(g_s/2) μB ψ†σψ in spherical
// components (M_r, M_θ); M_φ = 0 since ψ↑*ψ↓ ∝ e^{iφ}.
func (s State) Mag(r, th float64) (mr, mt float64) {
	u, d := s.Spinor(r, th)
	sz, srho := u*u-d*d, 2*u*d
	mz, mrho := -GS/2*sz, -GS/2*srho
	c, sn := math.Cos(th), math.Sin(th)
	return mrho*sn + mz*c, mrho*c - mz*sn
}

// JOrb is the orbital current density J_φ = −(eħ/m) Σ m_l|ψ_c|²/ρ
// = −2μB Σ m_l|ψ_c|²/ρ.
func (s State) JOrb(r, th float64) float64 {
	u, d := s.Spinor(r, th)
	rho := r * math.Sin(th)
	return -2 * (float64(s.mu)*u*u + float64(s.md)*d*d) / rho
}

// JSpin is the magnetization current (∇×M)_φ = (1/r)[∂r(r M_θ) − ∂θ M_r],
// by central differences on the analytic magnetization.
func (s State) JSpin(r, th float64) float64 {
	return curlPhi(s.Mag, r, th)
}

func curlPhi(m func(r, th float64) (float64, float64), r, th float64) float64 {
	hr, ht := 1e-5*(1+r), 1e-5
	_, a := m(r+hr, th)
	_, b := m(r-hr, th)
	dr := ((r+hr)*a - (r-hr)*b) / (2 * hr)
	c, _ := m(r, th+ht)
	d, _ := m(r, th-ht)
	return (dr - (c-d)/(2*ht)) / r
}

// J is the total azimuthal current density, orbital plus spin.
func (s State) J(r, th float64) float64 { return s.JOrb(r, th) + s.JSpin(r, th) }

// GJ is the Landé factor with the measured spin g-factor:
// g_j = 1 + (g_s − 1)[j(j+1) − l(l+1) + ¾]/(2j(j+1)).
func (s State) GJ() float64 {
	j, l := s.Jq, float64(s.L)
	return 1 + (GS-1)*(j*(j+1)-l*(l+1)+0.75)/(2*j*(j+1))
}
