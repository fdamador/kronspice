// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
)

// Rate-dependent materials: Maxwell's constitutive relation 𝔇 = (1/4π)K𝔈
// (Treatise Vol. II, Art. 619) with a permittivity K that depends on how
// fast the field changes. In the time domain the material carries its own
// polarization equations, which contain the first (Debye) and second
// (Lorentz) rates of change of the polarization; in the frequency domain
// they become ε(ω). Engineering sign convention e^{jωt}: ε = ε′ − jε″.
//
//	Debye term   τ dP/dt + P = ε0 Δε E                → Δε / (1 + jωτ)
//	Lorentz term d²P/dt² + γ dP/dt + ω0² P = ε0 Δε ω0² E → Δε ω0² / (ω0² − ω² + jγω)
//	D = ε0 ε∞ E + Σ P_k
//
// (P. Debye, Polar Molecules, 1929; H. A. Lorentz, The Theory of
// Electrons, 1909; time-domain forms as the auxiliary differential
// equations of Taflove & Hagness, Computational Electrodynamics, ch. 9.)
//
// Each term is also an exact circuit (Kron/Branin network analogy): a
// plate cell of area A and gap d is C∞ = ε0ε∞A/d in parallel with one
// branch per term — Debye: series R–C with C_k = ε0ΔεA/d, R_k = τ/C_k;
// Lorentz: series R–L–C with C_k = ε0ΔεA/d, L_k = 1/(ω0²C_k), R_k = γL_k.
// The cell's admittance is then exactly jω(ε0A/d)ε(ω).

// DebyeTerm is one relaxation: strength Δε and time τ [s].
type DebyeTerm struct{ DeltaEps, Tau float64 }

// LorentzTerm is one resonance: strength Δε, resonance ω0 [rad/s] and
// damping γ [1/s].
type LorentzTerm struct{ DeltaEps, Omega0, Gamma float64 }

// Dispersive is a material whose permittivity depends on the rate of
// change of the field.
type Dispersive struct {
	Name    string
	EpsInf  float64
	Debye   []DebyeTerm
	Lorentz []LorentzTerm
	Source  string
}

// WaterDebye25 is liquid water at 25 °C as a single Debye relaxation:
// ε(0) = 78.6, ε∞ = 5.4 (D. C. Elton, "The origin of the Debye relaxation
// in liquid water and fitting the high frequency excess response,"
// arXiv:1704.01667, Sec. II.B) and τ_D = 8.3 ps, the consensus of the
// dielectric measurements tabulated there (Table I: e.g. Barthel et al.,
// Chem. Phys. Lett. 165, 369, 1990, 8.3 ps; Buchner et al., Chem. Phys.
// Lett. 306, 57, 1999, 8.4 ps). Valid below about 100 GHz; faster
// relaxations (≈ 1 ps and below) and the infrared bands are not included.
var WaterDebye25 = Dispersive{
	Name:   "liquid water, 25 °C (single Debye)",
	EpsInf: 5.4,
	Debye:  []DebyeTerm{{DeltaEps: 78.6 - 5.4, Tau: 8.3e-12}},
	Source: "ε(0) 78.6, ε∞ 5.4: Elton, arXiv:1704.01667 §II.B; τ 8.3 ps: Barthel et al. 1990 / Buchner et al. 1999 via Elton Table I",
}

// Validate rejects non-physical parameters.
func (m Dispersive) Validate() error {
	if !(m.EpsInf >= 1) {
		return fmt.Errorf("%w: ε∞ = %g < 1", ErrOutsideModelValidity, m.EpsInf)
	}
	for _, d := range m.Debye {
		if !(d.DeltaEps > 0) || !(d.Tau > 0) {
			return fmt.Errorf("%w: Debye term needs Δε > 0 and τ > 0", ErrOutsideModelValidity)
		}
	}
	for _, l := range m.Lorentz {
		if !(l.DeltaEps > 0) || !(l.Omega0 > 0) || !(l.Gamma >= 0) {
			return fmt.Errorf("%w: Lorentz term needs Δε > 0, ω0 > 0, γ ≥ 0", ErrOutsideModelValidity)
		}
	}
	return nil
}

// Eps is the relative permittivity at angular frequency ω.
func (m Dispersive) Eps(omega float64) complex128 {
	e := complex(m.EpsInf, 0)
	for _, d := range m.Debye {
		e += complex(d.DeltaEps, 0) / complex(1, omega*d.Tau)
	}
	for _, l := range m.Lorentz {
		w2 := l.Omega0 * l.Omega0
		e += complex(l.DeltaEps*w2, 0) / complex(w2-omega*omega, l.Gamma*omega)
	}
	return e
}

// Static is ε(0).
func (m Dispersive) Static() float64 { return real(m.Eps(0)) }

// Response integrates the polarization equations for a field E(t) and
// returns D(t)/ε0 at the n+1 times 0, h, …, nh (the material starts
// unpolarized). See ResponseAt.
func (m Dispersive) Response(E func(t float64) float64, h float64, n int) []float64 {
	ts := make([]float64, n+1)
	for k := range ts {
		ts[k] = float64(k) * h
	}
	return m.ResponseAt(E, ts)
}

// ResponseAt is Response on any increasing list of times starting at 0.
// Debye terms use the exponential integrator for a field varying linearly
// over each step (exact for piecewise-linear E whose corners fall on the
// grid); Lorentz terms use classical RK4.
func (m Dispersive) ResponseAt(E func(t float64) float64, ts []float64) []float64 {
	out := make([]float64, len(ts))
	pd := make([]float64, len(m.Debye))
	pl := make([][2]float64, len(m.Lorentz))
	d := func(t float64) float64 {
		s := m.EpsInf * E(t)
		for _, p := range pd {
			s += p
		}
		for _, p := range pl {
			s += p[0]
		}
		return s
	}
	if len(ts) == 0 {
		return out
	}
	out[0] = d(ts[0])
	for k := 0; k+1 < len(ts); k++ {
		t0, h := ts[k], ts[k+1]-ts[k]
		if !(h > 0) {
			out[k+1] = d(ts[k+1])
			continue
		}
		e0, e1 := E(t0), E(t0+h)
		for i, dt := range m.Debye {
			// τ P′ + P = Δε E, E linear over the step: exact update.
			a := math.Exp(-h / dt.Tau)
			slope := (e1 - e0) / h
			// particular solution P_p(t) = Δε (E(t) − τ E′)
			pp0 := dt.DeltaEps * (e0 - dt.Tau*slope)
			pp1 := dt.DeltaEps * (e1 - dt.Tau*slope)
			pd[i] = pp1 + (pd[i]-pp0)*a
		}
		for i, l := range m.Lorentz {
			f := func(t float64, y [2]float64) [2]float64 {
				w2 := l.Omega0 * l.Omega0
				return [2]float64{y[1], l.DeltaEps*w2*E(t) - l.Gamma*y[1] - w2*y[0]}
			}
			y := pl[i]
			add := func(a [2]float64, s float64, b [2]float64) [2]float64 {
				return [2]float64{a[0] + s*b[0], a[1] + s*b[1]}
			}
			k1 := f(t0, y)
			k2 := f(t0+h/2, add(y, h/2, k1))
			k3 := f(t0+h/2, add(y, h/2, k2))
			k4 := f(t0+h, add(y, h, k3))
			pl[i] = [2]float64{y[0] + h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), y[1] + h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])}
		}
		out[k+1] = d(ts[k+1])
	}
	return out
}

// CellNetlist is the plate cell (area a, gap g) filled with the material,
// as a SPICE netlist driven by a 1 V AC source between "top" and ground:
// the source current is jω(ε0a/g)ε(ω).
func (m Dispersive) CellNetlist(a, g float64, acLine string) string {
	c0 := Eps0 * a / g
	var s strings.Builder
	fmt.Fprintf(&s, "plate cell filled with %s (A = %g m², d = %g m)\n", m.Name, a, g)
	s.WriteString("V1 top 0 DC 0 AC 1\n")
	fmt.Fprintf(&s, "CINF top 0 %.15g\n", c0*m.EpsInf)
	for i, d := range m.Debye {
		ck := c0 * d.DeltaEps
		fmt.Fprintf(&s, "RD%d top d%d %.15g\nCD%d d%d 0 %.15g\n", i, i, d.Tau/ck, i, i, ck)
	}
	for i, l := range m.Lorentz {
		ck := c0 * l.DeltaEps
		lk := 1 / (l.Omega0 * l.Omega0 * ck)
		fmt.Fprintf(&s, "RL%d top l%d %.15g\nLL%d l%d m%d %.15g\nCL%d m%d 0 %.15g\n", i, i, l.Gamma*lk, i, i, i, lk, i, i, ck)
	}
	s.WriteString(acLine + "\n.end\n")
	return s.String()
}

// CellAdmittance is the analytic admittance of that cell [S].
func (m Dispersive) CellAdmittance(a, g, omega float64) complex128 {
	return complex(0, omega) * complex(Eps0*a/g, 0) * m.Eps(omega)
}

// Loss tangent ε″/ε′.
func (m Dispersive) LossTangent(omega float64) float64 {
	e := m.Eps(omega)
	return -imag(e) / real(e)
}

// RefractiveIndex is √ε (non-magnetic medium), principal branch with
// non-positive imaginary part in this sign convention.
func (m Dispersive) RefractiveIndex(omega float64) complex128 { return cmplx.Sqrt(m.Eps(omega)) }
