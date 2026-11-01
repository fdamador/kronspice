// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"github.com/fdamador/kronspice/hyper/ga"
	"github.com/fdamador/kronspice/hyper/quat"
)

// Maxwell's equations as the single spacetime-algebra equation ∇F = J/(ε₀c)
// (doc.go; Hestenes 1966; Doran & Lasenby 2003, §7.1), on hyper/ga's
// Cl(1,3) with e₁ … e₄ = γ₀ … γ₃.

// STA is the spacetime algebra and the elements the field equations use.
type STA struct {
	Alg   *ga.Algebra
	Gamma [4]ga.Multivector // γ₀ … γ₃
	Sigma [3]ga.Multivector // σₖ = γₖγ₀
	ISig  [3]ga.Multivector // I·σₖ
	I     ga.Multivector
}

// NewSTA builds the spacetime algebra.
func NewSTA() *STA {
	a, _ := ga.New(1, 3, 0)
	s := &STA{Alg: a}
	for k := 0; k < 4; k++ {
		s.Gamma[k], _ = a.Basis(k + 1)
	}
	s.I = a.Mul(a.Mul(s.Gamma[0], s.Gamma[1]), a.Mul(s.Gamma[2], s.Gamma[3]))
	for k := 0; k < 3; k++ {
		s.Sigma[k] = a.Mul(s.Gamma[k+1], s.Gamma[0])
		s.ISig[k] = a.Mul(s.I, s.Sigma[k])
	}
	return s
}

// Field is the bivector F = E + I·c·B.
func (s *STA) Field(e, b quat.Q) ga.Multivector {
	out := ga.Multivector{}
	ev := [3]float64{e.X, e.Y, e.Z}
	bv := [3]float64{b.X, b.Y, b.Z}
	for k := 0; k < 3; k++ {
		out = ga.Add(out, s.Sigma[k].Scale(ev[k]))
		out = ga.Add(out, s.ISig[k].Scale(SpeedOfLight*bv[k]))
	}
	return out
}

// Current is the source vector J/(ε₀c) with J = cρ·γ₀ + Σ Jₖγₖ.
func (s *STA) Current(rho float64, j quat.Q) ga.Multivector {
	out := s.Gamma[0].Scale(SpeedOfLight * rho)
	for k, v := range [3]float64{j.X, j.Y, j.Z} {
		out = ga.Add(out, s.Gamma[k+1].Scale(v))
	}
	return out.Scale(1 / (Eps0 * SpeedOfLight))
}

// MaxwellSpacetime returns R = ∇F − J/(ε₀c) at every interior point
// (zero multivector at the edges, with ok[n] false): ∇F =
// γ₀·∂F/∂(ct) − Σₖ γₖ·∂F/∂xₖ, the derivatives of F built from those of E
// and B. Its grade-1 part holds Gauss's and the Ampère–Maxwell law, its
// grade-3 part Faraday's law and div B = 0 (SpacetimeParts splits them).
func (s *STA) MaxwellSpacetime(E, B, dEdt, dBdt QField, rho []float64, J *QField) (R []ga.Multivector, ok []bool) {
	g := E.G
	dE := [3]QField{Partial(E, 0), Partial(E, 1), Partial(E, 2)}
	dB := [3]QField{Partial(B, 0), Partial(B, 1), Partial(B, 2)}
	R = make([]ga.Multivector, g.Len())
	ok = make([]bool, g.Len())
	a := s.Alg
	for n := range R {
		if !g.Interior(n) {
			continue
		}
		ok[n] = true
		nf := a.Mul(s.Gamma[0], s.Field(dEdt.Q[n], dBdt.Q[n])).Scale(1 / SpeedOfLight)
		for k := 0; k < 3; k++ {
			nf = ga.Sub(nf, a.Mul(s.Gamma[k+1], s.Field(dE[k].Q[n], dB[k].Q[n])))
		}
		q := 0.0
		if rho != nil {
			q = rho[n]
		}
		j := quat.Q{}
		if J != nil {
			j = J.Q[n]
		}
		R[n] = ga.Sub(nf, s.Current(q, j))
	}
	return R, ok
}

// SpacetimeParts reads the four vector-calculus residuals out of a
// spacetime residual R, in VectorResidual's units: Gauss (div E − ρ/ε₀),
// Faraday (curl E + ∂B/∂t), div B and Ampère–Maxwell
// (curl B − μ₀J − ∂E/∂t/c²). With x⁰ = ct and F = E + IcB,
//
//	⟨R⟩₁ = Gauss·γ₀ + c·Σₖ Ampèreₖ·γₖ
//	⟨R⟩₃ = Σₖ Faradayₖ·Iγₖ − c·(div B)·Iγ₀
//
// read off by projection onto each basis element (scalar product divided
// by its square). The tests check these factors against VectorResidual on
// arbitrary (non-Maxwellian) fields and sources, where all four residuals
// are non-zero.
func (s *STA) SpacetimeParts(R ga.Multivector) (gauss float64, faraday [3]float64, divB float64, ampere [3]float64) {
	a := s.Alg
	g1 := R.Grade(1)
	g3 := R.Grade(3)
	c := SpeedOfLight
	proj := func(x, e ga.Multivector) float64 { return a.ScalarProduct(x, e) / a.ScalarProduct(e, e) }
	gauss = proj(g1, s.Gamma[0])
	for k := 0; k < 3; k++ {
		ampere[k] = proj(g1, s.Gamma[k+1]) / c
	}
	ig0 := a.Mul(s.I, s.Gamma[0])
	divB = -proj(g3, ig0) / c
	for k := 0; k < 3; k++ {
		faraday[k] = proj(g3, a.Mul(s.I, s.Gamma[k+1]))
	}
	return gauss, faraday, divB, ampere
}
