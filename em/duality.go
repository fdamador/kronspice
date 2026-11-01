// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"math"

	"github.com/fdamador/kronspice/hyper/quat"
)

// Dual-symmetric Maxwell equations with magnetic charge (CLAUDE.md §3f,
// hypothesis H2). HYPOTHETICAL: no magnetic monopole has been observed.
// Everything here is off unless Hypothetical.MagneticCharge is set, and
// its results are labelled hypothetical; their power belongs only in the
// ledger's hypothetical column.
//
// With magnetic charge density ρ_m [A·m/m³] and magnetic current density
// J_m [A·m/(m²·s)] (SI, q_m in A·m; J. D. Jackson, Classical Electrodynamics, 3rd ed., §6.11,
// eqs. 6.150), the equations become
//
//	div B = μ₀ρ_m,   curl E = −∂B/∂t − μ₀J_m,
//
// in spacetime-algebra form ∇F = J_e + I·J_m. They are invariant under the
// duality rotation (Jackson eq. 6.151) E′ = E cos ξ + cB sin ξ,
// cB′ = −E sin ξ + cB cos ξ, with the charges rotating the same way:
// q′ = q cos ξ + (q_m/c) sin ξ, q_m′/c = −q sin ξ + (q_m/c) cos ξ.

// ErrHypotheticalDisabled is returned when magnetic charge appears without
// the hypothetical mode switched on.
var ErrHypotheticalDisabled = errors.New("em: magnetic charge requires Hypothetical.MagneticCharge (hypothetical mode)")

// Hypothetical switches the hypothetical models on (CLAUDE.md §3f: off by
// default, labelled in all output).
type Hypothetical struct {
	MagneticCharge bool
}

// DualResidual is a QuaternionResidual computed with magnetic sources;
// Hypothetical is always true so a report cannot omit the label.
type DualResidual struct {
	QuaternionResidual
	Hypothetical bool
}

// MaxwellDual evaluates the dual-symmetric equations: R1 as in
// MaxwellQuaternion plus μ₀J_m on its vector part, R2 with μ₀ρ_m added to
// its scalar part (S.R2 = −div B + μ₀ρ_m).
func MaxwellDual(h Hypothetical, E, B, dEdt, dBdt QField, rho []float64, J *QField, rhoM []float64, JM *QField) (DualResidual, error) {
	if !h.MagneticCharge && (rhoM != nil || JM != nil) {
		return DualResidual{}, ErrHypotheticalDisabled
	}
	r := MaxwellQuaternion(E, B, dEdt, dBdt, rho, J)
	for n := range r.R1.Q {
		if JM != nil {
			r.R1.Q[n] = r.R1.Q[n].Add(JM.Q[n].Scale(Mu0))
		}
		if rhoM != nil {
			r.R2.Q[n] = r.R2.Q[n].Add(quat.Scalar(Mu0 * rhoM[n]))
		}
	}
	return DualResidual{QuaternionResidual: r, Hypothetical: true}, nil
}

// DualRotate applies the duality rotation by ξ to fields and charges.
func DualRotate(xi float64, e, b [3]float64, q, qm float64) ([3]float64, [3]float64, float64, float64) {
	c := SpeedOfLight
	cs, sn := cosSin(xi)
	var e2, b2 [3]float64
	for k := 0; k < 3; k++ {
		e2[k] = e[k]*cs + c*b[k]*sn
		b2[k] = (-e[k]*sn + c*b[k]*cs) / c
	}
	return e2, b2, q*cs + qm/c*sn, c * (-q*sn + qm/c*cs)
}

// DyonForce is the force and power on a particle with both charges:
// q(E + v×B) + q_m(B − v×E/c²) and qE·v + q_m B·v (hypothetical when
// q_m ≠ 0).
func DyonForce(h Hypothetical, q, qm float64, v, e, b [3]float64) ([3]float64, float64, error) {
	if qm != 0 && !h.MagneticCharge {
		return [3]float64{}, 0, ErrHypotheticalDisabled
	}
	fe, pe := LorentzVector(q, v, e, b)
	fm, pm := MagneticLorentzVector(qm, v, e, b)
	return [3]float64{fe[0] + fm[0], fe[1] + fm[1], fe[2] + fm[2]}, pe + pm, nil
}

func cosSin(x float64) (float64, float64) { s, c := math.Sincos(x); return c, s }
