// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "github.com/fdamador/kronspice/hyper/quat"

// SI field constants (CODATA 2018; c is exact, ε₀ = 1/(μ₀c²)).
const (
	Mu0  = 1.25663706212e-6 // vacuum permeability [H/m]
	Eps0 = 1 / (Mu0 * SpeedOfLight * SpeedOfLight)
)

// QuaternionResidual holds the residuals of Maxwell's equations as two
// quaternion equations (doc.go):
//
//	R1 = ∇E + ρ/ε₀ + ∂B/∂t    (S.R1 = −(div E − ρ/ε₀), V.R1 = curl E + ∂B/∂t)
//	R2 = ∇B − μ₀J − ∂E/∂t/c²  (S.R2 = −div B,          V.R2 = curl B − μ₀J − ∂E/∂t/c²)
type QuaternionResidual struct{ R1, R2 QField }

// MaxwellQuaternion evaluates Maxwell's equations in the quaternion form
// of the Treatise (Vol. II, art. 618, in SI; doc.go): E and B as pure
// quaternion fields, their time derivatives, charge density rho and current
// density J (nil for 0).
func MaxwellQuaternion(E, B, dEdt, dBdt QField, rho []float64, J *QField) QuaternionResidual {
	nE, nB := Nabla(E), Nabla(B)
	g := E.G
	r := QuaternionResidual{R1: QField{G: g, Q: make([]quat.Q, g.Len())}, R2: QField{G: g, Q: make([]quat.Q, g.Len())}}
	c2 := SpeedOfLight * SpeedOfLight
	for n := range r.R1.Q {
		q := 0.0
		if rho != nil {
			q = rho[n]
		}
		j := quat.Q{}
		if J != nil {
			j = J.Q[n]
		}
		r.R1.Q[n] = nE.Q[n].Add(quat.Scalar(q / Eps0)).Add(dBdt.Q[n])
		r.R2.Q[n] = nB.Q[n].Sub(j.Scale(Mu0)).Sub(dEdt.Q[n].Scale(1 / c2))
	}
	return r
}

// FieldsFromPotentials returns 𝔅 = V.∇𝔄 and 𝔈 = −∂𝔄/∂t − ∇ψ (Treatise,
// Vol. II, arts. 604 and 598, for conductors at rest; SI) from the vector
// potential A, its time derivative and the scalar potential ψ.
func FieldsFromPotentials(psi []float64, A, dAdt QField) (E, B QField) {
	g := A.G
	nA := Nabla(A)
	gp := Nabla(Scalars(g, psi))
	E = QField{G: g, Q: make([]quat.Q, g.Len())}
	B = QField{G: g, Q: make([]quat.Q, g.Len())}
	for n := range E.Q {
		B.Q[n] = quat.Vec(nA.Q[n].X, nA.Q[n].Y, nA.Q[n].Z) // V.∇𝔄
		E.Q[n] = dAdt.Q[n].Scale(-1).Sub(quat.Vec(gp.Q[n].X, gp.Q[n].Y, gp.Q[n].Z))
	}
	return E, B
}

// ToV converts a quaternion field's vector part to a VField.
func ToV(f QField) VField {
	out := VField{G: f.G, V: make([][3]float64, len(f.Q))}
	for n, q := range f.Q {
		out.V[n] = [3]float64{q.X, q.Y, q.Z}
	}
	return out
}

// ToQ converts a VField to a pure quaternion field.
func ToQ(v VField) QField {
	out := QField{G: v.G, Q: make([]quat.Q, len(v.V))}
	for n, x := range v.V {
		out.Q[n] = quat.Vec(x[0], x[1], x[2])
	}
	return out
}
