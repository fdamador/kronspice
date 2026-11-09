// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Step 2: the relaxed knot as a quantum rotor.
//
// Faddeev's Lagrangian, L = κ2(ṅ² − Σ(∂_i n)²) + κ4(Σ_i F_0i² − Σ_{i<j} F_ij²)
// with F_μν = n·(∂_μ n × ∂_ν n) (c = 1), has the static energy of
// ether_knot and the kinetic energy
//
//	T = κ2 ∫ |ṅ|² + κ4 ∫ Σ_i |ṅ × ∂_i n|².
//
// A rigid rotation with angular velocity ω about an axis gives ṅ = ω v,
// v = −(ê × x)·∇n, so T = ½ I ω² with I = 2[κ2 ∫|v|² + κ4 ∫ Σ_i |v × ∂_i n|²].
// About the symmetry axis v = ∂_φ n (a rotation there is the same as an
// isorotation about the vacuum, so the two share one inertia I∥); about
// a perpendicular axis the integrand depends on φ and is averaged over
// it here. Isorotations that move the vacuum have infinite inertia and
// are frozen. Units: energies in E0 = √(κ2κ4), lengths in ℓ = √(κ4/κ2),
// so I is in E0 ℓ²/c².
//
// Finkelstein–Rubinstein quantization (Krusch & Speight, Commun. Math.
// Phys. 264, 391, 2006): a Hopf soliton can be quantized as a fermion
// only when its charge Q is odd. For the symmetric top with J half-odd,
// the lowest state J = ½, |J3| = ½ has
//
//	E_rot = ħ²[J(J+1) − J3²]/(2I⊥) + ħ²J3²/(2I∥) = ħ²/(4I⊥) + ħ²/(8I∥).

import "math"

// Moments holds the energy-weighted radii and the inertias of a knot.
type Moments struct {
	E                 float64 // static energy κ2E2 + κ4E4 (node-local differences)
	Rho2, R2          float64 // energy-weighted ⟨ρ²⟩ and ⟨r²⟩
	IPar, IPerp       float64 // inertia about the axis and about a perpendicular axis
	IParVia, IParFull float64 // I∥ from the φ-average (check) and from v = ∂_φ n directly
}

// Inertia integrates the moments over the interior nodes (central
// differences, as ether_knot's Size), averaging the perpendicular
// rotation over nphi angles.
func (a *Axi) Inertia(nphi int) Moments {
	var mo Moments
	h := a.H
	w := 2 * math.Pi * h * h
	var se, sr, srr float64
	for j := 1; j < a.NZ-1; j++ {
		z := a.z(j)
		for i := 0; i < a.NR-1; i++ {
			r := a.rho(i)
			_, mr, mz, c := a.derivs(i, j)
			cr := scale(c, 1/r) // (1/ρ)∂_φ n in the unrotated frame
			x1, x2, x3 := cross(mr, mz), cross(mr, cr), cross(mz, cr)
			e := a.K2*(dot(mr, mr)+dot(mz, mz)+dot(cr, cr)) + a.K4*(dot(x1, x1)+dot(x2, x2)+dot(x3, x3))
			wv := w * r
			se += wv * e
			sr += wv * e * r * r
			srr += wv * e * (r*r + z*z)
			// about the axis: v = ∂_φ n = c
			ip := a.K2*dot(c, c) + a.K4*(dot(cross(c, mr), cross(c, mr))+dot(cross(c, mz), cross(c, mz)))
			mo.IParFull += 2 * wv * ip
			// about x, averaged over φ; all vectors in the frame rotated back
			// by R_z(−Mφ), which keeps every dot and cross product.
			var iperp, ipar float64
			for k := 0; k < nphi; k++ {
				ph := 2 * math.Pi * (float64(k) + 0.5) / float64(nphi)
				cs, sn := math.Cos(ph), math.Sin(ph)
				dx := sub(scale(mr, cs), scale(cr, sn))
				dy := add3(scale(mr, sn), scale(cr, cs))
				dz := mz
				x, y := r*cs, r*sn
				v := sub(scale(dz, y), scale(dy, z)) // (y∂_z − z∂_y) n
				q := a.K2*dot(v, v) + a.K4*(sq(cross(v, dx))+sq(cross(v, dy))+sq(cross(v, dz)))
				iperp += q
				u := sub(scale(dy, x), scale(dx, y)) // (x∂_y − y∂_x) n = ∂_φ n
				ipar += a.K2*dot(u, u) + a.K4*(sq(cross(u, dx))+sq(cross(u, dy))+sq(cross(u, dz)))
			}
			mo.IPerp += 2 * wv * iperp / float64(nphi)
			mo.IParVia += 2 * wv * ipar / float64(nphi)
		}
	}
	mo.E, mo.Rho2, mo.R2, mo.IPar = se, sr/se, srr/se, mo.IParFull
	return mo
}

func add3(p, q [3]float64) [3]float64 { return [3]float64{p[0] + q[0], p[1] + q[1], p[2] + q[2]} }
func sq(p [3]float64) float64         { return dot(p, p) }

// SpinCoefficient is c_r in E_rot = c_r (ħc)²/(E0 ℓ²) for the J = ½,
// |J3| = ½ state.
func (mo Moments) SpinCoefficient() float64 {
	return 1/(4*mo.IPerp) + 1/(8*mo.IPar)
}
