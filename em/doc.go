// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package em is the electromagnetic layer of the engine (CLAUDE.md §3c):
// the open-system energy ledger, the lumped-limit check, analytic coil
// formulas, and Maxwell's field equations on discrete fields in three
// equivalent forms.
//
// # Conventions (CLAUDE.md §3c rule 1)
//
// Quaternions follow Hamilton: i² = j² = k² = ijk = −1 (hyper/quat). A
// vector field is the pure quaternion v = vₓi + v_yj + v_zk, and Hamilton's
// operator ∇ = i∂ₓ + j∂_y + k∂_z acts from the left, so for a pure v
//
//	∇v = −div v + curl v      (S.∇v = −div v, V.∇v = curl v)
//
// and for a scalar field ψ, ∇ψ = grad ψ (Maxwell, Treatise, Vol. II,
// arts. 22–25 and 618; P. G. Tait, An Elementary Treatise on Quaternions,
// 3rd ed., Cambridge 1890, Ch. XII).
//
// Treatise symbols (Vol. II, art. 618, "General equations of the
// electromagnetic field") and their Go names, all in SI:
//
//	𝔄  electromagnetic momentum (vector potential)    A    [V·s/m]
//	ψ   electric potential                              Psi  [V]
//	𝔅  magnetic induction                              B    [T]
//	𝔈  electromotive intensity                         E    [V/m]
//	𝔇  electric displacement                           D = ε₀εr·E
//	ℌ  magnetic force                                  H = B/(μ₀μr)
//	ℭ  total current (conduction + displacement)       J + ∂D/∂t [A/m²]
//	𝔎  conduction current                              J
//	e   density of free electricity                     Rho  [C/m³]
//
// Gaussian (Treatise) ↔ SI: Maxwell writes 4π·ℭ = V.∇ℌ and 𝔇 = K𝔈/4π in
// the electromagnetic system of units. With the rationalized SI constants
// (CODATA 2018: ε₀ = 8.8541878128e−12 F/m, μ₀ = 1.25663706212e−6 H/m,
// c = 1/√(ε₀μ₀) = 299 792 458 m/s exactly) the same relations read
// curl H = J + ∂D/∂t and D = ε₀·E: the 4π is absorbed into the units
// (J. D. Jackson, Classical Electrodynamics, 3rd ed., Wiley 1999,
// Appendix on units and dimensions, Table 3).
//
// Maxwell's equations as two quaternion equations (Treatise relations
// combined; used by MaxwellQuaternion):
//
//	∇E = −ρ/ε₀ − ∂B/∂t          scalar part: div E = ρ/ε₀;  vector: curl E = −∂B/∂t
//	∇B =  μ₀J + (1/c²)·∂E/∂t    scalar part: div B = 0;     vector: curl B = μ₀J + ∂E/∂t/c²
//
// and as one equation of the spacetime algebra Cl(1,3) (D. Hestenes,
// Space-Time Algebra, 1966; C. Doran, A. Lasenby, Geometric Algebra for
// Physicists, Cambridge 2003, §7.1), used by MaxwellSpacetime: the basis
// e₁ … e₄ of hyper/ga's Cl(1,3) are γ₀ … γ₃ (γ₀² = +1, γₖ² = −1), the
// relative vectors σₖ = γₖγ₀, the pseudoscalar I = γ₀γ₁γ₂γ₃, the field
// bivector F = E + I·c·B with E = Σ Eₖσₖ, B = Σ Bₖσₖ, x⁰ = ct, and
// ∇ = γ^μ∂_μ (γ^0 = γ₀, γ^k = −γₖ):
//
//	∇F = J/(ε₀c),   J = cρ·γ₀ + Σ Jₖγₖ.
//
// The three forms (quaternion, spacetime algebra, plain vector calculus in
// vcalc.go) are equivalent; they are computed from the same discrete
// derivatives and must agree to round-off (triple-twin tests).
//
// # Discrete fields
//
// Fields live on a uniform Cartesian Grid. Spatial derivatives use the
// fourth-order central difference
//
//	f′(x) ≈ [f(x−2h) − 8f(x−h) + 8f(x+h) − f(x+2h)]/(12h),
//
// error h⁴·f⁽⁵⁾/30, defined at points at least two cells from every edge
// (the interior); edge points carry NaN. Time derivatives use the same
// stencil on snapshots at t−2δ, t−δ, t+δ, t+2δ.
package em
