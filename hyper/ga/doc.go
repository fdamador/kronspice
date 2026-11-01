// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package ga is geometric (Clifford) algebra Cl(p, q, r) for kronspice
// (CLAUDE.md §3b): the general sparse Multivector of any algebra up to 64
// basis vectors, and the fixed, allocation-free G2 and G3 for the plane
// and 3-space. It depends on the standard library and hyper/quat only.
//
// # Conventions (most GA bugs are convention bugs)
//
// Basis. An algebra Cl(p, q, r) has n = p + q + r basis vectors
// e₁ … eₙ. The first p square to +1, the next q to −1 and the last r to 0:
// Cl(3,0,0) is 3-space, Cl(1,3,0) spacetime with e₁ the time axis
// (signature + − − −, the convention of Hestenes, Space-Time Algebra,
// 1966, and of Doran & Lasenby, Geometric Algebra for Physicists, 2003,
// with their γ₀ … γ₃ as our e₁ … e₄), Cl(3,0,1) projective geometric
// algebra with e₄ the null vector.
//
// Blades. A basis blade is a bitmask: bit k−1 set means eₖ is a factor, in
// increasing index order, so e₁∧e₃ = e₁e₃ is 0b101 and the scalar is 0.
// The grade of a blade is its population count. Coefficients are always
// stored on the canonically ordered blade: e₃e₁ is −1·(0b101).
//
// Product of blades. eA·eB = s·eA⊕B, where ⊕ is XOR and the sign s is
// (−1)^(number of transpositions to sort the factors) times the product of
// eₖ² over the common factors (L. Dorst, D. Fontijne, S. Mann, Geometric
// Algebra for Computer Science, 2007, §19.1, "canonical reordering"). A
// common null factor makes the product 0.
//
// Products of multivectors (Dorst et al. Ch. 3 and 6):
//
//	geometric       AB
//	outer           A∧B = Σ ⟨⟨A⟩ᵣ⟨B⟩ₛ⟩ᵣ₊ₛ   (blades with no common factor)
//	left contraction  A⌋B = Σ ⟨⟨A⟩ᵣ⟨B⟩ₛ⟩ₛ₋ᵣ  (zero when r > s)
//	right contraction A⌊B = Σ ⟨⟨A⟩ᵣ⟨B⟩ₛ⟩ᵣ₋ₛ  (zero when s > r)
//	scalar product  A∗B = ⟨AB⟩₀
//
// Involutions. Reverse Ã reverses the order of the factors of each blade:
// sign (−1)^(k(k−1)/2) on grade k. Grade involution Â: (−1)^k. Clifford
// conjugate Ā = reverse of the grade involution: (−1)^(k(k+1)/2).
//
// Norm. A·Ã has scalar part Σ (±)coefficients²; Norm2 returns ⟨AÃ⟩₀, which
// is the squared Euclidean norm in a Euclidean algebra and can be negative
// or zero for other signatures; Norm returns √|⟨AÃ⟩₀|.
//
// Dual. A* = A·I⁻¹, with I = e₁e₂…eₙ the unit pseudoscalar; defined only
// when I is invertible (r = 0).
//
// Rotors. A rotor R = e^(−Bθ/2), B a unit bivector, rotates vectors by
// R·v·R̃ through θ in the plane of B. In G3, B = I·n for the axis n.
//
// Quaternions. The even subalgebra of G3 is Hamilton's quaternions under
// i ↔ −e₂₃, j ↔ −e₃₁ = e₁₃, k ↔ −e₁₂ (so ij = k holds, and the rotor
// cos(θ/2) − sin(θ/2)·I·n corresponds to the quaternion cos(θ/2) +
// sin(θ/2)·n, both rotating by θ about n; Doran & Lasenby §2.4.3). The
// even subalgebra of G2 is the complex numbers under i ↔ e₁₂.
package ga
