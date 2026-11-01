# M9d — quaternions and geometric algebra (`hyper/quat`, `hyper/ga`)

CLAUDE.md §3b. Both packages are leaves that use the standard library
only (`hyper/ga` also imports `hyper/quat` for the even-subalgebra
conversion). The plan named gonum's `num/quat`; the project has kept no
third-party dependencies since M6, and the quaternion algebra is a few
dozen lines, so it is written here instead.

## What exists

- **`hyper/quat`** — Hamilton's quaternions with Maxwell's sign
  convention (the scalar part of a product of vectors is minus the dot
  product): product, conjugate, norm, inverse, unit, exp/log, rotors
  (axis–angle), `r·v·r*`, the rotation matrix and Shepperd's inverse,
  slerp.
- **`hyper/ga`** — conventions written down in `hyper/ga/doc.go`.
  - `Algebra` = Cl(p, q, r) for up to 64 basis vectors, blades as
    bitmasks, `Multivector` sparse and sorted, products in a deterministic
    order.
  - Geometric, outer, left and right contraction, scalar product; reverse,
    grade involution, Clifford conjugate; grade projection, norm; inverse
    (versor formula, else the 2ⁿ linear system up to n = 10); dual (only
    for a non-degenerate metric); exp (closed form when x² is a scalar,
    otherwise scaling and squaring); rotors and the sandwich product.
  - `G2` and `G3`: fixed arrays, allocation-free products (18 ns per
    G3 product). `G2.Complex`, `G3.Quat` and `G3FromQuat` are the
    isomorphisms.

## Checks (`go test ./hyper/...`)

| Check | Result |
|---|---|
| Hamilton's relations i² = j² = k² = ijk = −1, ij = k … | exact |
| Associativity, distributivity, \|pq\| = \|p\|\|q\|, (pq)* = q*p*, uv = −u·v + u×v | 2000 random cases, ≤ 1e−14 |
| Rotations: matrix = `r·v·r*`, length kept, matrix ↔ quaternion round trip, exp∘log, slerp | 1000 cases, ≤ 1e−14 |
| Algebra laws in Cl(2,0,0), (3,0,0), (0,2,0), (1,3,0), (3,0,1), (2,2,1), (4,1,0): associativity, distributivity, reverse and conjugate anti-automorphisms, grade involution automorphism, eₖ² = metric, eᵢeⱼ = −eⱼeᵢ | 40 random full multivectors each |
| Outer product: grades add, u∧u = 0, antisymmetric and associative, a⌋(b∧c) = (a·b)c − (a·c)b, uv = u·v + u∧v, v₁∧…∧v₅ = det·I | 200 cases |
| **Twin: matrix representations** — Cl(3,0,0) on Pauli matrices, Cl(1,3,0) on Dirac γ matrices | every random product equal to the matrix product, 1e−12 |
| **Twin: quaternions** — Cl(0,2,0) ≅ ℍ; G3 even part ↔ `quat`: conversion round-trips exactly, products agree to 1e−15, `Rotor3` = `quat.Rotor`, both rotate vectors identically, R·R̃ = 1 | 1000 cases |
| Fixed G2/G3 equal the sparse algebra; G2 even part = complex product; no allocations | 500 cases |
| Inverse (x·x⁻¹ = 1 in four signatures, v⁻¹ = v/v²), a null vector refused, no dual in a degenerate algebra, e₁* = −e₂₃, exp closed form and series, rotor turns e₁ into e₂, a boost in Cl(1,3) keeps the interval | as stated |
| Determinism: identical products bit for bit | 20 repeats |

## Not done here

- Rotor logarithm in a general algebra (only the quaternion log).
- Phasor helpers for polyphase networks: they belong with M9f
  (`kron/polyphase`), which will use `hyper/quat`.
- Test vectors from ganja.js or the Python `clifford` package: the matrix
  representations serve as the independent reference instead.
