# M9h — Maxwell's equations in quaternion form (`em/`)

CLAUDE.md §3c. Conventions are written down in `em/doc.go`: Hamilton's
signs, ∇v = −div v + curl v, the Treatise symbols (𝔄, ψ, 𝔅, 𝔈, 𝔇, ℌ,
ℭ, 𝔎, e) mapped to Go names, and the Gaussian ↔ SI conversion.

## What exists

| File | Content |
|---|---|
| `em/nabla.go` | `Grid`, quaternion fields `QField`, fourth-order `Partial`, Hamilton's `Nabla` (∇ on the left), `TimeDerivative` from four snapshots |
| `em/maxwell.go` | SI constants (CODATA 2018); `MaxwellQuaternion`: residuals of ∇E = −ρ/ε₀ − ∂B/∂t and ∇B = μ₀J + ∂E/∂t/c²; `FieldsFromPotentials`: 𝔅 = V.∇𝔄, 𝔈 = −∂𝔄/∂t − ∇ψ |
| `em/biquat.go` | the spacetime algebra Cl(1,3) (`hyper/ga`): F = E + IcB, `MaxwellSpacetime`: ∇F − J/(ε₀c), `SpacetimeParts` reading the four laws out of grades 1 and 3 |
| `em/vcalc.go` | plain vector calculus written separately (div, curl, grad, `MaxwellVector`) — the third twin |
| `em/frames.go` | orthogonal frame changes (proper and improper) by the vector rule and by GA reflections; Lorentz boosts by the textbook formula and by the spacetime rotor R̃FR |

The derivative is the fourth-order central difference on interior points
(two cells from every edge; NaN at the edges).

## Checks (`go test ./em`)

| Check | Result |
|---|---|
| **Triple twin**: quaternion, spacetime-algebra and vector forms on arbitrary (non-Maxwellian) fields and sources, all four residuals non-zero | agree to 1e−9 of each residual's scale (round-off) |
| Hamilton's ∇ on polynomials of degree ≤ 4: scalar part −div, vector part curl; ∇ψ = grad ψ | exact to 1e−12 |
| Inside a uniformly charged sphere, E = ρr/3ε₀; inside a uniform current cylinder, B = μ₀J(−y, x, 0)/2 | residuals 1e−9 of ρ/ε₀ and μ₀J |
| Plane wave, general direction and polarization: ∇F = 0 | 3.7e−3, 1.5e−4, 5.9e−6 of k·E₀ at 10, 20, 40 points/λ (ratio ≈ 25 per halving: fourth order) |
| Fields from potentials (Lorenz gauge) and a gauge change A → A + ∇χ, ψ → ψ − ∂χ/∂t | fields reproduced; gauge change leaves E and B unchanged |
| **Parity and rotations**: a charge's E and a current segment's B (Biot–Savart, closed form) recomputed for transformed geometry vs E → ME, B → det(M)·MB vs GA (B as bivector I·B, reflections) — for inversion, a mirror and 50 random proper/improper rotations | all three equal to 1e−12 |
| **Moving observer**: spacetime rotor vs Jackson eq. 11.149, 200 random fields, speeds to 0.9c; E² − c²B² and E·B invariant; slow probe past a static charge sees B′ = −v×E/c² | 1e−12 |

Cost: ∇ of a vector field on 32³ points 2.9 ms; the spacetime-algebra
residual on 16³ points 10 ms (sparse multivector products; it is the
cross-check, not the working form).

## Not done here

- The time-domain field solver (`em/fdtd.go`, M9o) and quasi-static
  extraction of R, L, C, K (`em/extract.go`, M9i), which build on these
  operators.
- Non-uniform grids and material interfaces (ε, μ varying in space).
