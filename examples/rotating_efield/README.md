# Can a rotating electric field make a magnetic field?

The inverse of the homopolar dynamo of [homopolar](../homopolar): there a
conductor spinning in B gives a voltage; what does a *rotating E field* give?
The example keeps two questions apart (CLAUDE.md §3f, milestone M9n,
`doc/M9n.md`). First the real physics: a magnetic field comes from moving
charge (Rowland 1876) and from a changing electric field (the displacement
current of the Ampère–Maxwell law), checked with the M9h quaternion
operators. Then, opt-in and labelled, the hypothetical physics: the force a
magnetic charge would feel in these fields, the dual pair of Joot's
Theorem 3.16, kept out of the real momentum ledger.

```
go run ./examples/rotating_efield
go test ./examples/rotating_efield
```

## The five real cases

| Case | Set-up | Result | Expected |
|---|---|---|---|
| 1a. Rowland disk | σ = 10⁻⁵ C/m², R = 10 cm, 20,000 rpm; B_z on the axis by summing 400 current rings | 1.316×10⁻⁹ T at z = 0; 4.496×10⁻¹⁰ T at z = 5 cm | analytic (μ₀σω/2)[(R² + 2z²)/√(R² + z²) − 2\|z\|]: 1.316×10⁻⁹, 4.496×10⁻¹⁰ T |
| 1b. Rotating uniform E | 1 MV/m rotating at 1 MHz, as between quadrature-driven plate pairs | \|B\| = 1.75×10⁻⁶ T at r = 5 cm | (μ₀ε₀/2)·E·ω·r; Ampère–Maxwell residual 3.3×10⁻¹⁰, Faraday residual 4.2×10⁻⁷ (the neglected (ωr/c)² term) |
| 1c. Charged ring | 1 µC on 5 cm, at rest and spinning at 1000 rev/s | B at centre 0 T at rest, 1.26×10⁻⁸ T spinning | the E pattern is identical and static in both; B follows the moving charge (the mirror of the Faraday-disk paradox) |
| 1d. Rotating charge pair | ±1 nC on a 10 cm rod at 1000 rev/s | \|B\| = 8.99×10⁻¹² T on the axis 2.5 cm up | the E pattern really rotates; Ampère–Maxwell residual 1.3×10⁻⁴ of the displacement term (stencil truncation) |
| 1e. Röntgen–Eichenwald | εr = 3 disk, R = 10 cm, 1 mm thick, 50 rev/s in a fixed 1 MV/m | B_r = −3.46×10⁻¹⁰ T at r = 5 cm (ring sum) | thin-disk −(εr − 1)·E·v/c² = −3.5×10⁻¹⁰ T |

Case 1e is the electric mirror of the homopolar disk: a dielectric moving
in E is magnetized, M = P × v (Minkowski; Jackson §11.10). The ring sum
uses only the bound face charge ε₀(εr − 1)E carried round as two opposite
surface currents and Biot–Savart's μ₀, so the 1/c² = μ₀ε₀ comes out by
itself, not by hand. Its deficit against the thin-disk value falls
linearly in h/R (the test checks < 1.5 % at h/R = 0.01 and that halving h
halves it). That 1/c² is why the mirror effect is nanotesla where the
homopolar one is volts: 1 T at 100 m/s gives 100 V/m, while 1 MV/m at
100 m/s gives about 1 nT.

**Answer.** A rotating E pattern by itself makes no B. B comes from the
charge that moves (1a, 1c, 1d, 1e) or from ∂E/∂t where the field really
changes in time (1b, 1d), exactly as the Ampère–Maxwell law says; the
quaternion residuals confirm the fields are solutions.

## Hypothetical (opt-in, labelled)

With `em.Hypothetical{MagneticCharge: true}`, a magnetic charge
q_m = 1 A·m (none has ever been observed; CLAUDE.md §3f):

| | Hypothetical force | Real net force |
|---|---|---|
| At rest at the Rowland disk's centre | 1.32×10⁻⁹ N (q_m·B) | 0 (the momentum ledger's real column stays empty) |
| Moving at 100 m/s through 1 MV/m | 1.11×10⁻⁹ N (q_m·v×E/c²) | — |

These numbers go to the hypothetical column of `em.MomentumLedger` only.

## Files

| File | What it holds |
|---|---|
| `main.go` | the five real cases, the quaternion-operator residual checks and the hypothetical dyon force |
| `main_test.go` | each case against its analytic value; the Röntgen–Eichenwald h/R convergence |

## Sources

- H. A. Rowland, "On the magnetic effect of electric convection,"
  *Am. J. Sci.* 15, 30 (1878) — the rotating charged disk (1876).
- W. C. Röntgen, *Ann. Phys.* 35, 264 (1888); A. Eichenwald, *Ann. Phys.*
  11, 1 (1903) — magnetization of a dielectric moving in an electric
  field.
- J. D. Jackson, *Classical Electrodynamics*, 3rd ed., §11.10; W. K. H.
  Panofsky, M. Phillips, *Classical Electricity and Magnetism*, §18-3 —
  moving media, M = P × v.
- J. C. Maxwell, *Treatise*, Vol. II, Art. 619 — the displacement
  current in `4π ℭ = V.∇ℌ` (CLAUDE.md §3c).
- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16 —
  the dual (magnetic-charge) pair, reference only (§3f).
