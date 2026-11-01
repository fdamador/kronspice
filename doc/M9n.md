# M9n — Lorentz force, the homopolar dynamo, rotating E fields, momentum

CLAUDE.md §3f. Theorem 3.16 of P. Joot, *Geometric Algebra for
Electrical Engineers* (v0.3.6.4), is implemented from its equations
(no text or code taken; the repository has no license file).

## What exists

| File | Content |
|---|---|
| `em/lorentz.go` | F = E + IcB and the charge multivectors q(1 + v/c) and I(q_m/c)(1 + v/c) in G3; `LorentzGA`: force ⟨FQ†⟩₁ and power c⟨FQ†⟩₀; the vector forms; `MotionalEMF` — the EMF of a moving conductor from the Lorentz force on its charges |
| `em/duality.go` | **hypothetical**: dual-symmetric Maxwell (`MaxwellDual`: div B = μ₀ρ_m, curl E = −∂B/∂t − μ₀J_m), the duality rotation, dyon force; refused unless `Hypothetical.MagneticCharge` is set, results labelled |
| `em/biotsavart.go` | closed-form segment field, current loops, force on a loop in a field (Gauss–Legendre), circles |
| `em/momentum.go` | the **momentum ledger**: force rows (on, from, exterior, hypothetical); `Net`, `ByExterior`, `Hypothetical` (never in Net), `Check` (Newton's third law for internal pairs) |
| `examples/homopolar` | Tesla's two-disk dynamo (US 406,968) as two machine elements on a belt, vs a single disk with a rim brush; GA for efficiency at a power target |
| `examples/rotating_efield` | B from moving charge and from changing E, four cases; the hypothetical magnetic-charge force, labelled |
| `examples/tether` | electrodynamic tether in LEO: generator (drag) and thruster modes, momentum ledger |

Magnetic charge uses SI ampere-metres (force q_m·B at rest); Joot's
normalization is ours divided by c.

## Checks

| Check | Result |
|---|---|
| Grade selection = vector forms, electric and magnetic pairs | 2000 random cases, 1e−12 |
| Dyon force invariant under the duality rotation by any ξ; ξ = π/2 maps the electric pair onto the magnetic pair; magnetic charge refused outside hypothetical mode | 1e−10; as stated |
| Homopolar EMF from the Lorentz force along a radius = ½Bω(r_o² − r_i²) | 1e−12 |
| Dual Maxwell: uniformly magnetically charged ball, B = μ₀ρ_m r/3; ordinary Maxwell flags the same field (−div B = −μ₀ρ_m) | exact |
| Charging capacitor: B = μ₀ε₀Ė(−y, x)/2 satisfies Ampère–Maxwell with no conduction current | 1e−9 |
| Rowland disk: B on the axis by summing rings = (μ₀σω/2)[(R² + 2z²)/√(R² + z²) − 2\|z\|] | 1e−4 |
| Two current loops: equal and opposite forces (ledger `Check`); far apart → dipole–dipole 3μ₀m₁m₂/(2πz⁴) | 1e−6; 1 % |
| **Homopolar circuit**: two-disk open-circuit voltage = 2·EMF; shaft power from the ledger = Σ torque·ω; torque = K·I; ledger closes; two-disk with belt gives > 3× the single disk's power into the same load; GA meets 2 kW at η > 0.95 | as stated |
| **Rotating E**: rotating uniform E's B satisfies Ampère–Maxwell (3e−10); spinning vs resting charged ring — same static E, B = μ₀I/2a vs 0; rotating charge pair consistent with its displacement current (1.3e−4, stencil truncation) | as stated |
| **Tether**: EMF = vBL; generator F·v = −EMF·I, load + losses = EMF·I; thruster supply = losses + F·v; momentum reaction booked on the field source | 1e−9 |

## Results worth noting

- Tesla's arrangement: at 3000 rpm in 0.8 T, two 15 cm copper disks on a
  belt deliver 3.04 kW at η = 0.98 into 10 mΩ; one disk with a
  500 µΩ rim brush delivers 0.71 kW (contact resistances are illustrative
  inputs, not measurements).
- A uniform E of 1 MV/m rotating at 1 MHz makes about 1.75 µT at 5 cm —
  real, and small. A charged ring makes B only when its charge moves; an
  E pattern that is symmetric under the rotation does not change and makes
  none.
- Tether: 5 km in LEO gives 1155 V; as a generator 0.2 N of drag for
  884 W into a load; as a thruster 0.75 N from 15 kW (0.05 N/kW), the
  reaction going to the geomagnetic field and ionosphere.

## Not done here

- Contactor current–voltage models (orbit-motion-limited collection) and
  the charged-body buoyancy study of `doc/em-push-pull.md`.
- Armature reaction and brush heating in the homopolar model; the field
  from a field coil (here B is given per ampere).
