# Kelvin's 1899 charged-ring molecule

Lord Kelvin's "Magnetism and Molecular Rotation" (1899) pictures a molecule
as a rigid charged ring filled with an oppositely charged fluid that can
slide round it. This example builds that ring with standard electrodynamics
only, every charge under the Lorentz force, and runs the two things Kelvin
did with it: induction when a field is switched on, and his own Zeeman test,
which made him give the model up. It is the electron test bench that
[space_vortex](../space_vortex) (H13), [kron_electron](../kron_electron)
(H21) and [photon_electron](../photon_electron) (H23) are measured on
(`doc/kelvin-references.md`); the ring's hydrodynamic cousin is
[vortex_ring](../vortex_ring).

```
go run ./examples/kelvin_ring
go test ./examples/kelvin_ring          # -short skips the Zeeman ensemble
```

## The construction

The ring: a non-conducting tube of radius a carrying +e (mass w′), filled
with an incompressible fluid carrying −e (mass w) that slides without
friction. Twenty-four charges per half are integrated with RK4.

1. **Induction (Kelvin §§1–3).** A coaxial coil's current is ramped on,
   held and ramped off. Its field comes from Biot–Savart (`em.Loop`), the
   induced E from Faraday's law, and the force on each charge from
   `em.LorentzVector`. The flux is checked against Maxwell's
   elliptic-integral mutual inductance (`em.MutualCircles`).
2. **Zeeman test (Kelvin §§5–7).** The ring is held by isotropic angular
   springs (ω₀ = 1) in a steady field, kicked, and the spectrum of its
   motion and of its radiating magnetic moment is taken over twelve
   orientations. Three models are compared: Lorentz's bound charge (the
   standard classical Zeeman model), Kelvin's gyrostat as he wrote it
   (fixed internal spin along the axis) and the full ring (fluid free to
   slip, rigid-body rotation with `hyper/quat`, nothing linearised).

## Results

**Part 1**, ring a = 5 cm, ±1 µC, tube and fluid 10 g each; coil r = 10 cm,
5 cm below, 20 kA·turns, 1 s ramps:

| Quantity | Simulated | Expected |
|---|---|---|
| Flux through the ring | 6.98736×10⁻⁴ Wb (Biot–Savart) | 6.98732×10⁻⁴ Wb (Maxwell, elliptic) |
| Moment of momentum, tube / fluid | −1.11207×10⁻¹⁰ / +1.11207×10⁻¹⁰ kg·m²/s | ∓eΦ/2π = 1.11207×10⁻¹⁰ (Kelvin's corrected eAM/2π) |
| Net angular impulse on the ring | −5.5×10⁻³⁰ | 0: the two halves cancel, the coil gets none back |
| Fluid angular velocity | 4.44829×10⁻⁶ rad/s | Larmor's eB̄/2w = 4.44829×10⁻⁶ rad/s |
| Ring moment | −1.11×10⁻¹⁴ A·m² (diamagnetic), 1.44×10⁻¹⁴ N repulsion | Langevin diamagnetism |
| Ledger | source work −4.94683×10⁻¹⁶ J = kinetic 4.94683×10⁻¹⁶ J + mutual-field −9.89365×10⁻¹⁶ J | residual −8×10⁻¹⁶ |
| After switch-off | tube 1.5×10⁻²⁶, fluid −1.5×10⁻²⁶ kg·m²/s; net source work −5×10⁻³² J | everything at rest |

Kelvin's ring is a mechanical picture of Larmor's theorem. At atomic scale
(electron charge and mass, 1 T) the same fluid turns at 8.794×10¹⁰ rad/s.

**Part 2**, ω₀ = 1, ωL = eB/2w = 0.05, 12 orientations, resolution
2π/T = 0.0026:

| Model | Splitting vs cos θ | Ensemble upper line (FWHM) | Circular polarisation along B, radiated (upper / lower) | Across B: centre ∥B, sides ⊥B |
|---|---|---|---|---|
| Lorentz bound charge | 0.0513 at every θ | sharp, 0.0038 | −1.00 / +1.00 | 100 % / 100 % |
| Kelvin gyrostat as written | (ωL/2)·cos θ: 0 … 0.0243 | broad, 0.0188 | +0.72 / −0.71 | 67 % / 59 % |
| Kelvin ring, full Lorentz force | 0.0243 … 0.0260 at every θ | sharp, 0.0037 | −0.02 / −0.03 (motion −0.47 / +0.46; analytic 5/11) | 64 % / 71 % |

- **The gyrostat** reproduces Kelvin's reason for discarding the model
  (p. 357): with a fixed spin the splitting goes as cos θ, so randomly
  oriented molecules broaden the line instead of splitting it.
- **The full ring** is new. When the ring tilts the flux through it
  changes and the fluid's spin follows, since its canonical angular
  momentum w·a²·u − eΦ/2π is conserved. The gyroscopic vector becomes
  G = (ea²/2)(B − 2(n·B)n), the reflection of B in the ring's plane, whose
  length does not depend on orientation (`GyroVector`). Every ring splits
  by eB/2(w + w′) and the lines are sharp, which removes Kelvin's and
  Larmor's objection.
- **It still fails Zeeman's test.** The modes circle G, not B, so the
  light from the ring's own moment is nearly unpolarised along the field
  where Zeeman measured ±1, and the splitting is set by the total mass
  w + w′, not the electron's e/m. A ring of two opposite charges with
  different charge-to-mass ratios cannot give the normal Zeeman effect;
  Lorentz's single charge can, because all of its charge has one e/m
  (the condition of Larmor's theorem).

Any toroidal model of the electron must therefore give a splitting set by
e/m_e and the observed circular polarisations before the anomalous
effect, g and spin ½ are even asked for.

## Files

| File | What it holds |
|---|---|
| `main.go` | the ring, the coil, the induction run, the three Zeeman models and the ensemble spectra |
| `main_test.go` | induction identities and ledgers; the Zeeman test for all three models; \|G\| constant |

## Sources

- W. Thomson (Lord Kelvin), "Magnetism and Molecular Rotation," *Proc. R.
  Soc. Edinb.* 22, 631–635 and *Phil. Mag.* 48, 236–239 (1899);
  *Mathematical and Physical Papers* Vol. V, Art. 175, pp. 354–358
  (`reference/mathematicalphys05kelvuoft.pdf`); notes in
  `doc/kelvin-references.md`.
- H. A. Lorentz, the bound-charge model of the normal Zeeman effect;
  P. Zeeman, *Phil. Mag.* 43, 226 (1897).
- J. Larmor, Larmor's theorem (1897), quoted by Kelvin on p. 357.
- J. C. Maxwell, *Treatise*, Vol. II — mutual inductance of coaxial
  circles (`em.MutualCircles`).
