# Kaluza's fifth dimension, and the same idea in a Kron machine

Th. Kaluza, "On the Unification Problem in Physics," *Sitzungsberichte der
Preussischen Akademie der Wissenschaften* 1921, 966–972 (English: V. T.
Toth's revised translation). Kaluza added a fifth dimension to Einstein's
space-time and found that one 5-D geometry holds both gravity and
electromagnetism. A charge moving freely in 5-D is seen in 4-D as gravity
plus the Lorentz force, and **electric charge is momentum along the fifth
dimension**.

The link to kronspice is direct. Kaluza assumes nothing depends on the
fifth coordinate (his "cylinder condition"). Hoffmann (1949,
[hoffmann_motor](../hoffmann_motor)) uses the same words for Kron's
machines: the winding charges are coordinates that never appear in the
inductances.

```
go run ./examples/kaluza_kron -viz examples/kaluza_kron/viz   # about 1 second
python3 examples/kaluza_kron/viz/render.py
go test ./examples/kaluza_kron
```

## 1. Kaluza's equations check out

- **Eq. 3.** His closed forms for the new field quantities (Γ₀κλ = αF,
  Γκλ0 = −αΣ, Γ₀₀κ = −Γ₀κ₀ = g,κ) follow from his 5-D metric (eq. 2). On
  sample fields they agree to 2×10⁻¹³.
- **His coupling** α = √(κ/2) = 3.055×10⁻¹⁴ (he printed 3.06×10⁻¹⁴).
- **His electron e/m** = 1.77×10⁷ is e/m in electromagnetic units per gram
  (today 1.7588×10⁷).

## 2. A charge in a magnetic field as free motion in five dimensions

The exact 5-D geodesic (free-fall path, no approximations) for a charge
in a uniform magnetic field, using the modern form of Kaluza's metric with
the corner fixed. Units: c = G = 1, Gaussian charge, k = 2√G.

| Momentum along w | Charge ÷ mass | Orbit radius | Lorentz-force radius | Drift of p_w |
|---|---|---|---|---|
| 0.1 | 0.199 | 31.6054 | 31.6054 | 7×10⁻¹² |
| 1 | 1.414 | 4.4475 | 4.4475 | 2×10⁻¹⁴ |
| 10 | 1.990 | 3.1605 | 3.1605 | 7×10⁻¹⁵ |

The orbit frequencies also match the Lorentz force to six figures.
Nothing depends on the fifth coordinate, so momentum along it is
conserved. The 4-D path is exactly a Lorentz-force orbit with
q/m = 2·p_w/√(1 + p_w²). Kaluza's eq. 8a (q/m = 2α·u⁰) is the
small-momentum part of this curve.

## 3. Kaluza's difficulty with the electron, and Klein's answer

- **A ceiling on charge.** The curve can never exceed q/(m√G) = 2,
  however fast the motion along the fifth dimension. The electron needs
  2.04×10²¹ and the proton 1.11×10¹⁸. No real charged particle fits on a
  classical 5-D geodesic — the "serious difficulty" Kaluza reports on
  p. 971 (pointed out by Einstein).
- **Klein (1926).** The fifth dimension is a tiny circle. Quantum
  mechanics makes the momentum around it, and so the charge, come in
  whole steps n·e. The radius comes out at R = 2√G·ħ/(c·e) = 2 l_P/√α =
  **3.8×10⁻³⁴ m**, about 23 Planck lengths.
- **The price.** The first charged particle of this kind weighs
  5×10¹⁷ GeV, near the Planck mass — nothing like an electron.
- **Correction.** In my first written review I gave 1.3×10⁻³³ m. That
  paired one convention's coupling with another convention's unit of
  charge, a factor √(4π) too large.

## 4. The same construction in a Kron machine

A stator coil and a rotor coil (mutual inductance M·cos θ) with lossless,
shorted windings. Values (declared): L1 = L2 = 1 H, M = 0.5 H,
J = 0.01 kg·m²; it starts at rest with 1 A and 0.5 A in the coils.

- **The charges are cyclic coordinates.** The coordinates are the two
  winding charges and the rotor angle, and the inductances depend on the
  angle only. So the charges' momenta, the **flux linkages**, are
  conserved, just as Kaluza's momentum p_w is.
- **Eliminating them leaves the rotor alone** in a potential set by the
  conserved fluxes. This is the machine's version of the force the fifth
  dimension leaves behind.

Solved three ways:

| Time | Kron's geodesic in (q1, q2, θ) | Rotor alone, fluxes held | kronspice machine element |
|---|---|---|---|
| 0.5 s | −0.923503 rad/s | −0.923503 | −0.923217 |
| 1.0 s | 1.439010 | 1.439010 | 1.438928 |
| 2.0 s | 0.666266 | 0.666266 | 0.666265 |

- **Agreement.** The rotor-alone solution matches Kron's geodesic to
  10⁻¹⁰ rad/s, and kronspice to 4×10⁻⁴ rad/s (its windings are shorted
  through 10⁻⁷ Ω).
- **Conservation.** Flux linkages are conserved to 4×10⁻¹⁵ and energy to
  6×10⁻¹⁵. kronspice's ledger closes to 3×10⁻¹¹.

**Energy.** Nothing comes from the extra dimension: the 5-D path keeps
its norm and the machine keeps its energy.

Picture: `viz/kaluza_kron.png` (orbits, the charge ceiling, the machine
three ways).

## Where it stands

- **The mathematics lives on.** Kaluza–Klein reduction is a standard tool
  in string theory and supergravity.
- **No extra dimension has been observed.** Gravity follows Newton's
  inverse-square law down to 52 µm (Lee et al., *Phys. Rev. Lett.* 124,
  101101, 2020), and a circle of Klein's size is far beyond any
  experiment.
- **For kronspice.** Kaluza's geometry is exactly the structure of Kron's
  machine equations, and the machine shows it working.

## Files

| File | What it holds |
|---|---|
| `kaluza.go` | Kaluza's metric and three-index symbols; the KK metric for a uniform field; the 5-D geodesic |
| `checks.go` | eq. 3 check, Kaluza's α, charge-to-mass, Klein radius and mode mass |
| `machine.go` | the two-coil machine: Kron's geodesic, the rotor alone (Routh), kronspice |
| `main.go` | the report |
| `main_test.go` | eq. 3, geodesic = Lorentz, Klein numbers, machine three ways |
| `viz/render.py` | `kaluza_kron.png` |

## Sources

- Th. Kaluza, *Sitzungsber. Preuss. Akad. Wiss. Berlin (Math. Phys.)* 1921, 966–972; English translation by V. T. Toth (after T. Muta, HUPD-8401, 1984).
- O. Klein, "Quantentheorie und fünfdimensionale Relativitätstheorie," *Z. Phys.* 37, 895 (1926).
- B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Rev. Mod. Phys.* 21, 535 (1949) — the cylinder condition in Kron's machines.
- J. G. Lee et al., "New Test of the Gravitational 1/r² Law at Separations down to 52 µm," *Phys. Rev. Lett.* 124, 101101 (2020), arXiv:2002.11761.
