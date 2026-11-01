# Tewari's space-vortex electron on the test bench (H13)

P. Tewari's 2018 paper describes the electron as an irrotational vortex of
an incompressible, massless "space fluid" whose speed reaches c at the
surface of a central spherical void of radius r_e. From that void it
states the electron's charge, mass, angular momentum and magnetic moment,
and an electron orbit for hydrogen. This example puts those statements on
the electron test bench of [kelvin_ring](../kelvin_ring)
(`doc/kelvin-references.md`) and compares each with measurement. Register
entry H13 (CLAUDE.md §1); notes and the page checks in
`doc/tewari-references.md`; his generator test is
[space_power_generator](../space_power_generator).

```
go run ./examples/space_vortex
go test ./examples/space_vortex
```

## The construction

The paper's equations, in its own units (charge in cm³/s, mass in cm⁴/s):

```
q_e = π² r_e² c               (eq. 4)
m_e = (4π/3) r_e³ c           (eq. 6)
L   = (4/5) m_e c r_e         (eq. 15, angular momentum)
μ   = (3/4) q_e c r_e         (eq. 16, magnetic moment)
```

r_e is fixed from the measured charge by eq. 4, as in the paper, which
gives r_e = 4.029×10⁻¹¹ cm. The paper defines its gram from the measured
electron mass (eqs. 10–12), so charge and mass are both inputs. Eq. 16 as
printed has the dimensions of an SI moment; the Gaussian moment of the
same loop is that divided by c, and the g-factor is the same in either
reading. Measured values are CODATA 2018 and the NIST Atomic Spectra
Database (Gaussian units). The standard counterpart beside every row is
the Dirac electron with QED corrections, which gives the measured values.

## Results

| Quantity | Model | Measured | Model / measured | Note |
|---|---|---|---|---|
| Charge (eq. 4) | 4.803×10⁻¹⁰ statC | 4.803×10⁻¹⁰ | 1 | input: r_e fitted to e |
| Void radius vs classical radius | 4.029×10⁻¹¹ cm | 2.818×10⁻¹³ | 143 | no measured size; scattering bounds ≲ 10⁻¹⁸ cm |
| Angular momentum (eq. 15) vs ħ/2 | 8.802×10⁻²⁸ erg·s | 5.273×10⁻²⁸ | 1.669 | |
| Magnetic moment (eq. 16) | 1.451×10⁻²⁰ erg/G | 9.285×10⁻²¹ | 1.563 | \|μ_e\| = (g/2) μ_B |
| g-factor (eqs. 15–16) | 1.875 = 15/8 | 2.00232 | 0.9364 | independent of r_e |
| g with the ring moment ½qvρ | 0.9375 = 15/16 | 2.00232 | 0.4682 | p. 115 integrates dq·v·ρ; a ring current's moment is half that |
| Proton charge (eq. 30) | 12 e | e | 12 | atoms neutral to ~10⁻²¹ e |
| Fine-structure constant (UPSM eq. 5.2) | 15π²/16 = 9.253 | 7.297×10⁻³ | 1268 | book p. 130 blames experiment; α is measured to 1.5×10⁻¹⁰ |
| Density limit (UPSM §2.8) vs Sirius B | 3420 g/cm³ | ~2×10⁶ | 0.0017 | white dwarfs exceed it |
| Hydrogen orbit radius (eq. 35) vs a₀ | 9.67×10⁻¹⁰ cm | 5.292×10⁻⁹ | 0.183 | |
| Hydrogen orbit speed vs αc | 2.069×10¹⁰ cm/s (0.69 c) | 2.188×10⁸ | 94.6 | |
| Orbit frequency vs Lyman-α | 3.405×10¹⁸ Hz | 2.466×10¹⁵ | 1381 | a classical orbit radiates at v/2πr |

- **The g-factor is a real prediction.** g = 2m_e c μ/(e L) from eqs.
  15–16, and r_e, q and m all cancel, so 15/8 (or 15/16 with the ring
  moment) holds for any void radius. It is 6.4 % from the measured value,
  which is known to about 10⁻¹³.
- **The hydrogen orbit gives the wrong light.** A charge orbiting at
  0.69 c on 24 r_e radiates X-rays; hydrogen's lines are ultraviolet and
  its ground state has v = αc.
- **Zeeman effect: not treated in the paper.** The bench's discriminating
  test needs the model's prediction for the splitting and polarisation of
  a line in a field B (measured: Δω = eB/2m_e c per unit of m_l, σ±
  circularly polarised). Until that prediction exists H13 stays
  *untested*, with the rows above as the values it must also meet.

The test reproduces the paper's own numbers where it states them (r_e ≈
4×10⁻¹¹ cm, eq. 4 giving 4.7326×10⁻¹⁰ at that radius, 1 g = 8.8×10⁶ cm⁴/s
by eq. 11, void volume 2.67×10⁻³¹ cm³) before recording the comparisons.

## Files

| File | What it holds |
|---|---|
| `main.go` | the model's equations, the measured constants and the bench table |
| `main_test.go` | the paper's stated numbers, g = 15/8, α = 15π²/16 |

## Sources

- P. Tewari, "Structural relation between the vacuum space and the
  electron," *Physics Essays* 31(1), 108–129 (2018), DOI
  10.4006/0836-1398-31.1.108 (`reference/tewari_2018_physics_essays.pdf`).
- P. Tewari, *Universal Principles of Space and Matter*, Grain and Chaff,
  Toronto, 2017 (`reference/tewari_upsm_2017.pdf`) — eq. 5.2 (α), §2.8
  (density limit).
- CODATA 2018 (E. Tiesinga et al., *Rev. Mod. Phys.* 93, 025010, 2021);
  NIST Atomic Spectra Database (Lyman-α).
