# Toroidal coherence propulsion (Konstapel 2025): the Tier 1 bench

J. Konstapel, "Toroidal Coherence Engineering: From Hall's Photon Theory
to Practical Propulsion Systems" (Leiden, 1 Nov 2025, 17 pp.) proposes that
heated fiber-optic toroidal coils, phase-locked to the Earth's
electromagnetic field, invert a scalar "Sn field" (identified with the
scalar part of Maxwell's quaternion equations, after V. Robinson and
C. Hall) and reduce the inertial mass of a craft. Part 5 checks Hall's
explanation of lensed galaxy images as frustrated total internal
reflection at hydrogen clouds. Register entry H35
(CLAUDE.md §1).

```
go run ./examples/toroidal_coherence
go test ./examples/toroidal_coherence
```

The paper gives no equation with a value for its coupling σ, so the
hypothesis cannot be put into the simulator as a hypothetical model
(§3f). This example computes the **standard prediction** for each of the
paper's three Tier 1 measurements (§6.1), and the ordinary effects a bench
must remove before a result can be read.

## Part 1 — the (3n ± 1) "topological signature"

The paper predicts harmonics at (3n ± 1) × the fundamental. Three coils at
120°, each 1 mH + 2 Ω closed through 20 Ω in parallel with a diode (an
ordinary rectifying contact), driven at 10 kHz (inside the paper's
1–100 kHz), in kronspice transient. The array field is the space vector
(2/3)(i_a + a·i_b + a²·i_c), the field on the axis of three coils at 120°.

| n | one coil | array, forward | array, backward | linear control |
|---|---|---|---|---|
| 1 | 1 | 1 | 1.7e-9 | 1 |
| 2 | 3.25e-2 | 1.5e-10 | 3.25e-2 | 2.4e-15 |
| 3 | 1.06e-3 | 5.5e-10 | 5.3e-10 | 1.6e-15 |
| 4 | 3.12e-3 | 3.12e-3 | 1.0e-10 | 1.0e-15 |
| 5 | 2.58e-4 | 3.0e-10 | 2.58e-4 | 1.3e-15 |
| 6 | 8.50e-4 | 1.0e-10 | 1.2e-10 | 1.7e-16 |
| 7 | 1.01e-4 | 1.01e-4 | 2.1e-10 | 8.0e-16 |

Ledger: 11.98 mJ in − 8.65 mJ returned = 10.48 mJ dissipated − 7.16 mJ
released from storage (start-up), integration mismatch 1.6 nJ, worst
point residual 0.

**Standard prediction:** a balanced three-coil array with *any* ordinary
nonlinearity radiates only 3n ± 1 (3n+1 forward-rotating, 3n−1 backward),
because multiples of 3 cancel by symmetry; a single coil radiates every
harmonic. The paper's Tier 2 array would show 3n ± 1 with no new physics.
For the single Tier 1 coil the only unexpected result would be the
*absence* of multiples of 3, which the paper does not state.

## Part 2 — discrete infrared lines

A closed 10 m silica loop (n_g 1.46, declared) has longitudinal modes
every c/(n_g·ℓ) = 20.53 MHz; every closed loop shows discrete lines
(`examples/photon_electron` ring lines: AC peaks at n·v/ℓ). Thermal light
at 1000 K peaks at 2.90 µm with coherence length of order that wavelength,
3×10⁻⁷ of the loop: the loop filters it into modes but cannot phase-lock
it. A test needs a stated line spacing or linewidth that differs from
c/(n_g·ℓ).

## Part 3 — the balance

Rig (declared where the paper is silent): 50 kg on the balance, 100 kW
heater (paper), coil 1000 K (paper), 10 m (paper) × 50 mm bundle,
emissivity 0.8, silica c_p ≈ 1000 J/(kg·K), 250 A DC feed, 50 µT.

| Force on the balance | Equivalent mass |
|---|---|
| Claim: 0.001 %–0.1 % lighter | 500 mg – 50 g |
| Standard prediction: heat adds E/c² (heavier) | +0.39 µg |
| Radiation recoil, all emission one way (P/c) | ≤ 34 mg |
| Feed lead in the geomagnetic field | 1.3 g per metre |
| Air buoyancy lost round the hot coil (in air) | 16.7 g |

The coil radiates 71 kW at 1000 K, so the heater holds the temperature.
Ordinary effects reach the claimed band. A bench must run in vacuum,
cancel lead forces (twisted, reversed feed), compare against a dummy coil
without fiber at the same power, cycle the drive on and off, and fix the
"coherence protocol" in writing before the run (the paper attributes
null results to poor phase-lock, §5.3, §9.1). A well-shielded balance
can then bound the effect below 10⁻⁶ of the rig mass.

Existing bound: MICROSCOPE finds no composition-dependent free fall at
the 10⁻¹⁵ level, η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10⁻¹⁵
(Touboul et al., *PRL* 129, 121102, 2022), which leaves little room for
gravity as an electromagnetic coupling.

## Part 4 — the paper's figures

- 1 kN on a 1000 t craft gives 0.001 m/s², not 1 m/s² (§4.3); the 5 and
  100 m/s² rows inherit the ×1000 error.
- 6000 g on 1000 t with 1 kN needs a residual mass fraction 1.7×10⁻⁸,
  not the ">99.99 %" (10⁻⁴) the paper states.
- Coil length: 3–10 km per coil (§4.1), 3,000–5,000 km per coil (§4.3),
  ~2,000 km in total (§6.3).
- Aluminium cladding (melts at 933 K) is run at 1,500–2,000 K; fused
  silica softens near 1,940 K (paper: melts at 2,400 K).
- Lock to the 7.8 Hz Schumann mode within 0.1 %: that mode has Q ≈ 5
  (Sentman 1995; `examples/earth_machine`), a linewidth of about 20 %.
- E_eff = (1 − f)·γ·m0·c² for v > c (§8.1): γ is not real above c,
  whatever the mass.

## Part 5 — Hall's "FTIR lensing" of galaxy images

Hall's photon theory explains lensed galaxy images as total (and
frustrated total) internal reflection at thin hydrogen clouds with an
empty gap between them, not as gravitational lensing. Total reflection
from hydrogen to vacuum happens only for rays within the critical grazing
angle √(2(n − 1)) of a sharp boundary, and turns them by at most twice
that. With H₂'s measured refractivity (Peck & Huang, *JOSA* 67, 1550,
1977; n − 1 = 1.388×10⁻⁴ at STP, 589 nm), scaled to cloud density:

| H₂ density | n − 1 | Largest deflection |
|---|---|---|
| 1 cm⁻³ (diffuse) | 5.2×10⁻²⁴ | 1.3×10⁻⁶ arcsec |
| 10⁴ cm⁻³ (molecular cloud) | 5.2×10⁻²⁰ | 1.3×10⁻⁴ arcsec |
| 10⁶ cm⁻³ (dense core) | 5.2×10⁻¹⁸ | 1.3×10⁻³ arcsec |

Observed strong-lensing images are separated by about 1–30 arcsec,
10³–10⁷ times more. Real cloud edges are graded over parsecs, not the
sub-wavelength step that total reflection needs, and frustrated total
reflection needs a gap no wider than about a wavelength. The reflection is
also chromatic: blue light is turned 2.1 % more than 1 µm light, while
gravitational lensing turns every wavelength equally. **Test:** compare
the image positions of a multiply imaged quasar across bands; positions
that agree across bands rule out the reflection reading.

## Related register entries

- H23 (`examples/photon_electron`): the Williamson–van der Mark toroidal
  electron the paper builds on — charge 0.910 e, anomaly 17 % low, not
  self-confined under linear Maxwell.
- H18 (`examples/pais_cavity`): inertial-mass reduction of a craft —
  thrust ≤ P/c, and the gyroscope weight loss it cites was not reproduced.
- §3c: the quaternion form of Maxwell's equations is equivalent to the
  vector form; its scalar part adds no physics.
- H8, H10: extended electrodynamics with a scalar field, and Bearden's
  scalar claims, which the paper cites (untested).
