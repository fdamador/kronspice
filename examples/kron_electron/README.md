# A ring of charge as one of Kron's machines (H21)

Hypothesis H21 (CLAUDE.md §1): Kron's machine equations model the ring
electron, with voltage and current as the pressure and flow of energy at
that scale. This example puts it on the electron test bench of
[doc/kelvin-references.md](../../doc/kelvin-references.md), next to
Kelvin's 1899 ring ([kelvin_ring](../kelvin_ring), Lorentz force),
Tewari's space-vortex electron ([space_vortex](../space_vortex), H13) and
the toroidal-photon electron ([photon_electron](../photon_electron), H23).
The frames are those of [kron_machine](../kron_machine) (milestone
M9c.1), and the measurement that discriminates the result is magnetic
resonance, solved in the same frame in [kron_nmr](../kron_nmr).

```
go run ./examples/kron_electron -png viz/ring_field.png
go run ./examples/kron_electron -gif viz/ring_motion.gif
go test ./examples/kron_electron
KRONSPICE_EXPORT=1 go test ./examples/kron_electron -run TestExport   # netlists, KiCad schematics, results sheets in export/
```

## The idea, in Kron's terms

For a ring of charge q and mass m turning at θ̇ the current is
i = q·θ̇/2π, so the electrical and the mechanical coordinate are the same
axis. The kinetic energy ½·m·a²·θ̇² = ½·(m·a²·(2π/q)²)·i² is then a magnetic
energy on that axis, and the inertia enters Kron's **metric tensor** as an
inductance m·a²·(2π/q)², beside the ring's own self-inductance L. Kron
defined the metric exactly that way, as "all the self- and mutual
inductances and the rotor inertia" (*Equivalent Circuits of Electric
Machinery*, 1951, Epilogue p. 254, eqs. 3–4), and his eq. 11 (p. 259) says
that with e = R = 0 the machine equation is the motion of a charged
particle in an electromagnetic field. Everything below uses the N machine
element ([device/machine.go](../../device/machine.go)) unchanged.

The ring quantities are the thin-ring formulas: self-inductance
L = μ₀·a·(ln(8a/r) − 7/4) (Maxwell, *Treatise* Art. 706), capacitance from
Smythe §5.17, the inertial inductance m·a²·(2π/q)² above, and the mutual
inductance between coil and ring from `em.MutualCircles` (Maxwell's
elliptic form). For the twin test the printed metric is: coil
6.201×10⁻⁷ H, coil–ring mutual 3.494×10⁻⁸ H, tube and fluid inertial
inductance 9.87×10⁸ H each, and the ring's own L of 3.1×10⁻⁷ H dropped as
in Kelvin (or kept, which changes nothing at this scale).

## 1. Twin test: Kelvin's ring in Kron's frame (standard)

Kelvin's ring (tube +e, fluid −e, each 10 g, radius 5 cm; coil of 20
kA-turns, radius 10 cm, switched on over 1 s, 5 cm below) becomes a
three-winding machine: coil, tube, fluid, with Maxwell's mutual inductance
of coaxial circles between coil and ring and the two inertias on the charge
axes. Tube and fluid windings are closed through 0 V sources
(frictionless), so their flux linkages are conserved.

| Quantity | Kron machine | Lorentz-force bench |
|---|---|---|
| Flux through the ring | 6.98732×10⁻⁴ Wb (Maxwell elliptic) | 6.98736×10⁻⁴ Wb (Biot–Savart) |
| Moment of momentum, tube | −1.11207×10⁻¹⁰ kg m²/s | −1.11207×10⁻¹⁰ |
| Moment of momentum, fluid | +1.11207×10⁻¹⁰ | +1.11207×10⁻¹⁰ |
| Kelvin's eΦ/2π | 1.11207×10⁻¹⁰ | |
| Fluid angular velocity | 4.44827×10⁻⁶ rad/s | 4.44829×10⁻⁶ |
| Larmor eB̄/2w | 4.44827×10⁻⁶ | |
| Kinetic energy | 4.94677×10⁻¹⁶ J | 4.94683×10⁻¹⁶ J |
| Source work against the ring | = kinetic + mutual-field energy (−9.89×10⁻¹⁶ J), residual 1.5×10⁻¹³ | residual 8×10⁻¹⁶ |

The 1.1×10⁻⁵ difference in kinetic energy is the bench's Biot–Savart flux
against Maxwell's elliptic value used here. Adding the ring's own
self-inductance (3×10⁻⁷ H against 10⁹ H of inertial inductance) changes
nothing at this scale. The engine's ledger closes too (source 124.028 J,
machine's ½iᵀLi 124.02 J, mismatch 7.8×10⁻³ J over the ramp's
breakpoints), but it is dominated by the coil's own 124 J of field energy,
which is why the bench's own bookkeeping is reported.

## 2. Torque from the metric (standard)

The same ring, spinning at 1000 rad/s (current 0.159 mA, moment
1.25×10⁻⁶ A·m²), tilts about a diameter in 1 T: frame I with
M(θ) = B·π·a²·cos θ between a field winding and the ring, the tilt inertia
½·m·a² as a capacitor on the mechanical port. Kron's torque
½·iᵀ·(∂L/∂θ)·i is −μ·B·sin θ, and the oscillation released at 0.05 rad
has period 19.873 s against the pendulum value 4K(k)/ω₀ = 19.872 s
(2×10⁻⁵). The ring's current drifts by 6×10⁻¹¹ (flux conservation) and
the ledger mismatch is 6×10⁻¹⁰ J.

### Harmonics

The ring current itself has none: it is a steady current of a uniform
ring of charge, so its field is static. The one motion in the model that
generates harmonics is the tilt, a nonlinear pendulum J·θ̈ = −μB·sin θ,
and the solver's ω(t) carries exactly the pendulum's odd series (Jacobi
nome q = e^(−πK′/K), k = sin θ₀/2), measured by leak-free Fourier sums
over two periods:

| θ₀ | fundamental | ω₃/ω₁ run | exact | ω₅/ω₁ run | exact | ω₂/ω₁ |
|---|---|---|---|---|---|---|
| 0.05 rad | 0.05032 Hz | 3.906×10⁻⁵ | 3.907×10⁻⁵ | 2×10⁻⁸ | 1.5×10⁻⁹ | 2×10⁻⁸ |
| 1 rad | 0.0472 Hz | 0.01658 | 0.01658 | 2.70×10⁻⁴ | 2.71×10⁻⁴ | 2×10⁻⁸ |
| 2.5 rad | 0.03063 Hz | 0.1597 | 0.1598 | 0.0225 | 0.0225 | 8×10⁻⁷ |

Even harmonics are absent by the symmetry θ → −θ. The other harmonic
families are modal, not driven: the ring as a closed transmission line
resonates at every integer number of wavelengths round it, f_n = n·c/2πa
(for the Compton ring, n·m_e c²/h), and a ring whose charge was lumpy
rather than uniform would radiate at multiples of the rotation frequency,
which this uniform ring does not.

## 3. The electron as the ring (hypothesis)

Charge e; the total angular momentum of the metric is set to ħ/2 and the
rest follows (`RingModel.AtSpin`, core radius a/100).

| Ring | a [m] | v/c | μ/μ_B | g | KE/m_e c² | ½LI²/m_e c² | W_es/m_e c² |
|---|---|---|---|---|---|---|---|
| Compton ring, electron mass circulating | 3.86×10⁻¹³ | 0.494 | 0.494 | 0.989 | 0.122 | 0.0014 | 0.0078 |
| Half-Compton ring, mass circulating | 1.93×10⁻¹³ | 0.978 | 0.489 | 0.978 | 0.478 | 0.011 | 0.0155 |
| Compton ring, massless (field inertia) | 3.86×10⁻¹³ | 43.6 | 43.6 | 87.2 | 0 | 10.9 | 0.0078 |
| Classical-radius ring, massless | 2.82×10⁻¹⁵ | 43.6 | 0.318 | 0.637 | 0 | 1490 | 1.06 |
| Massless ring fitted to the measured moment | 8.86×10⁻¹⁵ | 43.6 | 1.001 | 2.002 | 0 | 475 | 0.338 |
| **Measured (Dirac + QED)** | < 10⁻¹⁸ (scattering) | | 1.00116 | 2.00232 | | | |

What the bench says:

- A ring whose mass circulates with its charge has
  g = 1/(1 + L·e²/4π²·m·a²) ≤ 1: its precession in a field is Larmor's
  eB/2m, half the electron's spin precession.
- Only a massless ring (inertia from its own field) reaches the measured
  moment at S = ħ/2, and it does so at 8.9 fm with v = 43.6 c and a
  magnetic energy of 475 m_e c². Every massless thin ring with S = ħ/2 moves
  at v/c = π/(2α(ln(8a/r) − 7/4)), fixed by the fine-structure constant;
  at v ≤ c (a ≥ ħ/2m_e c) the massless ring's moment is far below μ_B.
- Every ring radius lies far above the 10⁻¹⁸ m scattering bound.

### Natural frequencies of the rings (S = ħ/2, B = 1 T)

| Ring | Rotation of the charge | Precession (ESR) | Tilt libration | LC self-resonance | Wave once round at c |
|---|---|---|---|---|---|
| Compton ring, mass circulating | 6.1×10¹⁹ Hz | 13.8 GHz | 1.3×10¹⁵ Hz | 2.3×10¹⁹ Hz | 1.24×10²⁰ Hz |
| Half-Compton ring, mass circulating | 2.4×10²⁰ | 13.7 GHz | 2.6×10¹⁵ | 4.6×10¹⁹ | 2.47×10²⁰ |
| Compton ring, massless | 5.4×10²¹ | 1.22 THz | — | 2.3×10¹⁹ | 1.24×10²⁰ |
| Classical-radius ring, massless | 7.4×10²³ | 8.9 GHz | — | 3.1×10²¹ | 1.7×10²² |
| Massless ring fitted to the moment | 2.4×10²³ | 28.0 GHz | — | 1.0×10²¹ | 5.4×10²¹ |
| **Measured / reference** | | **28.025 GHz** (ESR) | | | **1.236×10²⁰ Hz** (m_e c²/h) |

- **Precession** is the spin's magnetic-resonance frequency, g·eB/4πm; it
  is what [kron_nmr](../kron_nmr) measures, and only the fitted massless
  ring reaches the measured value, by construction.
- **Rotation** is the charge's turning rate, the current's fundamental:
  about half the Compton frequency for the Compton ring (v ≈ c/2).
- **Tilt libration** is the part-2 pendulum, √(μB/J)/2π: optical
  (10¹⁵ Hz) for an electron-mass ring in 1 T; a spinning ring shows it as
  nutation.
- **LC self-resonance** uses the thin-ring inductance and capacitance; the
  wave once round the ring at c is the ring's transmission-line
  fundamental. For the Compton ring that is the Compton frequency m_e c²/h
  by the identity a = ħ/m_e c, and for the half-Compton ring it is 2m_e c²/h,
  the Zitterbewegung frequency. These are identities of the chosen radius,
  not predictions.
- No radial breathing mode: the model has no restoring force for the
  ring's radius (the electrostatic hoop stress is unbalanced, see step 3),
  so that frequency needs a declared stiffness.

Kelvin's 5 cm ring of part 2 (1 µC, 10 g, 1000 rad/s, 1 T): rotation
159 Hz, precession eB/4πm = 8×10⁻⁶ Hz, libration 0.0503 Hz (the 19.9 s
period), LC self-resonance 177 MHz.

So the status of H21 is *bounded*: the ring-as-machine reproduces the
classical ring exactly and inherits its obstacles (g ≤ 1 or v > c). The
Zeeman polarisation test stays with the 3-D bench of `kelvin_ring`; the
one-angle N element cannot precess. The standard prediction beside every
row is the Dirac electron with QED (g = 2.002 319 304 36,
μ = 1.001 16 μ_B); H13 on the same bench gives g = 15/8.

### Against existing measurements

| Quantity in the model | Measured value | Match? |
|---|---|---|
| Precession of a ring whose mass circulates with its charge, eB/4πm = 14.0 GHz/T (g = 1) | The normal Zeeman effect: orbital moments precess at Larmor's eB/4πm, g_L = 1 (Zeeman 1897; Larmor 1897; the Cd 643.8 nm singlet triplet in any teaching lab) | **Yes.** The ring is an orbital moment, and the measured orbital g-factor is 1. It does not match the electron's spin (g = 2.0023). |
| Precession of the fitted massless ring, 28.02 GHz/T | Electron spin resonance; g/2 = 1.001 159 652 180 59(13) (Fan et al. 2023; Hanneke et al. 2008) | By construction only: the ring was fitted to the measured moment, and needs v = 43.6 c. |
| Pendulum period 4K(k)/ω₀ and the odd-harmonic series of the tilt | Large-amplitude pendulum measurements agree with the elliptic-integral period to about 10⁻⁴ (Nelson & Olsson 1986) | **Yes**, for the classical pendulum; the model's tilt is that pendulum. |
| Wave once round the Compton ring, m_e c²/h = 1.236×10²⁰ Hz, and 2m_e c²/h for the half-Compton ring | An identity of the chosen radius. The Compton frequency itself is a CODATA constant, not a resonance observed in an electron. A channeling experiment reported a dip consistent with the de Broglie internal clock m_e c²/h (Gouanère et al. 2008); it has not been independently reproduced. Zitterbewegung at 2m_e c²/h has been simulated with a trapped ion obeying the Dirac equation (Gerritsma et al. 2010), not observed in an electron. | No measurement on the electron. |
| Charge rotation frequency, tilt libration (10¹⁵ Hz), LC self-resonance, the ring's modal series n·c/2πa | None exist: a uniform ring current radiates nothing at these frequencies, so there is nothing to detect | No measurement. |

So the two genuine matches are classical: the ring's precession is the
orbital Larmor frequency of the normal Zeeman effect, and its tilt is the
measured large-amplitude pendulum. Neither bears on the electron's spin,
which is exactly what the g ≤ 1 bound says.

## 4. The ring's field

`viz/ring_field.png`: field lines of the ring current from Biot–Savart, in
the meridional plane (left: the wire as two red dots, the axis vertical)
and in an oblique view at eight azimuths (right). Every line lies in a
meridional plane and threads the ring once: the far field is a dipole,
the near field circles the wire (`TestFieldLinesClose` checks that lines
started at 0.3, 0.6 and 0.88 of the radius return to their start to
10⁻⁹). Pictures of toroidal electron models with lines that wind around a
sphere need a current with a component along the ring's axis as well, a
second winding in Kron's frame; a single ring current cannot produce them.

## 5. Animation

`viz/ring_motion.gif` (90 frames, standard-library `image/gif`) shows the
two transients as they progress, each panel over its own run with a
progress bar along its bottom edge:

1. **Kelvin's switch-on from above (0–2 s).** The tube's marker (red) and
   the fluid's (blue) turn in opposite senses as the coil current (green
   bar) ramps up over the first second and then holds. The real angles are
   microradians, so they are scaled so that the fluid ends at half a turn.
   The red and blue bars are the moments of momentum, equal and opposite,
   growing with the flux and constant once it holds.
2. **The tilting ring, edge-on (three periods).** The ring (black bar,
   wire ends red) swings about its diameter in the vertical field (grey
   arrows); the red arrow is its magnetic moment, along the ring's normal.
   The 0.05 rad swing is drawn ×8.
3. **Its field.** The ring's dipole field lines tilt with it, so the
   whole pattern rocks at the pendulum period.

The animation is drawn from the solver's own time points (currents and
θ(t) integrated from ω(t)), not from an analytic motion.

## Tests

| Test | Check | Held to |
|---|---|---|
| `TestKelvinTwin` | moments of momentum ∓eΦ/2π, fluid at Larmor's rate, bench ledger closed, with and without the ring's self-inductance; against the Lorentz-force bench's 1.11207×10⁻¹⁰ kg m²/s, 4.44829×10⁻⁶ rad/s and 4.94683×10⁻¹⁶ J | 1e-6; bench 1e-5, kinetic energy 2e-5 |
| `TestTiltPeriod` | period against 4K(k)/ω₀; ring current conserved | 2e-4; drift 1e-9 |
| `TestTiltHarmonics` | 3rd and 5th harmonic ratios against the Jacobi series at θ₀ = 1 rad; no even harmonics; small-angle 3rd = θ₀²/64 | 2e-4; 1e-3 |
| `TestElectronRing` | g·(1 + Le²/4π²ma²) = 1; precession = g × Larmor; Compton ring below c; fitted ring carries the measured moment and is superluminal; every scenario above the size bound | 1e-12; 1e-9 |
| `TestMasslessSpeed` | v/c = π/(2α(ln(8a/r) − 7/4)) at a = 1 fm, ƛ_C and 0.1 nm | 1e-9 |
| `TestFrequencies` | Compton ring's c/2πa = m_e c²/h; precession = g·eB/4πm; fitted ring at the measured ESR; Kelvin ring libration period 19.869 s | 1e-9; 1e-12; 1e-4; 0.01 s |
| `TestFieldLinesClose` | field lines of the ring current close | 1e-9 |
| `TestExport` | Kelvin's ring (`.tran 5m 2 uic`) and the tilting ring (`.tran 50m 60 uic`) exported with the `N` element | |

## Files

| File | What it holds |
|---|---|
| `main.go` | thin-ring L and C, inertial inductance, `KelvinSpec`/`KelvinCircuit`/`RunKelvin`, `TiltSpec`/`TiltCircuit`/`RunTilt`, `TiltHarmonics` and `PendulumHarmonics`, `RingModel` with `AtSpin`, `Frequencies`, `FitMassless` and `Scenarios`, `FieldLine`, the PNG and GIF drawing, and the five-part report |
| `main_test.go` | the tests above |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |
| `viz/ring_field.png`, `viz/ring_motion.gif` | the field-line picture and the animation |

## Sources

- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley 1951,
  Epilogue pp. 253–265 (`reference/Equivalent Circuits of Electric
  Machinery Kron.pdf`): eqs. 3–4 (the metric), eq. 9 (frame II(b)),
  eq. 11 (charged particle in a field).
- W. Thomson (Kelvin), "Magnetism and Molecular Rotation," 1899, *Math.
  Phys. Papers* V, Art. 175 (`reference/mathematicalphys05kelvuoft.pdf`);
  the bench in [doc/kelvin-references.md](../../doc/kelvin-references.md).
- J. C. Maxwell, *Treatise*, Art. 706 (ring inductance) and the mutual
  inductance of coaxial circles (`em/inductance.go`); W. R. Smythe,
  *Static and Dynamic Electricity*, §5.17 (ring capacitance).
- CODATA 2018 (E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor, *Rev.
  Mod. Phys.* 93, 025010, 2021) for e, m_e, ħ, μ_B, g, α; D. Bourilkov,
  *Phys. Rev. D* 64, 071701 (2001) for the size bound.
- P. Zeeman, "On the influence of magnetism on the nature of the light
  emitted by a substance," *Phil. Mag.* 43, 226–239, 1897; J. Larmor, "On
  the theory of the magnetic influence on spectra," *Phil. Mag.* 44,
  503–512, 1897 — the normal Zeeman effect and Larmor's theorem, the
  orbital precession eB/4πm the mass ring reproduces.
- D. Hanneke, S. Fogwell, G. Gabrielse, "New measurement of the electron
  magnetic moment and the fine structure constant," *Phys. Rev. Lett.*
  100, 120801, 2008; X. Fan, T. G. Myers, B. A. D. Sukra, G. Gabrielse,
  "Measurement of the electron magnetic moment," *Phys. Rev. Lett.* 130,
  071801, 2023 — the spin g-factor the ring does not reach.
- R. A. Nelson, M. G. Olsson, "The pendulum — rich physics from a simple
  system," *Am. J. Phys.* 54, 112–121, 1986 — measured amplitude
  dependence of the pendulum period against the elliptic-integral
  formula, the basis of the tilt's period and harmonics.
- M. Gouanère, M. Spighel, N. Cue, M. J. Gaillard, R. Genre, R. Kirsch,
  J. C. Poizat, J. Remillieux, P. Catillon, "A search for the de Broglie
  particle internal clock by means of electron channeling," *Found. Phys.*
  38, 659–664, 2008 — a reported but unreproduced signal at m_e c²/h.
- R. Gerritsma et al., "Quantum simulation of the Dirac equation," *Nature*
  463, 68–71, 2010 — Zitterbewegung at 2mc²/h in a trapped-ion simulation,
  not in an electron.
- P. Tewari, *Physics Essays* 31, 108, 2018 — the space-vortex electron
  (H13) on the same bench, `examples/space_vortex`.
