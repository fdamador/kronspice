# Paramahamsa Tewari — space-vortex electron (hypothesis H13)

## Sources

**Read:** P. Tewari, "Structural relation between the vacuum space and the
electron," *Physics Essays* 31(1), 108–129 (2018),
DOI 10.4006/0836-1398-31.1.108. There is a free copy on tewari.org
(`uploads/3/9/2/2/39220475/16tewari.pdf`).

- The user attached the PDF (22 pp.), saved as
  `reference/tewari_2018_physics_essays.pdf`.
- Eqs. 9–16 and 29–35 were checked against the page images
  (pp. 114–115 and 118). The rest were read from the text layer, which
  drops minus signs in exponents.

**Read:** P. Tewari, *Universal Principles of Space and Matter: A Call
for Conceptual Reorientation*, Grain and Chaff, Toronto, 2017 (1st ed.
Crest, India, 2009), 257 pp., ISBN 978-0-9952058-0-2. The user attached
the PDF; saved as `reference/tewari_upsm_2017.pdf`. Page numbers below
are the printed ones.

**Still to obtain:**

- *Discovering Universal Reality* (2007), Natural Philosophy Alliance
  e-book. Its site blocks the web reader.
- Primary test reports of his "reactionless" generator (RLG). Only
  secondary accounts were found. A replica can be checked with the energy
  and momentum ledgers once a primary report gives dimensions, field,
  speed, and measured input and output.

## The model (paper §§II–V)

**Postulates.** Space is an eternal, incompressible, massless and
chargeless fluid whose flow speed is limited to c (§II).

**The electron.** An irrotational vortex, with u·r constant, breaks down
into a spherical void of radius r_e, where the flow reaches c (eqs. 1–3).
Integrating over the void's surface gives the following, in the paper's
own units (charge in cm³/s, mass in cm⁴/s):

| Quantity | Formula | Eq. |
|---|---|---|
| charge | q_e = π² r_e² c | 4 |
| mass | m_e = (4π/3) r_e³ c ("void volume × c") | 6–7 |
| energy | E = (4/5) m_e c² | 14 |
| angular momentum | L = (4/5) m_e c r_e | 15 |
| magnetic moment | μ = (3/4) q_e c r_e | 16 |
| permittivity | ε₀ = π/2c | 24 |
| permeability | μ₀ = 2/πc | 25 |

- Coulomb's law comes out with a coefficient of 0.73 instead of 1, which
  the paper calls negligible (§VI.A).
- The paper extends the model to the neutron and proton, hydrogen (an
  electron orbit of about 24 r_e at 0.69 c, eq. 35), gravity, light and
  the photoelectric effect, Planck's constant, solar-system and galactic
  "space vortices", and the size of the universe.

## The electron test bench (`examples/space_vortex`)

The same bench as Kelvin's rings (`doc/kelvin-references.md`). Measured
values are CODATA 2018 and NIST. r_e is fixed from the measured charge by
eq. 4, as in the paper: r_e = 4.03 × 10⁻¹¹ cm.

| Quantity | Model | Measured | Model / measured |
|---|---|---|---|
| charge (eq. 4) | 4.803e-10 statC | 4.803e-10 | 1 — **input**: r_e is fitted to e |
| radius | 4.03e-11 cm | none measured; scattering limits ≲ 10⁻¹⁸ cm | ≈ 4 × 10⁷ above the limit |
| angular momentum vs spin ħ/2 | 8.80e-28 erg s | 5.27e-28 | 1.67 |
| magnetic moment | 1.45e-20 erg/G | 9.28e-21 (1.00116 μ_B) | 1.56 |
| g-factor | **15/8 = 1.875** | 2.00232 | 0.936 |
| proton charge (eq. 30) | 12 e | e (atoms neutral to ~10⁻²¹ e) | 12 |
| hydrogen orbit radius vs Bohr radius | 9.67e-10 cm | 5.29e-9 | 0.18 |
| hydrogen orbit speed vs ground state αc | 0.69 c | 0.0073 c | 95 |
| radiation of that orbit vs Lyman-α | 3.4e18 Hz (X-ray) | 2.47e15 Hz (UV) | 1381 |

Notes on the table:

- **Mass is an input too.** The paper defines its gram from the measured
  electron mass (eqs. 10–12), so m_e is not predicted.
- **The g-factor needs no fitted value.** It is 2 m c μ/(e L) from
  eqs. 15–16, and r_e, q and m all cancel. So 15/8 is a genuine
  prediction of the model. It is 6.4 % from the measured
  g = 2.0023193 (known to about 10⁻¹³), which rules out eqs. 15–16 as
  written.
- **How eq. 16 is read.** As printed, μ = (3/4) q c r_e has the dimensions
  of an SI moment (current × area). The Gaussian moment of the same loop
  is divided by c. The g-factor is the same in either reading.
- **The hydrogen orbit gives the wrong light.** A charge orbiting at
  0.69 c on 24 r_e radiates at about 3 × 10¹⁸ Hz. Hydrogen's lines are
  at about 10¹⁵ Hz, and its ground state has v = αc.
- **The paper does not treat the Zeeman effect.** The model's prediction
  for the splitting and polarisation of a line in a field B is needed
  before the main bench test can be run. Measured: Δω = eB/2m_e c per
  unit of m_l, with σ± lines circularly polarised.

## Claims that do not discriminate

These parts of the paper agree with measurement, but they also agree with
standard physics for a reason built into the method:

- **Solar and planetary gravity (§XIII).** A velocity law v ∝ 1/√r fitted
  to one orbit (Mercury, the Moon) is Kepler's third law, so v²r = GM. The
  "surface acceleration" v²/R at the body's radius is then GM/R², which is
  Newton's surface gravity by construction. For the Sun:
  (47.9 km/s)² × 57.9 × 10⁶ km / R_sun² = 274 m/s².
- **Solar wind.** The "solar wind speed" of 436.7 km/s is the orbital speed
  √(GM/R) at the Sun's surface (436.8 km/s). Its closeness to the mean
  solar wind speed (≈ 440 km/s at 1 AU) is a coincidence. The wind
  accelerates over many solar radii, and its speed ranges from about 300
  to 800 km/s.

## Checked on the page

**Gram conversion.** Eqs. 10–12 (p. 114) give 1 g = 8.4, 8.8 and
8.6 × 10⁶ cm⁴/s. The exponent is +6; the text layer had lost the sign.
The test reproduces eq. 11.

**Proton charge.** Eq. 30 (p. 118) gives the proton charge as
q_p = 12 π² r_e² c = 12 q_e.

- The paper says so explicitly: "which is 12 times the electron charge".
- It makes hydrogen neutral by cancelling magnetic moments instead
  (eqs. 32–35). That is where the orbit speed of 0.69 c comes from.
- Measured, the proton's charge equals the electron's in magnitude, and
  atoms are neutral to about 10⁻²¹ e (gas-neutrality experiments, e.g.
  Bressi et al., Phys. Rev. A 83, 052101, 2011). A hydrogen atom with a
  net charge of 11 e would be ionised in every measurement.
- This is a second disagreement with measurement, as clear-cut as the
  g-factor.

**Magnetic moment (p. 115).** The paper integrates
dμ = dq·(ωr_e sin θ)·(r_e sin θ), that is, charge × speed × radius.

- The moment of a ring of charge dq moving at speed v on radius ρ is
  ½·dq·v·ρ (current dq·v/2πρ times area πρ²).
- With that standard factor, eq. 16 halves to (3/8) q_e c r_e, and the
  g-factor becomes 15/16 = 0.94 instead of 15/8.
- Both values are far from the measured 2.0023. The bench reports both.
- The angular-momentum integral, eq. 15 (dm × v × ρ), is correct as
  written.

**Orbit period.** The orbit period on p. 118 is 3 × 10⁻¹⁹ s, as the
paper's own numbers give.

## The book: what it adds

Chapters 2 and 5 repeat the electron formulas of the 2018 paper:

- the magnetic moment is eq. 2.16 (p. 54), with the same integrand
  dq·v·ρ;
- the spin is eq. 2.15.

The book also states two further predictions of its own.

**g-factor (§5.2, p. 131, eq. 5.4).** Tewari computes g = 1.875 himself
and writes: "Any other special meaning of the no. 2 for this constant is
most unlikely." The measured value is 2.00231930436 (CODATA 2018), known
to about 13 digits.

**Fine-structure constant (§5.1, p. 130, eq. 5.2).** In the model's own
constants, α = e²/(2ε₀hc) becomes (15/16)π² ≈ 9.25. The book says 1/137
"could be a cumulative experimental error". The measured value is
1/137.035999…, known to about 1.5 × 10⁻¹⁰. That is a factor of 1268, not
an error in the last digit. The bench computes 15π²/16 directly from
eqs. 4, 6 and 15.

**Density limit (§2.8, p. 54).** The book states that matter "cannot
possess density higher than 3.42 kg per c.c.", i.e. 3.42 × 10³ g/cm³.
White dwarfs, whose masses and radii are measured, average about
10⁶ g/cm³; Sirius B, at about 1.02 M☉ and 0.008 R☉, comes to roughly
2 × 10⁶ g/cm³.

## The Space Power Generator (§6.7, pp. 149–161; `examples/space_power_generator`)

**The machine.** The SPG is a homopolar machine (Fig. 6-10, p. 155):

- an electromagnet on a rotating iron core;
- a static magnetic yoke with three air gaps;
- output slip rings;
- a conductor whose U-shaped part sits in a radial hole, so that the
  fields of its two legs (C₁, C₂) cancel.

The book derives its EMF as V = ½ωBR² (eq. 6.16), the Faraday-disk
formula.

**The claim.** Most of the armature reaction is said to be avoided by
routing part of the conductor (L) close to the shaft, parallel to B.
Footnote 52 (p. 149) credits Bruce DePalma, who wrote to Tewari in 1978
about voltage induced in a co-rotating magnet and conductor, and Adam
Trombly's closed-path homopolar generator.

**The test (pp. 159–160).** The drive motor runs at 2700 r/min.

| Condition | Measured |
|---|---|
| No load | drive-motor input W₀ = 2943 W; excitation 650 W; EMF E = 3.41 V dc |
| Loaded, same speed | I = 3438 A; drive-motor input W_L = 8620 W |

From these the book computes:

- W_r = W_L − W₀ = 5677 W (eq. 6.20);
- W_E = E·I = 11 724 W (eq. 6.21);
- incremental efficiency 2.066 (eq. 6.22), or 2.49 for the SPG alone
  after dividing by the drive motor's efficiency of 0.83 (eq. 6.23).

The book calls this "a clear case of violation of the principle of
conservation of energy" (p. 161). It says a self-running system "does
exist" as a possibility if the output is fed back, but reports no such
test.

**The ledger.**

1. **Arithmetic.** It checks: 3.41 × 3438 = 11 724 W. Eq. 6.22 prints
   11 727 W, a typo.
2. **All inputs.** Drive motor plus excitation is 9270 W, against the
   claimed 11 724 W, a ratio of 1.26. The 650 W of excitation is left out
   of the book's ratio.
3. **The deciding assumption.** W_E multiplies the load current by the
   EMF measured at no load. Standard electrodynamics fixes the loaded EMF
   from the drive power: the power drawn from the shaft must equal
   E_load·I. The rotor's back torque is k·I with E = kω, wherever the
   reaction lands, here on the static yoke and the static conductor A–A₁.
   - The book's own numbers give 0.83 × 5677 W = 4712 W of extra shaft
     power, so E_load ≤ 4712 / 3438 = 1.37 V.
   - Even with a lossless drive motor, E_load ≤ 1.65 V.
   - So the standard prediction is that the EMF falls by at least 60 %
     at 3438 A. That is plausible: the load loop carries 3438 ampere-turns
     against the excitation, and the iron saturates.
   - The book reports no loaded-EMF measurement. It notes that the voltage
     drops "can also be" measured, but gives no numbers.

**The decisive bench measurement** is either of these at 3438 A:

- the sum of the voltage drops around the loaded loop (brushes, rings,
  internal conductor, shunt, load);
- the shaft torque, read with a torque transducer rather than inferred
  from the drive motor's input and an assumed efficiency.

| | Loaded EMF | Back torque |
|---|---|---|
| Standard prediction | ≤ 1.37 V | 16.7 N·m (as the drive motor delivered) |
| The claim | 3.41 V | 16.7 N·m, where 3.41 V at 3438 A requires k·I = 41.5 N·m |

A closed-loop, self-running test settles it outright.

## Status

H13 is *untested* until the model gives a Zeeman prediction. These
disagreements with measurement are recorded, and a revised model must
remove all of them:

- the g-factor: 15/8, which the book affirms (eq. 5.4), or 15/16 with the
  standard ring moment, against the measured 2.0023;
- the proton charge (eq. 30): 12 e, against e;
- the fine-structure constant (eq. 5.2): 15π²/16, against 1/137.036;
- the density limit (§2.8): 3.42 × 10³ g/cm³, exceeded by white dwarfs.

The SPG claim stands or falls on the loaded EMF (≤ 1.37 V predicted,
3.41 V assumed). It is a bench measurement still to be made.
