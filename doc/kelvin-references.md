# William Thomson (Lord Kelvin) — references added October 2026

Two volumes of Kelvin's *Mathematical and Physical Papers* (Cambridge
University Press) are in `reference/`. The user added them for Kelvin's
view of a charged particle as a ring-shaped (toroidal) structure of
electromagnetic energy, which is close to the project's working premise
(H1: toroidal fields carried by pulses; §3g, `doc/em-push-pull.md`).

| File | Volume | What it holds for kronspice |
|---|---|---|
| `mathematicaland01kelvgoog_text.pdf` (657 pp., Google scan) | Vol. I, 1882 (papers 1841–1853) | Art. XXVII, "On a Mechanical Representation of Electric, Magnetic, and Galvanic Forces" (*Camb. Dublin Math. J.* 2, 1847): electric and magnetic forces as displacements and rotations of an elastic solid — the start of Maxwell's mechanical models. Art. XLIV, potential of a closed galvanic circuit; Art. LII, magnetic induction in crystalline substances. Footnote on lines of magnetic force round an axis treated as vortex lines (text-layer line ≈ 2726). |
| `mathematicalphys05kelvuoft.pdf` (630 pp.) | Vol. V, 1911 (ed. Larmor): thermodynamics, cosmical physics, molecular theory, electrodynamics | Art. 175, "Magnetism and Molecular Rotation" (*Proc. R. Soc. Edinb.* 22, 631–635; *Phil. Mag.* 48, 236–239, 1899), pp. 354–358. "Æpinus Atomized" (1902), the electrion model of the atom. A remark on vortex atoms in the kinetic-theory papers. |

## Kelvin's ring model (Vol. V, Art. 175, pp. 354–358)

The model is concrete enough to simulate:

1. A ring of area A made of a rigid, non-conducting tube carrying a charge
   (vitreous, +e), filled with an incompressible fluid carrying −e.
2. When a magnetic field M normal to the ring is switched on, the induced
   EMF acts oppositely on the two charges. Tube and fluid turn in
   opposite directions, each taking a moment of momentum `eAM/2π`. That
   is the editor's corrected value: the 1899 text has `eAMr`, and a
   bracketed note on p. 354 corrects it.
3. If the tube is held by the ether, the tube's share leaves as a decaying
   train of transverse waves (`e^(−pt) sin qt`). The fluid keeps
   circulating for as long as the field is held: a persistent ring current
   set up by the field (p. 355).
4. Attached to the ether by springs, the ring is a radiating molecule with
   six vibration modes. Kelvin used it to explain Faraday's magneto-optic
   rotation by "inertia and pressure".

**Kelvin's own verdict (p. 357).** After Zeeman's discovery (1896) he
gave the model up: "I now … discard my old tempting gyrostatic
hypothesis for an irrefragable reason." Molecules carrying such rings,
oriented at random, broaden each spectral line in a magnetic field
instead of splitting it into Zeeman's sharp, circularly polarised pair or
triplet. He quotes Larmor's statement of the same objection. The editor
adds that "the effects are of course excessively minute" (p. 354).

Two points to keep in view:

- **The ring is a molecule, not the electron.** Kelvin's own electron
  ("electrion") model, in "Æpinus Atomized", is point-like. The toroidal
  *electron* belongs to later authors: A. L. Parson's magneton (1915) and
  more recent spinning-charged-ring models. Those are not in `reference/`
  yet.
- **Vortex atoms are in Vol. IV** (see the section on Vol. IV below).
  Vol. V only mentions them in passing.

## Simulation: `examples/kelvin_ring`

The 1899 ring is built with standard electrodynamics only: every charge
feels the Lorentz force through `em.LorentzVector`. Run it with
`go run ./examples/kelvin_ring`; tests are in `main_test.go`.

### Part 1 — induction (Kelvin §§1–3)

Set-up: a ring of 5 cm radius, ±1 µC, with tube and fluid of 10 g each. The
field comes from a coil of 10 cm radius, 5 cm below the ring, ramped to
20 kA-turns over 1 s, held for 1 s and ramped off. Every rotation figure
below is in kg·m²/s.

**Flux through the ring.**

- Biot–Savart integration: 6.98736 × 10⁻⁴ Wb.
- Maxwell's elliptic-integral mutual inductance: 6.98732 × 10⁻⁴ Wb.

**Rotation given to each half.** Kelvin's corrected value is
eAM/2π = eΦ/2π.

| | Simulated | eΦ/2π |
|---|---|---|
| Tube | −1.11207 × 10⁻¹⁰ | 1.11207 × 10⁻¹⁰ |
| Fluid | +1.11207 × 10⁻¹⁰ | 1.11207 × 10⁻¹⁰ |

- The two halves cancel, so the net angular impulse on the ring is
  −5.5 × 10⁻³⁰ (round-off). The coil gets no net reaction.
- The fluid turns at 4.448 × 10⁻⁶ rad/s. This is exactly Larmor's rate
  eB̄/2w. **Kelvin's ring is a mechanical picture of Larmor's theorem.**
  For an electron at 1 T the same rate is 8.79 × 10¹⁰ rad/s.

**Magnetic response.** The induced currents in the tube and the fluid
have the same sign, so the ring is diamagnetic: −1.1 × 10⁻¹⁴ A·m². This
is the classical Langevin diamagnetism. The coil repels the ring with a
force of 1.4 × 10⁻¹⁴ N.

**Energy ledger.**

- The coil's source does −4.95 × 10⁻¹⁶ J of work against the ring's
  back-EMF.
- The ring gains +4.95 × 10⁻¹⁶ J of kinetic energy.
- The mutual-field energy M·I_coil·I_ring is −9.89 × 10⁻¹⁶ J.
- Source work = kinetic energy + mutual-field energy, to a residual of
  8 × 10⁻¹⁶.

Switching the field off returns both halves to rest, and the net work by
the source is zero. The editor's remark holds: "the effects are of course
excessively minute".

### Part 2 — Kelvin's Zeeman argument (§§5–7), reproduced and extended

Set-up, in dimensionless units:

- The ring is held by isotropic angular springs, giving natural frequency
  ω₀ = 1.
- The field gives the fluid a Larmor frequency ωL = eB/2w = 0.05.
- The ring is tilted to 12 orientations, from cos θ = 0.04 to 0.96.
- The frequency resolution is 0.0026.

Each ring is kicked and its spectrum taken. Three models:

| Model | Splitting vs. orientation | Ensemble line | Polarisation along B (side lines) |
|---|---|---|---|
| Lorentz bound charge (standard classical Zeeman) | 0.0513 at every θ | sharp (FWHM 0.0038 = resolution) | ∓1.00: opposite circular polarisations, as Zeeman observed; across B: centre 100 % ∥B, sides 100 % ⊥B |
| Kelvin's gyrostat **as he wrote it** (fixed spin S along the axis) | (ωL/2)·cos θ: 0 … 0.0243 | **broad** (FWHM 0.0188) | mixed |
| Kelvin's ring under the **full Lorentz force** (fluid free to slip, nothing linearised) | 0.0243 … 0.0260 at every θ | **sharp** (FWHM 0.0037) | motion ∓0.47 (analytic 5/11 = 0.45); radiated moment about 0 |

**Kelvin's gyrostat.** The second row reproduces his reason for giving up
the model (p. 357). With a fixed internal spin, the splitting goes as
cos θ. Randomly oriented molecules therefore blur each line instead of
splitting it cleanly.

**The full ring.** The third row is new.

- *Why it differs.* When the ring tilts, the flux through it changes, so
  the fluid's spin changes with it. The fluid's canonical angular momentum
  `w a² u − eΦ/2π` is conserved, not its spin. Kelvin's gyrostat leaves
  this induction out.
- *What it gives.* With it included, the gyroscopic vector becomes
  `G = (e a²/2)(B − 2(n·B)n)`, the mirror image of B in the ring's plane.
  Its length does not depend on orientation, so every ring splits by the
  same amount: e B/2(w + w′) (`GyroVector`, `TestGyroVector`). The
  simulated splitting agrees within the resolution. That makes the lines
  **sharp**, which removes the objection Kelvin and Larmor raised.
- *Why it still fails Zeeman's test.* The modes circle G, not B. So along
  the field their circular polarisation averages to about ±0.45 in the
  ring's motion. In the light radiated by the ring's own magnetic moment
  it averages to about 0, whereas Zeeman measured ±1. Across the field the
  centre and side lines are only partly linearly polarised, not 100 %.
- *The scale is wrong too.* The splitting is set by the total ring mass,
  e B/2(w + w′), not the electron's e B/2m.

So a ring of two opposite charges with different charge-to-mass ratios
cannot reproduce the normal Zeeman effect, even with its full dynamics.
Lorentz's single charge does: all of its charge has the same
charge-to-mass ratio, which is the condition of Larmor's theorem.

**What this means for toroidal-particle ideas.** Any toroidal model of the
electron has to give two things:

- a splitting set by e/m_e;
- the observed circular polarisations.

Part 2 provides the test bench for that. It also shows that the anomalous
Zeeman effect, the electron's g-factor and spin ½ are further
requirements beyond these classical ones.

## Vol. IV — vortex atoms (`examples/vortex_ring`)

The user supplied Vol. IV (1910, *Hydrodynamics and General Dynamics*,
579 pp.): Internet Archive item `in.ernet.dli.2015.212310`,
<https://archive.org/download/in.ernet.dli.2015.212310/2015.212310.Mathematical-And_text.pdf>.
Neither the workspace nor the user's computer could download it, because
the network proxy refuses archive.org. So it is **not yet in
`reference/`**; save it there as `kelvin_mpp_vol4.pdf`. The notes below
come from the item's OCR text, read through the web reader.

**Contents relevant here** (article, year, page):

- Art. 1, "On Vortex Atoms" (1867), p. 1.
- Art. 2, "On Vortex Motion" (1868), p. 13.
- Art. 3, "The Translatory Velocity of a Circular Vortex Ring" (1867),
  p. 67.
- Art. 10, "Vortex Statics" (1875), p. 116.
- Art. 11, "On the Precessional Motion of a Liquid" (1877), p. 129.
- Art. 15, "Vibrations of a Columnar Vortex" (1880), p. 162.
- Art. 16, "On the Stability of Steady and of Periodic Fluid Motion"
  (1887), p. 166.
- Art. 21, "On the Stability and Small Oscillation of a Perfect Liquid
  full of Nearly Straight Coreless Vortices" (1889), p. 202.
- Ether papers: Arts. 72–74 (1890) and 82 (1900).

**What "On Vortex Atoms" claims** (pp. 1–12):

- Atoms are closed vortex rings in a frictionless, incompressible
  liquid. Helmholtz's theorem makes such a ring indivisible:
  "a closed line of vortex core is literally indivisible by any action
  resulting from vortex motion" (p. 3).
- Knotted and linked rings could give the variety of the chemical
  elements. Sodium, for example, is proposed as "two approximately equal
  vortex rings passing through one another like two links of a chain"
  (p. 5). Its two D lines would then come from two similar fundamental
  vibration modes.
- Tait's smoke rings, which bounce off one another, show elastic
  collisions, as the kinetic theory of gases needs (pp. 2, 11–12).
- Electricity and magnetism appear only as an analogy. The liquid's
  velocity round a vortex ring follows the lines of force of a circular
  current, and the ring's impulse "agrees in direction and magnitude with
  the magnetic moment" of the corresponding ring electromagnet (pp. 9–10).

**The atoms are not electrons.** These are atoms of matter, not
electrons, which were discovered thirty years later. The paper's
predictions are mostly qualitative. Its main quantitative statements are:

- the rotation rate (i−1)/i of deformation waves on a straight vortex
  (p. 4);
- a vibration period of less than 10⁻¹² s, matching sodium light (p. 5);
- the impulse–magnetic-moment equality (pp. 9–10).

**What the engine does with it.** Through Kelvin's own analogy, a vortex
ring's flow is the magnetic field of a current loop: circulation Γ plays
the current I, and velocity plays B/μ₀. The ring's kinetic energy is
½ρΓ²L/μ₀, where L is the loop's self-inductance. So the M9i inductance
code computes the vortex ring.

1. **Ring speed.** Maxwell's mutual inductance of coaxial circles is
   averaged over the core cross-section, and the ring's speed is taken as
   U = dE/dP (P = ρΓπR², at fixed circulation and core volume). The result
   approaches Kelvin's U = Γ/4πR·(ln 8R/a − ¼):

   | Core radius a/R | Difference from Kelvin |
   |---|---|
   | 0.2 | −0.97 % |
   | 0.1 | −0.26 % |
   | 0.05 | −0.07 % |

   The difference shrinks as (a/R)², the order of the terms the
   thin-ring formula leaves out. The formula is quoted from Saffman,
   *Vortex Dynamics*, §10.2; check it against Vol. IV p. 67 once the PDF
   is in `reference/`.
2. **Impulse ↔ magnetic moment.** The flow at the ring's centre is Γ/2R,
   just as the field at the centre of a loop is μ₀I/2R, and the impulse is
   ρΓπR², the counterpart of the magnetic moment IπR² (Kelvin p. 10).
3. **Leapfrogging** (Helmholtz 1858; Tait's rings). Two coaxial rings, 0.6R
   apart, pass through each other every 2.92 time units: 21 passes in 60.
   Their radii swing between 0.73R and 1.21R, and the pair advances 1.5
   times faster than a single ring. Impulse and energy are conserved to
   2 × 10⁻⁵. The energy here is the self terms plus the mutual-inductance
   interaction, so this also checks that the Biot–Savart velocities and
   the elliptic mutual inductance agree.

**Where the analogy stops.** The flow and the magnetic field have the
same shape, but different dynamics:

- A vortex ring travels because its vorticity is carried along by the
  flow it makes.
- A current loop's own magnetic field does not carry the loop. In
  Maxwell's equations the "ring electromagnet" stays where it is.

A self-propelled toroidal *electromagnetic* structure must therefore come
from the full field equations, not from the static analogy. The standard
example is Hellwarth and Nouchi's "flying doughnut", which travels at c.
It is already a validation case in §3g and is linked to H1.

**For the toroidal-electron premise.** A vortex ring in a neutral liquid
has no charge and no coupling to an external field. So the Zeeman test
bench of `examples/kelvin_ring` cannot be applied until a
charge–vortex coupling is specified. The model would then have to give the
electron's splitting e B/2m_e and Zeeman's polarisations, as set out
above.

## Not yet done

- **Radiated wave train.** Kelvin's tube, held by the ether, sheds its
  spin as a decaying wave train (p. 355). This would be a
  radiation-resistance edge (§3a.1, item 5), giving his `e^(−pt) sin qt`
  decay. The tube in Part 1 is free, and the springs in Part 2 are
  lossless.
- **Unequal masses and anisotropic springs.** Part 2 uses w = w′ and
  isotropic springs, the case most favourable to Kelvin. The general case
  only adds mode mixing.

Any toroidal-electron hypothesis that is added later goes in the §1
register beside H1. It must state its standard counterpart, and it must
meet the measured electron data:

- the Zeeman and anomalous-magnetic-moment values (g − 2);
- spin ½ with magnetic moment ≈ 1 Bohr magneton;
- no detected size down to about 10⁻¹⁸ m in high-energy scattering.

Sources of these limits are to be cited when the hypothesis is entered.
