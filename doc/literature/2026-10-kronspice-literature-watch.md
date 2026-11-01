# KronSpice literature watch — October 2026

Run: 2026-10-25 (scheduled task, unattended). First report in `doc/literature/`,
so the search window is **2026-08-01 → 2026-10-25** (≈ 2½ months), with a short
"Background" section for older items that bear on the register but are not yet
cited in CLAUDE.md §6.

Read before searching: CLAUDE.md §1 (register H1–H12), §3f (`em/duality.go`,
`em/hed.go`, `em/temporal.go`), §3i, §3j, §3k, §6.

**Nothing in CLAUDE.md or the code was changed.** Proposed edits are in the last
section, for review.

---

## Summary table

| # | Item (date) | Bears on | Effect on entry | Usable numbers |
|---|---|---|---|---|
| 1 | Pinchuk, Minotti, Modanese, Howard, Rentzber — *Scalar-longitudinal modes of extended electrodynamics in perfectly conducting waveguides* (arXiv 2610.28916, 20 Oct 2026) | H8; §3k tiers 1–2; §3a.1 Kron cavity; `em/fdtd.go` | **Extends / bounds** H8 | SL cutoffs = TM (Dirichlet) cutoffs; rectangular, circular, coaxial, L-shaped benchmarks |
| 2 | Rentzber — *Scalar-longitudinal radiation in extended electrodynamics with multipole theory and a compensated source model* (arXiv 2608.19405, 19 Aug 2026) | H8, H9; §3k tier 2 gates 2–3 | **Extends** H8, **bounds** H9's interpretation | Closed-form SL power; R = 1 cm, 1 GHz case; cos²θ pattern |
| 3 | Sautbekov — *Covariant electrodynamics with a scalar degree of freedom* (arXiv 2608.26839, 27 Aug 2026) | H8; §3b `hyper/ga`; van Vlaenderen S-field | **Extends** (GA formulation); critique recorded | none (formal) |
| 4 | Iadicicco, Modanese, Verolino — *Effective observer-split source terms in rotating frames and gravitomagnetic backgrounds in extended Aharonov–Bohm electrodynamics*, Gen. Rel. Grav. 58, 106 (published 10 Sep 2026) | H8; §3f rotating-field / homopolar examples; H12 | **Extends** (proposes a source for the scalar) | scaling laws only |
| 5 | Iadicicco, Verolino — *Optical-fiber activation of observer-split source terms in a rotating aluminum disc* (arXiv 2608.21403 / SSRN 7306298, 20 Aug 2026) | H8, H12; §3f `examples/homopolar/`, `examples/rotating_efield/` | **Extends** — bench protocol, no data | protocol only (10–20° gated sector, CW/CCW reversal) |

Searched with **no new item in the window**: Bearden / MEG replications;
Meyl scalar/vortex waves or critiques; Monstein & Wesley replications;
Biefeld–Brown thrust in vacuum; Giardino (new papers are quaternionic quantum
mechanics, not electrodynamics); Jack's temporal field; van Vlaenderen; Kron /
diakoptics / TAN (Maurice, Ravelo groups); Farnsworth-type multipactor;
Tesla bifilar / helical-line coils; spark-gap recovery in water mist or fog.
Unverified leads for those topics are listed at the end.

---

## New items (in window)

### 1. Scalar-longitudinal modes in perfectly conducting waveguides

A. Pinchuk, F. Minotti, G. Modanese, J. Howard, N. Rentzber, "Scalar-Longitudinal
Modes of Extended Electrodynamics in Perfectly Conducting Waveguides and
Divergence Cleaning Interpretation," arXiv:2610.28916 [physics.class-ph], v1,
20 Oct 2026. <https://arxiv.org/abs/2610.28916>

*Summary.* In extended electrodynamics (EED), where `C = ∇·A + (1/c²)∂Φ/∂t` is
dynamical rather than gauged away, a closed guide carries a fourth mode family
alongside TE/TM/TEM: scalar-longitudinal (SL) modes with B = 0 and E along ∇C.
With the wall condition C = 0 the SL family decouples, its cutoffs equal the
Dirichlet (TM) cutoffs of the cross-section, and there is no coaxial TEM
counterpart; charge-conserving sources with zero initial scalar data cannot
excite it. A real conductor, modeled through its constitutive law, gives
∂C/∂n = 0 instead and couples the scalar sector to the Maxwell fields. The
authors note that undamped hyperbolic divergence cleaning at speed c obeys the
same equations and wall condition, so the SL solutions double as exact
benchmarks for Maxwell solvers.

*Critique on record.* Only an automated review (Pith, AI) so far: maths judged
sound; main weakness that SL modes add no new cutoffs and do not exist as
independent modes in guides of ordinary metal; the power figures come from a
conserved quadratic form, not a physical flux. (Note: the arXiv abstract page
was rate-limited during this run; details are from the Pith summary page and
should be checked against the PDF.)

*Bears on.* H8; §3k tier 1 (TM cavity, Kron 1945 network) and tier 2; `em/fdtd.go`
/ `em/tlm.go` divergence handling.
*Effect.* **Extends** H8 with a closed-geometry prediction, and **bounds** it:
in a guide with ideal walls and conserved sources the SL sector is identically
zero, so a simulated or bench guide can show it only through the wall boundary
condition or an anomalous source.
*Validation case.* For each tier-2 model in `em/scalar/`: SL cutoff frequencies
of a rectangular guide must equal the TMₘₙ cutoffs `f_c = (c/2)√((m/a)²+(n/b)²)`
under C = 0 walls; switching to ∂C/∂n = 0 must produce sector coupling. The same
solutions test any divergence-cleaning option in `em/fdtd.go`.

### 2. Scalar-longitudinal radiation: multipoles and a compensated source

N. Rentzber (UCCS), "Scalar-Longitudinal Radiation in Extended Electrodynamics
with Multipole Theory and a Compensated Source Model," arXiv:2608.19405
[physics.class-ph], v1, 19 Aug 2026. <https://arxiv.org/abs/2608.19405>

*Summary.* Starting from Ohmura's (1956) Lagrangian with C dynamical, the paper
shows that the scalar sector is sourced only by the continuity anomaly
`Λ = ∂ρ/∂t + ∇·J`; for conserved sources with no initial C, EED reduces exactly
to Maxwell. For anomalous sources, the far field splits into the usual
transverse channel and a scalar-longitudinal channel (longitudinal E, B = 0,
co-propagating C, `E_L = cC r̂`) whose time-averaged powers add without cross
terms. It gives an all-orders multipole formula for the SL power, proves that
bound polarization and magnetization currents cannot excite it, and builds a
compensated polarized sphere with a neutral anomalous surface layer that
radiates only in the SL channel.

*Critique on record (in the paper itself).* The adopted energy density
`u = ε₀E²/2 + B²/2μ₀ + C²/2μ₀` and flux `S = (E×B + C E)/μ₀` are not derived from
a physical stress-energy tensor; the author notes that Minotti & Modanese's
tensor treatment (their ref. [10]) gives an *inward* flux of the same magnitude
for the same source. The results are stated as conditional predictions, not
evidence that charge conservation fails.

*Bears on.* H8 (tier-2 gates 2 and 3); H9 (Monstein & Wesley).
*Effect.* **Extends** H8 with quantitative predictions. **Bounds** the EED
reading of H9: a ball antenna fed by an ordinary (charge-conserving) circuit has
Λ = 0, so under this EED family it radiates no SL wave; an SL signal would need
an explicit anomaly source. The dipole pattern is a usable discriminator: SL
power ∝ **cos²θ** about the anomaly-dipole axis versus **sin²θ** for an ordinary
electric dipole.
*Validation numbers.*
- Dipole power `⟨P⟩ ≈ μ₀⟨|dp_Λ/dt|²⟩/(12πc)`, `p_Λ = ∫ r Λ d³r`.
- Shell model: `p_Λ,0 = (4π/3)R³λ_a`; exact `⟨P⟩ = (2πcμ₀/3) λ_a² R⁴ j₁²(kR)`
  (nulls at zeros of j₁); long-wavelength `⟨P⟩ = 2πμ₀ω²R⁶λ_a²/(27c)`.
- R = 1 cm, f = 1 GHz → kR ≈ 0.21, exact vs dipole formula differ by ≈ 0.9 %.
- With the paper's normalization `λ_a = ξ J_ref`, `J_ref = 10⁶ A/m²` (units as
  defined in the paper):
  `⟨P⟩ ≈ 3.9×10⁴ ξ² W`, on-axis peak `|E_L| ≈ 2.6×10³ ξ V/m` at 1 m;
  ξ ≈ 10⁻⁹ → a few µV/m and ≈ 40 fW.
These make a good unit test for `em/scalar/` (a source-free check, a j₁ null
check and the 0.9 % finite-kR check).

### 3. Clifford-algebra electrodynamics with a scalar S

S. Sautbekov, "Covariant Electrodynamics with a Scalar Degree of Freedom,"
arXiv:2608.26839 [physics.class-ph], v1, 27 Aug 2026.
<https://arxiv.org/abs/2608.26839>

*Summary.* Working in the spacetime Clifford algebra, the author keeps the
scalar part S of the geometric derivative of the four-potential next to the
field bivector, giving a generalized Maxwell system that reduces to the
standard one at S = 0. He derives an energy-momentum tensor, conservation laws
and a generalized Lorentz force, and finds free longitudinal "S-waves"
(longitudinal E, no B, speed c). The proposed source is again the effective
current divergence `∂ρ/∂t + ∇·j` of an open subsystem, e.g. a charge passing an
aperture.

*Critique on record.* Automated (Pith) review only: algebra sound, but for any
conserved total current the paper's own retarded solution gives S = 0, the
aperture example omits the screen's compensating response, and the
energy-momentum tensor is non-symmetric with an unusual momentum sign for pure
S-waves.

*Bears on.* H8 (this is essentially the van Vlaenderen S-field written in GA);
§3b `hyper/ga` (natural home for the grade-0 part of `∇A`).
*Effect.* **Extends** — a GA statement of the tier-2 equations that fits the
project's grade-selection style (§3f). Same conservation bound as item 2.

### 4. Observer-split source terms in rotating frames (published)

A. Iadicicco, G. Modanese, L. Verolino, "Effective observer-split source terms
in rotating frames and gravitomagnetic backgrounds in extended Aharonov–Bohm
electrodynamics," *Gen. Relativ. Gravit.* 58(9), 106 (2026), published
10 Sep 2026, open access, DOI 10.1007/s10714-026-03606-2
(<https://link.springer.com/article/10.1007/s10714-026-03606-2>; preprint
arXiv:2604.24787).

*Summary.* The authors first argue that rotation, a stationary gravitomagnetic
background or a coordinate change cannot create true microscopic charge
non-conservation for ordinary matter. They then show that splitting a conserved
current into quantities measured by a rotating observer gives an exact
observer-dependent term which, for weak rigid rotation, is proportional to
`Ω ∂_φ ρ` (angular velocity times the azimuthal gradient of charge density). They
propose using it as a phenomenological source for the AB scalar sector, with the
observer foliation fixed by the lab's synchronization and readout.

*Bears on.* H8; §3f (`examples/rotating_efield/`, `examples/homopolar/` — both
rotate charge distributions); H12 (gravitomagnetic/weight anomalies).
*Effect.* **Extends** — first concrete, falsifiable source mechanism for the
tier-2 scalar in a rotating apparatus. No absolute magnitudes.
*Predicted signatures (testable in simulation).* Sign reversal with rotation
reversal; zero for axisymmetric charge (the source integrates to zero around a
ring for any single-valued periodic density); for a Gaussian transient of
angular width σ, peak `|I_G| = |Ω| A e^{-1/2}/σ`; smoothing width λ reduces it
by `σ²/(σ²+λ²)`.

### 5. Optical-fiber gated rotating aluminum disc (bench proposal)

A. Iadicicco, L. Verolino, "Optical-Fiber Activation of Observer-Split Source
Terms in a Rotating Aluminum Disc," arXiv:2608.21403; SSRN 7306298 (posted
20 Aug 2026; earlier version 5 Aug 2026).
<https://arxiv.org/abs/2608.21403> · <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7306298>

*Summary.* Companion experiment to item 4: an illuminated fiber-optic ring
masked to a 10–20° active sector creates a localized angular density on a
rotating aluminum disc, gated ON/OFF and run CW and CCW, with material and dummy
controls. The paper gives a source coordinate `ΔG_opt = Ω χ_Al T_ev e^{-1/2}/σ`,
a response coefficient `Λ_opt = κ_opt χ_Al`, a sensitivity target, a
minimum-detectable-force expression, a null-bound recipe and an inverse
calibration to separate χ_Al and κ_opt. It reports no measurements and makes no
claim of an anomaly.

*Bears on.* H8, H12; §3f rotating examples.
*Effect.* **Extends** — a protocol whose controls (gating, rotation reversal,
dummy disc) match KronSpice's "standard prediction alongside" rule and could be
adopted for a §3f bench run. A null result would give a bound on |Λ_opt|.

---

## Background (before the window, not yet in §6)

- **F. Minotti, G. Modanese, "Scalar–tensor gravity and Aharonov–Bohm
  electrodynamics with bosons: applications to superconductors,"** *Eur. Phys.
  J. C* 86, 593 (1 Jun 2026), DOI 10.1140/epjc/s10052-026-15825-7.
  Couples two scalar fields of scalar–tensor gravity to AB electrodynamics; in
  superconductors the scalar S is generated semiclassically at
  normal/superconducting junctions. Fits threshold behaviour in the Poher and
  Podkletnov pulsed-discharge reports with one parameter, which the authors
  call phenomenological, not conclusive. Numbers: junction `S ≈ μ₀i_c/(8πλ_L)
  ≈ 2×10³ N A⁻¹ m⁻¹` (i_c ≈ 6 kA, λ_L ≈ 200 nm), saturation ≈ 7×10³; Poher
  case d ≈ 60 µm, 16 µs; Podkletnov d ≈ 1 cm, 200 ns. Superconducting RF
  cavities predicted far below detectability. **Bears on H12** (electrogravitic
  claims, via a different route than Biefeld–Brown) and H8. The underlying
  Poher/Podkletnov reports are themselves unreplicated; record that beside it.
- **K. Rębilas, "On the origin of 'longitudinal electrodynamic waves',"**
  *EPL* 83, 60007 (2008). Attributes the Monstein–Wesley signal to ordinary
  transverse radiation from currents the ball antenna induces in the ground
  (plasma model of the soil); argues no evidence for scalar theories.
  **Contradicts H9's interpretation.** Note: KronSpice's §3k already includes the
  ground image source; Rębilas implies the image must be a *lossy, conducting*
  ground with induced currents, not a perfect mirror.
- **E. J. Butterworth, C. B. Allison, D. Cavazos, F. M. Mullen, "Longitudinal
  electromagnetic waves? The Monstein–Wesley experiment reconstructed,"**
  *J. Sci. Explor.* 27(1) (2013),
  <https://journalofscientificexploration.org/index.php/jse/article/view/477>.
  A duplicate apparatus reproduced the qualitative pattern (closer to MW's
  theory than MW's own data) but did not support Rębilas' ground-current
  explanation; the authors leave the mechanism unidentified. **Partial
  replication of H9** (qualitative only; no published numbers seen).
- Other 2026 AB-electrodynamics preprints seen but not read in full:
  arXiv:2604.15913 (fluctuations in AB electrodynamics), 2605.10986 (Bohm
  potential and the AB scalar mode), 2605.08148 (Stueckelberg scalar in the
  electric AB effect — possibly a laboratory bound; worth reading next run).

---

## Cross-cutting observation

Items 1–4 (and the earlier Minotti–Modanese papers) converge on one point that
matters for `em/scalar/`: **in every current EED/AB variant the scalar is
sourced only by the continuity anomaly Λ** (or, in item 4, an
observer-dependent stand-in for it). Tier-2 gate 2 (charge conservation checked
every run) is therefore not just a sanity check — if it passes, every tier-2
model will predict exactly Maxwell for any lumped circuit or antenna KronSpice
simulates. To produce a non-trivial prediction a run must declare a Λ source,
just as `em/hed.go` declares β. Two further modelling choices are open in the
literature and must be stated per run: the conductor boundary condition
(C = 0 vs ∂C/∂n = 0, item 1) and the energy theorem (Rentzber's field-only
flux vs the Minotti–Modanese tensor, which disagree in sign, item 2).

---

## Proposed register and spec updates

None applied. For review:

1. **§1 register, H8 row** — keep *untested*; add to the "Modeled / tested in"
   cell: "Λ-sourced scalar; see §3k gate 5". Add a note under the table: "All
   published EED/AB variants (Ohmura 1956; Hively; Minotti–Modanese;
   Rentzber 2026; Sautbekov 2026) source the scalar only through
   Λ = ∂ρ/∂t + ∇·J; with charge conserved they reduce to Maxwell."
2. **§1 register, H9 row** — add "critiques: Rębilas 2008 (ground currents);
   partial replication: Butterworth et al. 2013". Consider a status note
   "EED reading bounded: Λ = 0 for a conventionally fed antenna (Rentzber
   2026)".
3. **§1 register, new row H13 (proposed)** — "Observer-split rotating-frame
   source for the AB scalar (Iadicicco, Modanese, Verolino 2026)", modeled in
   `examples/rotating_efield/` / §3f with the item-5 protocol; status
   *untested*. Alternatively fold into H8 as a sub-case.
4. **§1 register, H12 row** — add Minotti & Modanese (EPJC 2026) as a second
   theoretical route (superconducting junction discharges), flagged
   phenomenological, with Podkletnov/Poher listed as unreplicated sources.
5. **§3k tier 2** — add a fifth gate: "5. the Λ source, the conductor boundary
   condition for the scalar (C = 0 or ∂C/∂n = 0), and the energy theorem used
   (field-only flux or a named stress-energy tensor) are declared per run."
   Add Ohmura (1956) and the Aharonov–Bohm electrodynamics family
   (Minotti & Modanese) to the tier-2 list; note Sautbekov (2026) as a GA
   statement suitable for `hyper/ga`.
6. **§3k validation cases** — add: (a) SL waveguide cutoffs equal TM cutoffs
   under C = 0 walls, rectangular/circular/coax/L-shaped (Pinchuk et al. 2026);
   (b) Rentzber's compensated-shell power with j₁(kR) nulls and the 0.9 %
   finite-kR check at R = 1 cm, 1 GHz; (c) cos²θ vs sin²θ angular pattern as a
   reported output in the Monstein–Wesley simulation.
7. **§3k Monstein–Wesley bullet** — replace "Later Comment and Reply in EPL"
   with explicit citations (Rębilas 2008, Butterworth et al. 2013), and model
   the ground as lossy and conducting (induced-current radiation) in addition
   to the perfect image source.
8. **§3a.1 / `em/fdtd.go`** — note that undamped hyperbolic divergence cleaning
   shares the SL equations, so the item-1 solutions can validate a
   divergence-cleaning option.
9. **§6 References** — add items 1–5 and the three background references above,
   each marked *hypothesis source* or *critique* as appropriate, plus
   T. Ohmura, Prog. Theor. Phys. 16 (1956) (citation to be verified before
   adding).

---

## Unverified leads (could not open; check next run)

- "Adsorption-induced evolution of secondary electron emission characteristics
  and multipactor process," *J. Phys. D* (DOI 10.1088/1361-6463/ae8fe3) — §3i;
  date and content not verified (site blocks automated reading).
- "Data-Driven Analysis of Multipactor Dynamics via Dynamic Mode
  Decomposition," OSTI 3577308 — §3i; date not verified.
- "Humidity as a critical parameter for predicting breakdown voltage in
  submicrometer electrode gaps in air," *Phys. Rev. E* (DOI 10.1103/h4xn-4kvy)
  — §3j / H5 background (humidity, not mist; sub-µm gaps); date not verified.
- "Experimental study on DC breakdown characteristics of typical air gaps under
  high altitude and high humidity conditions," *Proc. SPIE* 14136, 141362J —
  §3j; content not verified.
- "Atmospheric EAD grid thruster performance: a theoretical–experimental
  approach on geometry and pressure effects," *J. Phys. D* (DOI
  10.1088/1361-6463/ae251d) — H12 standard (ion-wind) prediction vs pressure;
  not verified.

## Sources

- https://arxiv.org/abs/2610.28916 · https://pith.science/paper/2610.28916
- https://arxiv.org/abs/2608.19405 · https://arxiv.org/html/2608.19405
- https://arxiv.org/abs/2608.26839 · https://pith.science/paper/2608.26839
- https://link.springer.com/article/10.1007/s10714-026-03606-2
- https://arxiv.org/abs/2608.21403 · https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7306298
- https://link.springer.com/article/10.1140/epjc/s10052-026-15825-7
- https://www.researchgate.net/publication/258276804_On_the_origin_of_longitudinal_electrodynamic_waves
- https://journalofscientificexploration.org/index.php/jse/article/view/477
