# H47 (draft): the scalar part of Maxwell's quaternion equations couples to gravity

Reading notes for step 1, 8 November 2026 (user request). Equations are
taken from the papers' own text, with their equation numbers; nothing here
is built yet. Proposed register row at the end.

## 1. Where the scalar part is in Maxwell

Maxwell's quaternion equation (A), `𝔅 = V.∇𝔄`, carries his own condition
`S.∇𝔄 = 0` (Treatise Vol. II, Art. 619; CLAUDE.md §3c): the scalar part of
∇𝔄, the convergence of the vector potential, is set to zero by choice (the
Coulomb gauge). The Lorenz gauge does the same to `S = ∂_μA^μ = ∇·A +
(1/c²)∂φ/∂t`. So the scalar is present in the algebra and removed by
convention, not missing from the quaternion form. Extended electrodynamics
(H4 Jack's T, H8 van Vlaenderen, Hively) keeps it as a physical field.
H47 is the further step: that this scalar couples to gravity.

**The choice is a convention (user's note, 2026-11-08).** Nothing in
Maxwell's equations for E and B forces S = 0: Maxwell imposed it, and
Lorenz's condition (1867) is likewise a choice made by the theory's users.
What the choice does and does not change is the point to keep exact:

- *Dropping the condition alone.* In standard theory the Lagrangian
  (−ε₀c/4)F_μνF^μν is gauge-invariant, so S never enters E or B. Leaving S
  non-zero then changes the potentials only; the field equations for E and
  B are the same, and S has no source, energy or force.
- *Making S physical.* That needs a change to the theory itself: a term
  that contains S, such as the Aharonov–Bohm (Ohmura–Stueckelberg) term
  −(ε₀c/2)(∂_μA^μ)² in Minotti & Modanese's action (§2), or Jack's left and
  right derivatives (H4). Then S carries energy, has a source where local
  charge conservation fails, and can couple to other fields.

**Evans and the AIAS group.** M. W. Evans, T. E. Bearden and co-authors,
"Classical electrodynamics without the Lorentz condition: extracting energy
from the vacuum," *Physica Scripta* 61, 513–517 (2000) (author list to be
checked against the journal), discard the Lorenz condition and obtain
Lehnert's equations with a vacuum charge and current built from S, for
example j_A = −(1/μ₀) grad(div A + (1/c²)∂Φ/∂t), and propose devices that
extract energy from it. G. W. Bruhn's review argues that with the
condition dropped and nothing else changed, the extra terms cancel in the
Ampère and Gauss equations (curl B − (1/c²)∂E/∂t = 0 in vacuum), so no
vacuum current or charge appears: the first case above. On this reading the
AIAS result needs the second case, an added S term, which the paper does
not state. Bruhn's page gives no reply from the authors. F. Amador reports
(personal account, 2026-11-08) that when Evans and others evaluated Bruhn's
arguments they found that they did not stand; that account is the reference
here. He also reports that Evans wrote these same equations (S kept, the
modified Gauss and Ampère laws) in his UFT papers on the AIAS site; the UFT
numbers are to be added here. The step-2 solver lets both readings be computed: with S in the action (AB
term) and without it (Lorenz condition merely dropped). The AIAS papers already reviewed in this repo (`examples/AIAS/`,
UFT 1–2, and H15, H16) do not treat this point. Sources:
[Bruhn, "No Lorenz" review](https://www2.mathematik.tu-darmstadt.de/~bruhn/NoLorenz.html).

## 2. The theory with equations: Minotti & Modanese

F. Minotti, G. Modanese, "A new theory of tensor-scalar gravity coupled to
Aharonov-Bohm electrodynamics," arXiv:2408.05230 (v3), published in
*Mod. Phys. Lett. A* 40 (2025) (as cited in their 2026 paper).

**Fields.** Metric g_μν (signature +,−,−,−), Brans–Dicke scalar φ, a second
scalar ψ, potential A_ν with F_μν, and the Aharonov–Bohm (AB) scalar
`S = ∇^μ A_μ`.

**Action (eq. 1, SI).**

```
S = −(c³/16πG₀)∫√−g φR dΩ + (c³/16πG₀)∫√−g (ω(φ)/φ) ∇^νφ ∇_νφ dΩ
    + (c³/16πG₀)∫√−g φ[½∇^νψ ∇_νψ − U(ψ) − Jψ] dΩ + (1/c)∫ L_mat dΩ
    − (ε₀c/4)∫√−g [F_μν F^μν + 2(∇^μ A_μ)²] dΩ − (1/c)∫√−g j^ν A_ν dΩ
```

with the coupling (eq. 3)

```
J = β_mat(ψ,φ)(8πG₀/c⁴)(T^mat + T^AB) + β_AB(ψ,φ)(8πG₀ε₀/c²) S² + β_φ(ψ,φ) T^φ,
T^AB = −2ε₀c² ∇^μ(S A_μ)   (trace of the AB energy-momentum tensor)
```

The term `2(∇^μA_μ)²` is what keeps S physical (the AB, or Ohmura–Stueckelberg,
electrodynamics); `β_AB S²` is the coupling of S to the two scalars of gravity.

**Field equations.** Metric (eq. 4), Brans–Dicke scalar (eq. 5), ψ (eq. 6)
and the AB equation (eq. 7):

```
∇_μ F^μν = μ₀ j^ν − ∇^ν S − ∇^ν(β_AB φ ψ S)
```

**Weak field.** With β = C_φ φ + C_ψ ψ (perturbations about φ₀ = 1, ψ₀):

```
□S      = μ₀ ∂_ν j^ν                                    (21)  S has a source only where local charge conservation fails
□A^ν    = μ₀ j^ν − S β^ν                                (18)
□φ      = −K_φ S²,   □ψ = −K_ψ S²                       (19a, b; 2026 paper eqs. 13–14)
□χ      = ∂²φ/∂t² + (1/2c²) K_φ C_φ S²                  (2026 paper eq. 11)   force per mass f = −∇χ  (14)
K_φ = 8πG₀ε₀ C_φ /((2ω₀+3)c²),  K_ψ = 8πG₀ε₀ C_ψ /c²    (20a, b)
C_φ = φ₀ψ₀ ∂β_AB/∂φ,  C_ψ = φ₀ψ₀ ∂β_AB/∂ψ  at the vacuum values   (20c, d)
8πG₀ε₀/c² ≈ 1.6×10⁻³⁷ A² N⁻²
```

So a gravitational acceleration appears wherever S² is non-zero, through
the scalars φ, ψ. ω₀ is the Brans–Dicke parameter; ψ₀ and the β derivatives
are free (ψ₀ ≫ 1 is assumed). The paper's own estimates: an impulsive
S-pulse gives f ~ 10⁻²⁹ ψ₀ E₀² τ (SI, E₀ = cS₀); a source with charge
non-conservation I = −γ∇·j gives f ~ 10⁻⁷ γ² ψ₀ J² τ³ (eq. 45); a cavity
gives f_z ≈ γ² K_φ Q P/(c ε₀ h) (eq. 43). "There is presently much space
left for uncertainty."

## 3. Applied to superconductors: Minotti & Modanese 2026

F. Minotti, G. Modanese, "Scalar-tensor gravity and Aharonov-Bohm
electrodynamics with bosons: applications to superconductors,"
arXiv:2603.26986; *Eur. Phys. J. C* 86, 593 (2026).

**Source of S.** For Cooper pairs they use the Noether coupling A_μ j_c^μ
instead of local gauge coupling (eq. 18), giving (eq. 19)

```
□S = −(μ₀q/m) ∇·(|ψ_bos|² A)
```

so S is generated where the condensate density jumps, at a
normal/superconducting junction (pill-box condition eq. 23), with value
(eq. 24) `S ≈ −A·n/λ_L`, and for a uniform current `S ≈ μ₀ i_c /(8π λ_L)`.
Numbers: i_c ≈ 6 kA, λ_L ≈ 200 nm → S ≈ 2×10³ N A⁻¹ m⁻¹. In the bulk S
decays over λ_L; with the gravity coupling it saturates at
`S_sat² = 1/(Λ λ_L²)` (eq. 26), Λ = C_φK_φ + C_ψK_ψ (eq. 22).

**The fit (one parameter, C_φ).** Force per mass across a layer of
thickness d (eq. 27): `f ≈ c² d /(C_φ λ_L²)`.

| Experiment | d | pulse | reported f | C_φ |
|---|---|---|---|---|
| Poher | 60 µm | 16 µs | 1.5×10⁴ m s⁻² | 9×10²¹ |
| Podkletnov | 1 cm | 200 ns | 5×10⁶ m s⁻² | 4×10²¹ |

Derived: Λ ~ 5×10⁵ A² N⁻², S_c ~ 7×10³ N A⁻¹ m⁻¹, K_φ ~ 10⁻¹⁶ A² N⁻²;
threshold when eA/mc ~ 1. They note C_φ ~ √N_D (Dirac's large number,
N_D ≈ 4.1×10⁴²). A superconducting RF cavity (1959 MHz, Q 2063, 18 W) gives
f_z ≈ 6×10⁻⁹ m s⁻², far below detection.

**The authors' own status.** Phenomenological and exploratory; the
agreement of two experiments is "an indication of internal consistency …
rather than … conclusive evidence"; the large ψ₀ means an effective theory
in a restricted regime; solar-system tests need a strong suppression of the
scalar coupling (their Appendix B). No energy-conservation analysis; the
effective T⁰⁰ of the scalar sector is negative for K_φ > 0.

## 4. A source in rotating apparatus: Iadicicco, Modanese & Verolino

A. Iadicicco, G. Modanese, L. Verolino, "Effective observer-split source
terms in rotating frames and gravitomagnetic backgrounds in extended
Aharonov-Bohm electrodynamics," *Gen. Relativ. Gravit.* 58, 106 (2026),
doi:10.1007/s10714-026-03606-2; arXiv:2604.24787v2.

**First, a negative result.** Rotation, a stationary gravitomagnetic
background or a coordinate change cannot by themselves create microscopic
charge non-conservation for ordinary U(1)-invariant matter (eqs. 5–7):
∇_μJ^μ = 0 on shell. Gravitomagnetism enters only through the shift β^i.

**The proposal.** In a 3+1 split adapted to the apparatus (eq. 8,
t^μ = N n^μ + β^μ), the projected continuity equation is
`∂_t ρ + D_i(N j^i) = D_i(ρ β^i)` (eqs. 9–10), which defines an
observer-split term

```
I_split = (1/N) D_i(ρ β^i)          (11)
I_G     = Ω ∂_φ ρ                   (13, 23)  weak rigid rotation, β = Ω × r
□_g S   = μ₀ I_split                (19)      phenomenological closure, no new parameter
```

Signatures: sign reverses with Ω; zero for axisymmetric charge (the ring
integral vanishes, eq. 24); for a Gaussian charge transient of angular width
σ, `|I_G|max = |Ω| A e^(−1/2)/σ` (B3), reduced by smoothing over λ to
`σ²/(σ²+λ²)` of that (B7). Not a local, Lorentz-invariant source law: it
depends on the observer's foliation and vanishes in a small enough inertial
frame, so it is meaningful only coarse-grained over the apparatus.

**Bench protocol.** A. Iadicicco, L. Verolino, "Optical-fiber activation of
observer-split source terms in a rotating aluminum disc,"
arXiv:2608.21403v1 (7 Aug 2026): a 10–20° optically gated sector on a
rotating aluminium disc, I_opt ≈ Ω ∂_φ ρ_ind, response Λ_opt = κ_opt χ_Al;
controls: gating ON/OFF, CW/CCW, uniform 360° illumination, dummy and
insulating discs, independent electrometric calibration of ρ_ind. No
dimensions, speeds or measurements given; a null result bounds Λ_opt.

## 5. Standard prediction and what bounds it

- **Standard:** S = 0 for any conserved current (□S = μ₀∂_νj^ν = 0), so no
  scalar force; gravity is GR, weight unchanged by any circuit.
- **Brans–Dicke:** ω₀ > 4×10⁴ (Cassini; cited in the 2024 paper), unless
  ψ₀ ∂_φβ_mat ≈ 1 (their eq. 16).
- **Equivalence principle and fifth forces:** MICROSCOPE η < 10⁻¹⁵ (already
  cited in H35); the 2026 paper's Appendix B needs the matter coupling
  suppressed accordingly.
- **The data the fit uses:** Podkletnov and Poher reports are not
  independently replicated (H12 holds this caveat); the authors cite
  replication attempts without results.
- **Charge conservation:** every source of S here is effective (junction
  condensate coupling, observer split), not microscopic non-conservation;
  for a fully conserved current the retarded S vanishes (critique recorded
  in the October literature watch, item 3).
- **Energy:** neither paper gives an energy theorem for the coupled system;
  the scalar sector's effective energy density can be negative.

## 6. What kronspice can build (steps 2–4)

The §3k tier-2 gates apply: reduces to Maxwell at S = 0; charge
conservation checked every run; an energy theorem closing in its own ledger
column; every equation cited, free parameters declared per run.

1. **`em/scalar/` with AB electrodynamics.** Field equations with S
   (eq. 7 flat space, 18, 21) on the Yee grid beside `em/fdtd.go`;
   energy theorem with E = cS for S-waves (stated by Minotti & Modanese,
   to be derived and checked); checks from the literature watch: S-mode
   cutoffs = TM cutoffs (item 1), closed-form radiated power and the 0.9 %
   finite-kR case (item 2). S is a hypothetical column; S = 0 reproduces
   the standard field to round-off.
2. **The gravity coupling, weak field.** φ, ψ, χ from eqs. 11–14 (2026)
   with K_φ, C_φ, ω₀ declared; force per mass f = −∇χ on a test mass,
   booked in the momentum ledger (M9n) with the scalar sector's energy in
   its own column. First validation: reproduce the 2026 numbers (S at the
   junction, S_sat, C_φ from the Poher and Podkletnov inputs, the SRF
   cavity's 6×10⁻⁹ m s⁻²).
3. **The rotating disc.** Kronspice already rotates charge
   (`examples/rotating_efield`, `examples/homopolar`): drive a gated sector
   on a disc, compute I_G = Ω∂_φρ and S from eq. 19, and give the predicted
   signature (odd in Ω, zero for uniform illumination, ∝ |Ω|/σ) beside the
   standard prediction (S = 0, no force). Output: what a bench needs to
   resolve, and the bound a null result sets on Λ_opt and C_φ.

## 6a. Step 2 done: `em/scalar/` (8 November 2026)

Aharonov–Bohm electrodynamics on the Yee grid, `scalar.Region` with
`Mode: AB` (S physical) or `Mode: Standard` (S held at 0). Field equations,
derivation and the quaternion form are in the package comment and in
CLAUDE.md §3c ("With the scalar part kept"):

```
∂B/∂t = −curl E,   ∂E/∂t = c²(curl B − μ₀J − grad S),   ∂S/∂t = ρ/ε₀ − div E
⇒ □S = μ₀(∂ρ/∂t + div J);   quaternion form ∇G = μ₀J, G = ∇A = S + F
```

**Energy theorem (new; neither paper gives one).**
`∂u/∂t + div P = −J·E + c²ρS`, with `u = ½ε₀E² + B²/2μ₀ + S²/2μ₀` and
`P = (E×B + S E)/μ₀`. S carries energy, adds a longitudinal flux S·E (a
pure S-wave has B = 0), and does work c²ρS on charge, booked in its own
ledger column. The grid scheme satisfies it exactly step by step.

| Gate (§3k) | Test | Result |
|---|---|---|
| Reduces to Maxwell at S = 0 | `TestReducesToYee`: conserved dipole current, 160 steps | E, B = `em/fdtd.go` to 5e-16; S stays 1.9e-19 T (round-off) |
| Energy theorem closes | `TestEnergyTheorem`: dipole + created charge + initial S, 400 steps | residual 1.2e-14; c²ρS column 2.9e-4 J |
| Source law □S = μ₀I | `TestRetardedScalarWave`: charge created at a node | peak S = μ₀q̇/(4πr) to +0.23 % (axis), +0.70 % (diagonal); arrives at 1.000 r/c; B ≤ 6e-25 T (none) |
| Charge conservation in each mode | same test | AB: Gauss law in the extended form to 2e-19; standard mode: defect 1 |
| Cavity modes | `TestCavityScalarMode`: free S-mode in a conducting box | ω at the grid's Dirichlet value to 7e-7, continuum TM cutoff −0.05 %; B = 0; energy constant to 7e-15 |

**What it settles and what it does not.** The standard-mode run is the
"Lorenz condition dropped, nothing added" reading: created charge then has
no field to carry it and Gauss's law fails, so the AIAS vacuum current needs
the S term in the action (the AB reading) to exist at all. With it, S is a
consistent, energy-carrying field. Whether nature has it is not touched:
that needs a source of extra current (step 3's junctions, step 4's rotating
disc) and a measurement.

## 6b. Step 4 done: the rotating gated disc (8 November 2026)

`examples/scalar_disc`: the Iadicicco–Verolino disc with the observer-split
source I_G = Ω∂_φρ, S from em/scalar's equations (quasi-static), and the
S² gravity of Minotti & Modanese with their fitted K_φ = 10⁻¹⁶, C_φ = 6×10²¹.
Declared bench: 5 cm Al disc, lit ring at 4 cm, 15° sector, 1 nC, 3000 rpm,
sensor 2 cm above the track.

- Source = IMV eq. B3 to 3e-4, ring integral 0 (eq. 24), smoothing = eq. B7.
- S at the sensor 1.1×10⁻¹² T; exactly odd in Ω; exactly 0 with the whole
  ring lit. Standard prediction 0.
- **The papers disagree on symmetry.** IMV's linear response reverses with
  Ω (their decisive control); MM's coupling is through S², which does not.
  A CW/CCW run separates them.
- **Under the fit the force is unmeasurable.** 5×10⁻³⁷ m/s²; at the bench's
  limits (air breakdown 1.5 nC, Al burst 1.05×10⁴ rad/s) 1.3×10⁻³³ m/s².
  The junction S the fit relies on (≈ 2×10³ T) is 4×10¹³ times the disc's
  best. Reaching 10⁻⁹ m/s² needs C_φ ≈ 5×10³³ (9×10¹¹ × the fit), which a
  null at that sensitivity would bound.
- **The fit's Brans–Dicke parameter** is ω₀ ≈ 3.5 (2ω₀ + 3 =
  8πG₀ε₀C_φ/(c²K_φ) ≈ 9.9), against Cassini's ω > 4×10⁴, so the matter
  coupling must be suppressed as the authors' Appendix B says.
- **What the disc can measure** is IMV's linear coefficient Λ_opt, which no
  paper gives: S = 3.7×10⁻¹³ T per nC per 1000 rpm, so a response that
  reverses with Ω and vanishes for the whole ring, divided by S, is Λ_opt.
- **Systematic:** the gated charge's Coulomb field at the sensor is
  2×10⁴ V/m (2×10⁻⁵ m/s² on a 1 g sensor carrying 1 pC), 10³¹ times the S²
  force; it follows the charge, does not reverse and stays with the whole
  ring lit, so the controls separate it.

Step 3 (reproducing the 2026 junction numbers in the weak-field φ, ψ, χ
model) remains; step 4 used its quasi-static force law with their constants.

## 7. Open items before building

- The 2026 paper states ψ₀ ~ 10²¹ but does not give the conversion from
  C_φ to ψ₀; check against the published EPJC version.
- β_AB's form is free; only its derivatives at the vacuum enter the weak
  field. Declare C_φ, C_ψ, ω₀ per run.
- Energy theorem for AB electrodynamics coupled to φ, ψ: derive from the
  action (eq. 1) before the ledger column is trusted.
- The disc protocol gives no dimensions or speeds; take them from the PDF
  or choose and declare them.

## 8. Proposed register row

| # | Hypothesis | Modeled / tested in | Status |
|---|---|---|---|
| H47 | The scalar part of Maxwell's quaternion equations (the convergence S.∇𝔄 that Maxwell set to zero; S = ∂_μA^μ) is physical and couples to gravity, so circuits that generate S change the gravitational field (user hypothesis, 2026-11-08; equations from Minotti & Modanese, MPLA 40 (2025), arXiv:2408.05230, and EPJC 86, 593 (2026); source in rotating apparatus from Iadicicco, Modanese & Verolino, GRG 58, 106 (2026)) | `doc/h47-scalar-gravity.md`; planned `em/scalar/` (AB electrodynamics), weak-field φ, ψ, χ coupling, rotating gated disc after Iadicicco & Verolino (arXiv:2608.21403) | untested — theory is a Brans–Dicke + second-scalar gravity with β_AB S² coupling; S has a source only where local charge conservation fails (□S = μ₀∂_νj^ν), effectively at superconducting junctions (S ≈ μ₀i_c/8πλ_L ≈ 2×10³ N A⁻¹ m⁻¹) or as a rotating observer's split term Ω∂_φρ; one parameter C_φ ≈ 4–9×10²¹ fits the Podkletnov and Poher reports (unreplicated), authors call it phenomenological; standard prediction S = 0, GR weight; bounds: Brans–Dicke ω > 4×10⁴, MICROSCOPE η < 10⁻¹⁵, no energy theorem stated |

## Sources

- [Minotti & Modanese, arXiv:2408.05230](https://arxiv.org/html/2408.05230v3)
- [Minotti & Modanese, arXiv:2603.26986](https://arxiv.org/pdf/2603.26986) · [EPJC record](https://epjc.epj.org/articles/epjc/abs/2026/06/10052_2026_Article_15825/10052_2026_Article_15825.html)
- [Iadicicco, Modanese & Verolino, arXiv:2604.24787](https://arxiv.org/abs/2604.24787) · [PDF](https://arxiv.org/pdf/2604.24787v2)
- [Iadicicco & Verolino, arXiv:2608.21403](https://arxiv.org/abs/2608.21403)
- Maxwell, *Treatise* Vol. II, Art. 618–619 (1873), in `reference/`
