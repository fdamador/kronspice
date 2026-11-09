# CLAUDE.md — examples/ether_cavitation (can the nonlinear ether cavitate?)

Notes for Claude sessions working in this example. The repo's root CLAUDE.md
governs everything here: standard prediction beside every hypothesis, the
energy ledger closes (§3c), cite source and equation for every formula,
hypothetical results labelled, declared inputs stay declared. New example ⇒
a row in `examples/README.md`, its own `README.md`, and a citation in a
`book/` chapter (memory: examples-index-row). Every Go file starts with the
SPDX dual-license header.

## The question (user, 2026-11-08)

A liquid cavitates: under tension it tears into its own vapour, and the
cavity later collapses violently (`examples/bubble_collapse`). Can the
nonlinear ether of `examples/ether_knot` (H29) do the same? And if it can,
what does that tell us about the particle models that use this ether
(H13, H29, H31, H44, H28, H30)?

## Why the present ether cannot, and what has to change

`ether_knot` uses Faddeev's unit-vector ether, |n| = 1 everywhere. It has
no amplitude that could drop, so it **cannot cavitate by construction**.
That is also why its knot holds: a knot can untie only through a point
where the field vanishes (Kleckner, Kauffman & Irvine, *Nature Phys.* 12,
650, 2016, where every superfluid vortex knot unties at its cores).

So the work is to give the ether an amplitude a ≥ 0 and an equation of
state: φ = a·n, with

```
E = κ2 ∫ [ (∂a)² + a² (∂n)² ]  +  κ4 ∫ a^p |B_n|²  +  ∫ V(a)
B_n,k = ½ ε_kij n·(∂_i n × ∂_j n)          (as in ether_knot)
```

- **Cavity** = a region where a → 0. Inside it n is undefined, so the Hopf
  charge can change there.
- **λ → ∞ (stiff amplitude) must give back `ether_knot` exactly.** This is
  the twin test: E → 547.5 √(κ2κ4) as h → 0, E2/E4 → 1.
- **p is declared per run.** Run both p = 0 and p = 4.
  - p = 0: the Faddeev term does not depend on a. This is the
    Babaev–Faddeev–Niemi two-component condensate with the supercurrent
    frozen (*PRB* 65, 100512, 2002).
  - p = 4: the term scales as |φ|⁴. This is the "linear" σ-model form.
- **Origin of the model.** A charged two-component superfluid has exactly
  this structure. In it the density ρ = |ψ1|² + |ψ2|² can vanish, and a
  zero line carries a flux quantum (Babaev, *PRB* 79, 104506, 2009). The
  supercurrent coupling there was found to destabilize the knots (Jäykkä &
  Speight, *PRD* 82, 125030, 2010). Read both before step 2, check the
  citations against the papers, and run the supercurrent as an option.

### Two equations of state (both are genes)

| | V(a) | What a cavity is | Liquid analogue |
|---|---|---|---|
| **Q (quartic)** | (λ/4)(a² − 1)² | a = 0 is a maximum, not a phase. Cavities exist only as cores of defects, where the gradient energy pays for them (vortex cores, type II) | cavitation at a nucleus: no vapour without a seed |
| **S (sextic)** | λ a²(a² − 1)² + ε·(1 − a²)³-type tilt (state the exact form in code) | a = 0 is a second minimum, a true "vapour" phase, split from the ether by ε per volume | liquid–vapour coexistence; ε > 0 is the vapour pressure deficit; ε < 0 is a superheated liquid |

**Units:** κ2 = κ4 = 1, so lengths are in ℓ = √(κ4/κ2) and energies in
E0 = √(κ2κ4), as in `ether_knot`. The new numbers are:
- **β = λκ4/κ2²**, which is (ℓ/ξ)² with ξ = √(κ2/λ) the healing length;
- ε (sextic only);
- p.

The ether's "tensile strength" is reported in E0/ℓ³.

## Dictionary: liquid ↔ ether

| Liquid (`bubble_collapse`) | Ether (here) |
|---|---|
| pressure p, tension p < p_v | trace of the field stress tensor, from Derrick scaling (virial); tension = squeezing or stretching the knot |
| vapour | a = 0 phase (S), or a defect core (Q) |
| surface tension σ | wall tension σ_w of the 1-D a: 1 → 0 profile (step 1) |
| Blake threshold, critical radius R* = 2σ/Δp, barrier 16πσ³/3Δp² | same forms with σ_w and ε (S). Coleman's thin-wall bounce has R = 3σ/ε and B = 27π²σ⁴/(2ε³) (*PRD* 15, 2929, 1977) |
| nuclei (dust, gas pockets) | particles: knots as nucleation sites (step 4) |
| Rayleigh collapse, t_c ≈ 0.915 R√(ρ/Δp) | collapse of an a = 0 bubble in the ether (step 3), with the inward energy booked (H30) |
| flash and shock at collapse | outgoing ether waves (massless n waves, plus amplitude waves of mass ∝ √λ) |

## Standard counterparts (carried beside every result)

1. **Liquids:** classical nucleation and the Blake threshold. The
   `bubble_collapse` circuit is the working twin.
2. **Superfluids and superconductors:** vortex cores at the healing
   length; type I/II set by the Ginzburg–Landau parameter; reconnection
   unties knots (Kleckner 2016).
3. **The QED vacuum, which is the measured "ether":** Schwinger pair
   production, rate ∝ exp(−πE_s/E), E_s = 1.32×10¹⁸ V/m. This is the
   vacuum tearing under tension. It has not been observed: the strongest
   lasers reach ~10¹⁵ V/m. `em/eulerheisenberg.go` covers the field range
   below it and must not be used past its 0.1 E_s weak-field limit;
   `schwinger.go` gives the rate column.
4. **False-vacuum decay** (Coleman 1977; Callan & Coleman, *PRD* 16, 1762,
   1977): bubble nucleation of a lower vacuum. Whether our vacuum is
   metastable depends on the measured Higgs and top masses.

## Steps

**0. Liquid twin.** From `bubble_collapse`, extract:
- the Blake threshold;
- the critical radius;
- the Rayleigh collapse time.

These are written in the variables of the dictionary, so each ether result
has a liquid number beside it.

**1. The wall and the bubble in 1-D (`wall.go`, `radial.go`).**
- Planar profile a(x) for both equations of state. Wall tension σ_w(β, ε)
  is checked against the closed form; the Q profile is the tanh kink.
- Spherical critical bubble a(r) for S: R*, E* and the thin-wall limit
  checked against the formulas above.
- **Kron twin (`kron.go`):** the radial equation as a Kron 1945 ladder
  network (as `light_hydrogen` does for Schrödinger), with V′(a) as
  nonlinear shunt B elements (`device/behavioral.go`), solved by kronspice
  OP. It must equal the Go profile to discretization error.

**2. The knot in the soft ether (`soft.go`), quasi-static cavitation.**
- Copy `ether_knot/axi.go` here (package main cannot import it; keep a
  comment naming the source) and add the amplitude field a(ρ, z).
- Arrested Newton flow on (m, a).
- Twin: β → ∞ gives `ether_knot`.
- Scan β (both Q and S, both p) and record:
  - min a and where it sits;
  - E, E2/E4, core and rms radius;
  - the degree (Hopf charge) followed as in `claude_electron`.
- Find **β_c, where the core first cavitates**, and whether the charge
  survives below it.
- The same scan runs on the A_{M,K} family from `claude_electron` (A11,
  A12, A21).

**3. Driven cavitation (`drive.go`): the ether's tensile strength.**
Drives:
- **Squeeze:** Derrick rescaling x → λx imposed, or an external pressure
  term P∫(1 − a²).
- **Strike:** the `modes.go` strikes (squashing, breathing) at rising
  amplitude.
- **Spin:** rotation as in `claude_electron/spin.go`. Centrifugal
  stretching thins the core.
- **Collapse:** shrink the knot below its size, release it, and watch the
  core cavitate and the knot untie or survive.

For each drive, record:
- the threshold drive;
- the energy barrier to untying (a saddle/sphaleron search between Q = 1
  and Q = 0);
- what is radiated.

**Ledger (dynamic runs).** Field energy + energy out through an absorbing
boundary layer = work done by the drive, as in `knotted_light`. Report the
residual. A surplus is a bug.

**4. Nucleation and particles as nuclei (`nucleate.go`), S only.**
- **Free nucleation:** does the computed R* and barrier match thin-wall
  Coleman as ε → 0?
- **Heterogeneous nucleation:** a knot sits where the gradient energy is
  high, so it should lower the barrier, as dust does in a liquid. Compute
  the barrier with a knot at the centre against without one.
- If the ether is metastable (ε < 0), the knot's lifetime against
  catalysing decay is a bound on the model.

**5. Particle analysis (`particles.go`).** Each item has its standard
prediction beside it.

| Entry | What cavitation adds | Test | Standard counterpart |
|---|---|---|---|
| H29 knot | a stability limit: tensile strength and an untying barrier, so a lifetime | step 3 barrier against ambient (ZPF-like, declared) fluctuations | no stable knots are seen; the electron's stability > 6.6×10²⁸ yr (Borexino 2015) |
| H44 Claude electron | A12 untied by collapsing onto the axis; in a soft ether it can also untie through a cavitating core | barrier and route for A12 and A11; is the ℓ window (empty by 9×10⁶) moved by finite β? | Dirac + QED, pointlike |
| H13 Tewari | his electron is a vortex round a void, i.e. a cavity | build a vortex with an a = 0 core in the Q ether; size, energy, moment, g against his 15/8 and his density limit | g = 2.002 319… |
| H31 nucleons | Skyrme's |U| = 1 is the same stiff limit; relaxing it is the chiral-soliton / linear-σ route (Birse & Banerjee, *PLB* 136, 284, 1984) | does finite β lower the B = 4 binding (7.7 % in `nucleon_circuits/skyrme`) toward the measured 0.75 %? Compare with lightly bound Skyrme models (Gillard, Harland & Speight, *Nucl. Phys. B* 895, 272, 2015; check the citation) | ⁴He binds 28.3 MeV |
| H28/H39 neutrino | a neutral ring with an a = 0 core: a cavity vortex ring | its speed against size (Kelvin's formula, as in `vortex_ring`) and its energy against the E-ring of H28 | weak-only neutrino, m < 0.45 eV |
| H30 implosion and explosion | the ether bubble's collapse is an implosion with a real medium supplying the inward half | step 3 collapse with the ledger; compare the shape of the pulse with `bubble_collapse` | the inward energy is work by the drive, booked |

## Scales and bounds

κ2, κ4 and λ are not fixed by any measurement, so every result is a ratio,
as in H29:
- β_c;
- tensile strength in E0/ℓ³;
- barrier in E0.

Where a register entry fixes the scale (H44: E0 from m_e c² with a chosen
ℓ), convert and set the result against the measured bound of that entry.
Never present a dimensionless threshold as a prediction in SI without the
scale it assumed.

## First look (2026-11-08, before any code here)

**Local estimate.** Drop the healing term (∂a)². Each point then settles at

```
p = 0:  a² = 1 − 2κ2 (∂n)²/λ
p = 4:  a² = (λ − 2κ2 (∂n)²) / (λ + 4κ4 |B_n|²)
```

so a reaches 0 where (∂n)² ≥ λ/(2κ2). The quartic term does not move this
threshold, only the shape near it. In ether units this is
β ≤ 2·max((∂n)²)·ℓ².

**Measured on `ether_knot`'s relaxed Q = 1 knot** (box 8 ℓ, h = 0.1 ℓ,
E = 542.1, core 1.15 ℓ; one coarse run, no grid study):
- max (∂n)² = 16.7/ℓ²;
- max e4 = 69.4/ℓ⁴;
- both at the knot's centre on the axis (ρ = 0.15, z = 0), where m points
  to the antipode, not on the core ring.

So **β_c ≈ 33 (local estimate)**:
- above it the ether stays filled everywhere;
- below it the knot first cavitates at its centre, a point on the axis.

The healing term resists, so the true β_c is lower. Finding it is step 2's
first number. Note that the cavity forms on the axis, not on the ring:
the axis is where the double loop A12 collapsed in `claude_electron`.

## Step 1 (built, 2026-11-08)

Files: `eos.go`, `ladder.go`, `vec.go`, `wall.go`, `twist.go`,
`reversal.go`, `radial.go`, `kron.go`, `main.go`, `main_test.go`. Results
in README.md and `run_output.txt`. In brief:

- **Walls.** σ_Q = ⅔√(κ2λ) (tanh) and σ_S = ½√(κ2λ) (logistic), to h².
- **Tensile strength.** Eckhaus criterion ½V″(a0) > 3κ2k² ⇔ dj/dk > 0.
  - Q: k_c = 1/(√6 ξ) and j_c = 4√(κ2λ)/(3√6), exact.
  - S: binodal k_B ≈ √ε and spinodal k_c, tabulated.
  - A U(1) ring past k_c slips through a = 0.008, with the ledger closed.
- **S² ether.** A uniform twist escapes over the pole with no cavity
  (min a 0.96). In Faddeev's ether, tension exists only where topology
  forbids that escape. This decides where to look in step 2: knot cores,
  antipodal reversals.
- **Reversal slab.**
  - Q tears continuously at L = 4.934 ξ.
  - S tears with hysteresis: for ε = 0.2 λ, the bending ether tears below
    L_t = 4.12 ξ and the cavity sheet heals only above L_c = 12.65 ξ.
- **3-D critical bubble (S).**
  - A saddle with one unstable mode.
  - R*/R_thin → 0.995 and E*/E_thin → 0.970 at ε = −0.01 λ.
  - Under twist (local approximation, labelled), the barrier falls from
    1413 to 7.9 between binodal and spinodal.
- **Kron twin.** The kronspice OP equals the ladder to 1e-7 for the wall
  and 7e-6 for the bubble saddle. The ledger balances.
  - The PWL segments are centred on each node's expected value. Without
    this, Newton cycles between segments next to a saddle (`kron.go`
    comment).
- **Bridge to step 2.** The local criterion (∂n)²ξ² > ½ and the slab
  criterion (∂n)²ξ² > (π/4.934)² = 0.41 agree. With the knot's peak
  (∂n)² = 16.7/ℓ², they give **β_c ≈ 33–41**.
- **Not run:** `go test -race`, which needs cgo, not available on the
  Windows host. Run it in the Docker tools image (root §7). The package has
  no shared mutable state.

## Files (planned)

- `main.go`, flags `-part`, `-beta`, `-eos Q|S`, `-p`, `-viz`
- `wall.go`, `radial.go`, `kron.go`
- `soft.go`, `drive.go`, `nucleate.go`
- `particles.go`, `schwinger.go`
- `main_test.go`
- `viz/`, `run_output.txt`, `README.md`

**Tests:**
- β → ∞ twin = `ether_knot`;
- wall tension = closed form;
- thin-wall bubble = Coleman;
- Kron network = Go profile;
- gradient check on (m, a);
- dynamic ledger residual;
- liquid twin = `bubble_collapse`.

## Register entry (proposed H47; H46 is taken by harmonic_cavitation; add to root §1 when step 2 runs)

| # | Hypothesis | Modeled / tested in | Status |
|---|---|---|---|
| H47 | The nonlinear ether can cavitate: given an amplitude, it tears where its field bends faster than its healing length, the cavity lets a knot untie, and a collapsing cavity is an implosion; the threshold and barrier set the stability of the knotted-particle models (H13, H29, H31, H44, H28, H30) (user hypothesis, 2026-11-08) | `examples/ether_cavitation` | untested — first look: local estimate β_c ≈ 33 for the Q = 1 knot, cavity first at the knot's centre on the axis; standard prediction: liquids (Blake, Rayleigh), vortex cores at the healing length, Schwinger breakdown at E_s = 1.32×10¹⁸ V/m (not reached), false-vacuum decay (Coleman) |

## Next tasks

1. ~~Step 1~~ (built). Step 0, the `bubble_collapse` liquid numbers in
   the dictionary's variables, is still to do.
2. Step 2: the soft `axi.go`, the β → ∞ twin, the β scan, β_c with a grid
   study (h = 0.1, 0.07, 0.05).
3. Read Babaev 2009 and Jäykkä & Speight 2010 and check every citation in
   this file against its paper.
4. Step 3 drives with the ledger, then steps 4–5.
5. README, `examples/README.md` row, book citation (chapter 10 beside
   H29), register row H47.
