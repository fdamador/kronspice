# The Rodin coil: quadrature drive, field maps, resonance and output groups (H43)

The Rodin coil is a winding on a toroidal former whose wire goes over the
top, down through the hole and back out underneath, each pass advancing
150° round the axis. With 12 poles 30° apart the wire visits
1→6→11→4→9→2→7→12→5→10→3→8→1 and closes: a {12/5} star seen from above,
the aperture an inner dodecagon. One circuit goes **5 times round the axis
and 12 times through the hole**, a torus-knot path of the (5, 12) class.
Builds use two wires on mirrored paths, counter-wound. This example builds
the coil characterized by A. D. Ramírez Galindo, H. A. Pérez Olivas,
G. Basurto Islas, A. Hernández Rayas, F. González López and T. Córdova
Fraga, "Vortex magnetic induction: Mathematical, geometric and experimental
characterization", *Rev. Mex. Fís.* 70, 030901 (2024),
doi:10.31349/RevMexFis.70.030901 (outer diameter 9.95 cm, aperture 2.9 cm,
copper; centre field in their Table II), drives its two windings with a
**sine and a cosine**, maps the magnetic field they make, builds it up by
resonance, and adds **groups of output coils** round the toroid that charge
a **capacitor bank** through Schottky bridges, the drive frequency evolved
together with the groups. Every joule is booked against the two drive
sources.

Register entry **H43** (CLAUDE.md §1). The claims around the coil —
"vortex-based mathematics" (the 1-2-4-8-7-5 doubling sequence), free
energy, over-unity — come only from hobbyist and promotional sources
(M. Rodin's own material, forum posts) with no independent measurement.
They are registered as a hypothesis and tested here against the standard
prediction: an ordinary pair of coupled inductors (a transformer, and near
self-resonance a transmission line), with a ledger that closes.

```
go run ./examples/rodin_resonance            # full run, about 12 minutes (16 cores)
go run ./examples/rodin_resonance -short     # coarse grids, small searches, about 1 minute
go run ./examples/rodin_resonance -part 2    # one part: 1, 2, 2b or 3
python3 examples/rodin_resonance/viz/plot_all.py      # figures from the CSVs in viz/
python3 examples/rodin_resonance/viz/plot_cycle.py    # viz/cycle.gif
python3 examples/rodin_resonance/viz/make_report.py examples/rodin_resonance/rodin_resonance_report.pdf
go test -short ./examples/rodin_resonance    # about 20 s; without -short about 2 minutes
```

Every number below is in [run_output.txt](run_output.txt).

## How the Rodin coil relates to the photon electron (H23)

[photon_electron](../photon_electron) tests Williamson and van der Mark's
electron: one photon of the Compton wavelength on a **(2, 1)** torus-knot
path, twice round the ring and once through the hole, whose field comes
back to itself only after **720°**. It reproduces that property as a
circuit: a two-wire ring line with its conductors **crossed** once has
half-integer modes and returns a pulse inverted after one trip. Its README
section "The circuit mirror of the fdamador neutrino" adds the duality: the
fdamador neutrino (H28: E round the ring, B through the hole) is a **series
RLC loop**, and the photon electron's mirror is a **parallel RLC tank**
(`photon_electron/viz/circuit.png`).

The Rodin winding is the **same family of paths** with other winding
numbers: (5, 12) instead of (2, 1). The differences that matter:

- it carries a **current in a wire**, not light in free space, so nothing
  has to confine it: the wire does;
- at the drive frequencies used here (1 kHz, 100 kHz) the winding is far
  shorter than a wavelength (19.2 m of wire per winding; λ = 3 km at
  100 kHz), so it is a **lumped inductor**, and its field is the
  Biot–Savart field of the current, the same pattern at every instant;
- only near its **self-resonance** (1.966 MHz for one winding, far end
  open) does the winding become a **transmission line laid on a knot**,
  with standing waves along the wire (Part 2b);
- the counter-wound pair connected in series and closed on itself **is**
  the crossed ring of photon_electron: one loop through both windings. Its
  lowest modes come back inverted after one trip (trip parity −1.00), for
  the (5, 12) Rodin winding and for a (2, 1) winding on the same former
  alike. The 720° property belongs to the **connection** (one loop through
  two windings), not to the knot's winding numbers, which set only L, k and
  the frequencies.
- In the circuit-mirror language: a winding tuned by a series capacitor
  (Part 2) is the series RLC loop; the output coils tuned by series
  capacitors (Part 3) are too.

## Part 1 — geometry and validation

**Wire model ("taut wire").** Seen from above a pass is the straight chord
from pole P_k to pole P_k+5; its first half lies on the top surface of the
torus, its second half underneath, and it changes side through the hole at
the chord's midpoint. Former: R0 = 3.2125 cm, a = 1.7625 cm (from the
paper's outer diameter and aperture). One circuit of A is 1.6007 m of wire;
B, the mirror path (P_k → P_k−5) one bundle further out, 1.8529 m.

**Turns (declared, fitted).** The paper gives no turn count or wire gauge.
Wire 0.6 mm copper (declared); each winding lays its 12-pass path Nc times
in a round bundle (fill 0.7). Nc is fitted to the paper's measured centre
field by least squares, A and B in series in the **dipole** sense (the
paper's series connection, which puts a field in the middle):

| Nc | 8 | 9 | 10 | 11 | **12** | 13 | 14 |
|---|---|---|---|---|---|---|---|
| centre field [mT/A] | 1.7684 | 1.9521 | 2.1273 | 2.2995 | **2.4672** | 2.6250 | 2.7811 |

Nc = 12: bundle radius 1.242 mm, A to B bundle pitch 2.584 mm.

**Centre field (Biot–Savart, `em.SegmentB`) against Table II [mT]:**

| I [A] | paper theory | paper FEM | measured | kronspice | vs measured | vs paper theory |
|---|---|---|---|---|---|---|
| 1.60 | 3.90 | 4.08 | 4.13 | 3.948 | −4.4 % | +1.2 % |
| 1.25 | 3.01 | 3.00 | 3.05 | 3.084 | +1.1 % | +2.5 % |
| 0.60 | 1.30 | 1.10 | 1.16 | 1.480 | +27.6 % | +13.9 % |

With one fitted number the field is within the paper's own spread at 1.60
and 1.25 A. At 0.60 A it is 27.6 % above the measurement: Biot–Savart in
air is linear in the current, and the paper's own numbers are not (measured
2.58, 2.44, 1.93 mT/A at 1.60, 1.25, 0.60 A), so no linear model meets all
three. The discrepancy is the paper's, not an air-core effect.

**Inductance (segment method).** Every wire is a chain of straight
segments; L_ij = μ0/4π ∫∫ dl_i·dl_j/r (Neumann) with near pairs split
recursively and the straight-wire self term μ0ℓ/2π·(ln(2ℓ/r) − ¾) (Grover
1946; Ruehli 1972).

| check | result |
|---|---|
| 5 cm ring of 0.3 mm wire vs Kirchhoff μ0R(ln(8R/r) − 7/4) | 3.42106×10⁻⁷ vs 3.42147×10⁻⁷ H (−1.2×10⁻⁴) |
| L_A, L_B, M (path-sense currents), k | 263.44 µH, 310.53 µH, +111.97 µH, k = +0.3915 |
| series inductance, toroidal / dipole sense | 797.91 µH / 350.03 µH |
| reciprocity M_AB vs M_BA (sums in the two orders) | 111.9694754 µH both (2.4×10⁻¹⁶) |
| flux of A's Biot–Savart field through a 1 cm circle in the aperture vs Neumann | 516.443 vs 516.388 nH (1.1×10⁻⁴) |
| the same through a circle round the tube at 15° (3 mm from the bundles) | −1123.68 vs −1090.02 nH (3.1×10⁻²) |
| field energy ∫B²/2μ0 outside a 2.5 mm tube vs ½L_ext·I², (2, 1) circuit, 1 A | 125.02 vs 125.7 nJ (−0.54 %) |
| the same, one Rodin circuit | 695.14 vs 743.79 nJ (−6.54 %) |
| wire resistance (DC), A / B | 1.1414 Ω (19.209 m) / 1.3212 Ω (22.235 m) |

The Rodin path crosses itself on the torus surface, so the exclusion tubes
overlap at the crossings and the grid misses part of their volume: the
−6.54 % is that geometry, not the inductance (the smooth (2, 1) path agrees
to 0.54 %). The tube-circle flux check passes 3 mm from the bundles, where
the quadrature resolves the field less well.

## Part 2 — sine on A, cosine on B

Lumped two-winding circuit in kronspice (L_A, L_B and M from Part 1; valid
far below the first self-resonance), sources 3.383 V peak at 1 kHz through
1 Ω each, so that A alone carries 1.25 A peak (a current of Table II). The
cosine source starts from rest (a series EXP source cancels its value at
t = 0 and fades over two cycles). **Toroidal sense**: B connected so that
the through-the-hole parts of A and B add (field confined in the tube);
**dipole sense**: so that the round-the-axis parts add (a loop-like field
through the hole, the paper's connection).

| drive | \|I_A\| A | \|I_B\| A | P_A W | P_B W | A→B through M, W | ωM·I_A·I_B·sin φ (rms), W | W_m min–max µJ |
|---|---|---|---|---|---|---|---|
| A only | 1.2499 | 0 | 1.6728 | 0 | 0 | 0 | 0.022–205.8 |
| B only | 0 | 1.1156 | 0 | 1.4446 | 0 | 0 | 0.026–193.2 |
| in phase, toroidal | 1.0706 | 0.9521 | 1.2544 | 1.0253 | −0.0270 | −0.0268 | 0.26–405.3 |
| in phase, dipole | 1.4327 | 1.2944 | 2.1860 | 1.9569 | +0.0126 | +0.0122 | 0.07–322.9 |
| sin/cos, toroidal | 1.4624 | 0.8908 | 1.8500 | 1.3613 | **+0.4401** | **+0.4401** | 104.3–260 |
| sin/cos, dipole | 1.0297 | 1.3373 | 1.5904 | 1.6209 | **−0.4546** | **−0.4548** | 117.1–247.2 |

**The power the sources pass to each other.** With mutual inductance M,
the quadrature pair moves power from one source to the other through the
coupling, P = ω·M·I_A·I_B·sin φ (rms): +0.4401 W from A to B in the
toroidal sense, −0.4546 W (B to A) in the dipole sense, the numerical
integral of i_B·M·di_A/dt agreeing with the phasor formula to 4 figures.
In phase, φ ≈ 0 and the exchange is near zero. Every case closes its books
over the last cycle (sources = copper + source resistors to ≤ 1.6×10⁻⁵ W)
and over the run (kronspice ledger mismatch ≤ 9×10⁻⁷ of the energy in).

**The field** (|B| peak over a cycle and rms; ellipticity = minor/major
axis of the ellipse B(t) traces: 0 = it only swings, 1 = it turns at
constant size):

| drive | centre peak / rms mT | in the tube at a pole, peak mT (ellipticity) | largest peak, plane z = 0 / y = 0, mT | mean ⟨B²⟩/2μ0, z = 0 / y = 0, J/m³ | ⟨B²⟩-weighted ellipticity z = 0 / y = 0 |
|---|---|---|---|---|---|
| A only | 1.9339 / 1.3675 | 1.1292 (0.000) | 2.2098 / 2.1852 | 0.1137 / 0.06669 | 0 / 0 |
| B only | 1.0263 / 0.7257 | 1.0076 (0.000) | 2.0435 / 1.9149 | 0.08382 / 0.05338 | 0 / 0 |
| in phase, toroidal | 0.7860 / 0.5558 | 1.8229 (0.002) | 3.5089 / 3.2838 | 0.2606 / 0.08905 | 0.003 / 0.007 |
| in phase, dipole | 3.4075 / 2.4095 | 0.1892 (0.090) | 3.4925 / 3.6579 | 0.05043 / 0.1572 | 0.008 / 0.002 |
| sin/cos, toroidal | 2.6128 / 1.8475 | 1.3403 (0.065) | 2.8371 / 2.8258 | 0.1676 / 0.1249 | 0.090 / 0.094 |
| sin/cos, dipole | 2.3243 / 1.6435 | 1.2433 (0.078) | 2.4479 / 2.7285 | 0.1434 / 0.1214 | 0.108 / 0.098 |

(Planes ±7 cm; points within 3 mm of a bundle excluded.)

- **Connection sense** does what the geometry says: in phase and toroidal,
  the field sits in the tube (1.8229 mT at a pole) and the centre holds
  0.7860 mT; in phase and dipole, the centre holds 3.4075 mT and the tube
  0.1892 mT.
- **Quadrature turns the field vector but does not move the pattern.** The
  toroidal and the dipole parts of the field peak a quarter cycle apart, so
  the field at a point traces an ellipse (⟨B²⟩-weighted ellipticity
  0.090–0.108 against ≤ 0.008 in phase). The pattern itself stands still
  round the axis: the 12-fold Fourier component of |B|² on the circle
  ρ = R0 advances by at most 3×10⁻¹³° over a cycle, round-off. Each winding
  has the 12-fold symmetry of the poles, so their sum has it at every
  instant; nothing travels.
- **Accumulated field.** A linear coil's field is the superposition of
  the drives at each instant; the stored magnetic energy follows the
  currents (W_m above) and nothing accumulates from cycle to cycle. Real
  build-up comes only from resonance:

**Resonance at 100 kHz** (series capacitors; sources 2.677 V peak, 1.25 A
in A alone when tuned):

| case | \|I_A\| A | \|I_B\| A | φ_B − φ_A | mean W_m µJ | P_in W | τ µs | Q = π·f0·τ | ω0L/R |
|---|---|---|---|---|---|---|---|---|
| A only, tuned to L_A | 1.2480 | 0 | — | 102.6 | 1.6686 | 243 | **76.4** | **77.3** |
| A only, not tuned | 0.0162 | 0 | — | 0.0172 | 0.0003 | — | — | — |
| sin/cos, toroidal, each tuned to its own L | 0.0392 | 0.0368 | −90.2° | 0.2065 | 0.0032 | — | — | — |
| sin/cos, dipole, each tuned to its own L | 0.0367 | 0.0391 | +90.2° | 0.2075 | 0.0032 | — | — | — |
| sin/cos, toroidal, not tuned | 0.0204 | 0.0175 | 132.8° | 0.03753 | 0.0008 | — | — | — |
| in phase, toroidal, tuned to L + M | 1.1957 | 1.1958 | −0.1° | 285.2 | 3.1930 | 354 | 111 | — |
| sin/cos, toroidal, tuned to L + M | 0.8594 | 0.8327 | **−0.2°** | 142.5 | 1.5969 | 356 | 112 | — |

- One tuned winding builds up to its Q-limited level over τ = 2L/R
  (Q = 76.4 from the rise against ω0L/R = 77.3), the stored energy
  6,000× the untuned winding's at the same voltage; the ledger closes to
  4.4×10⁻⁵.
- **The quadrature pair cannot be tuned as a pair.** With k = 0.39 the
  coupled windings have two modes, L ± M; tuning each to its own L puts
  the drive between them (currents 0.04 A). Tuned to the L + M mode, the
  pair rings up (Q 111–112), but the currents then **lock in phase**
  (φ_B − φ_A = −0.2° under a sin/cos drive): a resonant pair carries its
  mode, and the quadrature character is gone. Resonance accumulates energy;
  it does not keep the rotating vector.

Figures: `viz/fields_xy.png`, `viz/fields_xz.png` (|B| peak with the
major-axis direction, rms, ⟨B²⟩/2μ0, ellipticity, for A only, in phase,
and sin/cos in both senses), `viz/energy.png` (W_m(t), B(t) in the tube,
the build-ups), `viz/cycle.gif` (one cycle of the sin/cos drive, both
senses, both planes).

## Part 2b — the winding as a line on a knot

Sectional model: every pass of every circuit is a section (12·Nc = 144 per
winding), the inductance matrix from the geometry (pass matrix of the
bundle path plus the leakage between the wires of a bundle, μ0ℓ/2π·
ln(r_b/r_w)), and declared capacitances: wire to wire in a bundle
(two-wire formula, enamel 0.03 mm, εr 3: 188.1 pF/m, circuit c beside c+1),
A–B where the bundles cross (×0.1), 5.54 pF to the surroundings shared by
all nodes. Lossless eigenmodes ω²·C·v = Dᵀ·L⁻¹·D·v:

| winding (12 circuits) | wire per winding | open (far end open) | series, toroidal | series, dipole | plain rings | crossed ring |
|---|---|---|---|---|---|---|
| Rodin (5, 12) | 19.209 m | 1.966, 13.66 MHz | 0.5314, 1.266 MHz | 0.7052, 1.578 MHz | 9.534 (−1.00), 9.814 (−1.00) MHz | **1.27 (−1.00), 6.258 (−1.00)** MHz |
| (2, 1) knot | 5.091 m | 8.039, 50.51 MHz | 3.578, 4.091 MHz | 2.279, 10.34 MHz | 38.85 (−1.00), 40 (−1.00) MHz | **4.144 (−1.00), 25.2 (−1.00)** MHz |

(In brackets: trip parity, the correlation of the node voltages along A
with those along B. For the crossed ring −1 means the wave comes back
inverted after one trip through one winding, so it needs two trips — 720° —
to come back the same way up; for the plain rings it says whether the two
closed windings swing together or opposite.) The crossed ring's lowest
modes are the 720° modes for both windings, at 0.133 (Rodin) and 0.107
((2, 1)) of the plain rings' lowest frequency. `viz/modes.png`.

Self-resonance and the winding's knot: below 1.966 MHz every pass of a
winding carries the same current and the field keeps the 30° symmetry of
the poles; above it the current changes along the wire (Part 3 shows what
that does round the toroid).

## Part 3 — groups of output coils round the toroid, the frequency evolved with them

Drive: sine on A, cosine on B, 3.383 V peak each through 1 Ω. Every output
coil is a round coil of 0.5 mm copper, tuned to the drive frequency by a
series capacitor (declared), its mutual inductance to every section of the
pair computed from the geometry (segment method), so the coils load the
drive in every solution. The sectional networks are solved by a complex
MNA solver written for this example (mutuals of either sign; kronspice's
K element takes 0 < k ≤ 1 and its N element, which takes any inductance
matrix, is transient-only and holds ≤ 64 windings); the solver is twinned
against a kronspice transient of the same network (test `TestACTwin`,
agreement ≤ 3×10⁻³, the trapezoidal step error).

### Where on the toroid, at which frequency: the (f, angle) maps

A reference coil at every angle θ round the axis, 90 frequencies from
1 kHz to 30 MHz × 180 angles (`viz/map.png`): (a) its open-circuit EMF from
the full sectional net (144 sections per winding), exact because an open
coil draws no current; (b) the power into its load when tuned and loaded by
its own wire resistance, the coil in the network. Two coils: **tube**
(20 turns round the tube, radius 25.7 mm, linking the toroidal flux) and
**top** (20 turns, radius 10 mm, lying flat on top of the tube).

| coil, sense | EMF max/min over θ at 1 kHz | at 100 kHz | at 1.86 MHz | at 5.28 MHz | at 18.9 MHz | largest load power (f, θ), share of the sources' power |
|---|---|---|---|---|---|---|
| tube, toroidal | 1.04 (m = 12) | 1.05 (m = 12) | 1.09 (m = 12) | 4.13 (m = 1) | 26.4 (m = 5) | 0.03586 W (5.68 kHz, 30°), 0.084 |
| tube, dipole | 1.05 (m = 12) | 1.05 (m = 12) | 1.43 (m = 2) | 112 (m = 5) | 5.01 (m = 1) | 0.03851 W (4.51 kHz, 210°), 0.064 |
| top, toroidal | 1.02 (m = 12) | 1.05 (m = 12) | 2.51 (m = 5) | 106 (m = 2) | 64.9 (m = 7) | 0.03318 W (18.9 MHz, 120°), 0.309 |
| top, dipole | 1.07 (m = 12) | 1.04 (m = 12) | 1.24 (m = 3) | 12.5 (m = 7) | 57 (m = 7) | 0.0331 W (18.9 MHz, 120°), 0.312 |

(m: the strongest angular harmonic of the EMF round the toroid.)

**Why.** Below the first self-resonance (1.966 MHz) every pass of a winding
carries the same current, so the field has the 12-fold symmetry of the
poles and a coil sees the same EMF within 2–7 % at every angle (harmonic
m = 12, the 30° pole spacing): **where** a group sits round the toroid does
not matter, and the sin/cos pair gives no angle-dependent phase because the
pattern does not turn (Part 2). Above it the current varies along the wire
as a standing wave, and since the wire goes round the axis 5 times per
circuit, advancing 150° per pass, the field round the toroid takes the low
harmonics of that wave (m = 1, 2, 5, 7): the EMF then changes by up to 112×
with angle, and the best angles move with frequency. That is where the
coupling between the frequency and the group positions lives. These modes
depend on the declared capacitances, and the winding's high input impedance
there leaves less power to take: for the tube coil the largest load power
is at 4.5–5.7 kHz, for the top coil at 18.9 MHz with comparable power
(0.033 W).

### The genetic search over frequency and groups

`evo` (value genome, seeded, 64 × 50): drive frequency 1 kHz–30 MHz;
number of groups 1–8; first angle θ0 and spacing Δθ; coil tilt β (0: linking
the toroidal flux, 90°: vertical axis); position ρ (0–9 cm), z (±7 cm);
coil radius 3–60 mm; turns 1–200; combination (each group into its own
bridge, all in series into one bridge, all in parallel into one bridge);
load R_L 0.1 Ω–10 kΩ; connection sense. Constraints: every coil ≥ 2 mm
clear of the bundles and outside the former, coils clear of each other,
each coil shorter than λ/20 at f and its tuning capacitor ≥ 30 pF (so it is
a lumped coil, not at its own self-resonance). Fitness: power into R_L (the
AC solution, coils loading the drive) times a bridge factor
(V/(V + 2·0.3 V))², which steers the search away from voltages too low to
rectify. Searches use the net with 2 sections per circuit, which the
kronspice transient can hold; each winner's power is the same on 4, 6 and
12 sections per circuit to 4 figures.

| | free search | distributed groups (≥ 2 coils, radius ≤ 20 mm) |
|---|---|---|
| best of generation 0 → evolved (fitness) | 1.664×10⁻⁷ → 0.181 W | 1.265×10⁻⁷ → 0.1172 W |
| drive frequency | **4,134 Hz** | **29.38 kHz** |
| connection sense | dipole | dipole |
| groups | **1**: one coil round the whole coil, coaxial (ρ = 0.3 mm), radius 53.3 mm, 200 turns, horizontal (β = 90°), 14.4 mm above the mid-plane | **2** at 79.7° and 251.1°, radius 20 mm, 200 turns, β = 50.1°, ρ = 8.2 mm, z = −34.8 mm (under the aperture), in series into one bridge |
| coil L, R, tuning C | 7677 µH, 5.733 Ω, 193.1 nF | 1894 µH, 2.15 Ω, 15.5 nF |
| load voltage amplitude, R_L | 3.322 V, 21.9 Ω | 1.92 V, 9.13 Ω |
| AC: sources = load + copper + source R | 0.4708 + 0.2793 = 0.2522 + 0.307 + 0.1909 W | 0.2061 + 0.1462 = 0.2018 + 0.126 + 0.02447 W |
| AC efficiency (load / sources) | 0.3363 | 0.5730 |
| reaction on the drives: sources without → with the groups | 0.5877 + 0.03219 W → 0.4708 + 0.2793 W | 0.058 − 0.04325 W → 0.2061 + 0.1462 W |
| the same groups combined separate / series / parallel | one coil: 0.2522 W each way | 0.1598 / **0.2018** / 0.1053 W |
| the layout at other frequencies (groups retuned), P_load | 1 kHz 0.0889, 2 kHz 0.194, 4 kHz 0.252, 10 kHz 0.17, 30 kHz 0.0355, 100 kHz 0.00354 W; ≥ 1 MHz infeasible | 1 kHz 0.0531, 2 kHz 0.12, 4 kHz 0.175, 10 kHz 0.198, 30 and 100 kHz 0.202 W; ≥ 1 MHz infeasible |

What the searches found:

- **Left free, the search builds a transformer**: one large coil coaxial
  with the toroid, linking the dipole flux through the hole, in the dipole
  sense, near the frequency where the windings' reactance is comparable to
  their resistance (4.1 kHz). Distributed groups round the tube never won:
  below the self-resonance the field is the same at every angle, so many
  small coils link less flux than one coil round the whole.
- **Made to use several small coils**, the search put two under the
  aperture (z = −34.8 mm), where the dipole flux leaves, nearly opposite
  each other, in series into one bridge. Their currents are in phase
  (138.3° both); in parallel they give 0.1053 W and on separate bridges
  0.1598 W against 0.2018 W in series, so the combination matters, and in
  the chosen series connection no group cancels another.
- Neither search went above the self-resonance: there the angle matters
  (maps above), but the winding's input impedance leaves less power to
  take, and the coils would need tuning capacitors below 30 pF.

### The winners charging the bank (kronspice transient)

The winners run in kronspice: the reduced net's 48 sections and the output
coils as the windings of one N element (full inductance matrix), node
capacitances, the soft-started sine and cosine sources, tuning capacitors,
Schottky bridges (is 3 µA, n 1.05, rs 0.03 Ω, cjo 100 pF) into one bank
(sized to fill within the run), chokes from each coil terminal to ground
as the DC reference (reactance 10⁴ × the load, no loss).

| | free search | distributed groups |
|---|---|---|
| run | 400 cycles, 96.77 ms | 400 cycles, 13.62 ms |
| bank | 2200 µF | 680 µF |
| bank voltage during the run | 1.097, 1.917, 2.662, **3.367 V** at 14.8, 29.2, 49.5, 96.8 ms | 0.311, 0.9056, 1.47, **1.86 V** at 1.8, 4.1, 7.39, 13.6 ms |
| bank energy at the end | 12.47 mJ | 1.176 mJ |
| energy from source A / B | 47.04 / 23.65 mJ | 2.393 / 1.266 mJ |
| **charging efficiency** (bank / sources) | **0.1764** | **0.3214** |
| copper (windings and output coils) | 45.08 % | 34.76 % |
| source resistances | 28.96 % | 8.10 % |
| diodes | 6.47 % | 21.65 % |
| left in coils and capacitors | 1.278 mJ | 0.1178 mJ |
| books: sources − (bank + losses + stored) | 35.6 µJ | 4.95 µJ |
| kronspice ledger integration mismatch | 35.8 µJ (4.0×10⁻⁴ of the energy in) | 4.97 µJ (1.1×10⁻³) |
| radiation of the output coil, 20π²N²(ka)⁴ | 3.59×10⁻¹⁵ Ω against 5.73 Ω of wire | 1.82×10⁻¹³ Ω against 2.15 Ω |

The books close to the transient's integration mismatch (the trapezoidal
error at 60 steps per cycle): every joule in the bank came from the two
sources, and the bank holds 18–32 % of what they gave. Radiation is
negligible. The output coils load the drive — the sources' power changes
when they are added — as a transformer's secondary loads its primary.
`viz/bank.png`, the coils in `viz/geometry.png`.

## Status of H43

| claim (hobbyist and promotional sources) | what the runs show | standard prediction |
|---|---|---|
| a "vortex" field | under sin/cos drive the field vector turns at each point (ellipticity 0.09–0.11), but the pattern stands still round the axis (12-fold advance ≤ 3×10⁻¹³° per cycle) | the superposition of two coils' Biot–Savart fields a quarter cycle apart |
| the field "accumulates" | only by resonance, to the Q-limited level (Q 76.4 against ω0L/R 77.3), paid by the source; the quadrature pair cannot be tuned as a pair (k = 0.39), and tuned to one mode its currents lock in phase | a series RLC tank |
| energy beyond the input (over-unity) | none: every ledger closes (lumped drives ≤ 9×10⁻⁷ at 1 kHz and ≤ 6.7×10⁻⁵ in the 100 kHz build-ups, AC books ≤ 5×10⁻¹³ W, transient ≤ 1.1×10⁻³ of the energy in, its integration error); the bank holds 0.18–0.32 of the source energy; output ≤ input in every case | a three-winding transformer: output ≤ input |
| a special role of the {12/5} (5, 12) path | it sets inductance, coupling and frequencies; the 720° modes appear for the series counter-wound loop with the (5, 12) and the (2, 1) path alike | a property of the connection, not of the knot |
| "vortex-based mathematics" (1-2-4-8-7-5) | no equation is given, so there is nothing to compute | — |

**Status: bounded.** Computed under the standard model with closed books,
the Rodin pair is a pair of coupled inductors (k = 0.39) and, near 2 MHz, a
line laid on a (5, 12) knot; the model reproduces the paper's centre field
with one fitted turn count at 1.60 and 1.25 A. What would move the entry:
a bench measurement in which the energy delivered to a load or bank exceeds
the energy measured at the drive terminals, with every input counted.
Bench test: the Part 3 winner (or the paper's coil with a coaxial
secondary), a power analyser on both drive channels, the bank energy ½CV²;
the standard prediction is 0.18–0.32 of the input in the bank with these
parts.

## Tests

| test | checks |
|---|---|
| `TestCentreFieldTableII` | the path closes and goes 5 times round the axis; centre field within 6 % of Table II at 1.60 and 1.25 A; least squares gives Nc = 12 |
| `TestInductanceChecks` | ring vs Kirchhoff (10⁻³); flux of Biot–Savart vs Neumann (10⁻³); reciprocity; 0 < k < 1 |
| `TestEnergyCheck` | field energy outside a tube vs ½L_ext·I² for the (2, 1) circuit (3 %) |
| `TestQuadratureDrive` | both senses: per-source books (10⁻³), exchange integral = ωM·I_A·I_B·sin φ (2 %), ledger, the 12-fold pattern does not advance |
| `TestResonantBuildUp` | tuned winding: Q from the rise vs ω0L/R (10 %), 1.25 A at resonance |
| `TestACTwin` | own complex MNA vs kronspice N-element transient of the same net (3×10⁻³); AC books (10⁻⁹) |
| `TestCrossedRing720` | the crossed ring's first mode has trip parity < −0.9, Rodin and (2, 1) |
| `TestGroupsLedger` | two output coils: AC books; kronspice transient with bridges and bank: books = ledger mismatch, 0 < bank < sources |

`go test -short ./examples/rodin_resonance` (about 20 s); in the
kronspice-tools image `go test -race -short` as well.

## Files

| file | contents |
|---|---|
| `main.go` | the four parts, the printed results, the CSV files in viz/ |
| `geom.go` | the taut-wire {12/5} Rodin path, the (2, 1) knot, bundles, the pair |
| `induct.go` | segment (partial-inductance) method; Biot–Savart field of a polyline |
| `model.go` | pass matrix, winding L and M, bundle leakage, reciprocity, flux and field-energy checks |
| `drive.go` | lumped two-winding kronspice circuit, soft-started sources, phasors, power exchange, field statistics, ring-travel test |
| `line.go` | sectional network (open, series, plain and crossed rings, driven pair), eigenmodes, trip parity |
| `linalg.go` | Cholesky, Jacobi eigenvalues, generalized symmetric eigenproblem (stdlib only) |
| `acsolve.go` | complex MNA of the sectional net with sources, output coils and loads |
| `groups.go` | output-coil geometry, mutuals, clearance, EMF and power maps |
| `part3.go` | genome, decoding, constraints, bridge factor, searches |
| `tranverify.go` | kronspice transient with the N element, Schottky bridges and bank; energy split per element |
| `rodin_resonance_test.go` | the tests above |
| `run_output.txt` | output of the full run |
| `viz/plot_all.py`, `viz/plot_cycle.py` | figures and the GIF from the CSV files |
| `viz/make_report.py` | the PDF report `rodin_resonance_report.pdf` |

## Approximations (declared)

- The paper gives no turn count or wire gauge: 0.6 mm copper, Nc = 12
  circuits per winding fitted to the measured centre field; each winding a
  round bundle on one path (fill 0.7); B one bundle pitch outside A.
- Taut-wire path (chords on the torus surface); the real winding's
  curvature and crossings are idealized (the field-energy check shows the
  crossings' effect, −6.5 %).
- Skin effect of an isolated round wire, no proximity effect; air core.
- Sectional model: sections equipotential inside; declared capacitances
  (two-wire formula between neighbouring wires of a bundle, A–B contact
  ×0.1, a sphere's capacitance to the surroundings); results above the
  first self-resonance depend on them.
- Searches on 2 sections per circuit (checked on 4, 6 and 12); fitness
  with a resistive load and a smooth bridge factor in place of the
  rectifier, replaced by real diodes in the kronspice transient.
- Radiation estimated (small-loop R_rad), not modelled.

## Sources

- A. D. Ramírez Galindo, H. A. Pérez Olivas, G. Basurto Islas,
  A. Hernández Rayas, F. González López, T. Córdova Fraga, "Vortex magnetic
  induction: Mathematical, geometric and experimental characterization",
  *Rev. Mex. Fís.* 70, 030901 (2024), doi:10.31349/RevMexFis.70.030901
  (geometry, Table II).
- M. Rodin, promotional material on the Rodin coil and "vortex-based
  mathematics" (academia.edu, forums): the hypothesis source, with no
  measurement.
- J. G. Williamson, M. B. van der Mark, *Ann. Fond. Louis de Broglie* 22,
  133 (1997), and [photon_electron](../photon_electron) (H23): the (2, 1)
  path, the crossed ring line, the circuit mirror.
- F. W. Grover, *Inductance Calculations*, 1946; A. E. Ruehli, *IBM J.
  Res. Dev.* 16, 470 (1972) (partial inductance); J. C. Maxwell,
  *Treatise* Vol. II art. 701; J. D. Jackson, *Classical Electrodynamics*,
  3rd ed., §5.2–5.3 (Biot–Savart).
- C.-W. Ho, A. Ruehli, P. Brennan, "The modified nodal approach to network
  analysis", *IEEE Trans. Circuits Syst.*, 1975; G. H. Golub,
  C. F. Van Loan, *Matrix Computations*, §8.4 (Jacobi).
- C. A. Balanis, *Antenna Theory* (small-loop radiation resistance);
  F. E. Terman, *Radio Engineers' Handbook*, 1943 (skin effect).
