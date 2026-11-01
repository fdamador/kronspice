# Vassilatos, *Secrets of Cold War Technology* — review for the pulse work (October 2026)

Reviewed at the user's request, with this premise in view: **a fast pulse
acts on space the way a fast slap acts on water — above some rate the
medium answers like a solid** — and the question whether Tesla's impulses
bring in second- and third-order rates of change that would call for new
Maxwell-equation concepts. This note records what the book actually says,
page by page, what established physics says beside it, and where each item
already lives (or should live) in kronspice.

## Files

| File | Pages | Text layer |
|---|---|---|
| `Secrets-of-Cold-War-Technology-by-Gerry-Vassilatos.pdf` (52.7 MB) | 324 | Yes (OCR, ~918 k characters; pp. 25–27 badly garbled). **Use this copy.** |
| `gerry-vassilatos-secrets-of-cold-war-technology.pdf` (44.3 MB) | 163 | None (image only) |

Page numbers below are the **printed** page numbers (the PDF page index is
one less: PDF page 34 carries printed page 35). All of the Tesla material
is Chapter 1 (pp. 22–90). Quotations were read from the OCR layer; check
against the page image before quoting them elsewhere.

## 1. The slap-on-water passage is in the book — as a gloss on line surges

The user's analogy is Vassilatos's own (p. 27), used for the switch-closure
surge on long dynamo lines:

> "the combined resistance of such systems seemed to influence the charge
> carriers before they were able to move away from the dynamo terminals.
> Like slapping water with a rapid hand, the [water] seemed solid. So also
> it was with the electrical force, charges running up against a seemingly
> solid wall. But the effect lasted only as long as the impact. Until
> current carriers had actually 'caught up' with the applied electrical
> field, the charges sprang from the line in all directions."

Vassilatos then has Tesla ask "why it was possible for electrostatic fields
to move more quickly than the actual charges themselves" (p. 27) and work
the effect up with "time, force, and resistance" as the variables (p. 28).

**Standard reading.** This is the transient on a transmission line. A step
applied to a line does not see the line's DC resistance; it sees the
natural (surge) impedance `Z0 = √(L/C)`, and the front carries `V = Z0·ΔI`
until the reflection returns after the round-trip time `2ℓ/v`. The
"super-concentration several orders of magnitude" above the dynamo voltage
(p. 27) is voltage doubling at open ends and build-up by reflections. The
field moving ahead of the charges is the wavefront moving at the line
velocity while the drift of the carriers is millimetres per second. All of
this is in Steinmetz (1911): natural impedance Lecture VI §25, travelling
waves and impulses Lecture VIII, step-up at a transition point Lecture IX
§41–42 — already the basis of §3e step 3.

The water analogy is exact, not loose, and it runs the other way too: a
hand slapping water feels a pressure `Δp = ρ·c·Δv` (Joukowsky's water-hammer
relation), the acoustic impedance `ρc = √(ρ·K)` playing the part of
`Z0 = √(L/C)`; the "solid" feel lasts until the pressure relief returns
from the nearest boundary, the water's `2ℓ/c`. Then the water flows and
the line conducts at its DC resistance. The dimensionless rule in both
cases is the Deborah number `De = τ / t_applied`: `De ≫ 1` solid-like,
`De ≪ 1` fluid-like.

So the premise holds, as standard physics, in two places that kronspice
already models:

| Where the "stiffening" happens | Characteristic time | kronspice |
|---|---|---|
| A conductor or coil long compared with the edge: lumped wire becomes a line, DC resistance becomes `Z0` | round trip `2ℓ/v` (and the `λ/10` lumped limit, `em/limits.go`) | `em/tline.go`, `em/helix.go`, Steinmetz line validation case (§3g), `examples/pulsed_bifilar` |
| A material whose polarization cannot follow: water `ε 78.6 → 5.4` | Debye `τ` = 8.3 ps (first order), Lorentz `ω0, γ` (second order) | `em/dispersion.go` (`Debye`, `LorentzTerm`), `examples/rate_of_change` (`Slap`) |

Maxwell's own viscoelastic model (1867: spring and dashpot in series, a
fluid that is elastic for stresses applied faster than its relaxation time
`τ = η/G`) is the mechanical twin of the Debye `R–C` branch; Maxwell's name
is already on the "water acts as a solid when struck fast" idea.

**Where the premise does not hold — the vacuum.** Empty space has
`Z0 = √(μ0/ε0) = 377 Ω` at every rate; it has no `τ` and no `ℓ`, so it
shows no threshold. Maxwell's vacuum is already the "solid" (his 1861
vortex-ether model was an elastic medium carrying transverse waves). A rate
threshold *in space itself* would be vacuum dispersion, which is bounded by
gamma-ray-burst timing: photons from GRB 090510 spanning MeV to tens of GeV
arrived within about a second after ~7×10⁹ light-years, bounding any linear
energy dependence of the vacuum light speed to above about the Planck
energy (Abdo et al., *Nature* 462, 331, 2009; value to be entered from the
paper). That is the standard prediction to set beside any claim that short
pulses make space answer differently.

## 2. What the book says about Tesla's impulses — claim by claim

| p. | Claim (Vassilatos) | Standard counterpart | kronspice |
|---|---|---|---|
| 25 | Maxwell considered longitudinal electric waves and "rejected" them; Hertz's waves were transverse | Maxwell's free-space equations have no longitudinal propagating mode; longitudinal **near** fields and TM guided modes exist within the standard theory (§3k tier 1, arXiv:1606.01738) | `em/dipole.go` near field; §3k tier 1 cases |
| 32 | The effect lives in "IMPULSE time"; charging must be a single impulse, no reversal; a capacitor between disrupter and dynamo strengthens it and protects the windings | Spectrum of a one-sign pulse has a DC component and a different harmonic set from the same-width alternating pulse; a one-sign field charges insulated bodies by rectification | §3g check 5 (polarity); `device/arc.go` |
| 32 | Stronger with higher voltage, faster make–break, shorter closure | Higher `dV/dt` ⇒ larger displacement current and radiated field (∝ second derivative, §3) | `em/fdtd.go` |
| 32 | Arc "blown out" by a crosswise permanent magnet; carbon–copper horns "rectify" | Magnetic blow-out gap: quench by stretching the arc; dissimilar electrodes give asymmetric breakdown | `device/arc.go` (quench parameter), §3g source |
| 35 | Duration thresholds: trains of pulses > 0.1 ms "pain and mechanical pressures", wires explode; ≤ 100 µs "pain and physical movements ceased"; 1.0 µs "strong physiological heat"; shorter "spontaneous illuminations… white light"; shorter still "cool room penetrating breezes" | A smooth function of the pulse spectrum (electrostatic force ∝ `V²`, ohmic heating ∝ `∫I²dt`, corona onset); no step is predicted at any duration | §3g check 4 — these are **the only numbers in the book**, they are the search targets, and no instrument, voltage, current or distance is given for any of them |
| 35 | 100 µs impulses "passed through all matter" | Skin depth: 100 µs ⇒ copper δ ≈ 1 mm, so low-frequency **magnetic** fields do pass thin copper; electric fields do not | §3g check 6 (shielding factor) |
| 36–37 | A 2-ft single-layer helix near the disrupter: discharge rises up the coil, 10,000 V per inch of axial length (10 kV → 240 kV over 24 in), **zero current measured at the crown**, "electrostatic induction laws" | Quarter-wave resonant secondary: open end is a voltage maximum and a **current node** (Tesla 593,138; §3e step 3). Step-up = resonant Q times the `√(Z₂/Z₁)` junction step-up (Steinmetz §41–42) | `examples/pulsed_bifilar`: 100 V → 7.4 kV when a pulse-train harmonic lands on `f₁`; `em/helix.go` |
| 37 | Output voltage rises with the "resistance of turns" | Thinner wire ⇒ higher `L` per length ⇒ higher `Z0 = √(L/C)`; the step-up follows `Z0`, not ohmic `R` (which lowers Q) | §3e genome: wire diameter |
| 37–38, 45 | Primary and secondary must have equal copper **mass**; each transformer "conducts a specific impulse duration" and is tuned by the arc distance | Equal mass at fixed wire sizes is equal wire length — a length match, i.e. matched electrical length. Arc distance sets breakdown voltage and so the repetition rate; §3e: pulses repeated at `2ℓ/v` build up, others do not | `examples/pulsed_bifilar` (build-up only on a harmonic of `f₁`) |
| 39 | Lamps across a heavy U-shaped copper bus-bar **shorted** by the bar light brilliantly: "another energetic current" | The bar's inductance (~1 µH per metre) drops `L·di/dt` across the "short": 1 µH × 100 A / 10 ns = 10 kV. A short circuit is a short only at DC | Runnable today: `.TRAN` of a PULSE edge into `L_bus ∥ lamp` |
| 41 | Dollard measured the surface-riding component travelling faster than light by "line interferometry" (1985) | No published data; the helical-line phase velocity measured **along the wire** exceeds `c` (up to ~2c in Drude's data, Knight §1.2) while energy moves along the axis below `c` — already in §3e | `em/helix.go` |
| 44–45 | Effects "grow with time", "additive and accumulative" because there are no reversals | Resonant build-up (amplitude grows until loss per cycle = input per cycle); DC charge accumulation on insulated bodies under one-sign pulses | §3e step 3 energy rule; ledger |
| 52–53 | Radiant matter is "incompressible": an impulse reaches all points of the path "instantly", 300,001 km as fast as 1 m | Retardation at `c`; GRB bound above | §3g check 1 (arrival times at spaced probes) |
| 54 | Impulses act as "tiny hammers"; electrode matter evaporates in vacuum globes | Field emission and sputtering under high-voltage pulses; cathode erosion | §3i multipactor surface models |
| 59 | With "carefully adjusted" impulse durations a spray becomes a needle that pierces plates | No numbers | — |
| 77, 81–82 | Colorado Springs / Wardenclyffe: 150,000 pulses per second, each ≈ 0.0062 ms, 100,000,000 V, 50 kV DC dynamo, rotary then magnetic disrupter, salt-water-bottle capacitor bank; the "rest period" is "difficult to ascertain" | 1/150,000 s = 6.67 µs, so a 6.2 µs pulse leaves a 0.47 µs rest: the "duration" is 93 % of the period. It reads as a half-cycle of a ~150 kHz oscillation back-computed, not as a measured impulse length. Tesla's Colorado Springs notes record oscillation frequencies in this range | Not a duration datum; treat as a frequency |
| 84 | Aerial terminal current "some 1100 amperes" | Plausible RF current in a resonant secondary of that size; no deviation implied | — |
| 85–86 | "Stationary waves" (incompressible-fluid fountains, unidirectional) vs "standing waves" (alternating electron waves); the sink-and-drop fountain analogy | Standing waves on the earth–ionosphere system and on the coil; a fountain is a free-surface impulse response, not a new wave type | `examples/ring_towers` (standing / travelling / cancel) |
| 88 | Teslian "frequency" = pulses per second; "resonance" = flow "with little or no resistance" | Repetition rate; resonance = reactive cancellation (`Z` real and minimum) | — |

**Verdict on the book as a source (unchanged from CLAUDE.md §3g).** It
gives no equations, no circuit values, no instrument readings and no
primary citation for any pulse measurement. Its one set of numbers (the
duration thresholds, p. 35) is uninstrumented; its Wardenclyffe duration is
a back-computed half-period. Every phenomenon it describes has a standard
quasi-static or transmission-line account that kronspice either already
reproduces or has as a planned validation case. It remains the source of
the H1 wording and of the §3g check-4 search targets, and nothing more.
Primary sources it points to that are **not yet in `reference/`**: Tesla,
"On the Dissipation of the Electrical Energy of the Hertz Resonator,"
*Electrical Engineer*, 21 Dec 1892 (p. 31); the Colorado Springs Notes
(1899–1900); L. Anderson's Wardenclyffe court-transcript treatise (p. 88).

## 3. Where second- and third-order rates of change already enter Maxwell's theory

Nothing in the book adds a derivative to the field equations. But the
question is a good one, because the standard theory already has places
where higher time derivatives are the physics, and kronspice covers only
some of them:

| Order | Where it enters | Equation | Standard or hypothetical | kronspice |
|---|---|---|---|---|
| 1st | Constitutive laws of C and L; displacement current | `i = C dv/dt`, `v = L di/dt`, `∂𝔇/∂t` | Standard (§3c table) | all of `device/`, `em/maxwell.go` |
| 1st | Debye relaxation (water's "slap" response) | `τ dP/dt + P = ε0 Δε E` | Standard | `em/dispersion.go` |
| 2nd | Lorentz oscillator (resonant polarization) | `d²P/dt² + γ dP/dt + ω0² P = ε0 Δε ω0² E` | Standard | `em/dispersion.go` (`LorentzTerm`) |
| 2nd | **Radiation**: the far field of a source goes with the second derivative of its dipole moment; Larmor power `P = μ0 \|d²p/dt²\|² / (6πc)` | from `∇F = J` | Standard | `em/dipole.go` (`LarmorPower`), `em/fdtd.go` |
| 2nd | Wave equation itself (`□A = μ0 J`) | second order in `t` and `x` | Standard | `em/fdtd.go`, `em/tlm.go` |
| **3rd** | **Radiation reaction** (self-force): the force a charge feels from its own radiation, `F = (q²/(6πε0 c³)) d³x/dt³` (Abraham–Lorentz); well-posed form by Landau & Lifshitz (*Classical Theory of Fields* §75–76), which replaces `d³x/dt³` by the derivative of the external Lorentz force | Standard, with known pathologies (runaway, pre-acceleration) in the raw form | **Not modeled.** It is the one legitimate "third order of rate of change" in Maxwell–Lorentz electrodynamics. It belongs in the energy ledger: the work done against the reaction force over a pulse must equal the Larmor energy radiated — a closed test with the existing `LarmorPower` |
| 4th | Bopp–Podolsky `(1 + a²□)□A = μ0 J` | hypothetical (H19) | `em/podolsky.go`, `examples/rate_of_change` |

Two points follow for the engine:

1. **No circuit element needs a second or third derivative in its own
   constitutive law.** Second-order behaviour (resonance) arises from
   first-order `L` and `C` in combination, exactly as in Maxwell. A
   "jerk" element would not be a new circuit law; the only place `d³/dt³`
   is physics is the radiation-reaction correction on a radiating branch,
   and that enters as a cited **radiation-resistance edge** (§3a.1 item
   5) whose value follows from the Larmor formula, not as a new device.
2. **Fast pulses raise the order that matters, they do not add equations.**
   A slow drive sees `R` (zeroth order); a faster one sees `L` and `C`
   (first); faster still, the wave equation and radiation (second); and at
   the extreme of a single charge, its own radiation reaction (third). Each
   is already in `∇F = J` plus the Lorentz force. The engine's job is to
   check which order a given pulse reaches (`em/limits.go`, `SourceBandwidth`)
   and switch model accordingly — lumped → line → field — which is the
   §3g pipeline.

### 3.1 Momentum: what the higher orders exchange with space

The derivative orders above are best read as a momentum story, because
the field momentum is already there at zeroth order and the higher orders
only move it.

- **Static fields hold momentum.** `g = ε0 E×B = S/c²` (Maxwell stress
  tensor; Poynting 1884). Switching off a solenoid around a charged disk
  sets the disk turning with the angular momentum the field had held
  (Feynman disk). Ledger: `em/momentum.go`.
- **Second order radiates it.** Energy `W` that leaves as radiation carries
  `W/c`. A symmetric dipole recoils not at all; an asymmetric source
  (one-sign pulse on an unbalanced antenna, dipole–quadrupole
  interference, a reflector) recoils with thrust ≤ `P/c` — the photon-rocket
  bound already stated against H18 (3.3 mN per MW). Lebedev 1901 and
  Nichols–Hull 1901 measured the receiving side (radiation pressure).
- **Third order is the recoil on the charge.** The Abraham–Lorentz
  self-force is the reaction to the momentum being radiated. Its two
  ledger tests are one: work against the self-force = Larmor energy
  radiated; impulse of the self-force = momentum radiated.
- **The reverse direction** is the Lorentz force density: incoming field
  momentum lands on charges (receiving antennas, radiation pressure).
  Locally, `∂g/∂t + ∇·T = −f` with `f = ρE + J×B`.

What this does not allow: net momentum on a closed source without
something leaving. By the local law a source gains momentum only by
radiating fields that carry the opposite momentum away (bounded by `P/c`)
or by pushing against an external field whose source takes the reaction
(§3f, `doc/em-push-pull.md`). A faster pulse raises the power a given
charge motion can radiate (Larmor ∝ acceleration²); it opens no third
channel. A one-sign current pulse radiates a far field whose time integral
is zero, because that field tracks `d²p/dt²` and its area is the net change
of current, which is zero once the pulse has ended.

### 3.2 Energy: momentum's rates of change are forces, not stores

Derivatives of momentum are rates — force, "yank", "tug" (N, N/s, N/s²) —
never joules; energy comes from force only through displacement,
`W = ∫F·dx`. The sense in which momentum in space does mean energy in
space is sharper than that:

- **One object, two axes.** Field energy density `u = ½(ε0E² + B²/μ0)` and
  momentum density `g = ε0 E×B` are the time–time and time–space parts of
  one stress–energy tensor; in the spacetime algebra `T(γ0) = ½ F γ0 F̃`
  has time component `u` and space component `c·g`. A boost mixes them;
  for radiation `u = c|g|`. Both exist at zeroth order, before any
  derivative is taken.
- **One conservation law, two components.** Its divergence gives
  Poynting's theorem `∂u/∂t + ∇·S = −J·E` (energy ledger) and the momentum
  law `∂g/∂t + ∇·T = −(ρE + J×B)` (momentum ledger). In empty space both
  right sides vanish: energy and momentum are transported, never created.
  The bridge between the two ledgers is velocity: power to matter is `f·v`.
- **Where momentum derivatives do put energy into space.** Larmor in terms
  of the force on a charge: `P = q²/(6πε0 c³ m²)·(dp/dt)²`. The
  third-order self-force is the term whose work `∫F_rad·v dt` equals minus
  that energy. The energy appears in space; the agent that produced
  `dp/dt` pays for it.

Test (one case closes both ledgers): a pulsed radiating element with the
Landau–Lifshitz self-force — work against the self-force = Larmor energy;
impulse of the self-force = radiated momentum (`W/c` for a beam); and
`u = c|g|` in the radiated pulse from `em/fdtd.go`.

## 4. The premise as a register entry, with its standard twin

To keep the user's idea in the project the way §1 of CLAUDE.md requires
(stated, tested, never dropped), it is proposed as:

> **H20 — Rate-dependent stiffening of the medium under fast one-sign
> pulses:** above some impulse rate (Vassilatos's thresholds, p. 35:
> 100 µs, 1 µs and shorter) space or the surrounding medium answers a pulse
> as a solid rather than a fluid, changing propagation speed, field ratio
> or energy transfer in a way the pulse spectrum does not explain.
> **Standard prediction:** in conductors and coils the stiffening is the
> surge impedance `Z0 = √(L/C)` with threshold `2ℓ/v`; in materials it is
> Debye/Lorentz dispersion with threshold `τ`; in vacuum there is none
> (`Z0 = 377 Ω` at all rates; GRB dispersion bound). **Tested in:**
> `examples/rate_of_change` (material), `examples/pulsed_bifilar` and the
> Steinmetz line case (line), `examples/pulse_duration` checks 1 and 4
> (space). **Status:** untested; the two non-vacuum parts are reproduced
> standard physics.

This keeps H1 as it is (the toroidal-field premise) and gives the
slap-on-water idea its own line, its own standard prediction and a bench
test that can actually be run: the same pulse train into (a) a long bus
bar, (b) a water-filled cell, (c) a vacuum gap with spaced E/B probes, each
duration from 100 µs down, looking for any change that `Z0`, `ε(ω)` and
`c` together do not predict.

## 5. Simulation cases kronspice can run now from this reading

1. **Shorted-bus lamps (p. 39).** `L_bus` (from `em/coil.go` straight-wire
   or `em/inductance.go`) in parallel with a lamp resistance, driven by a
   PULSE edge of 1–100 ns through a source impedance; report the lamp
   voltage and the ledger: all lamp energy from the source, none from the
   "short". A clean demonstration that the premise is `L·di/dt`.
2. **Zero current at the crown (p. 37).** Already in `examples/pulsed_bifilar`
   (`secondary.go` pickup at the middle); add the end-of-coil current
   probe and show the node.
3. **Duration sweep across the three thresholds (p. 35).** `examples/
   pulse_duration` (planned, M9o.1): the same source energy at 100 µs,
   1 µs, 100 ns, 10 ns into the three loads of §4, standard prediction at
   every probe, so a bench measurement can be read against it.
4. **Radiation reaction on the ledger (§3).** New standard model in `em/`:
   Landau–Lifshitz force on a pulsed current element; test that the work
   against it equals the Larmor energy radiated, to the integration
   tolerance. This is the only "third-order" item the review turns up and
   it is established physics.

## 6. Proposed CLAUDE.md changes (for review, not applied)

- §3g pulse-duration experiment, check 4: add the page: "Vassilatos
  p. 35" for the thresholds, and note that the Wardenclyffe 0.0062 ms
  (pp. 77, 81) is 93 % of the 150 kHz period and is not a duration datum.
- §6 references: give the Vassilatos entry the usable copy's filename
  (`Secrets-of-Cold-War-Technology-by-Gerry-Vassilatos.pdf`, OCR layer)
  and list the three primary sources above as still to be added.
- Hypotheses table: add H20 as in §4 of this note.
- §3c "How each circuit law traces back to Maxwell": add a row for the
  radiation-reaction (self-force) correction, third order, standard, with
  its ledger test — once `em/` has it.
