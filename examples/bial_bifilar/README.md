# Bi–Al capacitor element with a bifilar coil

The project's first experiment (milestone M5a; CLAUDE.md §3e, "As
implemented"). A bismuth–aluminium capacitor element and a bifilar coil
(Tesla, US 512,340) in parallel, self-driven: the element is the only
energy source, and the ledger must trace every joule stored or dissipated
to its declared port. Evolution trades the tank's quality factor Q at
resonance against the energy stored per volume. The coil model here is the
lumped first approximation; the pulse-driven, section-by-section coils are
in [pulsed_bifilar](../pulsed_bifilar), and the same coil in orbit is in
[coiled_tether](../coiled_tether). The element is also the bench for the
temporal-field hypothesis H4 (§3f), whose standard prediction is the
Seebeck EMF this example computes.

```bash
go run ./examples/bial_bifilar                 # ≈ 1 min: report in runs/<time>-bial_bifilar/
go run ./examples/bial_bifilar -quick          # a few seconds
go run ./examples/bial_bifilar -deltaT 300 -emf-factor 0.6 -line Al
go run ./examples/bial_bifilar -netlists doc/experiments   # SPICE netlists + schematics of the best designs
go test ./examples/bial_bifilar
KRONSPICE_EXPORT=1 go test ./examples/bial_bifilar -run TestExport   # base design: netlists, KiCad, results sheets in export/
```

Flags: `-quick` (small populations), `-out` (report directory, default
`runs`), `-line` (Cu, Al or Bi spectral line for the subharmonic report),
`-deltaT` (dry: plate temperature above the coil connections, K),
`-emf-factor` (wet: fraction of the standard cell EMF reached),
`-netlists <dir>` (read the `front_*.csv` files there and write, for the
highest-Q and highest-energy design of each variant and the highest-Q dry
design wound ordinarily, a plain SPICE netlist `.cir` that runs in
kronspice, ngspice, LTspice and QSpice, an LTspice `.asc` and an
equivalent-circuit `.svg`).

## The construction

**The element**, three variants, each with its declared energy port
(§3.1):

| Variant | Element model | Energy enters through |
|---|---|---|
| dry | Bi and Al plates with a solid dielectric; the Bi–Al leads form a thermocouple, EMF (S_Al − S_Bi)·ΔT = 75.5 µV/K·ΔT in series with the plate capacitance | thermal port |
| thermopile | the dry element charged by n Bi–Al junction pairs in series (n·75.5 µV/K·ΔT); genes add n (1–5000) and each leg's diameter and length; the legs add series resistance, volume, and a conductive heat leak n·(k_Bi + k_Al)·A/ℓ·ΔT that holding ΔT costs; legs shorter than 5 diameters are rejected (1-D conduction) | thermal port |
| wet | Al\|electrolyte\|Bi galvanic cell (≤ 1.970 V), electrolyte resistance, Randles electrodes (double layer ∥ charge transfer), geometric plate capacitance; plate capacity by Faraday's law | chemical port |

**The coil** (`em/coil.go`): Wheeler inductance; turn-to-turn capacitance
from the two-wire formula with V/N per gap (ordinary) or V/2 (bifilar),
which reproduces Tesla's 250,000× at 1,000 turns (`TestTeslaComparison`);
skin effect at the resonance; wire longer than λ/10 at resonance is
rejected as not quasi-static (`TestOutsideValidityIsRejected`).

**The circuit** (`Design.Build`): the coil as series R and L with its
self-capacitance across its terminals, in parallel with the element,
between node `a` and ground; a 1 A AC probe current into `a` makes v(a) the
tank impedance. For the switch-on transient the EMF steps on at t = 0 (an
idealized step; a real thermal or chemical source turns on over
milliseconds or longer).

**Resonance** (`Resonance`): the parallel resonance is where the impedance
phase falls through zero from inductive to capacitive; Q = (f0/2)·|dφ/df|
from the phase slope (exact for a parallel RLC and equal there to
f0/bandwidth, and defined when Q is too low for −3 dB points). A coarse
sweep brackets the crossing with the largest |Z|; bisection fixes f0 to
1e-12 (`TestResonanceOfKnownTank`). A tank with no crossing returns
`ErrNoResonance`, which the GA books as no-measurement.

**Genes.** Turns 10–2000, coil diameter 5–150 mm, wire diameter 0.1–2 mm,
pitch 1.05–3 wire diameters, insulation ε 2–6, plate side 5–150 mm, gap
10 µm–2 mm (wet: 50 µm–5 mm) and dielectric ε 2–100 (wet: the electrolyte
sets ε); the thermopile adds pairs 1–5000, leg diameter 50 µm–2 mm and leg
length 1–50 mm. **Objectives:** Q at resonance, and energy per volume as
the error log10(1 MJ/m³ / E_v), zero only at a reference far above anything
these elements reach. The trade-off is traced by a sweep of five
weighted-sum runs per variant (population 40, 30 generations, seed 1) and
the non-dominated designs are reported; NSGA-II arrived later in M9. Every
front design is re-simulated with an ordinary coil of the same geometry,
and its resonance is given as the nearest whole-number subharmonic of the
strongest Cu, Al or Bi emission line (NIST; adjacent subharmonics are ~1/n
apart with n ~ 10⁷–10¹¹, so any frequency is within 1/2n of one, and this
is reported, not a constraint).

**Assumptions**, stated in every report so each can be replaced by a bench
value:

| Assumption | Value | Source |
|---|---|---|
| Bi–Al Seebeck coefficient | 75.5 µV/K (S_Al − S_Bi relative to Pt: +3.5 and −72 µV/K) | Electronics Cooling (Nov 2006), Seebeck table, in `material/` |
| Dry: plates above coil connections | ΔT = 100 K, EMF 7.55 mV | assumption |
| Wet: cell EMF | 1.97 V (Al³⁺/Al −1.662 V, Bi³⁺/Bi +0.308 V), the thermodynamic upper bound | CRC electrochemical series; real cells are lower |
| Wet: electrolyte | ρ = 0.21 Ω·m (sea water, S = 35, 20 °C), ε = 80 | UNESCO 1981 practical salinity scale |
| Wet: double layer, charge transfer | 20 µF/cm², 10 Ω·cm² per electrode | Bard & Faulkner, ch. 13 (typical); the charge-transfer value must be measured |
| Plates | 1 mm Bi + 1 mm Al | assumption |
| Coil | copper; Wheeler inductance; two-wire turn-to-turn capacitance; skin effect at f0 | Wheeler 1928; Tesla US 512,340; Terman 1943 |
| Spectral line | Cu I 324.754 nm (air), 922.9 THz | NIST Handbook of Basic Atomic Spectroscopic Data |

## Results (run 2, `doc/experiments/bial_bifilar_run2.md`)

Each variant: 5,710 evaluations. Dry and wet are bit-identical to run 1
(same seed); run 2 adds the thermopile.

| Variant | Penalized | Front | Highest Q | Highest energy per volume |
|---|---|---|---|---|
| dry | 271 (coil too short for Wheeler, wire too thick, gap not small, no resonance, wire > λ/10) | 59 designs | 39 turns, 150 mm, 311.7 kHz, Q 1,050, 5.4 pJ/m³ | 15 turns, 5 mm, 10 µm gap, 78.9 kHz, Q 0.75, 1.254 µJ/m³ |
| wet | 599 (wire too thick, no resonance 216×, coil too short, gap) | 38 designs | 23 turns, 141 mm, 11.16 kHz, Q 6.4, 4.7 mJ/m³ | 10 turns, 5.0 mm, 15.63 kHz, Q 2.4, 89.3 J/m³ |
| thermopile | 423 (legs shorter than 5 diameters 168×, coil too short, wire > λ/10, gap, wire too thick) | 19 designs | 3 pairs, 11 turns, 21.15 MHz, Q 1,250, 86 pJ/m³ | 5,000 pairs, 19 turns, 27.75 MHz, Q 43.3, 31.25 J/m³ |

For the two thermopile designs the ordinary-coil re-simulation is reported
as outside the lumped model's validity (wire longer than λ/10 at the
ordinary coil's higher resonance), which the report prints in place of the
comparison.

What the fronts say:

- **Dry.** The only store is the field in the plate gap, ½·ε·(EMF/gap)²
  diluted by the plates: 1.255 µJ/m³ by that formula against 1.254 µJ/m³
  simulated for the best design. With a millivolt Seebeck EMF the energy
  per volume is necessarily microscopic. At the steady state no current
  flows (the plates block DC), so the EMF only charges the plates.
- **Wet.** The store is mostly ½·L·I² from the steady cell current, set by
  the assumed charge-transfer resistance: 4.079 A (96.9 mA/cm²) in the
  best design, which a real aluminium electrode would not sustain
  (passivation), so these are upper bounds until a cell is measured. The
  plate's charge lasts the same few hours (8 h 18 min) for every design
  because capacity and current both scale with plate area.
- **Thermopile.** 5,000 pairs raise the EMF to 37.75 V and the stored
  energy per volume from 1.25 µJ/m³ to 31 J/m³, the EMF² scaling
  (5,000² = 2.5×10⁷), but holding ΔT = 100 K across those legs conducts
  47.83 W of heat to store 1.4 mJ (outside the electrical ledger), and the
  plates charge through 16.86 kΩ with a 33 ms time constant.
- **Bifilar against ordinary.** At DC no voltage stands across the
  turn-to-turn capacitance, so the stored energy is the same for both
  windings; the bifilar winding raises the self-capacitance by about N²/4
  and lowers the resonance (on 13 of 59 dry designs by more than 0.1 %),
  and raises Q on none of them.
- **Ledger.** No evaluated design failed the operating-point check. Over
  the six switch-on transients the worst per-point residual is 1.6e-14 of
  peak power, the worst integration mismatch 7.9e-4 of the energy taken
  in, and model-supplied energy is 0 J: every joule comes from the
  element's declared port, heat or chemistry, and nothing else.

Outputs of a run: `report.md` (set-up with sources, Pareto fronts,
switch-on transients with ledgers, computed conclusions) and
`front_dry.csv`, `front_wet.csv`, `front_thermopile.csv`. The kept runs,
their fronts and the netlists and schematics of the best designs are in
[doc/experiments/](../../doc/experiments/).

**Superseded for prediction.** Knight (2016, §6, §13 point 8) shows that
turn-to-turn capacitance theories have no predictive power for
single-layer solenoids (Medhurst's 1947 measurements already showed the
model underestimates an ordinary coil's self-capacitance). The V/N model is
kept to reproduce Tesla's idealized comparison; predictions of ordinary-coil
self-capacitance and self-resonance belong to `em/helix.go`, and the
coupling between bifilar strands to the coupled-line model. The Q and
resonance values above used the superseded model and must be re-run
before they are relied on (CLAUDE.md §3e).

**Temporal-field test (H4, hypothetical).** Under rapid heating with an
imposed anisotropy, Jack's model predicts an EMF from the time variation
of temperature, with magnitude set by the Thomson coefficients and
conductivities of bismuth and aluminium (`em/temporal.go`). The standard
prediction, the Seebeck EMF from the instantaneous temperature difference
only, is what this example computes; the difference is the measurement
target.

## Tests

| Test | Check |
|---|---|
| `TestTeslaComparison` | bifilar stores (N/2)² = 250,000 × the ordinary turn-to-turn energy at 1,000 turns |
| `TestResonanceOfKnownTank` | parallel RLC: f0 = 1/2π√(LC), Q = R√(C/L) |
| `TestDrySteadyState` | stored ½·C·EMF², no steady current, energy through the thermal port |
| `TestWetSteadyState` | steady current EMF/(R_el + R_ct + R_coil); stored ½LI² + ½C_dl V² |
| `TestThermopileScaling` | n pairs: EMF × n, stored × n², leg resistance, volume and heat leak |
| `TestOutsideValidityIsRejected` | wire longer than λ/10 at resonance is refused |
| `TestQuickExperiment` | the whole experiment small: fronts written, ledger closes |
| `TestExport` | base design (200-turn bifilar, 40 mm plates) as dry, wet and dry-ordinary, sheets at the dry tank's f0 (the wet tank has no resonance) |

## Files

| File | What it holds |
|---|---|
| `main.go` | variants, gene space, objectives, the GA sweep and the report |
| `model.go` | assumptions, the element and coil models, `Build`, `Resonance`, `Measure` |
| `netlist.go` | plain SPICE3 netlists of measured designs (PWL switch-on, `.ac` two decades either side of f0) |
| `asc.go`, `schematic.go` | LTspice `.asc` and SVG equivalent-circuit writers |
| `model_test.go`, `asc_test.go`, `export_test.go` | the tests above |
| `export/` | base design: `bial_bifilar_dry`, `bial_bifilar_wet`, `bial_bifilar_dry_ordinary` |

## Sources

- N. Tesla, "Coil for Electro-Magnets," U.S. Patent 512,340, 1894
  (`reference/`) — the bifilar winding and the 250,000× comparison.
- H. A. Wheeler, "Simple inductance formulas for radio coils," *Proc. IRE*
  16, 1398, 1928; F. E. Terman, *Radio Engineers' Handbook*, 1943 —
  inductance and turn-to-turn capacitance.
- R. G. Medhurst, "H.F. resistance and self-capacitance of single-layer
  solenoids," *Wireless Engineer*, 1947 — measured self-capacitance.
- D. W. Knight, "The self-resonance and self-capacitance of solenoid
  coils," 2016, DOI 10.13140/RG.2.1.1472.0887 (`reference/self-res.pdf`)
  — why turn-to-turn models fail; the helical-line model that replaces
  them.
- Seebeck coefficients relative to Pt: *Electronics Cooling*, Nov 2006
  (Bi −72, Al +3.5, Cu +6.5 µV/K), as recorded in `material/material.go`
  with the resistivity, density and electrode-potential sources.
- CRC *Handbook of Chemistry and Physics*, electrochemical series —
  standard potentials of Al and Bi.
- UNESCO, "Background papers and supporting data on the Practical Salinity
  Scale 1978," 1981 — sea-water conductivity.
- A. J. Bard, L. R. Faulkner, *Electrochemical Methods*, ch. 13 —
  double-layer capacitance and the Randles circuit.
- NIST Handbook of Basic Atomic Spectroscopic Data — Cu I, Al I, Bi I
  lines.
- P. M. Jack, arXiv:math-ph/0307038, 2003 — the temporal field (H4);
  hypothesis source only.
