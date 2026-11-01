# Bi–Al capacitor element with a bifilar coil

The project's first experiment. A bismuth–aluminium capacitor element and a
bifilar coil (Tesla, US 512,340) in parallel, self-driven: the element is the
only energy source. Evolution trades the tank's quality factor Q at
resonance against the energy stored per volume.

| Variant | Element model | Energy enters through |
|---|---|---|
| dry | Bi and Al plates with a solid dielectric; the Bi–Al leads form a thermocouple, EMF (S_Al − S_Bi)·ΔT = 75.5 µV/K·ΔT in series with the plate capacitance | thermal port |
| thermopile | the dry element charged by n Bi–Al junction pairs in series (n·75.5 µV/K·ΔT); genes add n (1–5000) and each leg's diameter and length; the legs add series resistance, volume, and a conductive heat leak n·(k_Bi + k_Al)·A/ℓ·ΔT that holding ΔT costs | thermal port |
| wet | Al\|electrolyte\|Bi galvanic cell (≤ 1.970 V), electrolyte resistance, Randles electrodes (double layer ∥ charge transfer), geometric plate capacitance | chemical port |

Coil: Wheeler inductance; turn-to-turn capacitance with V/N per gap (ordinary)
or V/2 (bifilar), which reproduces Tesla's 250,000× at 1,000 turns; skin
effect at the resonance; wire longer than λ/10 at resonance is rejected as
not quasi-static. Every non-dominated design is re-simulated with an
ordinary coil of the same geometry, and its resonance is given as the
nearest whole-number subharmonic of the strongest Cu, Al or Bi emission line
(NIST).

```bash
go run ./examples/bial_bifilar                 # ≈ 1 min: report in runs/<time>-bial_bifilar/
go run ./examples/bial_bifilar -quick          # a few seconds
go run ./examples/bial_bifilar -deltaT 300 -emf-factor 0.6 -line Al
go run ./examples/bial_bifilar -netlists doc/experiments   # SPICE netlists + schematics of the best designs
```

`-netlists <dir>` reads the `front_*.csv` files in a directory and writes,
for the highest-Q and highest-energy design of each variant (and the
highest-Q dry design wound ordinarily), a plain SPICE netlist (`.cir`, runs
in kronspice, ngspice, LTspice and QSpice) and its equivalent-circuit
schematic (`.svg`).

Outputs: `report.md` (set-up with sources, Pareto fronts, switch-on
transients with their energy ledgers, computed conclusions) and
`front_dry.csv`, `front_wet.csv`, `front_thermopile.csv`. Kept runs are in
`doc/experiments/`: run 1 (dry, wet) and run 2 (adds the thermopile; dry
and wet are bit-identical to run 1, same seed).

Run 2 result: 5000 pairs raise the EMF to 37.75 V and the stored energy per
volume from 1.25 µJ/m³ to 31 J/m³ — the EMF² scaling (5000² = 2.5·10⁷) — but
holding ΔT = 100 K across those legs conducts about 48 W of heat to store
1.4 mJ, and the plates charge through 16.9 kΩ with a 33 ms time constant.

All material constants and modelling assumptions are listed, with sources,
at the top of every report; the wet-cell electrode values are typical
figures that must be replaced by measurements of a built cell.
