# 2N3904 common-emitter amplifier: evolve for gain and bandwidth

The M7 example (CLAUDE.md §9; [doc/M7.md](../../doc/M7.md)), the first
run of the value-genome GA (§4) on a nonlinear circuit. A common-emitter
stage around a 2N3904 (Fairchild datasheet Gummel–Poon model), driven from
a 600 Ω source into a 10 kΩ load on 12 V. Evolution chooses eight standard
values: the bias divider R1/R2, the collector load RC, the unbypassed and
bypassed emitter resistors RE1/RE2, the Miller capacitor CF, the input
coupling CIN and the emitter bypass CE (resistors E24, capacitors E12).
Its multi-objective successor is [mos_amplifier_pareto](../mos_amplifier_pareto);
the same engine on a linear fixture is [rc_lowpass_values](../rc_lowpass_values).

```
go run ./cmd/evolve -config examples/amplifier_gain/config.json   # ≈ 25 s
go test ./evo -run TestAmplifierExample                            # the M7 acceptance test
```

Outputs go to `runs/<timestamp>-amplifier_gain/`: a copy of the config,
`stats.csv` (best, mean and worst fitness per generation), `best.cir` per
generation and a final `checkpoint.json` that `-resume` continues from.

## The circuit

`circuit.cir` is a plain SPICE deck: VCC 12 V, VIN with RS 600 Ω coupled
through CIN into the base, the divider R1/R2, RC to the collector, RE1 in
series with the bypassed RE2 (so RE1 alone sets the mid-band gain
≈ RC ∥ RL / RE1), CF from collector to base (the Miller capacitor that
sets the upper edge), COUT 10 µF into RL 10 kΩ, and `.ac dec 20 10 100meg`.
The transistor is the datasheet's SPICE model (IS 6.734 fA, BF 416.4,
VAF 74.03 V, CJC 3.638 pF, CJE 4.493 pF, TF 301 ps, RB 10 Ω); the
datasheet print "Ise=6.734" is a known misprint for 6.734f, corrected in
the netlist. The engine solves it with the SPICE3 Gummel–Poon equations of
`device/bjt.go` (forward and reverse Early effect, high injection,
non-ideal base currents, depletion and diffusion charges, `pnjlim`
limiting), the operating point by Newton with the M4 fallback chain, and
the AC response at a polished operating point (two extra Newton iterations
so the small-signal transconductances are exact to round-off, doc/M7.md).

### Genes and objectives

| Gene | Range | Series |
|---|---|---|
| R1 | 10 kΩ – 1 MΩ | E24 |
| R2 | 1 kΩ – 100 kΩ | E24 |
| RC | 470 Ω – 47 kΩ | E24 |
| RE1 | 10 Ω – 2.2 kΩ | E24 |
| RE2 | 100 Ω – 10 kΩ | E24 |
| CF | 1 pF – 1 nF | E12 |
| CIN | 100 nF – 100 µF | E12 |
| CE | 1 µF – 1 mF | E12 |

| Objective | Target | Weight |
|---|---|---|
| gain at 10 kHz, V(o)/V(in) | 30 dB | 1 |
| upper −3 dB edge | 200 kHz | 1 |
| lower −3 dB edge | 100 Hz | 1 |
| collector voltage | 6 V | 0.2 |
| emitter voltage | 1.5 V (so V_CE ≈ 4.5 V) | 0.2 |

GA settings: population 100, 300 generations, tournament 5, mutation rate
0.25 with σ = 0.1 decade, seed 1, 2 s timeout per evaluation, stop at a
weighted error of 0.02.

## Results

Seed 1 (the config's seed) after 300 generations:

| Part | Value |
|---|---|
| R1 / R2 | 30 kΩ / 6.8 kΩ |
| RC | 5.1 kΩ |
| RE1 / RE2 | 68 Ω / 1.3 kΩ |
| CF | 33 pF |
| CIN / CE | 4.7 µF / 18 µF |

Re-simulated from the written `best.cir`: 29.92 dB at 10 kHz, edges at
100.0 Hz and 199.8 kHz, V_CE 4.9 V. The spec the acceptance test holds
(`evo/amplifier_test.go`) is 30 ± 0.5 dB, edges within 5 % and
V_CE ≥ 3 V, the tolerance of the 5 % parts the design is made of; it also
requires the AC and OP energy ledgers to balance to 1e-9. 10 of
seeds 1–12 meet it; the other two miss the gain by less than 0.1 dB beyond
0.5 dB (seeds 3 and 10, at 29.44 and 29.41 dB; full table in doc/M7.md).

| Seed | Gain | Lower edge | Upper edge | V_CE |
|---|---|---|---|---|
| 1 | 29.92 dB | 100.0 Hz | 199.8 kHz | 4.91 V |
| 2 | 30.37 dB | 100.1 Hz | 201.3 kHz | 4.62 V |
| 5 | 29.92 dB | 100.0 Hz | 199.9 kHz | 4.35 V |
| 9 | 30.24 dB | 99.9 Hz | 200.0 kHz | 4.27 V |

`best.cir` opens in ngspice, LTspice and QSpice as it is.

## How the search behaves

- **The E-series floor.** With E24/E12 values in ~4 % steps and five
  coupled targets, the weighted error settles near 0.02: no single-part
  change of the best design improves it, so that is the floor of the grid
  rather than of the search. A local search that moves several parts at
  once (the bias resistors together) would likely close the remaining
  gap; it is noted under M8 as future work.
- **Operators** (CLAUDE.md §4, as implemented in M5): tournament selection,
  BLX-0.5 crossover on log magnitudes (Eshelman & Schaffer 1993), Gaussian
  mutation in decades, snapping to the E-series, elitism.
- **Penalties, never a crash.** A candidate that fails to converge, gives a
  singular matrix, fails validation, breaks the ledger or times out gets a
  penalty fitness from the M5 classes (no-measurement < no-convergence <
  singular < invalid < ledger < timeout < other < panic), so the run keeps
  going and learns from it.
- **Determinism.** One PCG stream per generation derived from the seed and
  results stored by individual index, so a run is bit-identical for any
  worker count and resumes exactly from its checkpoint (§1 goal 5).

## Files

| File | What it holds |
|---|---|
| `circuit.cir` | the stage, the 2N3904 model and the `.ac` sweep; starting values |
| `config.json` | genes, objectives and GA settings |
| `export/` | the starting circuit as `.cir`, `.kicad_sch` and a results sheet (MNA only: the torn solve is linear-only) |
| `../../evo/amplifier_test.go` | `TestAmplifierExample`, the M7 acceptance test |
| `../../doc/M7.md` | milestone notes: BJT model, golden references, the seed table |

## Sources

- Fairchild Semiconductor, "2N3904 NPN General Purpose Amplifier,"
  datasheet (1997), SPICE model section.
- H. K. Gummel, H. C. Poon, "An integral charge control model of bipolar
  transistors," *Bell Syst. Tech. J.* 49, 827 (1970).
- L. W. Nagel, *SPICE2: A Computer Program to Simulate Semiconductor
  Circuits*, UC Berkeley ERL Memo M520 (1975), ch. 4.
- G. Massobrio, P. Antognetti, *Semiconductor Device Modeling with SPICE*,
  2nd ed., ch. 2 — the Gummel–Poon equations as implemented in
  `device/bjt.go`.
- SPICE3f3 User's Manual §4.3 — BJT model parameters and defaults.
- L. J. Eshelman, J. D. Schaffer, "Real-coded genetic algorithms and
  interval-schemata," *Foundations of Genetic Algorithms* 2 (1993) — BLX-α
  crossover.
