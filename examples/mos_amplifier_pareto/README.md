# NMOS common-source amplifier: a gain / bandwidth / power Pareto front

The M9 example (CLAUDE.md §9; [doc/M9.md](../../doc/M9.md)): the MOSFET
levels 1–3 and NSGA-II together. A common-source NMOS stage on 5 V, driven
from a 1 kΩ source into a 100 kΩ ∥ 5 pF load, with a Level 3
(semi-empirical) 1 µm-class process that includes weak inversion (`NFS`),
so low-power designs biased near or below threshold follow the exponential
subthreshold law. Evolution chooses the gate bias divider R1/R2, the drain
load RD, the bypassed source resistor RSS (all E24) and the channel width
and length W and L (continuous, 2–500 µm and 1–10 µm, through
instance-parameter genes `M1.w`, `M1.l`). The single-objective BJT
predecessor is [amplifier_gain](../amplifier_gain).

```
go run ./cmd/evolve -config examples/mos_amplifier_pareto/config.json   # ≈ 6 s
go test ./evo -run 'TestMOSParetoExample|TestNSGA2'                       # acceptance and ZDT benchmarks
```

Besides the usual outputs the run writes `front.csv` (every front design:
its three errors, the constraint and the gene values, sorted by the first
objective) and `front/NNN.cir`, one netlist per design, each of which opens
in ngspice, LTspice and QSpice as it is.

## The circuit

`circuit.cir`: VDD 5 V; VIN through RSIG 1 kΩ and CIN 1 µF into the gate;
the divider R1/R2; RD to the drain; RSS bypassed by CSS 100 µF at the
source; `M1 d g src 0 nproc`; COUT 1 µF into RL 100 kΩ ∥ CL 5 pF;
`.ac dec 10 100 1g`. The process model is `LEVEL=3` with parameters derived
from TOX, NSUB and UO as in SPICE3 (TOX 20 nm, NSUB 3e16, UO 550, VMAX
1.5e5, THETA 0.06, ETA 0.05, KAPPA 0.4, NFS 4e11, overlap and junction
capacitances), the same n-channel process as the `mos3_*` golden circuits.
`device/mosfet_level23.go` writes the Level 3 current in forward-mode dual
numbers, so gm, gds and gmbs are the exact derivatives of the stamped
current (doc/M9.md); the gate capacitances are Meyer's.

### Genes and objectives

| Gene | Range | Series |
|---|---|---|
| R1, R2 | 100 kΩ – 10 MΩ | E24 |
| RD | 100 Ω – 1 MΩ | E24 |
| RSS | 10 Ω – 100 kΩ | E24 |
| M1 W | 2 – 500 µm | continuous |
| M1 L | 1 – 10 µm | continuous |

The GA is NSGA-II (`"algorithm": "nsga2"`; population 100, 100
generations, mutation 0.2, σ 0.15, seed 1): no weights are traded against
each other; the run returns the set of designs none of which is beaten on
every objective by another.

| Objective | Kind | Meaning |
|---|---|---|
| gain at 10 kHz, V(o)/V(in) | `gain_max`, ceiling 40 dB | more is better (error = shortfall/20) |
| upper −3 dB edge | `bandwidth_max`, ceiling 1 GHz | more is better (error = decades short) |
| power from VDD at the operating point | `supply_power`, weight 1000 (mW) | less is better |
| gain at 10 kHz ≥ 20 dB | `gain_max` with `"constraint": true` | a constraint, not an objective: a design below 20 dB is dominated by every design above it |

## Results

Seed 1 (the config's seed): in generation 0 one random design of 100 meets
the 20 dB constraint; after 100 generations the front holds 99 designs.
Five of them:

| Design | Gain | Bandwidth | Power | RD | RSS | W / L |
|---|---|---|---|---|---|---|
| highest gain | 34.45 dB | 1.10 MHz | 404 µW | 62 kΩ | 180 Ω | 500 µm / 1.84 µm |
| widest band | 20.17 dB | 33.8 MHz | 19.4 mW | 1 kΩ | 130 Ω | 499 µm / 1.22 µm |
| lowest power | 24.10 dB | 344 kHz | 22.2 µW | 1 MΩ | 27 kΩ | 301 µm / 6.47 µm |
| ≥ 30 dB at least power | 31.41 dB | 392 kHz | 56.7 µW | 430 kΩ | 910 Ω | 427 µm / 3.43 µm |
| ≥ 10 MHz at least power | 23.49 dB | 11.7 MHz | 2.88 mW | 3.3 kΩ | 750 Ω | 500 µm / 1.32 µm |

The shape is the textbook one: bandwidth costs power (a small RD against
the 5 pF load needs current to keep the gain), and the gain tops out near
35 dB on every seed because the 100 kΩ load caps the drain resistance, so
the gain is gm·(RD ∥ RL ∥ ro) with gm/I_D at most its weak-inversion value
1/(n·V_T): the highest-gain designs use the widest device the genes allow,
in moderate inversion. Seeds 1–6 reach 34.5–35.2 dB, 34–36 MHz at
17–21 mW, and 20–34 µW, all at gain ≥ 20 dB.

`evo/mos_pareto_test.go` re-simulates every front netlist from scratch and
checks the constraint, the recorded errors, the energy ledger and mutual
non-dominance, and that the front spans the trade-off (a design above
34 dB, one above 30 MHz, one below 50 µW).

## How the search behaves

NSGA-II as in Deb et al. 2002, in `evo/nsga2.go`:

- **Ranking:** fast non-dominated sort (§III-A) and crowding distance
  (§III-B); ties are broken by index, so runs are bit-identical across
  worker counts and resume exactly from a checkpoint.
- **Constrained dominance (§VI):** an objective marked `"constraint": true`
  is a violation (0 = met). Feasible designs dominate violators, violators
  are ordered by total weighted violation, and failed simulations
  (singular matrix, no convergence, timeout) come last, ordered by their
  M5 penalty class. A failure never crashes the run.
- **Loop (§III-C):** binary crowded tournament, the same BLX-α and Gaussian
  operators as the weighted GA, then the best `Population` of parents ∪
  offspring by rank and crowding. Selection is elitist by construction.
- **Benchmarks** (`TestNSGA2ZDT`): on ZDT1 the 60 front members lie within
  0.0065 of the true front; on ZDT2 (non-convex, where a weighted sum can
  only find the two ends) within 0.0040 with 31 interior points.

## Files

| File | What it holds |
|---|---|
| `circuit.cir` | the stage, the Level 3 process model and the `.ac` sweep |
| `config.json` | genes (including the `M1.w`, `M1.l` instance parameters), the three objectives and the constraint, NSGA-II settings |
| `export/` | the starting circuit as `.cir`, `.kicad_sch` and a results sheet (MNA only: the torn solve is linear-only) |
| `../../evo/nsga2.go` | NSGA-II |
| `../../evo/mos_pareto_test.go` | `TestMOSParetoExample`, the M9 acceptance test |
| `../../doc/M9.md` | milestone notes: MOSFET levels, golden references, NSGA-II |

## Sources

- K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A fast and elitist
  multiobjective genetic algorithm: NSGA-II," *IEEE Trans. Evol. Comput.*
  6(2), 182–197 (2002).
- A. Vladimirescu, S. Liu, "The simulation of MOS integrated circuits
  using SPICE2," UCB/ERL M80/7 (1980), §2.4 — the Level 3 drain current,
  process parameters and junctions.
- H. Shichman, D. A. Hodges, "Modeling and simulation of insulated-gate
  field-effect transistor switching circuits," *IEEE J. Solid-State
  Circuits* 3(3), 285–289 (1968) — Level 1, the base the levels share.
- J. E. Meyer, "MOS models and circuit simulation," *RCA Review* 32, 42–63
  (1971) — gate capacitances.
- SPICE3f3 User's Manual §4.4 and §4.10 — MOSFET parameters and the
  process-parameter equations.
- L. J. Eshelman, J. D. Schaffer, "Real-coded genetic algorithms and
  interval-schemata," *Foundations of Genetic Algorithms* 2 (1993) — BLX-α
  crossover.
