# NMOS common-source amplifier: a gain / bandwidth / power Pareto front

The M9 example. A common-source NMOS stage on 5 V, driven from a 1 kΩ
source into a 100 kΩ ∥ 5 pF load, with a Level 3 (semi-empirical) 1 µm-class
process that includes weak inversion (`NFS`), so low-power designs biased
near or below threshold follow the exponential subthreshold law. Evolution
chooses the gate bias divider R1/R2, the drain load RD, the bypassed source
resistor RSS (all E24) and the channel width and length W and L (continuous,
2–500 µm and 1–10 µm, through instance-parameter genes `M1.w`, `M1.l`).

The GA is NSGA-II (`"algorithm": "nsga2"`): no weights are traded against
each other; the run returns the set of designs none of which is beaten on
every objective by another.

| Objective | Kind | Meaning |
|---|---|---|
| gain at 10 kHz, V(o)/V(in) | `gain_max`, ceiling 40 dB | more is better (error = shortfall/20) |
| upper −3 dB edge | `bandwidth_max`, ceiling 1 GHz | more is better (error = decades short) |
| power from VDD at the operating point | `supply_power`, weight 1000 (mW) | less is better |
| gain at 10 kHz ≥ 20 dB | `gain_max` with `"constraint": true` | a constraint, not an objective: a design below 20 dB is dominated by every design above it |

```bash
go run ./cmd/evolve -config examples/mos_amplifier_pareto/config.json   # ≈ 6 s
```

Besides the usual outputs the run writes `front.csv` (every front design:
its three errors, the constraint and the gene values, sorted by the first
objective) and `front/NNN.cir`, one netlist per design, each of which opens
in ngspice, LTspice and QSpice as it is.

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
in moderate inversion. Seeds 1–6 reach 34.5–35.2 dB, 34–36 MHz and
20–34 µW. `evo/mos_pareto_test.go` re-simulates every front netlist from
scratch and checks the constraint, the recorded errors, the energy ledger
and mutual non-dominance.
