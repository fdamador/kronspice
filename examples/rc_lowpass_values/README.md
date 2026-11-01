# RC low-pass: evolve two standard values for a 1 kHz cutoff

The M5 example, the first and simplest genetic-algorithm experiment
(CLAUDE.md §4, value genome). A first-order RC low-pass driven by a 1 V AC
source; evolution chooses the resistor from the E24 series and the
capacitor from the E12 series so that the −3 dB frequency of V(out)/V(in)
lands at 1 kHz. Everything is in the config: the netlist is a plain SPICE
file and `cmd/evolve` reads the genes, the objective and the GA settings
from `config.json`.

```bash
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json   # ≈ 1 s
go test ./evo -run TestRCLowpassExample                              # runs this config end to end
```

## Circuit, genes and objective

```
V1 in  0   AC 1
R1 in  out 10k        gene: 100 Ω … 1 MΩ, E24
C1 out 0   10n        gene: 100 pF … 10 µF, E12
.ac dec 100 10 100k
```

| Objective | Target | Error |
|---|---|---|
| upper −3 dB edge of V(out) | 1000 Hz | \|f − 1000\| / 1000 |

GA settings: population 40, up to 40 generations, seed 1, a 2 s deadline
per evaluation, and `target_fitness` 0.006, which stops the run as soon as
a design is within 0.6 % of the target.

## Result

Seed 1 stops after generation 0 (40 evaluations): the random initial
population already contains a design under the target.

| Part | Value |
|---|---|
| R1 | 160 Ω |
| C1 | 1 µF |

fc = 1/(2πRC) = 994.7 Hz, an error of 0.53 %. That is the closest any
E24 × E12 pair gets to 1 kHz (the other such pair is 1.6 kΩ with 100 nF,
the same product), so the result is the floor of the value grid, not of
the search (doc/M5.md). The test `evo.TestRCLowpassExample` runs this
config exactly as `cmd/evolve` does, checks the cutoff is within 10 Hz,
re-parses the written `best.cir` and measures the same cutoff from it.

Outputs go to `runs/<timestamp>-rc_lowpass_values/`: a copy of the config,
`stats.csv` (best, mean and worst fitness per generation), `best/genNNNN.cir`,
`best.cir` and `checkpoint.json`. `best.cir` opens in ngspice, LTspice and
QSpice as it is.

## Files

| File | What it holds |
|---|---|
| `circuit.cir` | the template netlist with the starting values |
| `config.json` | genes, objective and GA settings |
| `export/` | the template's netlist, KiCad schematic and results sheet (MNA against Kron's torn solve), written by `KRONSPICE_EXPORT=1 go test ./internal/export` |

## See also

- [doc/M5.md](../../doc/M5.md) — the value-genome engine, run directories, penalties.
- [negative_resistance_gain](../negative_resistance_gain) and
  [rectifier_filter](../rectifier_filter) — the other two config-driven M5 examples.
