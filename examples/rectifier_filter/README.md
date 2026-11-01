# Half-wave rectifier: evolve the filter for mean output and ripple

An M5 example (CLAUDE.md §4, value genome) and the first config-driven run
that is nonlinear and transient: a 10 V peak, 1 kHz sine through a diode
into a capacitor and load. Evolution chooses the filter capacitor, the load
resistor and the diode's series resistance so that, once the start-up
transient has settled, the output averages 9 V with no more than 0.2 V of
peak-to-peak ripple.

```bash
go run ./cmd/evolve -config examples/rectifier_filter/config.json   # ≈ 20 s
```

## Circuit, genes and objectives

```
V1 in  0   SIN(0 10 1k)
D1 in  out dr
R1 out 0   1k          gene: 100 Ω … 100 kΩ, E24
C1 out 0   1u          gene: 100 nF … 1 mF, E12
.model dr D(IS=1e-14 N=1 RS=0.5)     gene: RS 0.05 … 50 Ω (a .model parameter)
.tran 10u 8m
```

| Objective | Window | Target | Error |
|---|---|---|---|
| mean of V(out) | 4 ms – 8 ms | 9.0 V | \|mean − 9\| / 9 |
| peak-to-peak ripple of V(out) | 4 ms – 8 ms | ≤ 0.2 V | max(0, ripple − 0.2) / 0.2 |

Both objectives weigh 1; the ripple is a constraint, costing nothing once
it is met. GA settings: population 30, up to 25 generations, seed 11, a 5 s
deadline per evaluation and `target_fitness` 0.002.

The third gene shows that a `.model` parameter evolves like a part value:
the diode's RS sets the charging current's peak and the drop at the crest.

## Result

Seed 11 reaches the target after generation 13 (394 evaluations):

| Part | Value |
|---|---|
| C1 | 820 µF |
| R1 | 270 Ω |
| dr RS | 0.119 Ω |

Re-simulated from the written `best.cir` over 4–8 ms: mean 9.007 V
(0.08 % from the target) and 0.154 V of ripple, under the 0.2 V limit
(doc/M5.md: "within 0.08 %, ripple met"). The time constant R1·C1 =
221 ms against the 1 ms period is what holds the ripple down; the mean is
then set by the crest voltage, 10 V less the diode drop at the charging
peak.

Outputs go to `runs/<timestamp>-rectifier_filter/` (`stats.csv`,
`best/genNNNN.cir`, `best.cir`, `checkpoint.json`); `best.cir` opens in
ngspice, LTspice and QSpice as it is. Every transient is checked by the
energy ledger as it runs: a candidate whose books do not balance is
penalized, never reported.

## Files

| File | What it holds |
|---|---|
| `circuit.cir` | the template netlist with the starting values |
| `config.json` | genes, objectives and GA settings |
| `export/` | the template's netlist, KiCad schematic and operating-point sheet, written by `KRONSPICE_EXPORT=1 go test ./internal/export` |

## See also

- [doc/M5.md](../../doc/M5.md) — the value-genome engine; [doc/M4.md](../../doc/M4.md) — the diode and the Newton fallback chain.
- [rc_lowpass_values](../rc_lowpass_values) and
  [negative_resistance_gain](../negative_resistance_gain) — the other config-driven M5 examples.
