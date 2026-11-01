# 2N3904 common-emitter amplifier: evolve for gain and bandwidth

The M7 example. A common-emitter stage around a 2N3904 (Fairchild datasheet
Gummel–Poon model), driven from a 600 Ω source into a 10 kΩ load on 12 V.
Evolution chooses eight standard values: the bias divider R1/R2, the
collector load RC, the unbypassed and bypassed emitter resistors RE1/RE2,
the Miller capacitor CF, the input coupling CIN and the emitter bypass CE
(resistors E24, capacitors E12).

| Objective | Target | Weight |
|---|---|---|
| gain at 10 kHz, V(o)/V(in) | 30 dB | 1 |
| upper −3 dB edge | 200 kHz | 1 |
| lower −3 dB edge | 100 Hz | 1 |
| collector voltage | 6 V | 0.2 |
| emitter voltage | 1.5 V (so V_CE ≈ 4.5 V) | 0.2 |

```bash
go run ./cmd/evolve -config examples/amplifier_gain/config.json   # ≈ 25 s
```

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
V_CE ≥ 3 V, the tolerance of the 5 % parts the design is made of. 10 of
seeds 1–12 meet it; the other two miss the gain by less than 0.1 dB beyond
0.5 dB. With E-series values in ~4 % steps and five coupled targets, the
weighted error settles near 0.02: no single-part change of the best design
improves it, so that is the floor of the grid rather than of the search.

`best.cir` opens in ngspice, LTspice and QSpice as it is.
