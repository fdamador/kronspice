# An amplifier grown from an embryo

Topology evolution with active parts (CLAUDE.md §4; doc/M8.md). The embryo
(`embryo.cir`) is a 12 V rail on `vcc`, a 1 V AC source behind 1 kΩ on
`in`, a 10 kΩ load on `out` and the 2N3904 model of
[amplifier_gain](../amplifier_gain); nothing connects them. Evolution
chooses resistors, capacitors and transistors, how they are wired and their
values, for 20 dB of gain at 1 kHz. Two runs grow it two ways:

```
go run ./cmd/evolve -config examples/amplifier_topology/config.json             # direct, weighted sum, ≈ 50 s
go run ./cmd/evolve -config examples/amplifier_topology/config_map_elites.json  # Koza's developmental + MAP-Elites, ≈ 65 s
go test ./evo -run TestAmplifierTopology                                         # both, re-simulated from best.cir
```

## Why the first attempts failed

With only a gain target, no run of either encoding found an amplifier
(population 1000, 200 generations, three seeds each): every passive
circuit scores the same (at most −0.83 dB into this load), and a randomly
placed transistor is almost always cut off or saturated, so it amplifies
no better than a wire and the search has no gradient toward one.

The `active_bias` objective supplies it: error 0 when some transistor sits
in its forward-active region at the operating point, rising with the
distance outside it. Its first form asked only for |V_CE| ≥ v_min and
|I_C| ≥ i_min — and five of six runs met it with a transistor whose base
and collector were on the 12 V rail and emitter on ground, a diode across
the supply drawing amperes and amplifying nothing. The region now also
holds 0.4 V ≤ V_BE ≤ 1.0 V, V_BC ≤ 0.4 V and I_C ≤ i_max (100 mA), and the
same runs grow real amplifiers. Evolution finds every gap in a fitness
function; this one is recorded in `evo.ActiveBias`'s comment.

## Run 1: direct encoding, weighted sum (`config.json`)

Parts R (10 Ω–1 MΩ, E24), C (1 nF–100 µF, E12), Q (2N3904); objectives gain
(weight 1), active bias (v_min 2 V, i_min 100 µA, weight 0.5), parts
(0.002 each). Seed 2:

```
Q_g4  collector g_n1, emitter 0, base in      V_CE 8.93 V, I_C 1.42 mA
R_g2  vcc  g_n1  1.5k     collector load
R_g1  g_n1 in    13k      collector-to-base feedback bias
R_g3  vcc  g_n5  2.7k  ┐  second bias path to the base
R_g6  g_n5 in    240k  ┘
C_g5  g_n1 out   330n     output coupling
```

20.05 dB at 1 kHz on the written netlist. A common-emitter stage with
collector feedback, the textbook self-stabilizing bias. Seeds 1 and 3 also
reach 20 dB (8 and 11 parts).

## Run 2: Koza's developmental encoding with MAP-Elites (`config_map_elites.json`)

The embryo gets two modifiable wires, `in`→`out` and `vcc`→`out`, each
developed by a circuit-constructing program (`evo/encode_develop.go`);
`(Q value cc cc cc)` puts a transistor in place of the component under the
writing head with three lead wires that grow on. With the weighted sum this
encoding learns to bias a transistor but never cuts the `in`→`out` wire
that shorts the signal past it: that takes several changes, each worthless
alone.

MAP-Elites (`ga.algorithm: "map_elites"`; Mouret & Clune 2015) keeps
instead the best design in every cell of a four-dimensional grid —
transistor count (4 bins), gain error (12), bias error (4), program size
(8) — and breeds from the occupied cells uniformly, so a poor design that
is the only one of its kind survives as a stepping stone. Seed 2:

```
Q_g4  collector g_n3, emitter 0,   base in    V_CE 8.12 V, I_C 0.18 mA
R_g3  g_n3 vcc 20k       collector load
Q_g1  collector vcc,  emitter out, base g_n3  V_CE 4.55 V, I_C 1.36 mA
R_g2  in   out 11k       feedback from output to input
```

20.04 dB with four parts: a common-emitter stage driving an emitter
follower, with the whole amplifier biased and its gain set by one feedback
resistor — the kind of compact solution Koza reports. Seed 3 reaches
17 dB, seed 1 stays at the wire.

| Seeds 1, 2, 3 | gain error after 200 × 1000 | parts |
|---|---|---|
| direct, gain only | 1.043, 0.838, 1.043 | — |
| developmental, gain only | 1.041 ×3 | — |
| direct, + active bias | 0.0024, 0.0025, 0.0018 | 11, 6, 8 |
| developmental, + active bias | 1.041 ×3 | — |
| developmental, + bias, MAP-Elites | 1.041, 0.0019, 0.150 | —, 4, 12 |
| direct, + bias, MAP-Elites | 0.0016, 0.0008, 1.041 | 4, 4, — |

(gain error = |gain − 20 dB|/20; 1.04 is a passive circuit.)

## Files

| File | What it holds |
|---|---|
| `embryo.cir` | rail, source, load, 2N3904 model, sweep |
| `config.json` | direct encoding, weighted sum with the bias term |
| `config_map_elites.json` | developmental encoding, MAP-Elites grid |
| `export/` | the bare embryo as `.cir`, `.kicad_sch` and a results sheet |
| `../../evo/fitness.go` | `ActiveBias` |
| `../../evo/mapelites.go` | MAP-Elites archive and descriptors |
| `../../evo/amplifier_topology_test.go` | both runs, re-simulated |

## Sources

- J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, *Genetic
  Programming III* (1999) — transistor-inserting construction functions,
  evolved amplifiers.
- J.-B. Mouret, J. Clune, "Illuminating search spaces by mapping elites,"
  arXiv:1504.04909 (2015) — MAP-Elites.
- The 2N3904 model as in [amplifier_gain](../amplifier_gain).
