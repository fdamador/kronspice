# Harmonic filters for a six-pulse rectifier, scored by geometric-algebra power

The M9f example (CLAUDE.md §3b item 3): a genetic algorithm whose fitness is
a power-quality measure computed with geometric algebra. A 60 Hz, 120 V
(phase) three-phase supply with 0.05 Ω and 0.5 mH per line feeds a diode
bridge into a 2200 µF, 10 Ω DC load. At the point of common coupling each
phase carries two series L–C branches to a grounded star point, nominally
for the 5th and 7th harmonics; the four values are the genome. The report
shows the IEEE 1459 quantities before and after, names the pairs of
harmonics that carry the most distortion power, and checks the quaternion
instantaneous power of the three phases against ΣP.

```bash
go run ./examples/power_filter                 # ≈ 2 min; -gens 10 -pop 24 -seed 1 are the defaults
go test ./examples/power_filter                # GA and IEEE 1459 agree; a short GA improves PF (-short skips the GA)
```

## Fitness

For each phase the last two of 14 simulated cycles are resolved into
harmonics up to the 25th (`analysis.HarmonicsOf`), and the geometric
product M = V·I in Cl(2H) (`analysis.GeometricPower`) gives the scalar
part P, the same-harmonic bivectors Q and the cross-harmonic bivectors D,
grouped by pair of harmonics. Fitness is 1 − PF with PF = ΣP/ΣS over the
three phases. The genes are L5, L7 (0.1–10 mH) and C5, C7 (10–1000 µF),
all E12; a branch resistance of 0.05 Ω is fixed.

## Result (defaults: 10 generations × 24, 222 evaluations)

| | no filter | evolved |
|---|---|---|
| PF | 0.9033 | 0.9559 |
| P | 7389 W | 7673 W |
| phase a: S / P / Q1 | 2731 VA / 2463 W / +515 var | 2701 VA / 2559 W / −73 var |
| D_I | 1018 VA | 784 VA |
| THD_I / THD_V | 40.4 % / 9.9 % | 30.6 % / 12.8 % |
| largest distortion pairs | 1–5 961 VA, 1–7 333 VA, 1–11 192 VA | 1–5 578 VA, 1–7 461 VA, 1–11 156 VA |

Evolved branches: L5 3.9 mH with C5 82 µF (281 Hz) and L7 0.15 mH with
C7 10 µF (4109 Hz). The GA and IEEE 1459 columns agree per phase (S to
0.3 %, the commutation notches above the 25th harmonic), and the quaternion
mean power equals ΣP (7389 W and 7673 W).

The branches do not sit at 300 Hz and 420 Hz. Maximizing PF trades
harmonic absorption against the leading fundamental reactive power the
capacitors add: the 5th-harmonic branch is tuned a little low, which
removes most of the 1–5 distortion and swings Q1 from +515 to −73 var, and
the second branch here is mostly a capacitor. A longer run (30 × 20, 562
evaluations) in doc/M9f.md reached PF 0.957 with branches at 284 Hz and
673 Hz. The program exits non-zero if the evolved filter does not improve
the power factor.

## Files

| File | What it holds |
|---|---|
| `main.go` | the rectifier, the analysis at the PCC, the GA problem and the report |
| `main_test.go` | GA against IEEE 1459, largest pair 1–5, quaternion mean = ΣP, a short GA run |
| `export_test.go` | the unfiltered rectifier and one with 1 mH, 270 µF / 150 µF branches (300 Hz, 420 Hz) as netlists, KiCad schematics and results sheets in `export/` (`KRONSPICE_EXPORT=1`) |

## Sources

- IEEE Std 1459-2010, *Definitions for the Measurement of Electric Power
  Quantities Under Sinusoidal, Nonsinusoidal, Balanced, or Unbalanced
  Conditions*.
- M. Castro-Núñez, R. Castro-Puche, "The IEEE Standard 1459, the CPC power
  theory, and geometric algebra in circuits with nonsinusoidal sources and
  linear loads," *IEEE Trans. Circuits Syst. I*, 2012.
- [doc/M9f.md](../../doc/M9f.md) — quaternion polyphase, GA power and its checks.
