# A low-pass filter grown from an empty embryo

The M8 example: topology evolution. The embryo (`embryo.cir`) is only the
test fixture, after Koza et al., *Genetic Programming III* (1999), ch. 25:
a 1 V source with 1 kΩ source resistance driving node `in`, and a 1 kΩ load
on node `out`. Nothing connects `in` to `out`. Evolution decides which
parts exist (inductors, capacitors, resistors, E24 values), how they are
wired among ground, `in`, `out` and eight internal nodes, and their values.

| Mask | Requirement |
|---|---|
| passband 10 Hz–1 kHz | −7.5 dB ≤ gain ≤ −5.5 dB (matched: −6.02 dB) |
| stopband 3 kHz–100 kHz | gain ≤ −46 dB (40 dB below the passband) |

Fitness is the mask error (0 when the response stays inside the mask,
otherwise at least 0.1, growing with how far and how often it leaves) plus
0.002 per part, so once the mask is met evolution keeps removing parts.

```bash
go run ./cmd/evolve -config examples/filter_topology/config.json   # ≈ 40 s
```

Seed 1, after 200 generations (60,000 evaluations), four parts:

```
L_g1 in  out 0.2       series inductor
C_g4 in  out 11n       across it: a notch at 1/(2π√(LC)) = 3.39 kHz
C_g2 0   in  180n      shunt capacitors at both ends
C_g3 out 0   220n
```

That is a third-order elliptic (Cauer) low-pass, π-section with the trap
across the series arm. On a sweep ten times denser than evolution used
(200 points per decade) the passband stays between −6.02 and −7.12 dB and
the stopband reaches at most −46.77 dB.

All of seeds 1–8 meet the mask; seven do it with four parts and one of the
two textbook elliptic forms — the trap across the series inductor (seeds
1, 3, 6) or a series-LC trap to ground (seeds 4, 7, 8) — and seeds 2 and 5
with six parts (`doc/M8.md`).

`best.cir` opens in ngspice, LTspice and QSpice as it is.
