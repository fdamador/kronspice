# Bridges of water diodes and spark gaps, for different jobs

[sparkgap_diode](../sparkgap_diode) built two kinds of one-way device from
water: the **water cell** (aluminium against lead or iron in a phosphate
solution, the 1890s rectifier of Pollak and Nodon: one-way up to about
140 V) and the **polarity spark gap** (a point–plane gap in water that
breaks down at a lower voltage with the point positive). This example
combines them, with plain air spark gaps, into circuits for six jobs.

```
go run ./examples/sparkgap_diodebridge -viz examples/sparkgap_diodebridge/viz   # about 20 seconds
python3 examples/sparkgap_diodebridge/viz/render.py
go test ./examples/sparkgap_diodebridge
```

| Job | Circuit | Result |
|---|---|---|
| A. DC supply from the mains | 4 lead–aluminium cells | 88 V DC from 100 V peak, 85 % |
| B. Surge protection | air-gap arrester + aluminium–aluminium clamp in front of A | bridge cells stay at 117 V (unprotected: 221 V, sparking over) |
| C. Spark pulser | water spark-gap bridge charges a capacitor; an air gap fired at each mains zero crossing dumps it | 120 pulses/s, 505 A peak, 50 % of mains power to the load |
| D. Voltage doubler | Greinacher circuit of cells, 2 per position | 174 V DC from 100 V peak, 85 % |
| E. Varactor with memory | the aluminium oxide film tunes a 1 mH tank | 643 Hz to 2,405 Hz, set by the film's forming voltage |
| F. Frequency doubler | any bridge: 60 Hz in, 120 Hz out | 95 % of the output's AC power at 120 Hz (cells) |

Every run has an energy ledger: the cells and gaps add no energy; every
joule at a load comes from the mains or the surge source.

## A. Mains DC supply

Four lead–aluminium cells, 100 V peak at 60 Hz, 0.5 Ω source, 2,200 µF,
50 Ω load: **87.9 V DC**, 5.2 V ripple, 85 % efficient. The same circuit in
kronspice gives 87.96 V (0.01 % apart), which checks this example's own
solver (below).

## B. Surge-protected supply

A 4 kV, 1.2/50 µs surge from a 2 Ω source (the IEC 61000-4-5 test-wave
shape and impedance; its 8/20 µs current shape is not modelled) arrives on
the mains line at the mains peak. Protection, in the usual two stages:

1. a **0.2 mm air gap** across the line (breaks down at 1,423 V, Paschen's
   law for air) takes the bulk of the surge;
2. after 10 µH and 0.5 Ω of line, an **aluminium–aluminium water cell**
   across the bridge input. Both of its plates are valve metal, so it
   blocks both ways up to 140 V and then sparks over — a two-way clamp,
   like a pair of Zener diodes.

| | Largest voltage on a bridge cell | Extra energy into the bridge cells | Arrester | Clamp |
|---|---|---|---|---|
| Unprotected | **221 V** (sparks over at 140 V) | 15.4 J | — | — |
| Protected | **117 V** | 0.64 J | 22.9 J (1 strike) | 0.006 J |

## C. Spark pulser

A 6 kV, 60 Hz transformer (1 kΩ) feeds a bridge of four water polarity gaps
(polarity ratio 3, the ratio that rectified in sparkgap_diode), which
charges 0.1 µF. An air gap discharges the capacitor through 1 µH into a
load. Three designs:

| Output switch and load | Pulses/s | Peak current | Load power ÷ mains power | What happens |
|---|---|---|---|---|
| 1 mm air gap fires by itself, 1 Ω (rings) | 120 | 1,200 A | 0.6 % | the capacitor swings negative and the bridge gaps re-strike backwards (15 J): the transformer is shorted |
| 1 mm air gap fires by itself, 8 Ω (overdamped) | 120 | 482 A | 1.8 % | the mains keeps the output arc lit after the pulse (power follow) |
| **2 mm air gap fired at each mains zero crossing, 8 Ω** | **120** | **505 A** | **50 %** | the transformer cannot follow at its zero; 1.1 J per pulse |

Firing at a fixed point of the mains cycle is what Tesla-coil builders do
with a synchronous rotary spark gap. Here the trigger is ideal (declared).
In the working design the bridge gaps still make tiny backward sparks, about
0.02 J in total, because the transformer and the full capacitor nearly
balance.

## D. Voltage doubler

Greinacher's circuit (a pump capacitor and two diodes; 2,200 µF, 200 Ω) with
lead–aluminium cells, 100 V peak:

| Cells per position | Output | Efficiency | Reverse peak |
|---|---|---|---|
| 1 | 151 V | 17 % — sparks over | 160 V |
| 2 | **174 V** | **85 %** | 182 V |

Each position sees twice the peak in reverse, so it needs two cells, as in
the historical "two or more valves in series for higher pressures".

## E. A varactor with memory

The aluminium plate's oxide film is a capacitor whose thickness is set when
it is formed: about 1.3 nm per volt (declared, typical barrier anodic
alumina). Formed at a higher voltage, the film is thicker and its
capacitance smaller. A 1 mH tank across a cell biased in its blocking
direction, solved by kronspice AC:

| Formed at | Film | Resonance (kronspice) | 1/(2π√LC) | Q |
|---|---|---|---|---|
| 10 V | 61.3 µF | 643 Hz | 643 Hz | 24 |
| 20 V | 30.6 µF | 909 Hz | 909 Hz | 33 |
| 40 V | 15.3 µF | 1,285 Hz | 1,286 Hz | 48 |
| 70 V | 8.8 µF | 1,701 Hz | 1,701 Hz | 62 |
| 140 V | 4.4 µF | 2,405 Hz | 2,405 Hz | 87 |

Biased the other way (−5 V) the cell conducts and the tank has no
resonance. **This is not a semiconductor varactor**: its capacitance does
not follow the bias, and below the forming voltage it stays fixed. It is set
by the highest voltage the film has seen — a capacitor with memory, tuned by
forming. Re-forming higher thickens the film at once; thinning it again
takes the slow dissolution of the oxide, which is not modelled. Q is low at
low forming voltage because the film then leaks more (declared leakage).

## F. Frequency doubler

A full bridge's output repeats twice per cycle: 60 Hz in, 120 Hz out.
Amplitudes divided by the source peak, resistive load:

| Bridge | DC | 120 Hz | 240 Hz | Share of the AC power at 120 Hz |
|---|---|---|---|---|
| perfect diode bridge | 0.637 | 0.424 | 0.085 | 95 % |
| lead–aluminium cells, 100 V | 0.580 | 0.411 | 0.081 | 95 % |
| water spark gaps, 6 kV | 0.495 | 0.534 | 0.108 | 89 % |

A filter tuned to 120 Hz picks out the doubled frequency. The spark-gap
bridge's sudden jumps put more into 120 Hz and the higher harmonics.

Picture: `viz/diodebridge.png` (surge, pulser, varactor, doubler).

## How it is solved

Circuits with only cells run in kronspice (D, E, and the twin of A). A spark
gap is a switch with memory — it strikes at its breakdown voltage and stays
on while its arc carries current — which kronspice's static elements
cannot hold, so circuits with gaps run in this example's small nodal solver
(`sim.go`): backward Euler for capacitors and inductors, Newton on the
cells' curves, gap states updated after each step. It is checked against an
exact RC discharge and against kronspice on the supply. Its energy ledger
splits what capacitors and inductors take in into stored energy and
backward Euler's numerical damping, and closes to a millionth of the
energy supplied or better.

Spark gaps (declared unless cited):

- water polarity gap: 3 mm point–plane; breakdown with the point positive
  from Stygar's water relation (extrapolated to 60 Hz), negative ×3; arc
  50 V + 0.5 Ω, out below 0.1 A after 10 µs;
- after an arc the breakdown recovers as 1 − e^(−t/1 ms), and never below
  the Paschen minimum of water vapour, 480 V (Škoro et al. 2011);
- air gaps: Paschen's law for air (Raizer), arc 30 V + 0.2 Ω, never below
  air's Paschen minimum, about 330 V.

## Files

| File | What it holds |
|---|---|
| `cell.go`, `gap.go` | the water cell and polarity gap (copied from sparkgap_diode; keep in step) |
| `sim.go` | the nodal solver for circuits with spark gaps, with its energy ledger |
| `apps.go` | the six circuits |
| `main.go` | the report |
| `main_test.go` | model match, solver checks, the twin, and each job's result |
| `viz/render.py` | `diodebridge.png` |

## Sources

- *Encyclopaedia Britannica*, 11th ed. (1911), vol. 27, "Transformers" — Pollak and Nodon aluminium rectifiers.
- H. Greinacher, *Z. Phys.* 4, 195 (1921) — the voltage doubler.
- IEC 61000-4-5 — surge test waveform (1.2/50 µs, 2 Ω).
- Yu. P. Raizer, *Gas Discharge Physics* (Springer, 1991) — Paschen's law for air.
- W. A. Stygar et al., *Phys. Rev. ST Accel. Beams* 9, 070401 (2006) — water breakdown.
- M. Škoro et al., *Phys. Rev. E* 84, 055401 (2011) — Paschen curve of water vapour.
