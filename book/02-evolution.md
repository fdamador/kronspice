# 2. Evolution

kronspice exists to sit inside a loop that designs circuits by trial: a
genetic algorithm proposes thousands of candidates, the simulator scores
each one, and the best survive to be varied again. The engine was built
for that job (CLAUDE.md §4): it simulates in memory with no files or
processes per candidate, returns a typed error for every way a generated
circuit can be broken instead of stalling, and runs candidates in
parallel with nothing shared. This chapter follows that loop through four
examples: a transistor stage whose eight standard values are evolved
([amplifier_gain](../examples/amplifier_gain)), a filter grown from an
empty fixture ([filter_topology](../examples/filter_topology)), a MOS
stage evolved on three objectives at once
([mos_amplifier_pareto](../examples/mos_amplifier_pareto)), and the same
loop run with Kron's diakoptics so that only the mutated part of a
circuit is re-solved ([torn_ga](../examples/torn_ga)).

## 2.1 Run it: the amplifier before and after

The starting point of the first example is an ordinary common-emitter
stage, a 2N3904 on 12 V, driven from a 600 Ω source into a 10 kΩ load,
with the values any first guess might use:

```
VCC vcc 0 12
VIN in 0 DC 0 AC 1
RS in s 600
CIN s b 1u
R1 vcc b 47k
R2 b 0 10k
RC vcc c 4.7k
RE1 e e2 100
RE2 e2 0 1k
CE e2 0 100u
CF c b 10p
COUT c o 10u
RL o 0 10k
Q1 c b e Q2N3904
.model Q2N3904 NPN(IS=6.734f BF=416.4 VAF=74.03 ... RB=10)
.ac dec 20 10 100meg
```

The transistor model is the Fairchild datasheet's Gummel–Poon parameter
set, solved with the SPICE3 equations in `device/bjt.go`. Run through
the command-line tool, the sweep's two ends read:

```
Sweep  141 points, 10 Hz to 100 MHz

  node    at 10 Hz              at 100 MHz
  v(in)   0 dB ∠ 0°             0 dB ∠ 0°
  v(c)    14.98 dB ∠ -78.79°    -11.8 dB ∠ 38.2°
  v(o)    14.87 dB ∠ -69.74°    -11.8 dB ∠ 38.2°

Average-power ledger: worst relative residual 0 over the sweep
```

The guess gives about 15 dB. The experiment asks for 30 dB at 10 kHz, a
passband from 100 Hz to 200 kHz, and a collector near 6 V with the
emitter near 1.5 V. After evolution the written `best.cir` has

```
R1 vcc b 30000
R2 b 0 6800
RC vcc c 5100
RE1 e e2 68
RE2 e2 0 1300
CE e2 0 1.8e-05
CF c b 3.3e-11
CIN s b 4.7e-06
```

and, re-simulated through `kronspice -json` and read from the document:
29.917 dB at 10 kHz, the −3 dB edges at 100.1 Hz and 199.6 kHz, and an
AC ledger residual of zero at every one of the 141 frequencies. Every
value is a standard part: resistors from the E24 series, capacitors from
E12.

## 2.2 The value genome

The simplest genome fixes the topology and evolves values. Each gene
names an element (or a `.model` parameter, as the rectifier of chapter 1
did with its diode's series resistance), a range, and optionally an
E-series to snap to. Values live on a log scale, so a mutation of
σ = 0.1 decade means "about 25 % either way" at any magnitude, and
crossover blends the parents' exponents (BLX-α, Eshelman and Schaffer
1993). Snapping happens at decode, so the circuit the simulator sees is
always buildable.

| Gene | Range | Series |
|---|---|---|
| R1, R2 | 10 kΩ – 1 MΩ, 1 kΩ – 100 kΩ | E24 |
| RC, RE1, RE2 | 470 Ω – 47 kΩ, 10 Ω – 2.2 kΩ, 100 Ω – 10 kΩ | E24 |
| CF, CIN, CE | 1 pF – 1 nF, 100 nF – 100 µF, 1 µF – 1 mF | E12 |

The five objectives are measurements made on the simulation (gain at a
frequency, the two −3 dB edges, two operating-point voltages), each
turned into a relative error and summed with weights of 1, 1, 1, 0.2 and
0.2. The run on this machine, population 100 for 300 generations, made
29,402 evaluations in 15.8 s and stopped with a weighted error of 0.0202:
gain 0.0041, upper edge 0.0011, lower edge 0.00003, collector voltage
0.070 and emitter voltage 0.005 (the last two at weight 0.2). That total
is the floor of the value grid, not of the search: with parts in 4 %
steps and five coupled targets, no single-part change of the best design
improves it (doc/M7.md). Ten of twelve seeds meet the acceptance
tolerance of 30 ± 0.5 dB, edges within 5 % and V_CE above 3 V; the other
two miss the gain by under 0.1 dB beyond it.

## 2.3 Failure is a penalty, never a crash

A random or mutated circuit is often broken: a floating node, a loop of
voltage sources, a diode biased so hard that Newton diverges, a time step
that collapses, or a design so large it never finishes. A simulator that
stalls or panics on any of these stops the whole generation. kronspice
returns a typed error for each, and `evo.Penalty` maps the error to a
fitness worse than any working design, graded so that evolution can
still tell "almost works" from "broken":

| Class | Penalty | Errors |
|---|---|---|
| no-measurement | 1e6 | the objective's quantity is undefined (no −3 dB edge in the sweep) |
| no-convergence | 2e6 | Newton, a too-small time step, too many steps |
| singular | 3e6 | singular or non-finite matrix |
| invalid | 4e6 | floating node, bad value or model, self-connected element, not quasi-static |
| ledger | 5e6 | the energy ledger does not balance |
| timeout | 6e6 | the per-evaluation deadline passed |
| error, panic | 1e7, 2e7 | anything else; a panic is caught and reported as a bug |

Two of these deserve a word. The ledger class is deliberately severe: a
candidate whose energy does not add up is a numerical artefact, and
evolution exploits artefacts eagerly if they score well. The timeout is a
`context` deadline that reaches into every Newton iteration, time step and
frequency point, so one pathological candidate costs at most its
allotted seconds. In the filter run below, between 7 and 53 of each
generation's 300 candidates took a penalty, and none stopped the run.

Determinism is the other guarantee. Every random draw comes from a PCG
stream derived from the seed and the generation number, results are
stored by individual index before selection, and so a run is
bit-identical for any number of workers and resumes exactly from its
`checkpoint.json` (CLAUDE.md §1 goal 5; `TestTopologyResume`).

## 2.4 Growing a topology

The second example starts with no circuit at all. Its embryo, after Koza
et al. (1999, ch. 25), is only the test fixture:

```
VIN vs 0 DC 0 AC 1
RS vs in 1k
RL out 0 1k
.ac dec 20 10 100k
```

Nothing connects `in` to `out`, and the tool says so:

```
  node    at 10 Hz      at 100 kHz
  v(vs)   0 dB ∠ 0°     0 dB ∠ 0°
  v(in)   0 dB ∠ 0°     0 dB ∠ 0°
  v(out)  -Inf dB ∠ 0°  -Inf dB ∠ 0°
```

The genome (`evo/encode_topology.go`) is a variable-length list of parts,
each with a kind (L, C or R), two nodes and a log value, decoded into a
copy of the embryo. Parts may use ground, the two terminals and eight
internal nodes. Mutation can change a value, add or delete a part,
rewire one end, change the kind, split a part in series through a new
node or add one in parallel; crossover is one-point on the lists. At
decode, parts that carry no current (self-loops, parts with an end on an
internal node nothing else touches) are pruned, because they change no
voltage and would otherwise turn every bit of mutation debris into a
floating-node penalty. The objective is a mask: between −7.5 and −5.5 dB
from 10 Hz to 1 kHz (a matched passband sits at −6.02 dB), and at most
−46 dB from 3 kHz to 100 kHz. A design inside the mask scores 0; outside
it scores at least 0.1, so "under 0.1" reads "spec met". A parsimony term
of 0.002 per part then keeps removing parts that are not needed.

With population 300 and seed 1, the run on this machine reached its
final design by generation 122 after 36,656 evaluations (the README's
full 200-generation run reaches the same one): fitness 0.008, which is a
mask error of 0 plus four parts.

```
L_g1 in  out 0.2
C_g2 0   in  1.8e-07
C_g3 out 0   2.2e-07
C_g4 in  out 1.1e-08
```

That is a third-order elliptic low-pass: a series inductor with a
capacitor across it, which makes a notch at 1/(2π√(LC)) = 3.39 kHz, and
shunt capacitors at both ends. Run through the tool and read from the
JSON document:

| f | \|V(out)\| |
|---|---|
| 10 Hz | −6.02 dB |
| 1 kHz | −7.12 dB |
| 2.0 kHz | −27.2 dB |
| 3.2 kHz | −53.9 dB |
| 3.5 kHz | −60.0 dB |
| 10 kHz | −49.0 dB |
| 100 kHz | −68.0 dB |

Over the whole sweep the passband stays between −7.12 and −6.02 dB and
the stopband never rises above −46.78 dB. Third order is the least that
meets about 1 dB of ripple and 40 dB of attenuation at three times the
passband edge (a Butterworth needs five reactive parts), and seven of
eight seeds find a four-part elliptic form, the trap either across the
series arm or as a series-resonant shunt. The parsimony pressure needs
time: with an early stop at "spec met", seed 1 had 15 parts; run on,
it pared them to four.

## 2.5 More than one objective

The third example asks for three things that pull against each other:
gain, bandwidth and supply power, for an NMOS common-source stage on 5 V
with a Level 3 process model, with the gate divider, drain and source
resistors (E24) and the channel width and length (continuous) as genes.
A weighted sum would answer with one design and hide the trade-off.
NSGA-II (Deb et al. 2002, `evo/nsga2.go`) keeps instead the whole
**Pareto front**: the designs none of which is beaten on every objective
by another. A fourth objective, gain of at least 20 dB, is marked as a
constraint; a design that violates it is dominated by every design that
meets it, and failed simulations rank after all of them by their penalty
class.

The run on this machine, population 100 for 100 generations, made
10,000 evaluations in 3.3 s and returned a front of 99 designs, every one
above 20 dB, written as `front.csv` and one netlist each in `front/`.
Three of them, each re-simulated through the tool:

| Design | Gain at 10 kHz | −3 dB edge | Supply power | RD | RSS | W / L |
|---|---|---|---|---|---|---|
| highest gain | 34.71 dB | 629 kHz | 166 µW | 130 kΩ | 13 kΩ | 500 µm / 2.30 µm |
| widest band | 20.22 dB | 32.4 MHz | 17.4 mW | 1 kΩ | 240 Ω | 497 µm / 1.27 µm |
| lowest power | 20.68 dB | 355 kHz | 16.1 µW | 1 MΩ | 36 kΩ | 74 µm / 2.41 µm |

The shape is the textbook one. Bandwidth costs power, because a small
drain resistor against the 5 pF load needs current to keep its gain. Gain
tops out near 35 dB on every seed, because the 100 kΩ load caps the
effective drain resistance while gm/I_D cannot exceed its weak-inversion
value, so the highest-gain designs use the widest device the genes allow
in moderate inversion. (The example's README quotes the milestone's
front, whose extreme designs differ in detail from this run's; the
envelope is the same.) On the ZDT1 and ZDT2 benchmarks the same code
places its front within 0.0065 and 0.0040 of the true one, ZDT2 being the
non-convex case where a weighted sum can only find the ends.

## 2.6 Tearing inside the loop

A mutation usually changes one part of a circuit, yet an ordinary solver
factors the whole matrix again for every candidate. Kron's diakoptics
(chapter 3; CLAUDE.md §3a) tears the circuit into pieces, factors each
on its own and joins them through a small interconnection system, so the
pieces a mutation did not touch can be reused from a cache. The fourth
example makes the comparison honest by running the same GA three times
from the same seed: a two-pole RC pre-filter (four genes) in front of a
fixed plant, a 40-section LC line of 125 unknowns, joined by one resistor
RT, evolved for a 5 kHz cutoff at the plant's output.

The exported netlist runs through the tool in 0.18 s for its 201
frequencies, and its results sheet shows the tearing: two pieces of 4 and
121 unknowns, one link, an interconnection system of order 1, the link
current through RT 0.4546 mA ∠ −20.73° at 10 kHz, torn and untorn node
voltages equal to 4.5e-15, and the power taken by the R, L and C branches
identical in both, which is Roth's invariance. The three GA runs on this
machine, population 40 for 30 generations:

| Solver | Time | Best design | Fitness |
|---|---|---|---|
| MNA (sparse LU) | 1.19 s | R1 13 Ω, C1 47 nF, R2 3 kΩ, C2 390 nF | 5.645e-6 |
| Kron node form, untorn | 26.51 s | same | same |
| Kron diakoptics, torn at RT | 2.65 s | same | same |

Tearing makes the Kron path 10.0× faster than not tearing, and the plant's
121-unknown piece came from the cache about as often as the pre-filter's
4-unknown piece was factored (233,263 hits against 225,821
factorizations). It is still 2.2× slower than the sparse MNA engine, whose
Gilbert–Peierls LU reuses its pattern and is cheap on a chain; what
remains per frequency in the torn path is building the primitive network
and the node-pair rows, which doc/M9b.md lists as the next thing to
cache. So MNA is the default and tearing is opt-in, chosen when the
fixed part of a problem is large and its factorization fills in.

## 2.7 What to keep

- A gene is an element value or a model parameter on a log scale,
  snapped to a standard series at decode; the E-series grid, not the
  search, usually sets the floor of the error.
- Every failure has a class and a penalty. The ledger class is severe on
  purpose, and a deadline bounds every candidate.
- A run is reproducible from its seed for any number of workers, and
  resumes exactly from its checkpoint.
- A topology genome needs pruning of dead parts and a parsimony term;
  given both, four standard parts reproduce a textbook elliptic filter
  from nothing.
- When objectives conflict, keep the front, not a weighted winner, and
  express requirements as constraints.
- Tearing wins against not tearing by an order of magnitude; against
  sparse MNA it does not yet, and the measurement says exactly why.

## Run it yourself

```
go run ./cmd/kronspice examples/amplifier_gain/circuit.cir
go run ./cmd/evolve -config examples/amplifier_gain/config.json          # ≈ 16 s
go run ./cmd/evolve -config examples/filter_topology/config.json         # ≈ 10 s to the four-part design
go run ./cmd/evolve -config examples/mos_amplifier_pareto/config.json    # ≈ 3 s, writes front.csv and front/
go run ./examples/torn_ga                                                # ≈ 30 s, three solvers, same best
go run ./cmd/kronspice -json runs/<run>/best.cir                         # re-simulate any result
```

Each run writes `runs/<timestamp>-<name>/` with the config, `stats.csv`,
the best netlist per generation and `checkpoint.json`; every `best.cir`
opens in ngspice, LTspice and QSpice as it is.

## Sources

- K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A fast and elitist
  multiobjective genetic algorithm: NSGA-II," *IEEE Trans. Evol. Comput.*
  6(2), 182–197, 2002.
- J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, *Genetic
  Programming III: Darwinian Invention and Problem Solving*, Morgan
  Kaufmann, 1999, ch. 25.
- L. J. Eshelman, J. D. Schaffer, "Real-coded genetic algorithms and
  interval-schemata," *Foundations of Genetic Algorithms* 2, 1993.
- H. K. Gummel, H. C. Poon, "An integral charge control model of bipolar
  transistors," *Bell Syst. Tech. J.* 49, 827, 1970; Fairchild 2N3904
  datasheet, SPICE model section.
- A. Vladimirescu, S. Liu, "The simulation of MOS integrated circuits
  using SPICE2," UCB/ERL M80/7, 1980 (the Level 3 model).
- G. Kron, *Diakoptics: The Piecewise Solution of Large-Scale Systems*,
  Macdonald, 1963; P. W. Aitchison, "Diakoptics or tearing: a
  mathematical approach," *Quart. Appl. Math.* 41, 265–272, 1983;
  J. P. Roth, "An application of algebraic topology: Kron's method of
  tearing," *Quart. Appl. Math.* 17, 1959.
- J. R. Gilbert, T. Peierls, "Sparse partial pivoting in time proportional
  to arithmetic operations," *SIAM J. Sci. Stat. Comput.* 9, 862, 1988.
- [doc/M5.md](../doc/M5.md), [doc/M7.md](../doc/M7.md),
  [doc/M8.md](../doc/M8.md), [doc/M9.md](../doc/M9.md),
  [doc/M9b.md](../doc/M9b.md): the milestone reports behind each example.
