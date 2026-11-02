# 1. Circuits and the ledger

kronspice is, before anything else, a circuit simulator: a netlist goes
in and node voltages, element currents and an energy ledger come out. This
chapter runs four small circuits through it and reads every line the
engine prints, because the habits formed here, especially the habit of
reading the ledger, are what the later chapters rely on. The three
examples it draws on are the config-driven experiments of milestone M5:
[rc_lowpass_values](../examples/rc_lowpass_values),
[negative_resistance_gain](../examples/negative_resistance_gain) and
[rectifier_filter](../examples/rectifier_filter). Each of them evolves
component values; here we look only at the circuits and what the engine
says about them, and leave how the values were found to chapter 2.

## 1.1 A circuit as the engine sees it

The smallest netlist worth running is a resistive divider:

```
Divider
V1 in 0 DC 10
R1 in out 1k
R2 out 0 1k
.op
.end
```

The first line is the title. Each element line is a name whose first
letter gives the kind (V for a voltage source, R for a resistor), two
nodes, and a value with the usual SPICE suffixes (`1k` is 1000 Ω). Node
`0` is ground. `.op` asks for the operating point, the DC solution, and
`.end` closes the deck. Run it:

```
go run ./cmd/kronspice divider.cir
```

```
Divider
Analysis: op

Node voltages
  v(in)   10 V
  v(out)  5 V

Elements
  name  kind            nodes   voltage  current  power absorbed
  V1    voltage source  in→0    10 V     -5 mA    -50 mW
  R1    resistor        in→out  5 V      5 mA     25 mW
  R2    resistor        out→0   5 V      5 mA     25 mW
```

The engine formulates the circuit by modified nodal analysis (Ho, Ruehli
and Brennan 1975): the unknowns are the voltages of every node except
ground plus the current through every voltage source, and each element
stamps its contribution into one matrix. For a linear circuit one LU
solve gives the answer. Currents are reported from the element's first
node to its second, so V1's −5 mA means 5 mA flows out of its positive
terminal, and its absorbed power is negative: it supplies 50 mW, which the
two resistors absorb as 25 mW each.

With `-json` the same run gives one document that a program can read
without parsing the table: `analyses[0].op.nodes.out` is 5, and each
element appears with `voltage`, `current`, `power_absorbed`,
`stored_energy` and `model_supplied`. The fields are listed in
[llms.txt](../llms.txt).

## 1.2 The energy ledger

The report does not stop at the elements. It ends with the part of
kronspice that makes it different from a conventional simulator:

```
Energy ledger (+ into the circuit from the exterior, − out)
  element  channel     status   power   note
  V1       electrical  modeled  50 mW   impressed source: exchanges energy with the exterior
  R1       electrical  modeled  0 W     terminals connect inside the network
  R1       thermal     modeled  -25 mW  Joule heating i²R to surroundings
  R2       electrical  modeled  0 W     terminals connect inside the network
  R2       thermal     modeled  -25 mW  Joule heating i²R to surroundings
  (every other channel of each element is not modeled: 0 W each; use -all-channels to list them)

  inflow        50 mW
  outflow       0 W
  dissipation   50 mW
  storage rate  0 W
  residual      0 W (relative 0)
```

A circuit is an open system. Energy comes in from outside through
sources, and leaves as heat, as radiation, or through a load; what
Maxwell's equations require is not that the energy inside stays constant
but that every joule crossing the boundary is accounted for (Poynting's
theorem, in circuit form Tellegen's theorem; CLAUDE.md §3c). So kronspice
keeps a ledger with one row per element and per **channel**. The
electrical channel is the circuit's own terminals; a source's electrical
row is its exchange with the exterior, and a resistor's electrical row is
zero because its terminals connect inside the network. The thermal
channel carries Joule heat out to the surroundings.

Every element has a port on every channel in the registry, whether or not
a coupling is modeled for it: light, magnetic and electric field,
mechanical, acoustic, ionizing radiation, chemical, neutrino,
gravitational, and a hypothetical column (CLAUDE.md §3.1). The line in
parentheses says what happened to all of those here: they are present,
report exactly 0 W, and are marked "not modeled". `-all-channels` lists
them. The point of carrying them is that a report always shows which
energy paths were considered and which were not, and that a model added
later for one of them plugs into an existing port.

The totals close the books:

```
residual = inflow + model-supplied − outflow − dissipation − storage rate
```

Here inflow is V1's 50 mW, dissipation is the resistors' 50 mW, and the
residual is zero. The relative residual, residual over the sum of all row
powers, is what the tests check in every circuit and every analysis in
the repository, and what the genetic algorithm checks for every candidate
it evaluates. A residual above round-off is a bug in the numerics or the
model. It is never a result, and chapter 10 is largely about why that
rule matters.

## 1.3 Frequency: the RC low-pass

[rc_lowpass_values](../examples/rc_lowpass_values) is a first-order RC
filter driven by a 1 V AC source. With the values the genetic algorithm
settled on, 160 Ω and 1 µF, the netlist is:

```
RC low-pass, evolved values
V1 in 0 AC 1
R1 in out 160
C1 out 0 1u
.ac dec 10 10 100k
.end
```

`.ac dec 10 10 100k` sweeps ten points per decade from 10 Hz to 100 kHz.
The engine linearizes the circuit at its operating point and solves the
complex system once per frequency. Read from the JSON document,
`ac.v.out` holds a magnitude and a phase at each `ac.freq`:

| f | \|V(out)\| | phase |
|---|---|---|
| 10 Hz | 0.000 dB | −0.6° |
| 100 Hz | −0.044 dB | −5.7° |
| 1 kHz | −3.033 dB | −45.2° |
| 10 kHz | −20.09 dB | −84.3° |
| 100 kHz | −40.05 dB | −89.4° |

The corner, where the response is 3 dB down and the phase is −45°, sits
at 1/(2πRC) = 994.7 Hz; the sweep point at 1 kHz is just past it. That
0.53 % miss from the 1 kHz target is not the search's fault: 160 Ω × 1 µF
is the closest any product of an E24 resistor and an E12 capacitor comes
to 1 kHz, so the result is the floor of the value grid (doc/M5.md).
The AC ledger closes at every frequency as an average-power balance.

## 1.4 Negative resistance and model-supplied power

Can a two-resistor divider have gain? [negative_resistance_gain](../examples/negative_resistance_gain)
asks the genetic algorithm for a gain of two at 1 kHz from

```
V1 in 0 DC 1 AC 1
R1 in out 1k
R2 out 0 ?
```

and the search, allowed to flip the sign of R2, settles at −2000.86 Ω,
within 0.04 % of the exact answer. With R2 = −2 kΩ the transfer function
is R2/(R1 + R2) = −2000/(1000 − 2000) = 2. kronspice accepts negative
resistors, as ngspice does, and the operating point of that circuit
shows what they mean:

```
Elements
  name  kind            nodes   voltage  current  power absorbed
  V1    voltage source  in→0    1 V      1 mA     1 mW
  R1    resistor        in→out  -1 V     -1 mA    1 mW
  R2    resistor        out→0   2 V      -1 mA    -2 mW

Energy ledger
  V1  electrical  modeled  -1 mW  impressed source: exchanges energy with the exterior
  R1  thermal     modeled  -1 mW  Joule heating i²R to surroundings
  R2  electrical  modeled  0 W    idealized negative resistance: delivers v²/|R| as model-supplied power

  inflow        0 W
  outflow       1 mW
  dissipation   1 mW
  model-supplied 2 mW (idealized active elements; not inflow from the exterior)
    R2  2 mW
  residual      0 W (relative 0)
```

The output node sits at 2 V with the input at 1 V. R2 delivers
v²/|R| = 4/2000 = 2 mW; 1 mW heats R1 and 1 mW is pushed back into the
source, which now absorbs power (its row is −1 mW, an outflow). The
books balance, but look where the 2 mW is booked. It is not inflow. An
ideal negative resistance stands for a real device, a tunnel diode, a
Gunn diode or an op-amp negative-impedance converter, that takes this
energy from a bias supply the model leaves out. So the ledger has a
separate column, **model-supplied**, listed per element, and the inflow
total stays at zero.

The distinction is between two transformations that are easy to confuse
(CLAUDE.md §3c). Duality, exchanging the roles of voltage and current as
a gyrator does, maps a capacitor to an inductor and a resistance R to a
conductance 1/R, still a positive resistor. Sign inversion, reversing the
current, maps R to −R, dissipation to generation, and leaves C and L
alone. A negative resistor is the second kind: an active element in
disguise. A fitness function that measures harvested energy therefore
has to count the model-supplied column as an input, and the objective
`evo.ModelSuppliedLimit` exists to bound it, so that evolution cannot
"discover" free energy by evolving a sign.

## 1.5 Nonlinear and in time: the rectifier

[rectifier_filter](../examples/rectifier_filter) is the first circuit in
this chapter that a single matrix solve cannot handle: a 10 V peak, 1 kHz
sine through a diode into a capacitor and a load. With the evolved
values:

```
Half-wave rectifier, evolved values
V1 in 0 SIN(0 10 1k)
D1 in out dr
R1 out 0 270
C1 out 0 820u
.model dr D(IS=1e-14 N=1 RS=0.119)
.tran 10u 8m
.end
```

The diode is the SPICE3 level-1 model named by the `.model` card, with
its saturation current, emission coefficient and series resistance. At
every time step the engine solves the nonlinear equations by Newton's
method with SPICE's junction voltage limiting, and when plain Newton
fails it falls back through gmin stepping, source stepping and
pseudo-transient continuation (milestone M4). `.tran 10u 8m` asks for
8 ms of transient with a 10 µs print step; the integration is
trapezoidal with backward Euler at the start and at source breakpoints,
and the step size is chosen from the local truncation error. The report
begins with how that went:

```
Time points  843 recorded, 842 steps accepted, 7 rejected
Step size    201.806 ns to 10 µs
Interval     0 s to 8 ms

Node voltages at the end
  v(out)  9.03021 V
```

Read from the JSON document, which resamples the waveforms on the 10 µs
grid, the output over the last half of the run, 4 ms to 8 ms, averages
9.008 V with 0.153 V of peak-to-peak ripple. Those two numbers were the
experiment's objectives: a 9 V mean with no more than 0.2 V of ripple.
The time constant R1·C1 = 221 ms against the 1 ms period is what holds
the ripple down, and the mean is then set by the crest, 10 V less the
diode drop at the charging peak, which the evolved series resistance of
0.119 Ω keeps small.

The ledger of a transient is a balance at every time point plus energies
integrated over the run:

```
Energy ledger over the interval
  energy in            61.7434 mJ
  energy out           31.6088 fJ (returned to sources and ports)
  energy dissipated    28.3105 mJ
  stored: start → end  2.86867e-63 J → 33.4333 mJ
  integration mismatch -471.007 nJ (time-stepping error, not a physical channel)
  worst point residual 0 (relative)
```

Every line can be checked. The capacitor ends the run at 9.03 V, so it
holds ½CV² = ½ × 820 µF × 9.03² = 33.43 mJ, the stored energy reported.
The source delivered 61.74 mJ; 28.31 mJ left as heat in the load, the
diode's junction and its series resistance; the rest is in the
capacitor. The books close to −471 nJ, which is 7.6 parts per million of
the input. That mismatch is not a missing channel: the per-point residual
is zero, so the power balance holds at every accepted step, and the
difference comes from integrating those powers with the trapezoidal
rule over the time grid. The report labels it as time-stepping error
for exactly that reason, so that nobody reads it as physics.

One more thing the experiment shows is that a `.model` parameter can be
a gene like any part value: the diode's RS evolved alongside C1 and R1,
and the engine re-parses the written netlist, so the result opens
unchanged in ngspice, LTspice and QSpice.

## 1.6 What to keep

- A netlist in, and node voltages, element currents and a ledger out,
  in a text report for people or a JSON document for programs.
- The ledger has a row for every element on every channel. Modeled rows
  carry power; unmodeled rows carry exactly zero and say so.
- The residual must be at round-off in every analysis. A residual is a
  bug to find, never a finding.
- Power from an idealized active element is model-supplied. It balances
  the books and is never inflow from outside.
- A transient's integration mismatch measures the time stepping, and is
  reported apart from the physics.

## Run it yourself

```
go run ./cmd/kronspice divider.cir                 # 1.1 and 1.2
go run ./cmd/kronspice -json rc.cir                # 1.3, read the ac section
go run ./cmd/kronspice neg.cir                     # 1.4
go run ./cmd/kronspice -json rect.cir              # 1.5, read the tran section
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json
go run ./cmd/evolve -config examples/negative_resistance_gain/config.json
go run ./cmd/evolve -config examples/rectifier_filter/config.json
```

The three `evolve` runs reproduce the values used above; each example's
README gives its genes, objectives and the run's statistics.

## Sources

- C.-W. Ho, A. E. Ruehli, P. A. Brennan, "The modified nodal approach to
  network analysis," *IEEE Trans. Circuits and Systems* 22, 504–509,
  1975.
- L. W. Nagel, *SPICE2: A Computer Program to Simulate Semiconductor
  Circuits*, UC Berkeley ERL Memo M520, 1975.
- B. D. H. Tellegen, "A general network theorem, with applications,"
  *Philips Research Reports* 7, 259–269, 1952.
- J. C. Maxwell, *A Treatise on Electricity and Magnetism*, 1873, Vol. II
  (the field equations behind the ledger; CLAUDE.md §3c).
- doc/M4.md (diode, Newton and the fallback chain), doc/M5.md (the
  value-genome engine, penalties, negative resistance).
