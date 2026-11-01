# A low-pass filter grown from an empty embryo

The M8 example (CLAUDE.md §9; [doc/M8.md](../../doc/M8.md)): topology
evolution (§4, "Topology genome"). The embryo (`embryo.cir`) is only the
test fixture, after Koza et al., *Genetic Programming III* (1999), ch. 25:
a 1 V source with 1 kΩ source resistance driving node `in`, and a 1 kΩ load
on node `out`. Nothing connects `in` to `out`. Evolution decides which
parts exist (inductors, capacitors, resistors, E24 values), how they are
wired among ground, `in`, `out` and eight internal nodes, and their values.
The value-genome examples it builds on are [rc_lowpass_values](../rc_lowpass_values)
and [amplifier_gain](../amplifier_gain).

```
go run ./cmd/evolve -config examples/filter_topology/config.json   # ≈ 40 s
go test ./evo -run 'TestFilterTopologyExample|TestTopologyResume'   # the M8 acceptance tests
```

## The construction

**Embryo.** Three elements and the sweep: `VIN vs 0 AC 1`, `RS vs in 1k`,
`RL out 0 1k`, `.ac dec 20 10 100k`. A matched passband therefore sits at
−6.02 dB. The export sheet shows the engine's view of the bare fixture:
`out` at 0 V, and no automatic tearing possible because the circuit has no
path to tear.

**Genome** (`evo/encode_topology.go`). A variable-length list of parts,
each with a kind, two nodes and a log₁₀ value, grown into a copy of the
embryo as `L_g1`, `C_g2`, … The embryo is never modified. This is a direct
netlist encoding, not Koza's developmental program (doc/M8.md says why the
direct one was enough). At decode, self-loops and parts with an end on an
internal node nothing else touches are pruned, repeatedly: they carry no
current, so every node voltage is unchanged, and without the pruning,
debris from mutation would turn into floating-node penalties. What remains
invalid (capacitor-only islands, inductor loops) is penalized.

| Setting | Value |
|---|---|
| terminals | `in`, `out` (plus ground) |
| internal nodes | 8 |
| parts | 3 at start, at most 20 |
| L | 10 µH – 10 H, E24 |
| C | 100 pF – 100 µF, E24 |
| R | 1 Ω – 1 MΩ, E24 |
| GA | population 300, 200 generations, seed 1, 2 s timeout per evaluation |

**Operators.** Change value (Gaussian, in decades), add, delete, rewire one
end, change kind (keeping the relative position in the range), split a
part in series through a fresh node, add a part in parallel; one operator
per child plus a geometric number more (mean 1.43). Crossover is one-point
on the part lists. New parts attach to nodes already in use plus one fresh
internal node, so they mostly join the circuit.

**Objectives.**

| Mask | Requirement |
|---|---|
| passband 10 Hz–1 kHz | −7.5 dB ≤ gain ≤ −5.5 dB (matched: −6.02 dB) |
| stopband 3 kHz–100 kHz | gain ≤ −46 dB (40 dB below the passband) |

Fitness is the mask error (`evo.FilterMask`: 0 when the response stays
inside the mask at every sweep point and band edge, otherwise
0.1 + the fraction of points outside + the mean excess in dB/20, after
Koza's filter scoring) plus 0.002 per part (`evo.PartCount`), so once the
mask is met evolution keeps removing parts. The 0.1 floor means any design
that leaves the mask scores worse than any design that meets it, so
"fitness < 0.1" reads "spec met".

## Results

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
the stopband reaches at most −46.77 dB; that re-simulation of the written
`best.cir` is the acceptance test.

| Seed | Parts | Form |
|---|---|---|
| 1 | 4 | series L with C across it (notch 3.4 kHz), shunt C at both ends |
| 2 | 6 | same idea, with two small resistors |
| 3 | 4 | series L with C across it, shunt C both ends |
| 4 | 4 | series L, shunt C, series-LC trap to ground at the input (notch 3.9 kHz) |
| 5 | 6 | two series L, trap across, shunt C both ends, one small resistor |
| 6 | 4 | series L with C across it, shunt C both ends |
| 7 | 4 | series L, shunt C, series-LC trap to ground |
| 8 | 4 | series L, shunt C, series-LC trap to ground |

All of seeds 1–8 meet the mask; seven do it with four parts and one of the
two textbook forms of the third-order elliptic low-pass, the trap across
the series inductor or a series-resonant trap in a shunt arm. Third order
is the lowest that can meet about 1 dB ripple and 40 dB attenuation from
3× the passband edge (a Butterworth needs five reactive parts, a Chebyshev
four), reached here from an empty fixture with standard values.

`best.cir` opens in ngspice, LTspice and QSpice as it is.

## How the search behaves

- **Parsimony needs time.** With an early stop at "spec met", seed 1
  stopped at generation 50 with 15 parts; run to 200 generations, the
  0.002-per-part penalty pared it to 4.
- **Penalties, never a crash.** Every decoded circuit passes
  `circuit.Validate` before simulation; in this example 15–25 % of each
  generation fails validation or measurement and takes a penalty from the
  M5 classes (no-measurement < no-convergence < singular < invalid <
  ledger < timeout < other < panic). None panics.
- **Deterministic and resumable.** Operators draw from the per-generation
  PCG stream; a topology run checkpointed at generation 2 and resumed
  equals the uninterrupted run (`TestTopologyResume`).
- **Not done here** (doc/M8.md): Koza's developmental encoding,
  transistors and sources as evolvable parts, multi-part local search.

## Files

| File | What it holds |
|---|---|
| `embryo.cir` | the fixture: source, source resistance, load, sweep |
| `config.json` | the topology genome settings, mask and parsimony objectives, GA settings |
| `export/` | the bare embryo as `.cir`, `.kicad_sch` and a results sheet |
| `../../evo/encode_topology.go` | the genome, decode and pruning, operators |
| `../../evo/topology_test.go` | `TestFilterTopologyExample`, `TestTopologyResume` |
| `../../doc/M8.md` | milestone notes and the seed table |

## Sources

- J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, *Genetic
  Programming III: Darwinian Invention and Problem Solving*, Morgan
  Kaufmann (1999), ch. 25 — the embryo fixture, the filter mask scoring and
  the evolved low-pass filters this example repeats with a direct encoding.
