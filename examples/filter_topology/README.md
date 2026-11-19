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
go run ./cmd/evolve -config examples/filter_topology/config.json                 # direct, ≈ 40 s
go run ./cmd/evolve -config examples/filter_topology/config_developmental.json   # Koza's developmental, ≈ 12 s
go test ./evo -run 'TestFilterTopologyExample|TestTopologyResume'               # the M8 acceptance tests
go test ./evo -run 'TestFilterDevelopmentalExample|TestDevelop'                 # the developmental ones
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
netlist encoding; Koza's developmental program is the second run below
("Koza's developmental encoding"). At decode, self-loops and parts with an end on an
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
- **Not done here** (doc/M8.md): transistors and sources as evolvable
  parts, multi-part local search.

## Koza's developmental encoding

`config_developmental.json` evolves the same filter against the same mask
with Koza's encoding (`evo/encode_develop.go`; doc/M8.md). The genome is a
program, not a netlist: the embryo gets one modifiable wire from `in` to
`out`, and a program tree rewrites it through a writing head. R, L and C
turn the component at the head into a part (an arithmetic subtree of +, −
and constants sets its value, in decades from the middle of the kind's
range); SERIES, PARALLEL and VIA (a tap to ground or a terminal) divide it
into new components, each with its own head; END stops a head and
SAFE_CUT removes its component. Heads run breadth-first. A wire still in
place at the end joins its nodes (a 0 V source `V_w1` if it joins two
embryo nodes). Crossover swaps subtrees of the same type; mutation
replaces a subtree or moves a constant. FLIP reverses a component; it
matters only for parts with per-terminal ECE terms (`kappa_a ≠ kappa_b`
in a part kind, hypothetical, H48) and is drawn only then, so this
Maxwell run never uses it (doc/M8.md).

The evolved program is printed with the result; seed 1's starts
`(PARALLEL (PARALLEL SAFE_CUT (L 0.804 (SERIES END (VIA0 …`. Most of its
ops build nothing that survives (introns, as Koza reports); the circuit
they grow, after pruning:

```
L_g1 in   g_n1 62m
C_g5 g_n1 in   39n      across L_g1: notch at 3.24 kHz
L_g3 g_n1 out  120m
C_g4 g_n1 0    160n
C_g2 0    out  270n
```

| Seed | Parts | Form |
|---|---|---|
| 1 | 5 | L–C–L–C ladder with a C across the first L (trap) |
| 2 | 5 | C–L–C–L ladder, the first series L split in two |
| 3 | 4 | L–C–L–C ladder (fourth-order all-pole) |
| 4 | 4 | C–L–C–L ladder |
| 5 | 4 | shunt C, series L, series-LC trap to ground at the output (third-order elliptic) |
| 6 | 6 | L–C–L–C ladder, both series L split in two |
| 7 | 4 | C–L–C–L ladder |
| 8 | 4 | series-LC trap to ground at the input, series L, shunt C (third-order elliptic) |

All eight meet the mask. Where the direct encoding nearly always finds
the four-part third-order elliptic, this one more often finds a four-part
fourth-order all-pole ladder, the other minimal answer: SERIES and VIA
build ladders naturally, while a trap needs a PARALLEL whose copy becomes
the other kind.

## Files

| File | What it holds |
|---|---|
| `embryo.cir` | the fixture: source, source resistance, load, sweep |
| `config.json` | the topology genome settings, mask and parsimony objectives, GA settings |
| `config_developmental.json` | the same experiment with Koza's developmental encoding |
| `export/` | the bare embryo as `.cir`, `.kicad_sch` and a results sheet |
| `../../evo/encode_topology.go` | the genome, decode and pruning, operators |
| `../../evo/encode_develop.go` | the developmental encoding: program trees, development, operators |
| `../../evo/topology_test.go` | `TestFilterTopologyExample`, `TestFilterDevelopmentalExample`, `TestTopologyResume` |
| `../../evo/develop_test.go` | known programs → known circuits, operator invariants, repair, resume |
| `../../doc/M8.md` | milestone notes and the seed table |

## Sources

- J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, *Genetic
  Programming III: Darwinian Invention and Problem Solving*, Morgan
  Kaufmann (1999), ch. 25 — the embryo fixture, the filter mask scoring,
  the circuit-constructing functions, and the evolved low-pass filters this
  example repeats with a direct and a developmental encoding.
