# Torn GA: diakoptics in a running evolution

A two-pole RC pre-filter (R1, C1, R2, C2 in E24/E12) is evolved in front of
a fixed plant, joined to it by the resistor RT, for a 5 kHz cutoff at the
plant's output. The same run (same seed) is made three times: AC by MNA, by
Kron's node-pair form untorn, and by Kron's diakoptics torn at RT with the
plant's factored pieces cached and shared by all workers. This is the
acceptance example of milestone M9b (CLAUDE.md §3a, "Diakoptics" and "Why
this matters for the genetic algorithm"; [doc/M9b.md](../../doc/M9b.md)):
a mutation changes only the pre-filter, so only its piece is re-factored
and the plant's factors are reused from one candidate to the next.

```
go run ./examples/torn_ga                      # 40-section LC line
go run ./examples/torn_ga -plant grid -n 12    # 12×12 RC mesh
go run ./examples/torn_ga -gens 30 -pop 40 -seed 1
go test ./examples/torn_ga                     # same best design, torn faster than untorn
go test ./examples/torn_ga -bench Generation   # one generation, fresh seeds
```

The program prints each run's time and best design and exits non-zero if
the three best designs differ.

## The construction

**Circuit** (`Netlist`): `V1 in 0 AC 1`, the pre-filter `R1 in m 1k`,
`C1 m 0 10n`, `R2 m f 1k`, `C2 f 0 10n`, the coupling `RT f p0 10`, and one
of two plants from `p0` to `out`:

- **line:** n sections of 1 Ω, 1 µH, 1 nF (√(L/C) ≈ 31.6 Ω), matched by
  `RL out 0 31.6`; 40 sections by default, 125 unknowns;
- **grid:** an n×n RC mesh, 10 Ω between neighbours and 1 nF to ground at
  every node, driven at one corner and loaded by 100 Ω at the other; the
  kind of plant (a resistive plane, a substrate) whose LU factorization
  fills in; 12×12 by default, 144 nodes.

**Problem** (`Problem`): genes R1, R2 (10 Ω – 100 kΩ, E24) and C1, C2
(100 pF – 10 µF, E12); one objective, `evo.Cutoff` at `out`, target 5 kHz,
on `.ac dec 50 100 1meg` (201 frequencies per candidate).

**Three solvers.** The MNA run is the ordinary `analysis.AC` path. The
untorn Kron run sets `evo.Problem.Tear` with no torn branches, so every
candidate goes through the node-pair form `Y' = Aᵀ·Y·A` (M9a). The torn run
names RT: `kron.TornAC` cuts it out, so the node-pair matrix is
block-diagonal (the pre-filter in one block, the plant in the other) plus
RT's rank-one term, and solves

```
z = B⁻¹·c                       each piece on its own
(I + D·Uᵀ·B⁻¹·U)·w = D·Uᵀ·z     interconnection system, order 1 here
x = z − B⁻¹·U·w
```

with `w` the current through RT, Kron's intersection current through the
open path tearing created. The torn solve refines against the full matrix
and checks its componentwise backward error, reporting a nearly floating
piece instead of returning a wrong answer. A `kron.TornCache` keeps
factored pieces by their exact content (every entry's bit pattern, confirmed
entry by entry on a hit), bounded by count and by 128 MiB of storage; a
cached and an uncached solve are bit-identical.

## Results

Population 40, 30 generations (1,240 evaluations), seed 1. The milestone's
run, on two CPUs:

| Plant | MNA | Kron node form, untorn | Kron diakoptics, torn at RT | Same best design |
|---|---|---|---|---|
| 40-section LC line (125 unknowns) | 6.4 s | 138.8 s | 12.5 s (11.1× faster than untorn) | yes |
| 12×12 RC grid (144 nodes) | 21.3 s | 177.0 s | 31.5 s (5.6× faster than untorn) | yes |

The same line run on a 16-CPU machine (2026-11-01), best design in all
three R1 13 Ω, C1 47 nF, R2 3 kΩ, C2 390 nF, fitness 5.6e-6:

| Solver | Time | Pieces cached / factored |
|---|---|---|
| MNA | 1.18 s | |
| Kron node form, untorn | 21.80 s | 9 / 229,533 |
| Kron diakoptics, torn at RT | 2.67 s | 233,046 / 226,038 |

Torn is 8.2× faster than untorn Kron here. On the line the plant's pieces
are taken from the cache about as often as the pre-filter's piece is
factored: every candidate factors its own 4-unknown piece and reuses the
121-unknown line.

**Against MNA.** The milestone compares tearing with not tearing, and the
torn run wins by an order of magnitude. Against the MNA engine (sparse
Gilbert–Peierls LU with a reused pattern, M6) the torn Kron solve is still
slower, 2.0× on the line and 1.5× on the grid in the two-CPU run and 2.3×
on the line above: once the fixed piece is cached, what remains per
frequency is building the primitive network, the branch relations and the
node-pair rows, which costs about as much as MNA's whole sparse
factorization of a chain. Reusing the per-network structure across the
frequencies of one sweep, and building only the changed piece of a mutated
candidate, are the steps listed in doc/M9b.md that would let the torn path
overtake sparse MNA; until then MNA is the default and tearing is opt-in.

The export sheet (`export/torn_ga_line40_results.md`) shows the problem
circuit torn at RT into pieces of 4 and 121 unknowns joined by one link,
the torn and MNA node voltages agreeing to 4.5e-15, and the power taken by
the R, L, C branches identical in both, Roth's invariance.

**Tests.** `TestSameBestFasterTorn` runs a 20-section line (population 16,
6 generations, seed 3) and requires the same best values from all three
solvers and a torn time below the untorn one. `BenchmarkGeneration` times
one generation of 8 candidates with a fresh seed each iteration (new
candidates, as in a running GA: the plant's pieces stay cached, the
pre-filter's do not) for line/40 and grid/12, by MNA and torn.

## Files

| File | What it holds |
|---|---|
| `main.go` | `Netlist` (line and grid plants), `Problem`, `Run` (the three solvers), the report |
| `main_test.go` | `TestSameBestFasterTorn`, `BenchmarkGeneration` |
| `export_test.go` | the 40-section line problem torn at RT, exported with `KRONSPICE_EXPORT=1` |
| `export/` | `torn_ga_line40.cir`, `.kicad_sch`, `_results.md/.csv` |
| `../../kron/diakoptics.go` | the torn solve and `TornCache` |
| `../../doc/M9b.md` | milestone notes: algorithm, accuracy, cache, the two-CPU timings |

## Sources

- G. Kron, *Diakoptics: The Piecewise Solution of Large-Scale Systems*,
  Macdonald (1963).
- G. Kron, *Tensor Analysis of Networks*, Wiley (1939), p. 133 — a system
  divided by opening branches, the parts analysed separately and recombined
  by a connection tensor.
- H. H. Happ, *Diakoptics and Networks*, Academic Press (1971), Ch. XII.
- P. W. Aitchison, "Diakoptics or tearing: a mathematical approach,"
  *Quarterly of Applied Mathematics* 41, 265–272 (1983) — the algorithm in
  matrix form; Lemma 1 (the interconnection matrix is non-singular when the
  whole is) and Lemma 2 (every torn piece needs its own reference).
- J. P. Roth, "An application of algebraic topology: Kron's method of
  tearing," *Quarterly of Applied Mathematics* 17 (1959), Sec. 5 — power is
  invariant under tearing.
- J. R. Gilbert, T. Peierls, "Sparse partial pivoting in time proportional
  to arithmetic operations," *SIAM J. Sci. Stat. Comput.* 9, 862 (1988) —
  the MNA sparse LU the torn path is measured against.
