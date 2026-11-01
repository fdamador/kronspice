# Milestone M9b — status

M9b from CLAUDE.md §9: `kron` diakoptics with manual tearing, with cached
torn solutions in `evo`. Done when the value-genome GA runs faster with
tearing than without, with the same results.

## Done-when checklist

| Test | Status |
|---|---|
| Torn solve equals the untorn solve (and MNA) to 1e−9 on the linear golden circuits | ✅ `TestTornTwin`: the 13 linear golden circuits and 60 random circuits, every single-branch cut, every second branch, and all branches at once: 2 289 torn solves, OP and AC; worst 5.1e−11 relative to MNA |
| Power is invariant under tearing (Roth 1959 Sec. 5) | ✅ same test: impressed power `−J·E` of the torn solution equals the untorn one; Tellegen's sum is zero |
| `m = b − N + S` after tearing | ✅ same test: the torn network's counts satisfy it, and meshes drop by `t − ΔS` for `t` torn branches |
| A cut that leaves a singular block is reported, not solved | ✅ `TestTornSingularPiece` (a floating piece named in the error, Aitchison Lemma 2), and 187 such cuts in `TestTornTwin`, each reported |
| Cached torn solutions in `evo` | ✅ `evo.Problem.Tear`; `TestTearSameResults`, `TestTearFallback`; a shared, concurrency-safe cache (`-race` clean) |
| Value-genome GA faster with tearing than without, same results | ✅ `examples/torn_ga` (below): same best design, torn Kron 11× faster than untorn Kron |

## What was built

**Algorithm (`kron/diakoptics.go`).** The named R, C and L branches are cut
out. In the node-pair frame the remaining matrix B is block-diagonal, one
block per piece, and the torn branches add back a low-rank term
`A = B + U·D·Uᵀ` (`D` their admittances, each column of `U` one torn
branch's incidence). Then

```
z = B⁻¹·c                       each piece on its own
(I + D·Uᵀ·B⁻¹·U)·w = D·Uᵀ·z     interconnection system, order t
x = z − B⁻¹·U·w
```

and `w` is exactly the vector of torn-branch currents — Kron's
intersection currents through the open paths that tearing created. This is
the matrix form of Aitchison's algorithm (J. Eng. Math. 1983, Sec. 2,
quoted in CLAUDE.md §3a): the interconnection matrix is non-singular
exactly when A is (his Lemma 1), and every block must be non-singular
(Lemma 2), which is checked and reported per piece with its node names.

**Accuracy.** In exact arithmetic the torn operator is exact; in floating
point it loses accuracy when a piece is *nearly* floating (a small
capacitor its only reference at a low frequency), because B is then
ill-conditioned although A is not. The torn solve therefore refines
against the whole A and checks the componentwise backward error
`max |A·x − c|ᵢ / (|A|·|x| + |c|)ᵢ`: it stops at 1e−14 (the level of a
backward-stable LU), and a solve that cannot reach 1e−12 in 12 passes is
reported (`ErrSingular`, "nearly floating") instead of returned. In the
random-circuit test this rejection caught one cut whose answer would
otherwise have been wrong by a factor of 3 000.

**Pieces and cache.** The pieces are the connected components of B's
pattern, numbered deterministically. Each piece is equilibrated and
factored on its own (`solver.CDense.Factor`), and its columns `B⁻¹·e` at the
torn-branch terminals are computed once. A `kron.TornCache` keeps factored
pieces by their exact content (every entry's bit pattern and the port
positions; a hit is confirmed entry by entry, so a hash collision cannot
return a wrong factor). It is bounded by count and by storage (128 MiB);
when full, the pieces reused since the last sweep are kept — the fixed
parts of a GA problem stay, one-off candidate pieces go. A cached and an
uncached solve are bit-identical (`TestTornCache`).

`solver.CDense.SolveFactored` now only reads the factors (the solution
vector is its working storage), so one factored piece is solved by any
number of goroutines at once.

**Node-form assembly.** The node-pair system is now assembled from the
relation terms in sparse rows (two passes: count, then fill one backing
array) instead of the dense `Aᵀ·Y·A` product; the dense product remains
the check (`TestNodeAssemblyEqualsTransform`). Each sum is added in the
same order as before, so the untorn results are unchanged.

**API.** `kron.Tear(c, Tearing{Torn: names}, cache)`, `Diakoptics.Solve`
(returns `TearStats`: pieces, sizes, interconnection order, cache hits,
refinements, backward error), `Diakoptics.Counts`, `kron.TornOP`,
`kron.TornAC`.

**GA (`evo`).** `Problem.Tear = &evo.Tear{Torn: names, Cache: …}` solves
every AC analysis of a value-genome candidate by `kron.TornAC`, with one
cache shared by all workers (Run makes one if none is given). Tearing never
changes which candidates are accepted: nonlinear circuits, circuits with
idealized active elements (whose model-supplied power must come from the
MNA ledger), and any sweep the torn path rejects are solved by MNA.

## Results (`examples/torn_ga`)

A two-pole RC pre-filter (four E-series genes) drives a fixed plant through
the torn resistor RT; fitness is the −3 dB point at the plant's output,
target 5 kHz, 201 frequencies per candidate. Population 40, 30 generations
(1 240 evaluations), seed 1, two CPUs:

| Plant | MNA | Kron node form, untorn | Kron diakoptics, torn at RT | Same best design |
|---|---|---|---|---|
| 40-section LC line (125 unknowns) | 6.4 s | 138.8 s | 12.5 s (11.1× faster than untorn) | yes |
| 12×12 RC grid (144 nodes) | 21.3 s | 177.0 s | 31.5 s (5.6× faster than untorn) | yes |

On the line the plant's pieces were taken from the cache 235 394 times and
the pre-filter's piece factored 223 690 times.

**Against MNA.** The milestone compares tearing with not tearing, and the
torn run wins by an order of magnitude. Against the MNA engine — sparse
Gilbert–Peierls LU with a reused pattern — the torn Kron solve is still
slower on the line: once the fixed piece is cached, what remains per
frequency is building the primitive network, the branch relations and the
node-pair rows, which costs about as much as MNA's whole sparse
factorization of a chain (2.0× slower than MNA). On the grid, where the
LU fills in, the gap narrows to 1.5×. The remaining per-frequency
overhead is listed below as the next step; until then MNA stays the
default and tearing is opt-in.

## Not done here

- Tearing for transient analysis and the nonlinear Newton Jacobian (M9c).
- Automatic cut selection (a graph partitioner); cuts are named by hand
  (`Tearing`, or `Problem.Tear`). A `.tear` netlist directive.
- Reusing the per-network structure (pieces, row patterns) across the
  frequencies of one sweep, and building only the changed piece of a
  mutated candidate instead of the whole network: the steps that would let
  the torn path overtake sparse MNA on chain-like plants.
- Tearing along the exterior open paths (§3a.1 item 7) to solve channel
  networks separately.
