# M9e — Kron's network in geometric algebra (`kron/gaframe.go`)

CLAUDE.md §3b item 1. A network of separate resistive coils has a
diagonal primitive impedance, which is a metric on the branch space:
e_k·e_k = z_k (or y_k). Kron's mesh axes are the columns of the
connection tensor C and his node-pair axes the columns of A; as vectors of
the branch space, their outer product is a blade whose square is the Gram
determinant,

    ⟨B·B̃⟩ = det(Cᵀ·Z·C) = det Z′      (and det(Aᵀ·Y·A) = det Y′).

Building the blade one axis at a time, |B_k| = |B_{k−1}|·|rej(c_k)|, so
the rejection of each axis from the ones before tells how independent it
is, and a zero names the axis that makes the network singular.

## What exists

- `Network.GAFrames()` — the mesh and node frames at DC for networks of
  R, L, C and independent sources: the mesh frame on the branches that
  pass DC (metric R; 0 for inductors and voltage sources; capacitors and
  current sources are junction pairs), the node frame on the resistors
  with inductors and voltage sources merging their nodes (metric 1/R).
- `GAFrame.Blade()` — the explicit blade in the branch algebra
  (`hyper/ga`), the metric carried by scaling each axis by √|w| (positive,
  negative and zero weights become the p, q and r parts of the
  signature); up to 24 branches.
- `GAFrame.Independence()` — the rejection ratios and their product
  (twice-iterated Gram–Schmidt in the metric), any size.
- `GramDet` — the conventional determinant the blade is checked against.
- `kron.Precheck(c, opt)` — rejects before any matrix is built a circuit
  whose mesh frame has a zero-impedance axis: a loop of ideal voltage
  outputs (V, E, H, voltage-output B), and with `DC` also loops closed by
  inductors and machine windings. The error wraps `kron.ErrSingular` and
  names the loop. It works on nonlinear circuits too (it is structural).
  With `Conditioning` it also reports, for a linear circuit, the least
  independent mesh and node axes.
- `evo.Problem.Precheck` — runs it on every decoded candidate; a
  rejection is penalized as "singular" (`evo.Penalty` now maps
  `kron.ErrSingular`).

Floating nodes (a node-pair axis with no admittance) were already caught
by `circuit.Validate`, so the precheck does not repeat them.

## Checks (`go test ./kron -run 'GA|Precheck'`, `./evo -run Precheck`)

| Check | Result |
|---|---|
| Blade² = det(Gram) = product of rejections, mesh and node frames | 300 random passive networks (L, C, sources), 1e−9; 600 explicit blades |
| Duality of the frames for resistor networks: det(Cᵀ·R·C) = det(Aᵀ·G·A)·Π R (matrix-tree theorem) | 100 random networks, 1e−9 |
| Spanning trees of K₄ with unit resistors: node blade² = 16 (Cayley) | exact to 1e−12 |
| Linear golden circuits: blade² = det; every golden circuit passes the precheck | 9 linear, 40 total |
| Precheck names V1 ∥ V2, a loop with an E output, and V–L (DC only) | as stated; MNA also fails on each |
| No false positives: random circuits with extra ideal sources | every precheck rejection is also singular for MNA |
| Ill-conditioning (1 µΩ next to 1 GΩ) is reported (least independent axis 4.5e−5), not rejected | as stated |
| GA: precheck leaves a sound problem's best fitness unchanged; a V1 ∥ V2 template is penalized "singular" with the loop named | as stated |

## Cost (12 × 12 resistor grid, 265 branches)

| | time |
|---|---|
| structural precheck | 5.8 µs |
| MNA operating point | 1.7 ms |
| rejection pass, node frame (143 axes) | 11 ms |
| Gram matrix + dense LU of the same frame | 6.5 ms |
| precheck with conditioning (both frames) | 26 ms |

So the structural check is about 300× cheaper than the solve it saves
and is what the GA uses. The rejection pass is a diagnostic: about as
costly as a dense factorization of Z′ (1.7× here), and 15× an MNA
operating point, so it is off by default. The explicit blade is for
checking and teaching; its size grows as C(b, m).

## Not done here

- AC frames: the impedance metric is complex there, which a real Clifford
  algebra does not carry; the AC singularity cases (series LC at
  resonance with no resistance) are left to the solver.
- Controlled sources, transducers and mutual inductance are not a
  diagonal metric; `GAFrames` returns `ErrUnsupported` for them (the
  structural precheck still covers their voltage outputs).
