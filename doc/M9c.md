# Milestone M9c — status

M9c from CLAUDE.md §9: `kron` tearing of the Newton Jacobian for
nonlinear circuits. Done when the diode and BJT golden circuits match via
the torn solve.

## Done-when checklist

| Test | Status |
|---|---|
| Diode / BJT golden circuits match via the torn solve | ✅ `TestJacobianTearingGolden`: **all 36 golden circuits** — diodes (forward, clipper, zener, charge storage), BJTs (CE amplifier OP and AC, PNP, output curves, switching transient), MOSFETs levels 1–3 (OP, DC curves, VTCs, AC, inverter transients) and the linear ones — under 78 cuts each torn at a resistor, capacitor or nonlinear device, and at all of them at once. Every operating point, DC sweep point, AC frequency and transient Newton iteration is solved torn, and every run matches its ngspice golden values within the golden tolerances. 35 of the 36 circuits split into several pieces under at least one cut. |
| Torn equals untorn | ✅ same test: torn and untorn runs agree to 2.8e−12 of scale (78 993 torn systems, 0 untorn fallbacks) |
| Factorizations of unchanged pieces reused between Newton iterations (CLAUDE.md §3a rule 3) | ✅ `TestJacobianReuse`: in a rectifier driving a torn RC ladder, the ladder's piece is factored once per time-step size and reused at every other Newton iteration (959 reuses, 1 093 factorizations in 1 026 systems); waveform equal to MNA to 1.7e−13 V |
| Pieces whose devices "did not change operating region beyond tolerance" | ✅ `JacobianTearing.ReuseTol`: at 1 % the rectifier's piece also keeps its factors while the diode moves little (244 factorizations instead of 1 093), refinement makes up the difference, waveform equal to 1.5e−11 V |
| Singular systems stay singular; concurrency; allocation | ✅ `TestJacobianErrors` (a voltage-source loop is reported as by MNA), `TestJacobianConcurrent` (one tearing shared by goroutines, `-race`), `TestJacobianSolveAllocs` (a torn solve with reused pieces does not allocate) |

## What was built

**A torn linear-solver backend (`kron/jacobian.go`).** Newton–Raphson
linearizes the devices at each iterate and solves `J·x = b`, where `J` is
the MNA matrix the devices stamp. Tearing applies to `J` unchanged, so
`kron.JacobianTearing` is a solver backend: `analysis.Options.Solver`
(new, a `SolverFactory` of real and complex systems) makes every system of
an analysis torn. Junction limiting, gmin and source stepping,
pseudo-transient continuation and transient integration all run on torn
systems with no change to the analysis code.

**The cut.** Elements are named; the pairs of their non-ground terminals
are removed from `J`'s pattern, and the connected components of what
remains are the pieces (an unknown of a device — a branch current, an
internal node — joins the piece it couples to). Resistors, capacitors,
diodes, transistors, G and voltage-controlled B elements can be torn;
elements whose coupling runs through a branch-current unknown (V, L, E, F,
H, K) are refused with the reason. Tearing *at a nonlinear device* puts
its linearized conductances into the coupling entries: the device sits
between pieces, Kron's intersection network.

**The algorithm.** Every entry of `J` that joins two pieces is a coupling
entry. With `B` the block-diagonal part and `E = P·M·Q` the coupling
entries (`P` selects their rows, `Q` their columns):

```
z = B⁻¹·b                        each piece on its own
(I + M·Q·B⁻¹·P)·w = M·Q·z        interconnection, order = number of coupling rows
x = z − B⁻¹·P·w
```

This is Aitchison's algorithm for a matrix that is block-diagonal apart
from a few coupling entries (J. Eng. Math. 1983, Sec. 2; CLAUDE.md §3a).
The diagonal stamps of a torn resistor stay in the blocks of its two
pieces, so each piece sees it as a conductance to the datum and no piece
floats because of the cut (Aitchison's Lemma 2 holds wherever the untorn
circuit has a DC path to ground). The interconnection system is solved by
an equilibrated dense LU; so is each piece.

**Accuracy and fallback.** Every solve is refined against the whole `J`
and accepted at a componentwise backward error of 1e−12 or better (it
stops at 1e−14). A singular piece, or a solve that does not reach the
bound, is solved by an ordinary LU of the whole `J` and counted as a
fallback, so tearing never changes which circuits converge or what they
converge to. No golden run needed a fallback.

**Reuse.** Each piece remembers the block it was factored from. A block
that is bit-identical keeps its factors (the linear pieces of a nonlinear
circuit: factored once per analysis, and once per time-step size in a
transient). With `ReuseTol` > 0 a block that changed by less than that
fraction of its largest entry also keeps its factors — Kron's "devices
that did not change operating region beyond tolerance" — and refinement
against the exact `J` removes the difference; if it cannot, the piece is
refactored and the solve retried before any fallback.

## Speed

`BenchmarkJacobianTran` (transient, 300 steps, two CPUs):

| Circuit | MNA | torn | torn, ReuseTol 1 % |
|---|---|---|---|
| rectifier + 3-section RC ladder (8 unknowns) | 1.8 ms | 3.3 ms | 4.1 ms |
| rectifier + 12×12 RC mesh (150 unknowns) | 104 ms | 114 ms | 292 ms |

The torn solve gives the same answers but is not yet faster: on small
systems the bookkeeping outweighs an 8×8 LU, and in a transient the fixed
mesh's right-hand side changes every step, so its piece is solved (though
not refactored) at every iteration with a dense O(n²) substitution against
MNA's sparse one. Stale reuse trades factorizations for refinement steps
and loses here. The ways forward are listed below; until then
`Options.Solver` is opt-in.

## Reference added

F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits: Principles,
Theory and Applications*, Ellis Horwood / Wiley / SNTL, 1988, 413 pp.
(ISBN 0-85312-606-2) is now in CLAUDE.md §6 as the reference for the
nonlinear and parametric circuit work. It is not yet in `reference/`, so
nothing here cites its pages; once a copy is added, its treatment of
nonlinear and time-varying (parametric) elements will be read and cited
where it applies — the next parametric step (pumped capacitances and
inductances, the parametric amplifier) in particular.

## Not done here

- Sparse factors for large pieces, and a sparse substitution, so a large
  fixed piece costs less per iteration than MNA's whole solve.
- Tearing a branch-current coupling (inductors, sources) by naming the
  branch unknown.
- Kron's own frame for nonlinear elements: the node-pair (open-path) form
  takes voltage-controlled (N-type) elements directly and current-
  controlled (S-type) ones with their current as an axis (CLAUDE.md §3c);
  `kron.OP` and `kron.AC` still refuse nonlinear circuits, and the torn
  Jacobian works in the MNA frame.
- Using the torn Jacobian in `evo` (the AC tearing of M9b is there;
  nonlinear candidates still use MNA).
- Automatic cut selection.
