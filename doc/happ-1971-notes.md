# Happ, *Diakoptics and Networks* (1971) — notes for kronspice

H. H. Happ, *Diakoptics and Networks*, Mathematics in Science and
Engineering vol. 69, Academic Press, New York, 1971, 331 pp. Dedicated to
Gabriel Kron. Copy: `reference/happ_1971_diakoptics_and_networks.pdf`
(scan with a text layer; equations and matrices are images, so read
numbers from the page images). PDF page = printed page + 17.

## Chapter map

| Happ | Content | kronspice |
|---|---|---|
| I | Direct, index and modified index notation; tensor concepts | `kron/transform.go` (typed frames) |
| II | Primitive reference frame, primitive network, transformation equations, topology of a dead network, live contours | `kron/primitive.go`, `kron/connection.go` |
| III | The orthogonal network: open- and closed-path contours, current and voltage transformations, equations of the orthogonal frame | `kron/orthogonal.go`; CLAUDE.md §3a.1 |
| IV | Additional cases; equivalent sources; validity of power invariance (§5) | ledger twin test, §3a.1 item 6 |
| V–VII | Restricted primitive network; constrained contours | — (not used yet) |
| VIII–IX | Transformations between contour sets; non-singularity of C; current and voltage parts | `kron/orthogonal.go` frame changes |
| X | Boundary cases: all-open-path and all-closed-path networks | node and mesh forms (`kron/node.go`, `kron/mesh.go`) |
| XI | Ignoring or eliminating open paths, or closed paths | `kron/happ.go` `ClosedPathsZ`, `OpenPathsY`; unmodeled ports as inactive open paths (§3a.1 item 4) |
| XII | Diakoptics: contours unconstrained and constrained, interconnecting by groups of branches, dual case | `kron/happ.go` `PiecewiseZ` (Table 12.1), `TwoLevel` (§4, Fig. 12.12); `kron/diakoptics.go` (Aitchison form) |
| App. I | Numerical illustration of Ch. III (node-to-datum case Fig. I.2, not-node-to-datum Fig. I.3) | `kron/happ_test.go` |
| App. II | (B + N − 1) current sources present: netting method | — |

## Appendix I, reproduced (`TestHappAppendixI`)

Six branches with mutual impedances (eq. 1.3, p. 276), voltage sources
e_b = (1, 2, 3, 1, 4, 2), impressed junction currents I° = (1, 1, 1), the
contours of Fig. I.2 (connection tensor C^b_s and inverse A^s_b, Table I.1,
p. 278). kronspice's `Transform` gives Happ's Z_ss (p. 279) exactly and
e_s = (−4, −5, 1; 3, 9, 4) (p. 277). The closed-path currents
(eq. 3.76) are exactly i^c = (1, 82, 15)/61; Happ's printed
(0.016390, 1.344262, 0.245904) come from his six-digit (Z_cc)⁻¹ and differ
by up to 6×10⁻⁶. Branch currents J^b (p. 280), i^b (p. 283), the zero
open-path parts of i^s, A·C = I, and the power invariant e·J between the
primitive and orthogonal frames all check.

## Table 12.1, the twelve-step piecewise algorithm (`PiecewiseZ`, `TestHappTable12_1`)

Happ's factorized solution of the orthogonal network (p. 237): steps 1–2
solve with every closed path open, 3–6 compute the closed-path currents
from Z_cc = C_cᵀ Z_bb C_c (factored, never inverted), 7–9 add their
contribution, 10–12 give the open-path voltages. On Appendix I it gives
i^c = (1, 82, 15)/61 exactly, V_b and E_o = (8.491808, 10.540966,
8.672140) as printed on p. 284 (within Happ's rounding, 5×10⁻⁵), and both
of Kron's invariants: e·J the same in primitive and orthogonal frames, and
E_o·I° equal to the power of the junction-pair sources.

## Appendix I §3, not node-to-datum (`TestHappAppendixI3`)

Fig. I.3 and Table I.2 (p. 291): same closed paths, open paths that do
not share a datum. The bus currents of Fig. I.1 become I° = (1, 1, 2)
(p. 290); the branch excitation C_o I° = (0, 0, −1, 1, −1, −1) is
unchanged, so J^b and i^c are unchanged. E_o becomes (8.491808, 1.868826,
8.672140) (eq. 3.33) and E_b matches eq. 3.41. The junction-pair power
E_o·I° = 27.704914 is the same in both open-path frames (Ch. VIII).

## Chapter XI, ignoring and eliminating paths (`TestHappChapterXI`)

- Open paths ignored (§2, p. 213): current sources become equivalent
  voltage sources −Z_bb I^b; only C_c is needed (`ClosedPathsZ`).
- Closed paths ignored (§3, p. 218): voltage sources become equivalent
  current sources −Y_bb e_b; only the junction-pair rows A_o are needed
  (`OpenPathsY`). This is the Happ form of the rule that inactive open
  paths may be skipped (§3a.1 item 4).
- Both give the Table 12.1 branch currents to 10⁻¹², and the nodal route
  gives Table 12.1's E_o in either open-path frame.
- Eliminating instead of ignoring: (Z_cc)⁻¹ = Y^cc − Y^co (Y^oo)⁻¹ Y^oc
  (eq. 11.10) and (Y^oo)⁻¹ = Z_oo − Z_oc (Z_cc)⁻¹ Z_co (eq. 11.28) hold to
  10⁻¹².

## Two-level diakoptics, single common bus (`TwoLevel`)

Ch. XII §4 (pp. 250–264), algorithm Fig. 12.12 (Table 12.5 rightmost
column, e_T = 0, I^L = 0). Each subdivision is hinged at the common bus;
its radial equivalent Z_k = Y_k⁻¹ is applied by factoring Y_k. Z_TT is
block-diagonal (eq. 12.104); Z_cc = Uᵀ Z_TT U + Z_LL has the order of the
number of links. `SetSubdivision` refactors only the changed block and
reuses the other blocks' columns of Z_TT·U, the caching that matters for
the GA.

**Happ's example reproduced (`TestHappSingleCommonBus`).** Fig. 12.9,
every branch 1 Ω except the 2A–2B link (0 Ω), 0.5 A into all eight buses
(p. 261). From the Fig. 12.13 Z_TT blocks the subdivision trees are
A: G–3A–1A–2A, B: 2B–{G, 1B, 3B}, C: G–2C–1C. Every step of Fig. 12.14
matches: V_T⁰ = (2.5, 3, 1.5; 2, 1.5, 2; 1.5, 1), e_c′ = (−0.5, 1.5, −1, 0),
i^c = (1/6, 2/3, −13/30, −7/30) exactly, the injected I^T′, V_T¹, and
V_T = (1.733333, 1.566666, 1.233333; 1.9, 1.566666, 1.633333; 1.466666,
1.2); so does Z_cc⁻¹ of Fig. 12.13. `TestHappTwinMNA` builds the same
network as a netlist and gets these voltages from the MNA operating point
and from `TornOP` (Aitchison form, links b1, b3, b4 cut) to 10⁻¹².

**Two forms, one answer.** Happ's interconnection unknowns are the link
(closed-path) currents with Z_cc of order "links"; `kron/diakoptics.go`
solves (I + D Uᵀ B⁻¹ U) w = D Uᵀ z, whose unknowns w are the same torn-branch
currents. The matrices are related by I + D Uᵀ B⁻¹ U = D·Z_cc with
Z_cc = Uᵀ B⁻¹ U + D⁻¹ (B⁻¹ = Z_TT, D⁻¹ = Z_LL) — so `diakoptics.go` can tear a branch of zero
admittance (dropped) while Happ's form handles a zero-impedance link
(Z_b2 = 0) directly. Cost is the same order.

## Still to do

- Multilevel diakoptics (Ch. XII, successive closing of interconnections,
  pp. 264 ff.) and the dual algorithm of Table 12.6.
- Constrained (tree/link) contours, Table 12.3, and the T/L algorithms of
  Table 12.5's other columns (link sources e_L and currents I^L present).
- Appendix II netting of (B + N − 1) current sources.
- Use `TwoLevel` from `evo` when a genome mutates one subdivision.
