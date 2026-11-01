# Kron's nomenclature in modern electrical engineering terms

Gabriel Kron wrote from the 1930s to the 1960s, in a vocabulary of his own drawn
from tensor analysis and differential geometry. Much of it describes things every
circuit engineer now knows under other names. This file gives the translation.
For the history behind the terms, see [HISTORICAL.md](HISTORICAL.md).

**Match** rates each translation:

- **Same:** the same thing under a different name.
- **Close:** the modern term is the nearest equivalent, but the two differ, as the
  note explains.
- **None:** no modern equivalent; the entry only describes the term.

The modern terms and the ratings are our own reading, not Kron's. The definitions
are checked against the sources listed [at the end](#sources). "1939" means
*Tensor Analysis of Networks* and its page numbers, and "1962" means "Camouflaging
Electrical 1-Networks as Graphs".

## Words that mislead

These look familiar but do not mean what a modern reader expects.

| Kron's word | What he means | What it does not mean |
|---|---|---|
| **Coil** | Any two-terminal element or branch of the network | A wound inductor |
| **Branch** | Parts of a network carrying the same current: one coil, or several in series (1939 p. 76). In his later work, and in Happ, it means an individual element, as it does today | One element, in the 1939 sense |
| **Mesh** | Any independent closed loop, in any network | A window of a planar circuit |
| **Tensor** | A matrix or array together with its rule for changing basis | A bare array, as in machine learning |
| **Orthogonal network** | A network described by both closed and open paths | A network with perpendicular parts |
| **Primitive** | The elements before they are connected | Simple or crude |
| **Apparent coil** | A zero-impedance or zero-admittance branch holding an ideal source | Apparent power |
| **Reference frame** | Any choice of independent variables, even with nothing rotating | Only a rotating d-q frame |
| **Graph** | A structure of nodes and branches, which Kron rejected as a model of an excited network (1962) | The standard tool it is today |

## Network structure

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Coil | Branch, element | Same | Coils = meshes + junction-pairs (1939 p. 75) |
| Junction | Node | Same | |
| Junction-pair | Node pair; a voltage across two nodes | Same | Count = junctions − sub-networks (1939 p. 75) |
| Datum | Reference node, ground | Same | |
| Node-to-datum case | Node voltages referred to ground | Same | Happ's name for a set of junction pairs that all share one node (1971 p. 23) |
| Sub-network | Connected component | Same | The `S` in `m = b − N + S` |
| Mesh, closed path, closed mesh | Independent loop, fundamental loop | Same | Modern "mesh" is the narrower planar case |
| Open mesh, open path | Path between the two nodes of a node pair; cut-set counterpart of a loop | Close | "Any circuit through the coils" from one junction of a pair to the other (1939 p. 76). The modern literature has no common word for the path itself |
| Contour | Loop or path chosen as a basis element | Close | Happ's umbrella term for closed and open paths |
| Dead network, dead contour | Topology alone, with no sources or values | Same | Happ 1971 pp. 25–26 |
| Live contour, excited network | Network with currents and voltages applied | Same | |
| Stretching | Node splitting with a zero-impedance branch | Same | 1939 pp. 79–80 |
| Physical subdivision | Partitioning a network into parts | Same | 1939 Ch. XIX |
| Compound coil, three-phase apparatus | Three-phase element treated as one block | Same | 1939 Ch. XIX |

## Currents, voltages and sources

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Impressed quantity | Independent source, excitation | Same | "Originate outside the network" (1939 p. 76) |
| Response quantity | Solution unknown, network response | Same | |
| `e` (impressed voltage) | Independent voltage source in series with a branch | Same | Summed around a loop to give the right side of loop analysis |
| `I` (impressed current) | Independent current source across a branch; node injection | Same | "Currents withdrawn by outside loads" also count (1939 p. 407) |
| `i` (mesh current) | Loop current | Same | |
| `E` (junction-pair voltage) | Node-pair voltage; node voltage in the datum case | Same | |
| Lowercase and capital letters | Quantities on closed paths (`e`, `i`) and open paths (`E`, `I`) | Close | Lowercase and capitals do not otherwise follow modern usage |
| Solenoidal `e`, `i` | Loop quantities, cycle-space component | Close | Borrowed from field theory for divergence-free flow (1962 p. 165) |
| Lamellar `E`, `I` | Cut-set quantities, cut-space component | Close | Borrowed from field theory for potential flow (1962 p. 165) |
| Active junction-pair | Port or driven node | Close | Classed as input or output (1939 pp. 566–567) |
| Inactive junction-pair | Internal node with no external injection | Same | No current `I` flows; can be eliminated (1939 p. 485) |
| Apparent coil | Ideal source treated as a branch of zero impedance or admittance | Close | Fills the role of the extra current unknown in modified nodal analysis, but Kron removes it from the analysis instead (1939 pp. 408, 412) |

## Matrices and transformations

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| n-matrix | Multi-dimensional array | Same | 1-matrix is a vector, 2-matrix a matrix |
| Geometric object | Quantity with a basis-change rule | Close | One symbol standing for its form in every frame (1939 p. 67) |
| Impedance tensor `z`, `Z` | Impedance matrix | Same | Primitive `z` is the branch impedance matrix, with mutuals off the diagonal |
| Admittance tensor `y`, `Y` | Admittance matrix | Same | `Y′` of a datum frame is the node admittance matrix, Y-bus |
| Connection tensor, transformation tensor `C` | Branch-to-loop matrix; in general, a change of basis | Close | See [the connection matrices](#the-connection-matrices) |
| `Cₜ` | `Cᵀ`, the transpose | Same | |
| `A = Cₜ⁻¹` | Inverse transpose of `C` | Same | Node-pair counterpart of `C` |
| Singular transformation | Rectangular, full-column-rank `C` | Close | Kron calls any non-square `C` singular; the "missing" axes complete it to a square one (1939 p. 429) |
| Equation of constraint | Kirchhoff current law written as `B·i = 0` | Close | `C` spans the null space of `B` (*Tensors for Circuits* pp. 46–50) |
| Transformation law `Z′ = Cₜ·Z·C` | Congruence transformation; loop-impedance matrix | Same | |
| Invariance of power | Power-invariant transformation; Tellegen's theorem | Close | Kron states it as a consequence of `i = C·i′`, and `e′ = Cₜ·e` then follows (1939 pp. 101–102). Tellegen's theorem (1952) is a related result from the Kirchhoff laws alone |
| Group | Set of transformations closed under composition | Same | Standard mathematics |
| Compound tensor, subdivision | Block matrix, partitioning | Same | 1939 Ch. IX |
| Doubly compound tensor | Nested block matrix | Same | 1939 Ch. X |
| Junction tensor, coil tensor | Connection matrix of a multi-phase element | Close | Changes an apparatus's own frame, as with a transformer connection (1939 p. 501) |
| Spinor transformation | Complex basis change with conjugate transpose | Close | Used for symmetrical components with complex power (1939 Ch. XIII). Not a physics spinor in the usual sense |

### The connection matrices

Kron's central object is a square, non-singular `C`, split by closed and open paths
(1962 pp. 164–165). Happ (1971 p. 25) ties the parts to linear graph theory.

| Kron | Modern term |
|---|---|
| `C = [C_o  C_c]`, `A = [A_o  A_c]` | Square connection matrices; no common modern name |
| `C_c` | Transpose of the tie-set (circuit) matrix |
| `A_o` | Cut-set matrix; augmented incidence matrix in the node-to-datum case |
| `C_cᵀ·A_o = 0` | Cut space orthogonal to cycle space |

kronspice's `Square` in [kron/orthogonal.go](kron/orthogonal.go) is `[C | P]`, where
`C` is Kron's `C_c` and `P` plays the part of `A_o`. It is not Kron's `[C_o  C_c]`.

## Analysis methods

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Mesh network method | Loop analysis (mesh analysis) | Same | `Z′·i′ = e′` |
| Junction network method | Node-pair or cut-set analysis; nodal analysis in the datum case | Same | `Y′·E′ = I′` (1939 Ch. XIV) |
| Orthogonal network method | Mixed or hybrid analysis, with loop and cut-set unknowns together | Close | Not modified nodal analysis, which mixes node voltages with selected branch currents. Mesh and nodal analysis follow by ignoring one set of paths (Happ 1971 Ch. XI) |
| All-open-path, all-closed-path network | Boundary cases of that formulation | Same | Happ 1971 Ch. X |
| Equation of voltage, of current | Loop equations, node equations | Same | |
| Dualism | Duality | Same | 1939 Ch. XIV |
| Network analyzer | Analog AC network computer | Same | Hardware of the era, replaced by digital load-flow programs |
| Equivalent circuit (of a field) | Discretized field model; lumped-element or transmission-line model | Close | The 1944 paper is cited as the origin of the transmission-line matrix method (Russer) |

## Reduction and tearing

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Reduction formula, elimination of variables | Schur complement; block Gaussian elimination | Same | 1939 Ch. X, p. 242 onward |
| (Named after him) Kron reduction | Schur complement of the node admittance matrix | Same | A modern name, not Kron's |
| Generalization of Thévenin's theorem | Multi-port Thévenin equivalent | Close | 1939 Ch. X |
| Star-mesh transformation | Star-mesh transformation | Same | Generalizes Y-Δ; a special case of the Schur complement |
| Tearing, diakoptics | Substructuring, domain decomposition; bordered block-diagonal solve | Close | Modern terms cover the numerics. Kron's method also has a physical reading, with the cut as an interconnecting network |
| Torn model | Primitive or partitioned network | Close | Happ calls the primitive a "torn model" of a network yet to be specified (1971 p. 15) |
| Interconnecting, interconnection system | Schur complement system on the cut variables | Close | `L·B⁻¹·K` in Aitchison's form, order = number of torn couplings |
| Piecewise solution | Divide-and-conquer direct solve | Close | The result is exact, not an iterative approximation |
| Multilevel diakoptics | Nested dissection, hierarchical substructuring | Close | Happ 1971 |

## Power-system terms

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Phase axes | abc frame | Same | |
| Sequence axes | Zero, positive and negative sequence frame | Same | 1939 Ch. XX |
| Sequence tensor | Fortescue transformation matrix | Same | |
| Sequence impedance tensor | Sequence impedance matrix | Same | Reduced to diagonal form in 1939 Ch. XIII |
| Integrated transmission system | Interconnected power system | Close | Inferred from the title of his 1951–53 series, which we have not read |

## Rotating machines

These follow the 1951 Epilogue as summarized in
[kron/machine.go](kron/machine.go).

| Kron | Modern term | Match | Notes |
|---|---|---|---|
| Primitive machine, generalized machine | Generalized machine model | Same | All machine types are obtained from one model by transformation (Hoffmann, *Tensors for Circuits* p. ix) |
| Metric tensor `g` | Inductance matrix with rotor inertia added | Close | Same matrix, read as the geometry of a coordinate space |
| Coordinates `q`, `θ` | Winding charges and rotor angle | Same | |
| Christoffel symbols `[βγ, α]` | Terms from `∂L/∂θ`: motional voltage and electromagnetic torque | Close | Also the source of the speed voltages in the winding rows |
| Speed voltage | Rotational or motional EMF | Same | |
| Speed-voltage tensor `G` | Coefficient matrix of `ω` in the dq voltage equations | Same | `G` from `∂(L·C)/∂θ` |
| Holonomic frame | Axes fixed to the iron; physical winding variables | Close | Time-varying `L(θ)` |
| Non-holonomic frame | Transformed frame, as with Park's dq | Close | `C` depends on `θ`, and `L` becomes constant |
| Two-reaction frame | Direct-quadrature (d-q) axes | Same | After Blondel and Park |
| Non-Riemannian dynamics | Machine equations in rotating coordinates | Close | The title's geometric reading has no place in modern practice |
| Operator `p` | `d/dt` | Same | Heaviside notation, `s` after Laplace transform (1939) |
| Hunting | Rotor angle oscillation, electromechanical mode | Same | Standard older term, still in use |
| Curvature tensor | Stability indicator | None | Kron read curvature as dynamic stability (1939 p. 196), and a 1955 paper's title ties it to damping and synchronizing torques. Modern practice uses eigenvalues |

## Notation

| Kron | Modern term | Match |
|---|---|---|
| Direct notation | Matrix notation | Same |
| Index notation, dummy index | Einstein summation | Same |
| Covariant index (lower) `e_a` | Dual vector, cochain | Same |
| Contravariant index (upper) `i^a` | Vector, chain | Same |
| First generalization postulate | Scalar laws keep their form when written in matrices | None |
| Second generalization postulate | The same equation holds in every basis | None |

The two postulates are heuristics, not theorems. The first says a set of equations
written as one matrix equation has the same form as for one degree of freedom
(1939 p. 62). The second says one symbol stands for its many forms across reference
frames (1939 p. 67). Hoffmann criticized his appeal to such postulates in place of proofs
(*Tensors for Circuits* p. viii).

## Polyhedra and multi-dimensional networks

Kron's late work has no accepted translation. The right-hand column gives the
nearest modern framework, and Kron himself rejected some of these identifications.
For the argument, see [HISTORICAL.md](HISTORICAL.md#polyhedra-and-n-networks-195868).

| Kron | Nearest modern idea | Match | Notes |
|---|---|---|---|
| p-simplex | p-simplex | Same | Node, branch, triangle, tetrahedron |
| p-network | Set of p-cells on their own | Close | A 1-network has branches but no nodes or planes, which standard theory does not treat separately (1962 pp. 161, 163) |
| 0-network | Set of nodes on their own | Close | Its closed paths are node pairs, and its single open path is the collection of all the nodes (1962 p. 169) |
| 1-network | Set of branches, without their nodes | Close | Kron's model of an excited circuit |
| 2-network | Set of planes (2-cells) on their own | Close | Its open paths are the planes covering the closed paths of the 1-network (1962 p. 170) |
| n-network | Set of n-cells on their own; the top-dimensional part of an n-dimensional complex | Close | The p-network with p = n, the last in the sequence 0-network, 1-network, … n-network (1962 p. 171). His obituary's "n-dimensional networks" refers to this late work |
| (n-1)-network | Set of (n-1)-cells; in Kron's 1965 title, a dual induced by a 1-network | Close | "Invisible Dual (n-1)-Networks Induced by Electric 1-Networks". Known to us only from the title and the abstract's record |
| Polyhedron | Simplicial or cell complex | Close | n + 1 connected p-networks, from the 0-network to the n-network (1962 p. 171) |
| Planar polyhedron | 2-complex | Same | |
| Dual polyhedron | Dual cell complex | Same | 1962 p. 167 |
| Incidence matrix `M` | Boundary operator | Same | Kron's relation to exterior forms and Stokes' theorem matches discrete exterior calculus, which came later |
| Polyhedral wave, magnetohydrodynamic wave | None | None | Kron's own usage. Do not read it as modern magnetohydrodynamics without checking |
| Crystal computer | None | None | A proposal; our sources document no working device |
| Self-organizing automata | None | None | Not the modern sense of self-organizing networks, as far as we can tell from 1962; the original papers have not been read |

## Related later terms

These come from later authors working in Kron's framework.

| Term | Modern term | Match | Notes |
|---|---|---|---|
| TAN, tensorial analysis of networks | Kron's method applied to electromagnetic compatibility | Same | Maurice and colleagues |
| Chord | Off-diagonal coupling term in `Z`; current-controlled voltage source | Close | Mutual inductance is the classic chord, and it may be any function, including delayed. See [kron/chords.go](kron/chords.go) |

## Key equations side by side

| Kron | Modern | Name |
|---|---|---|
| `e = z·i` | `v = Z·i` | Branch law, mesh case |
| `E + e = z·(i + I)` | `v_b + v_s = Z·(i_b + i_s)` | Generalized branch equation |
| `i = C·i′` | `i_b = Bᵀ·i_loop` | Branch currents from loop currents, with `B` the tie-set matrix |
| `e′ = Cₜ·e` | `e_loop = B·v_s` | Kirchhoff voltage law, loop sources |
| `Z′ = Cₜ·z·C` | `Z_loop = B·Z·Bᵀ` | Loop analysis |
| `Y′ = Aₜ·y·A`, or `Y′ = C⁻¹·y·Cₜ⁻¹` | `Y_node = 𝐀·Y·𝐀ᵀ` | Nodal analysis, with `𝐀` the reduced incidence matrix (rows = nodes) |
| Reduction formula, 1939 Ch. X | `Z_red = Z₁₁ − Z₁₂·Z₂₂⁻¹·Z₂₁` | Schur complement, eliminating block 2 (our block labels) |
| `m = b − N + S`, `p = N − S` | Same | Loop and cut-set counts |

The last column of each row uses the common modern convention. Kron's `C` has one
column per loop, so it is the transpose of the tie-set matrix `B`. Sign conventions
for `v_s` and `i_s` vary between textbooks.

## Sources

Definitions were checked against these, in `reference/` unless linked.

- G. Kron, *Tensor Analysis of Networks*, Wiley, 1939. Terms were read in the
  chapters cited, not the whole book. See also
  [doc/kron-3a-verification.md](doc/kron-3a-verification.md).
- G. Kron, "Camouflaging electrical 1-networks as graphs," *Quarterly of Applied
  Mathematics* 20(2), 161–174, 1962.
- G. Kron, *Tensors for Circuits*, Dover, 1959: Hoffmann's introduction, and pp.
  43–50 as checked in the verification file.
- H. H. Happ, *Diakoptics and Networks*, Academic Press, 1971: preface and Ch. II.
- P. W. Aitchison, "Diakoptics or tearing — a mathematical approach," *Quarterly of
  Applied Mathematics* 41, 265–272, 1983, as used in CLAUDE.md §3a.

The rotating-machine terms come from [kron/machine.go](kron/machine.go) and Kron's
*Equivalent Circuits of Electric Machinery* (1951) as summarized there; the 1951 scan
has no text layer, and we did not read it again for this file. Kron's *Diakoptics*
(1963) has not been read, so its own terms, such as "intersection network", are
left out.
