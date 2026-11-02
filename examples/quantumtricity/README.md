# quantumtricity — H22: the quark model as circuit electricity at the quantum level

Register entry H22 (CLAUDE.md §1) holds the user's hypothesis that the
Ohm's-law wheel (P, V, I, R) and the Standard Model wheel describe one
structure, the quark model being the principles of electricity at the
quantum level. As stated it makes no prediction, so the register holds it
in a testable form with the standard prediction beside it. This example
runs that test. The name is the user's and is used for this directory only.

```
go run ./examples/quantumtricity
go test ./examples/quantumtricity
```

## The one real link: both are gauge theories on a graph

A lattice gauge theory (Wilson 1974) has a group element on every link of a
graph, a gauge transformation at every node, and gauge-invariant quantities
on closed loops. Kron's network has the same three objects: the branch, the
junction and the mesh. `kron/gauge.go` writes the correspondence out:

| Lattice gauge theory | Kron's network (`kron/`) |
|---|---|
| node, link, plaquette | junction, branch, mesh |
| `Grad χ` (node → link) | `E = A·v`, the incidence tensor |
| `Curl a` (link → plaquette), holonomy | `e′ = Cᵀ·e`, the loop sum of the connection tensor |
| `Curl∘Grad = 0` | `Aᵀ·C = 0` (connection.go) |
| `Div = Gradᵀ` | `Aᵀ·i`, Kirchhoff's current law |
| `Curlᵀ` | `i = C·i′`, branch currents from mesh currents |
| Gauss generator `Div Q − ρ` | node charge balance; its time derivative is KCL |
| pure gauge `a = Grad χ`, zero holonomy | the branch voltage, `E = Grad(v_node)`: KVL |
| Hamiltonian `Σ Q²/2C + Σ Φ_p²/2L_p` | the ledger's `½Cv²` and `½Li²` (mesh current `i_p = Φ_p/L_p`) |
| Wilson loop `exp(iqΦ_p/ħ)` | mesh flux linkage in flux quanta |

For the group U(1) this is Maxwell's electrodynamics on Kron's 1944 field
network (Tonti 2001 is the modern statement), and in the lumped limit it is
circuit theory. For SU(3) it is the lattice form of the quark model (Kogut &
Susskind 1975). The question H22 asks is whether anything in the circuit
behaves as the SU(3) version does.

## Step 1 — abelian (standard): circuit theory is U(1) on Kron's network

On a two-mesh RLC circuit (`V1 R1 L1 L2 R2 C1`, three fundamental meshes)
solved in Kron's mesh frame at 1 kHz, read off the solution:

| Identity | Residual |
|---|---|
| `Curl∘Grad = 0` (`Aᵀ·C = 0`) | 0 |
| KVL as zero holonomy of the branch voltage | 2e-16 |
| the voltage is pure gauge, `E = Grad(v_node)` | 2e-16 |
| mesh EMF = `−jω` × flux holonomy (Faraday with the flux in the inductors) | 1e-16 |
| KCL as the rate of the Gauss generator, `Div I = 0` | 1e-17 |
| continuity, `Div I_C = jω·Div Q = jω·ρ_node` | 1e-18 |

Mesh 0 links 6.1e-5 Wb at the phasor peak, 2.9e10 flux quanta: a Cooper
pair carried round it would have its Wilson loop turned through 2π that many
times. The Hamiltonian's two terms are the ledger's energies element by
element (`C1: ½C|v|²`, `L1, L2: ½L|i|²`).

Then the same `U1Field` on Kron's 1944 field network, a 21×17×6 Yee grid
closed by conducting walls, is stepped beside `em.FDTD` on the TM110 mode for
4000 steps:

| Lattice vs FDTD | |
|---|---|
| E, edge by edge | 6.6e-13 of max E |
| energy: lattice Hamiltonian vs FDTD ledger | 8.6e-15 of W0 |
| exact discrete energy drift | 2.0e-14 |
| Gauss generator at interior nodes | 0 |
| ω against Yee's own dispersion / the continuum | 1.1e-7 / 5.2e-4 |

The two agree to round-off because Hamilton's equations on the lattice,
`Q̇ = Curlᵀ(Φ/L)` and `ȧ = −Q/C`, are Yee's leapfrog. The field network and
the circuit are one construction: a branch capacitor is an edge's `ε·(dual
face)/length`, a mesh inductance is a face's `μ·(face)/(dual length)`.

## Step 2 — non-abelian (the hypothesis): SU(3) on the same graph

`kron.SU3Field` puts an SU(3) element on every branch of the same two-mesh
graph. Gauge transformations `U_k → g_from U_k g_to†` leave the Wilson-loop
traces unchanged (6e-15) and turn the Gauss generator covariantly
(`G_n → g_n G_n g_n†`, 1e-14) but not invariantly: it changes by order 1,
because it is a matrix, not a number. Two meshes meeting at a node compose
by matrix multiplication:

| Field on the links | composition defect `‖W_pW_q − W_qW_p‖/3` |
|---|---|
| random SU(3), 200 fields | mean 0.75, min 0.37 |
| diagonal (Cartan) SU(3), the abelian U(1)×U(1) part | 0 |

A circuit quantity that "composes as SU(3) currents do" would show such a
defect. The search over every quantity the simulator carries on a branch:

| Circuit | Quantities | Largest composition defect |
|---|---|---|
| two-mesh RLC, Kron mesh frame | voltage, current, charge, flux, complex power | 3e-16 (the round-off of addition) |
| diode pair, operating point (D1 0.669 V 1.73 mA, D2 0.652 V 0.87 mA) | voltage, current; KCL 9e-15, KVL 0 | 0 |

None found. Every Kirchhoff sum commutes; a nonlinear constitutive law
(the diodes) changes the numbers, not the composition law. The difference
between the two wheels is exactly the difference between U(1) and SU(3): the
circuit is the abelian, classical limit.

## The quantum level of V, I and R

Since the 2019 SI the three units of the Ohm's-law wheel are realized from
the exact `e` and `h` (computed in the run from CODATA 2018):

| | | Source |
|---|---|---|
| volt | `K_J = 2e/h = 483 597.848 4 GHz/V` | Josephson 1962 |
| ohm | `R_K = h/e² = 25 812.807 459 Ω` | von Klitzing, Dorda & Pepper 1980 |
| ampere | `I = e·f`; 1 nA = 6.24e9 electrons/s | Pekola et al. 2013 |

One lepton and one gauge boson of the Standard Model wheel; no quark. The
only quark content a circuit sees is the lattice charge Z of its conductors
(copper: 29 × uud = 29 e). Conduction current is electrons. The fractional
charges e/3 and e/5 measured in the fractional quantum Hall effect
(de-Picciotto et al. 1997; Saminadayar et al. 1997) are collective
excitations of electrons, not quarks.

## Status: work in progress

Step 1 passes; step 2 is a null result so far, which records a bound: H22 is
*bounded*. The entry stays in the register and the example keeps growing.
What would move it: a circuit quantity with a non-zero composition defect,
i.e. one that carries an index transported along branches and whose node
sums do not commute.

Still open:

- **The matter sector.** The dictionary above matches the gauge fields
  only (photon ↔ gluon). A quark is a charged matter field at the nodes; its
  circuit counterpart is the electron, entering the Gauss law as the node
  charge ρ. Step 2 carried the quark side only as abstract colour charges
  ρ_n. A charge field on the nodes (electrons for U(1), colour triplets for
  SU(3)) with its hopping current along the links is the next construction.
- **Where the wheel's terms land.** V is the electric field on a link, I the
  matter current along it, P = V·I the work term J·E between matter and
  field, and R has no counterpart: a gauge theory is lossless, resistance
  is coupling to a heat bath. To be written into the dictionary with a test.
- **The TLM twin.** The register asks for agreement with `em/tlm.go` as well
  as `em/fdtd.go`; TLM is not yet written.

## Sources

- K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* 10, 2445, 1974.
- J. Kogut, L. Susskind, "Hamiltonian formulation of Wilson's lattice gauge
  theories," *Phys. Rev. D* 11, 395, 1975.
- G. Kron, "Equivalent circuit of the field equations of Maxwell — I,"
  *Proc. IRE* 32, 289, 1944; E. Tonti, *PIER* 32, 1, 2001.
- K. S. Yee, *IEEE Trans. Antennas Propag.* 14, 302, 1966.
- U. Vool, M. Devoret, "Introduction to quantum electromagnetic circuits,"
  *Int. J. Circuit Theory Appl.* 45, 897, 2017 (branch flux and charge).
- M. Gell-Mann, *Phys. Rev.* 125, 1067, 1962 (the λ matrices).
- D. Marcos, P. Rabl, E. Rico, P. Zoller, *Phys. Rev. Lett.* 111, 110504,
  2013 (circuits that emulate gauge fields).
- B. D. Josephson, *Phys. Lett.* 1, 251, 1962; K. v. Klitzing, G. Dorda,
  M. Pepper, *Phys. Rev. Lett.* 45, 494, 1980; J. P. Pekola et al., *Rev.
  Mod. Phys.* 85, 1421, 2013; BIPM, *SI Brochure*, 9th ed., 2019.
- R. de-Picciotto et al., *Nature* 389, 162, 1997; L. Saminadayar et al.,
  *Phys. Rev. Lett.* 79, 2526, 1997.
