# Kron's field network and his 1962 polyhedron

The standard-physics half of register entry H20 (CLAUDE.md §1): Kron's
equivalent circuit of Maxwell's equations, and the curve-fitting polyhedron
of his 1962 paper, solved by the circuit engine and checked against every
number he printed. The drawings and tables are transcribed in
[doc/kron-figure-catalogue.md](../../doc/kron-figure-catalogue.md).

```
go run ./examples/kron_polyhedron
go test ./examples/kron_polyhedron ./kron
```

## What was checked against what

| Source | Kronspice | Result |
|---|---|---|
| Kron 1945 Fig. 14, cavity potentials read on the network analyzer | The lattice solved exactly (`kron.Kron1945Cavity`) | Worst 5.9 % of the largest potential, mean 2 %; his own unbalance was 8 % |
| Kron 1945 Fig. 15, unbalanced currents | Recomputed from his Fig. 14 by his eq. (1) | 23 of 25 junctions within 0.33 A; the other two look like his slips (below) |
| Balasubramanian 1968 Fig. 15, waveguide matching error | The TE lattice (`kron.Waveguide`), with a closed-form match leaving ripple 1e-16 | +0.69, +1.50, +3.47, +6.43 % against about 1.0, 1.6, 3.4, 6.2 % read from his curve |
| TE₁₁ box resonance | Lattice dispersion, the spatial part of Yee's | Equal to 2e-14; continuum to 0.36 %, second order; Kron's solve equals MNA to 4e-10 |
| Kron 1962 Tables I–VI, given columns | `CurveFit.Characteristics` on his printed z | Worst 1.9e-6 |
| Tables I–VI, first planar frame | The integer frame below | Worst 8.9e-6 |
| Tables I–VI, oscillatory column | The closed-form frame below | Worst 2.8e-5 |
| Kron 1962 p. 67, power input is the error-square | One-mesh circuit, MNA ledger and Kron's mesh form | Dissipation equals the summed error-square; ledger residual 4e-13 W |

## The two findings

**A planar reference frame is a vector of whole numbers.** Kron never
states his frames. With 13 rows and 12 coefficients, a fit has one residual
direction `w`, and `ŷ = y − w(wᵀy)/(wᵀw)` whichever 12 columns span the model.
The errors in the first planar column have rank one across the six tables, and
that direction is

```
w = ( −3, 1, 1, 1 |  4, −1, −3, −4, −4, 2 over the branch lengths |  2, −10, −4 over 9 )
```

The genetic algorithm then recovers these nine integers from the printed
column alone, with every candidate solved as a circuit (17,702 evaluations).

**The oscillatory column returns its input.** It equals the given column for
all six functions, to rounding. All 13 characteristics are linear in the four
values at the points, so they span only a 4-dimensional subspace. A frame with
`w` orthogonal to that subspace fits every function exactly, with zero power,
and `kron.CurveFit.OscillatoryFrame` writes one down:

```
N = (9, M₂/2, M₁/2, −13, −10, M₃/2), any even M
```

That is why the oscillatory polyhedron "can estimate satisfactorily even a set
of random numbers", of which Kron wrote "the significance of this is not yet
understood" (p. 69): 1000 random sets here fit to 7e-15. It is also why six
"widely different" functions fit without adjustment. A model with 12
coefficients is fitting 4 independent numbers, so this agreement shows no
estimating power. It is a property of the frame, and it does not by itself
support the physical reading of H20.

This is an inference from the tables, with three limits. Any member of the
family matches the printed column, so which frame Kron used cannot be
determined. He wrote that "no exact fit is possible", and a residual below
about 1e-5 cannot be excluded from the printed figures. Kron 1959 and the 1960
General Electric report have not been read and might state otherwise.

## Where the GA did not succeed

On Kron's objective, minimum power input, the search over integer frames did
not reach zero. In a sweep of 200 × 300 runs, none of 17 completed runs (three
mutation steps, up to six seeds each) got there. The best was 3.7e-4 W, against
90 W for frame 1. The sweep was cut off by its time limit before the 0.2-decade
step ran. Zero power is a linear condition, so it is solved by the closed form,
and the GA adds nothing there. `TestGAMinimumPower` records the shortfall.

## What is not reproduced

- **The second planar frame.** Its function rows follow frame 1's residual
  times a constant 0.3032, but its divided differences are not one linear map
  of z. That means several misprints or a different frame per table; six tables
  cannot separate the two.
- **Least-square + linear.** With two corrections it is linear to 2e-5. Two
  outer branches come out with the sign reversed and the third exact, but the
  inner rows have a rank-3 error map, not a single residual direction, so no
  frame of this form gives it.
- **The full 3-D cell with its ideal transformers** (Kron 1962 Fig. 1). Its
  values are in `kron.FieldCell`, but only the 2-D TE reduction is built and
  tested. Stage B of the catalogue's plan is still open.
- **q-networks in general.** `kron.CurveFit` handles this one polyhedron. The
  second-difference formula `S = −(Δz/L)(A/L)` is read off the tables, not
  derived from Kron's construction.

## Apparent errors in the printed sources

These are inferences from internal consistency, not published errata.

- Table II, least-square + linear, second function row: −12.339281 printed,
  −12.330283 needed.
- Table I, same column, fourth divided difference: −6.4341227 printed,
  −0.64341227 needed (a slipped decimal point).
- 1945 Fig. 15, row 2 column 3: 2.15 printed, −2.16 computed, and the figure's
  own e_u is −0.51 (a dropped minus sign).
- 1945 Fig. 15, row 1 column 0: 2.12 printed, −5.94 from the printed potentials.
  The cause is unresolved.

Kron's p. 62 claim that a polyhedron "improves always" the least-square
estimate fails in his own Table I: z₁ is linear, the least-square fit is exact,
and frame 1 moves it from 9.000 to 8.879.

## Sources

- G. Kron, *J. Math. Anal. Appl.* 5, 46–69, 1962; *J. Appl. Phys.* 16, 172–186,
  1945 (`reference/`).
- N. V. Balasubramanian, PhD thesis, Liverpool, 1968 (`reference/642931.pdf`).
- H. Nicholson, ACSE Research Report 26, Sheffield, 1974, on the least-squares
  reading of the polyhedron (`reference/`).
- E. Tonti, *PIER* 32, 1–44, 2001 (`reference/`).
