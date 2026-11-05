# Knotted light: does topology alone hold a knot together? (H29, step 1)

Hypothesis H29 (CLAUDE.md §1, user hypothesis 2026-11-05): particles are
knots of field energy, held together by their topology, since a knot
cannot untie without being cut. Kelvin proposed the same for atoms in 1867
(vortex knots in the ether; see [vortex_ring](../vortex_ring)). Step 1
asks the question that every ring model in this repository has run into
([photon_electron](../photon_electron), [fdamador_neutrino](../fdamador_neutrino)):
in ordinary, linear Maxwell electrodynamics, does a knotted field stay
together?

```
go run ./examples/knotted_light      # about 30 s
go test ./examples/knotted_light     # -short skips the 3-D solve
```

## The field

Rañada's knotted light, the hopfion: an exact, finite-energy solution of
Maxwell's equations in vacuum whose field lines are closed loops, every
pair linked once (the Hopf fibration). In the Bateman form of Kedia et al.
(2013), lengths in L0 and times in L0/c:

```
α = (r² − t² − 1 + 2iz) / (r² − (t − i)²)
β = 2(x − iy) / (r² − (t − i)²)
F = E + icB = ∇α × ∇β = ∇ × (α∇β)
```

It is set in `em/fdtd.go` as discrete curls of Re(α∇β) and Im(α∇β), so
div E = div B = 0 exactly on the grid. The helicity sign of the formula is
checked numerically against Maxwell's equations rather than assumed. The
solver also carries the vector potential A (temporal gauge: A ← A − Δt·E
after each step, since the solver advances E and then B ← B − Δt·curl E),
so curl A = B holds at every step and the magnetic helicity ∫A·B dV, the
integral measure of field-line linking, can be followed. For a null field
like the hopfion (E·B = 0 and |E| = c|B| everywhere), magnetic helicity is
conserved exactly.

## Results

10 cells per L0, region ±6 L0, absorber 12 cells; every row against the
exact solution at the same time:

| t (L0/c) | Energy within 2 L0, solver | Exact | RMS radius (L0) | Helicity, solver | Exact | Linking of two B lines, exact |
|---|---|---|---|---|---|---|
| 0 | 0.966 | 0.966 | 1.00 | 1.0000 | 1.0000 | 1.0000 |
| 0.5 | 0.955 | 0.956 | 1.11 | 0.9999 | 0.9999 | 1.0000 |
| 1 | 0.905 | 0.910 | 1.41 | 0.9992 | 0.9995 | 1.0000 |
| 1.5 | 0.730 | 0.751 | 1.80 | 0.9982 | 0.9988 | 1.0000 |
| 2 | 0.399 | 0.392 | 2.23 | 0.9971 | 0.9976 | 1.0000 |
| 3 | 0.019 | 0.020 | 3.16 | 0.9935 | 0.9933 | 1.0000 |
| 4 | 0.0012 | 0.0012 | 4.12 | 0.9873 | 0.9824 | 1.0000 |

- **The topology is kept.** The exact field's traced lines link once at
  every time, and the solver's helicity follows the exact helicity to
  within 0.005. The slow fall in both is the field's tail leaving the
  finite box, not the knot coming undone.
- **The energy is not.** By t = 3 only 2 % of the energy is within 2 L0;
  the knot grows as √(1 + t²) and in the end expands at the speed of light.
- **Books:** energy ledger residual 4.8×10⁻¹³ of the starting energy, div B
  at 5×10⁻¹⁵, curl A − B at 3×10⁻¹⁴ of max B. The solver's starting
  helicity is 0.992 of the exact value (grid discretization at 10 cells per
  L0).

Field lines traced on the solver's own field link once up to t = 1, then
stop closing: a closed field line is a fragile object, and the grid's
small errors turn near-closed loops into open ones. That is why the
helicity, an integral that does not depend on following one line, is the
solver's topology measure.

**Conclusion of step 1.** Linear Maxwell carries the knot's topology along
perfectly and does nothing to hold its energy together: the knot stays
linked while it flies apart. This is the standard prediction (Rañada 1989;
Irvine & Bouwmeester 2008). Derrick's theorem (1964) says why: no static,
finite lump of an ordinary field exists in three dimensions. A knot that
holds needs a nonlinear term, as in Skyrme's model of baryons (1961) and
Faddeev and Niemi's knot solitons (1997). Step 2 is to add such a term,
declared and labelled, and look for the size at which a knot stops
spreading.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestLinkingKnownCurves` | Gauss linking integral: Hopf link = ±1, separated circles = 0 | 1e-3 |
| `TestHopfionExact` | the field is null (F·F = 0); traced B lines close and link once at t = 0 and t = 3 | 1e-12; 1e-3 |
| `TestKnotInSolver` | coarse grid: energy spreads as the exact solution says; helicity kept as exact; curl A = B; ledger | 0.02; 0.02; 1e-12; 1e-10 |

This example has no circuit, so it has no `export_test.go`.

## Files

| File | What it holds |
|---|---|
| `main.go` | the hopfion (Bateman form, helicity sign check), the vector-potential and helicity bookkeeping beside `em/fdtd.go`, field-line tracing and the Gauss linking integral, the run and report |
| `main_test.go` | the tests above |
| `run_output.txt` | the run above |

## Sources

- A. F. Rañada, "A topological theory of the electromagnetic field," *Lett. Math. Phys.* 18, 97–106 (1989).
- W. T. M. Irvine, D. Bouwmeester, "Linked and knotted beams of light," *Nature Physics* 4, 716–720 (2008).
- H. Kedia, I. Bialynicki-Birula, D. Peralta-Salas, W. T. M. Irvine, "Tying knots in light fields," *PRL* 111, 150404 (2013) — the Bateman construction used here.
- H. K. Moffatt, "The degree of knottedness of tangled vortex lines," *J. Fluid Mech.* 35, 117–129 (1969) — helicity as linking.
- G. H. Derrick, *J. Math. Phys.* 5, 1252 (1964); T. H. R. Skyrme, *Proc. R. Soc. A* 260, 127 (1961); L. Faddeev, A. J. Niemi, "Stable knot-like structures in classical field theory," *Nature* 387, 58–61 (1997).
- W. Thomson (Kelvin), "On vortex atoms," *Phil. Mag.* 34, 15–24 (1867).
