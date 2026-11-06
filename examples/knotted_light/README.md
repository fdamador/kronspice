# Knotted light: does topology alone hold a knot together? (H29 step 1, H30)

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
go run ./examples/knotted_light -viz examples/knotted_light/viz   # also the 3-D field lines
python3 examples/knotted_light/viz/render.py                      # knot_3d.png, knot_time.png
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

## H30: implosion and explosion

Hypothesis H30 (CLAUDE.md §1, user hypothesis 2026-11-05, after Walter
Russell's "rhythmic balanced interchange" of compression and expansion):
a particle is a cycle of implosion and explosion, and its energy books
balance only when the inward half is counted as well as the outward one.

The exact knot already has both halves. Maxwell's equations are
time-symmetric, and the hopfion before t = 0 is the expanding knot run
backwards: it comes in from far away, is tightest at t = 0, and leaves
(energy within 2 L0, exact: 0.021 at t = −3, 0.963 at 0, 0.021 at +3).
Step 1 above watched only the second half.

**How the incoming half is fed in.** A closed Huygens surface, a cube of
half-width 5.75 L0, surrounds the knot (total-field/scattered-field;
Taflove & Hagness, *Computational Electrodynamics*, 3rd ed., ch. 5).
Inside it the solver carries the whole field; outside, only what departs
from the declared incoming wave, which is zero up to discretization
error. The surface is a sheet of equivalent sources: 160,080 edge currents
and 158,700 magnetic face currents (`em.FaceSource`, added for this, with
its own exact ledger term and test). So the solver's ledger books every
joule the surface puts in and takes out. The declared incoming wave
stands in for whatever an ether, an absorber or any other medium would
have to supply. The run reads its power off the ledger; where it would
come from is the open question.

```
go run ./examples/knotted_light               # step 1 and H30, about 3 min
go run ./examples/knotted_light -h29only      # step 1 only
python3 examples/knotted_light/viz/render_h30.py   # viz/h30_cycle.png from run_output.txt
```

**Part 1: one implosion, then one explosion** (t = −10 … 10 L0/c; energies
as fractions of one knot's energy):

| t (L0/c) | Within 2 L0, solver | Exact | Inside the surface | Put in by the surface | Taken out |
|---|---|---|---|---|---|
| −10 | 0.0000 | 0.0000 | 0.004 | 0 | 0 |
| −4 | 0.0012 | 0.0012 | 0.988 | 0.984 | 0 |
| −2 | 0.377 | 0.393 | 0.990 | 0.986 | 0 |
| 0 | 0.957 | 0.966 | 0.990 | 0.986 | 0 |
| 2 | 0.413 | 0.389 | 0.990 | 0.986 | 0.0001 |
| 6 | 0.0000 | 0.0000 | 0.839 | 0.986 | 0.150 |
| 10 | 0.0000 | 0.0000 | 0.004 | 0.986 | 0.980 |

- **Every joule is accounted for.** The surface put in 0.986 of a knot and
  took out 0.980; the other 0.0066 leaked past it into the absorber (the
  grid's dispersion makes the solver's wave differ slightly from the
  declared one). Ledger residual 1.8×10⁻¹³ of the knot's energy.
- **The knot survives the implosion.** At t = 0 two magnetic field lines
  traced on the solver's own field come back to within 0.034 L0 of their
  starting points, and their linking number is 1.0000. (Lines traced on
  an inexact field do not close exactly. Each is closed at its nearest
  return, a chord much shorter than the loop, which keeps the Gauss
  integral at a whole number.) The solver's field is 4.3 % from the exact
  one within 2 L0, after the wave has crossed 10 L0 of grid.

**Part 2: a pulsing knot.** Knots tightest at t = 0, 8 and 16 L0/c are
fed through the same surface, a train of implosions and explosions with
period 8 L0/c. The first starts inside the surface at t = −4; the next two
come in through it.

| Peak at t (L0/c) | Within 2 L0, solver | Exact | Linking of two solver lines (gap) |
|---|---|---|---|
| 0 | 0.9571 | 0.9668 | 1.0000 (0.020 L0) |
| 8 | 0.9579 | 0.9674 | 1.0000 (0.035 L0) |
| 16 | 0.9573 | 0.9668 | 1.0000 (0.037 L0) |

- The knot is rebuilt at the centre every period, to the same energy and
  with the same linking.
- Over the middle cycle (t = 4 … 12) the surface's net work is
  +0.0046 knot energies, which equals the energy that leaked past it: in
  equals out. Ledger residual 1.5×10⁻¹².
- The surface's "put in" and "taken out" columns count each step's net
  flow. When an arriving and a departing knot overlap, they show only the
  difference, so the gross flow each way (one knot per period) cannot be
  read from them.

**Conclusion.** Under linear Maxwell, an imploding and exploding knot
holds its topology through both halves, and a knot can be kept pulsing
indefinitely. The cost is an inward wave, supplied from outside, of one
knot's energy per period (here U·c/(8 L0)), which the ledger now books as
inflow. That is the standard prediction too: the equations are
time-symmetric, so the contraction is the expansion reversed. What the
run does not supply is a source for the inward wave. Candidates to test
next, each as a declared medium with its own ledger column: the absorber
response of distant matter (Wheeler & Feynman, *Rev. Mod. Phys.* 17, 157,
1945), the zero-point field (inward, but random; compare H28 Part 4), and
a vacuum or ether that is itself nonlinear (Born–Infeld 1934;
Euler–Heisenberg 1936, the measured QED vacuum), which could also hold the
knot in place (H29 step 2).

**Part 3: the measured QED vacuum as the source (`vacuum.go`).** The one
nonlinear vacuum with experimental support is Euler–Heisenberg's
(light-by-light scattering: ATLAS, *Nature Phys.* 13, 852, 2017; CMS,
*Phys. Lett. B* 797, 134826, 2019). In it waves act on one another only
where they cross, through E² − c²B² and E·B; a lone plane wave has both
zero and never acts on itself. Both measurements run on Kron's plane-wave
network with the Euler–Heisenberg vacuum as its nonlinear primitive
(`kron.EHLine`, validated on light-by-light phase shifts to 0.4 %), each
beside Maxwell's linear vacuum, at 0.08–0.09 of the Schwinger field
(1.3×10¹⁸ V/m), the top of the range the theory's weak-field form covers:

| Measurement | Euler–Heisenberg vacuum | Needed for a knot rebuilt every period |
|---|---|---|
| **A.** A compressed field at rest explodes into two halves; nonlinear energy back between them after they separate | 6×10⁻²⁸ of the explosion (round-off; the linear lattice's own ripple there is 3×10⁻¹⁷) | 1 |
| Nonlinear part of the field after the halves crossed (phase shifts, all outgoing) | 3×10⁻¹³ | — |
| **B.** A probe crossing a sharp counter-propagating field of 0.08 E_s; energy turned round | 2.4×10⁻¹⁷ (Fresnel bound for its index step 16a: 1.7×10⁻¹²) | 1 |
| Nonlinear share of the field energy at the peak, the most the vacuum could redirect | 1.3×10⁻⁶ (≈ 6a) | 1 |

- **Nothing comes back.** The halves interact while they cross, but what
  the vacuum gives them is a phase shift that travels outward with them;
  once apart, each is a plane wave the vacuum no longer touches.
- **In three dimensions** the knot's expanding shell is locally a plane
  wave within a few L0, so the same holds; whatever the vacuum could ever
  turn round is bounded by the nonlinear part of the field energy, about
  10⁻⁶ of the knot at the edge of the theory's range, against the whole
  knot per period the pulsing knot needs.

So the measured QED vacuum is not the source of the inward wave: H30's
pulsing knot still needs it supplied from outside, and this route is
bounded at ≲ 10⁻⁶ of a knot per period (computed in 1-D: below 10⁻¹⁶).
The other declared media (the Wheeler–Feynman absorber, the zero-point
field, Born–Infeld) are untested.

`viz/h30_cycle.png` plots both parts: energy within 2 L0 (solver and
exact) and the surface's cumulative work in and out.

## Pictures

`viz/knot_3d.png`: 24 magnetic field lines of the exact hopfion at t = 0,
seeded on three nested rings (light: near the axis, dark: near the core);
all close, and the first pair links once. `viz/knot_time.png`: the same
family at t = 0, 1, 2 and 3 L0/c on one scale, seeds scaled by the knot's
size √(1 + t²), with a grey sphere at the energy's rms radius (1.00, 1.41,
2.23, 3.16 L0); all 24 lines close at every time and the linking stays
1.000. Seeds sit half a step off the y-axis, whose field line in this
solution passes through infinity. `viz/knot_lines.json` holds the traced
lines for the pictures and for an interactive viewer.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestLinkingKnownCurves` | Gauss linking integral: Hopf link = ±1, separated circles = 0 | 1e-3 |
| `TestHopfionExact` | the field is null (F·F = 0); traced B lines close and link once at t = 0 and t = 3 | 1e-12; 1e-3 |
| `TestKnotInSolver` | coarse grid: energy spreads as the exact solution says; helicity kept as exact; curl A = B; ledger | 0.02; 0.02; 1e-12; 1e-10 |
| `TestKnotEnergy` | the knot's energy density at t = 0 is spherically symmetric; radial total = 3-D sum | 1e-12; 2e-3 |
| `TestImplosionThroughSurface` | H30, coarse grid: ≥ 0.9 of a knot enters through the Huygens surface, the centre reaches the exact energy, linking 1 at t = 0, in = out, ledger | 0.05; 0.02; 0.05; 1e-10 |
| `TestVacuumInward` | H30 Part 3: in the Euler–Heisenberg vacuum the halves of an explosion interact while crossing but return nothing between them (below the linear lattice's own ripple), a sharp counter-propagating field reflects less than its Fresnel bound, the nonlinear energy share is 6a | 1e-24; Fresnel; 1 % |
| `em.TestFDTDMagneticSourceLedger` | the solver's magnetic face sources: ledger closes with them, sources in the absorber refused | 1e-10 |

This example has no circuit, so it has no `export_test.go`.

## Files

| File | What it holds |
|---|---|
| `cycle.go` | H30: the Huygens surface (edge and face sources from the declared incoming train), part 1 and part 2 runs |
| `vacuum.go` | H30 Part 3: the Euler–Heisenberg vacuum (`kron.EHLine`) as the source of the inward wave — explosion and reflection, beside the linear vacuum |
| `viz/render_h30.py` | `viz/h30_cycle.png` from `run_output.txt` |
| `main.go` | the hopfion (Bateman form, helicity sign check), the vector-potential and helicity bookkeeping beside `em/fdtd.go`, field-line tracing and the Gauss linking integral, the run and report |
| `main_test.go` | the tests above |
| `viz.go`, `viz/render.py` | `-viz`: the traced 3-D lines (`viz/knot_lines.json`) and the two pictures |
| `run_output.txt` | the run above |

## Sources

- A. F. Rañada, "A topological theory of the electromagnetic field," *Lett. Math. Phys.* 18, 97–106 (1989).
- W. T. M. Irvine, D. Bouwmeester, "Linked and knotted beams of light," *Nature Physics* 4, 716–720 (2008).
- A. Taflove, S. C. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*, 3rd ed. (Artech House, 2005), ch. 5 — the total-field/scattered-field surface.
- J. A. Wheeler, R. P. Feynman, "Interaction with the absorber as the mechanism of radiation," *Rev. Mod. Phys.* 17, 157–181 (1945).
- W. Russell, *The Universal One* (1926); *The Secret of Light* (1947) — compression and expansion as "rhythmic balanced interchange" (qualitative; no equations).
- M. Born, L. Infeld, *Proc. R. Soc. A* 144, 425 (1934); W. Heisenberg, H. Euler, *Z. Phys.* 98, 714 (1936) — nonlinear vacuum electrodynamics.
- H. Kedia, I. Bialynicki-Birula, D. Peralta-Salas, W. T. M. Irvine, "Tying knots in light fields," *PRL* 111, 150404 (2013) — the Bateman construction used here.
- H. K. Moffatt, "The degree of knottedness of tangled vortex lines," *J. Fluid Mech.* 35, 117–129 (1969) — helicity as linking.
- G. H. Derrick, *J. Math. Phys.* 5, 1252 (1964); T. H. R. Skyrme, *Proc. R. Soc. A* 260, 127 (1961); L. Faddeev, A. J. Niemi, "Stable knot-like structures in classical field theory," *Nature* 387, 58–61 (1997).
- W. Thomson (Kelvin), "On vortex atoms," *Phil. Mag.* 34, 15–24 (1867).
