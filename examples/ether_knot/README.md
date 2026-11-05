# Ether knot: can a nonlinear ether hold a knot at a fixed size? (H29, step 2)

Hypothesis H29 (CLAUDE.md §1): particles are knots of field energy, held
together by their topology. Step 1 ([knotted_light](../knotted_light))
found that linear Maxwell keeps a knot's topology but not its energy: the
knot links once forever while flying apart at the speed of light. H30
found that it can be kept pulsing only by an inward wave supplied from
outside. Step 2 asks whether a nonlinear, fluid-like ether can hold a knot
by itself, at a size of its own.

```
go test -short ./examples/ether_knot     # field link, gradients, 2-D vs 3-D (< 1 s)
go test ./examples/ether_knot            # + the relaxation and controls (~15 s)
go run ./examples/ether_knot             # grid study, controls, 3-D check (~30 min on 2 cores)
```

## Research first: which ether?

| Candidate | What the literature says | Can it hold a knot at rest? |
|---|---|---|
| Light itself made nonlinear: Born–Infeld (1934), Euler–Heisenberg (1936, the QED vacuum, seen as light-by-light scattering by ATLAS, *Nature Phys.* 13, 852, 2017) | The only finite-energy, source-free, stationary solution of Maxwell–Born–Infeld is zero (Kiessling, arXiv:1107.2333, Theorem 1), so no such soliton moves slower than light | No |
| An ordinary superfluid (Gross–Pitaevskii), Kelvin's vortex atoms | All 1,458 simulated vortex knots untie by reconnection, the vortex cores being places where the field vanishes and can be cut (Kleckner, Kauffman & Irvine, *Nature Phys.* 12, 650, 2016) | No |
| A fluid whose order parameter never vanishes: a unit vector n(x). A charged two-component superfluid reduces to it (Babaev, Faddeev & Niemi, *Phys. Rev. B* 65, 100512, 2002) | Faddeev's knots (Faddeev & Niemi, *Nature* 387, 58, 1997); stable, computed by Battye & Sutcliffe (*PRL* 81, 4798, 1998): Q = 1 energy 504 on their 100³ lattice | Yes: tested here |

**The link to knotted light.** The ether's energy is

```
E = κ2 ∫ (∂_i n)² d³x  +  κ4 ∫ |B_n|² d³x,      B_n,k = ½ ε_kij n·(∂_i n × ∂_j n)
```

For Rañada's Hopf map, the starting field here, B_n is exactly Rañada's
knotted magnetic field times −16, at every point
(`TestEtherCurvatureIsKnottedLight`). So the second term is the knotted
light's magnetic energy, the term step 1 already had, and the first is
the new ingredient: the ether's stiffness. Under a change of size
x → λx the stiffness grows as λ and the magnetic term falls as 1/λ
(Derrick 1964): alone, the magnetic term spreads the knot and the
stiffness shrinks it; together they fix a size, where the two are equal
(the virial theorem). Units: κ2 = κ4 = 1, so lengths are in
ℓ = √(κ4/κ2) and energies in √(κ2κ4).

## How it is computed

The Q = 1 knot is axially symmetric (Battye & Sutcliffe), with the field
turning once round the axis: n = R_z(φ) m(ρ, z). Rañada's map has exactly
this form. So the energy is minimized on the (ρ, z) half-plane, a
thousand times cheaper than in 3-D (`axi.go`), with m kept a unit vector
(no angle singularities at the core or on the axis). Relaxation is
arrested Newton flow, as Battye & Sutcliffe: the field accelerates
downhill and is stopped whenever the energy rises. The start is Rañada's
map pulled in smoothly to reach the vacuum inside the box (the map itself
approaches the vacuum only as 1/r; holding a box edge at the vacuum
would cut it). The relaxed knot is then rebuilt on a 3-D grid to check
its energy and its topology (`main.go`).

## Results

**1. Stiffness + magnetic: the knot settles at a fixed size.**

| Box (ℓ) | Grid h (ℓ) | Energy | Stiffness E2 | Magnetic E4 | E2/E4 | Core radius (ℓ) |
|---|---|---|---|---|---|---|
| 10 | 0.10 | 541.08 | 261.15 | 279.94 | 0.933 | 1.15 |
| 16 | 0.10 | 540.24 | 262.11 | 278.14 | 0.942 | 1.15 |
| 16 | 0.07 | 544.00 | 268.06 | 275.94 | 0.971 | 1.23 |
| 16 | 0.05 | 545.70 | 270.66 | 275.05 | 0.984 | 1.23 |
| h → 0 | (extrapolated in h²) | **547.5** | | | **0.997** | |

- The box edge hardly matters (1 part in 500 from 10 ℓ to 16 ℓ); the grid
  does, and as h → 0 the stiffness and magnetic terms become equal to
  0.3 %: Derrick's balance, reached by the relaxation without being
  imposed.
- **Energy 547.5·√(κ2κ4), core radius ≈ 1.2 ℓ.** Battye & Sutcliffe's
  100³ lattice gave 504 (their spacing is not stated). Coarse grids give
  lower energies here too (540 at h = 0.1), so the 8 % gap is recorded,
  not explained away; a continuum value from the literature in the same
  normalization is still to be checked.

**2. Controls (Derrick's theorem, in the same code).**

- *Magnetic energy alone* (κ2 = 0: knotted light with no ether): the
  knot spreads, rms radius 1.0 → 5.1 ℓ, energy 312 → 0.09. Step 1's
  result, now as an energy minimization.
- *Stiffness alone* (κ4 = 0): the knot shrinks, core radius
  0.95 → 0.75 → 0.35 ℓ in 200 steps; then its core reaches the grid
  spacing and it unties through the lattice (the core's m_z goes from
  −1 to +0.99, the energy falls from 271 to 4.9 by step 1000). In the continuum it shrinks to a point. A
  grid cannot hold topology below its own spacing.

**3. The relaxed knot in 3-D.** The finest knot (box 16 ℓ, h = 0.05)
rebuilt on 101³ over ±8 ℓ: two pairs of preimage curves of n (the field
lines of B_n), each closed at its nearest return, link at 1.0001 and
1.0002 (gaps 0.019 and 0.011 ℓ): Hopf charge 1, the knot intact. Its 3-D
energy reads 297.3 + 255.4 against 270.7 + 275.0: the ±8 ℓ box cuts the
tail of a knot relaxed in a 16 ℓ box (stiffness high) and the 3-D grid is
three times coarser (magnetic term low by O(h²)). A knot relaxed in an
8 ℓ box and rebuilt the same way matched to 0.5 % in stiffness (258.7
against 260.0), with the magnetic term 6 % low.

## Vibrations: the knot's own tone

Struck, the knot rings at one clear frequency: squashing (taller and
thinner, then back) at **ω ≈ 1.36–1.42 c/ℓ with Q ≈ 6–7**; its other
motions die within a cycle or two.

```
go run ./examples/ether_knot -ring examples/ether_knot/vib     # ~8 min on 2 cores
python3 examples/ether_knot/ring.py examples/ether_knot/vib/ring_*.json
python3 examples/ether_knot/modes.py examples/ether_knot/vib/ops_b16_h0.15.txt
```

**How.** Small motions about the relaxed knot obey 2K ẍ = −H x, from
Faddeev's Lagrangian (kinetic energy ∫|ṅ|² + Σ|ṅ × ∂_i n|²): H is the
energy's Hessian on the sphere, K the kinetic form, both built in
`modes.go` for motions that keep the axial symmetry. Two motions cost
nothing and come out at ω² ≈ 0 against the real modes' ~1 (sliding along
z, about 2×10⁻⁴; turning the field about the vacuum direction, 10⁻⁷).

**Why striking, not just eigenvalues.** Small waves of this ether travel
at c with no lowest frequency, as light does, so a vibrating knot leaks
energy as waves and has no sealed tone. In a box, the low eigenmodes are
standing waves of the box: the lowest falls 0.428 → 0.324 → 0.252 c/ℓ as
the box grows 8 → 12 → 16 ℓ, and its share within 3 ℓ of the knot falls
42 % → 9 % → 2 %. So the knot was struck instead, with three localized
pushes (breathing, squashing, twisting, under exp(−r²/8ℓ²), the zero modes
removed), left to ring, and its motion fitted to damped cosines (matrix
pencil, `ring.py`) up to the time the first echo from the wall returns:

| Push | Box 16 ℓ, h 0.15 | Box 24 ℓ, h 0.15 | Box 16 ℓ, h 0.10 |
|---|---|---|---|
| **Squashing** | ω 1.416, Q 7.0 | ω 1.415, Q 7.0 | ω 1.385, Q 6.0 |
| Breathing | ω 0.43, Q 1.9 | ω 0.45, Q 1.9 | ω 0.50, Q 2.0 |
| Twisting | ω 0.93, Q 1.7 | ω 0.86, Q 1.1 | ω 0.83, Q 1.3 |

- **Squashing is the knot's tone.** A bigger box leaves it unchanged
  (1.416 → 1.415), the finer grid lowers it 2 %; extrapolated in h², about
  1.36 c/ℓ. Q ≈ 6–7: the energy falls by e in about one cycle, so it rings
  a handful of times.
- **Breathing,** the size change Derrick's balance fixes, is a broad
  resonance near 0.45–0.5 c/ℓ, gone in about two cycles (Q ≈ 2).
- **Twisting** hardly rings (Q ≈ 1): no tone of its own.
- **What it means.** The knot has a frequency of its own, set by its size:
  f ≈ 1.4 c/(2πℓ), with ω·(core radius)/c ≈ 1.7. Since no measurement fixes
  ℓ, it is a ratio, not yet a number: for ℓ = 1 m it would be 67 MHz,
  for ℓ = 1 fm, 6.7×10²² Hz. The squashing-to-breathing ratio, about 3,
  needs no constant at all.

## What this means

A nonlinear ether of this kind does what linear Maxwell cannot: it holds
a knot of field energy, whose magnetic part is knotted light, at a fixed
size and energy, with its linking number intact and no inflow needed.
That is the standard result (Faddeev & Niemi 1997), now reproduced in
kronspice and tied to the knotted light of step 1.

What it does not do:

- **It does not fix κ2 and κ4.** The knot's size is ℓ = √(κ4/κ2) and its
  energy 547.5·√(κ2κ4); no measured "ether stiffness" exists to set them.
  Matching one particle's mass and size uses both numbers up, so the
  model's first prediction would be a ratio between two particles (for
  example the Q = 2 knot against Q = 1: Battye & Sutcliffe 835 against
  504).
- **Its frequency is a ratio.** The knot rings at ω ≈ 1.4 c/ℓ (above),
  which becomes a number only once ℓ is fixed.
- **It is not electromagnetism.** n is an ordinary field with its own
  energy. B_n has the form of a magnetic field, but how n would couple to
  electric charges, and whether such a knot carries charge or spin ½, is
  open (Rañada's charge quantization and the Skyrme model's fermion
  quantization are the leads).

Next steps: the Q = 2 knot (a prediction that uses no new constant), the
vibrations that break the axial symmetry, and a Faddeev ether coupled to
the electromagnetic field of `em/fdtd.go`.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestEtherCurvatureIsKnottedLight` | B_n of Rañada's map = −16 × his knotted B, at four points | 1e-6 |
| `TestGradient` | analytic gradients of the 3-D and axisymmetric energies against finite differences | 1e-6; 1e-5 |
| `TestAxiEqualsThreeD` | axisymmetric and 3-D energies of the same field (E4 converges as h²) | 1 %; 6 % |
| `TestVibrations` | coarse grid: the two zero modes at ω² ≈ 0; a squashing push rings at ω ≈ 1.45 c/ℓ (zero crossings) | 0.02; 0.15 |
| `TestEtherHoldsKnot` | coarse grid: both terms settle near the virial balance with a core ~ℓ; magnetic only spreads; stiffness only shrinks, then unties through the lattice | 0.2; ×1.1; ×0.9, 5 % |

## Files

| File | What it holds |
|---|---|
| `main.go` | the research summary, the 3-D grid (energy, gradient, preimage linking), the compact Hopf start, the run |
| `axi.go` | the axisymmetric energy, gradient, relaxation, size measures, and the rebuild in 3-D |
| `modes.go` | the linearized motion: Hessian H and kinetic form K, zero-mode checks, the three strikes and the ringing (`Ring`) |
| `ring.py`, `modes.py` | damped-cosine fit of the ringing (matrix pencil); lowest eigenmodes and their share near the knot (scipy) |
| `vib/`, `vib_output.txt` | the ringing records of the `-ring` run and its analysis (`ops_*.txt` operators are left out of git) |
| `main_test.go` | the tests above |
| `run_output.txt` | the run |

## Sources

- L. Faddeev, A. J. Niemi, "Stable knot-like structures in classical field theory," *Nature* 387, 58–61 (1997).
- R. A. Battye, P. M. Sutcliffe, "Knots as stable soliton solutions in a three-dimensional classical field theory," *PRL* 81, 4798 (1998), [hep-th/9808129](https://arxiv.org/abs/hep-th/9808129).
- E. Babaev, L. D. Faddeev, A. J. Niemi, "Hidden symmetry and knot solitons in a charged two-condensate Bose system," *Phys. Rev. B* 65, 100512 (2002), [cond-mat/0106152](https://arxiv.org/abs/cond-mat/0106152).
- D. Kleckner, L. H. Kauffman, W. T. M. Irvine, "How superfluid vortex knots untie," *Nature Phys.* 12, 650 (2016), [arXiv:1507.07579](https://arxiv.org/abs/1507.07579).
- M. K.-H. Kiessling, uniqueness of stationary Maxwell–Born–Infeld solutions, [arXiv:1107.2333](https://arxiv.org/abs/1107.2333).
- ATLAS Collaboration, "Evidence for light-by-light scattering in heavy-ion collisions," *Nature Phys.* 13, 852 (2017).
- G. H. Derrick, *J. Math. Phys.* 5, 1252 (1964); A. F. Rañada, *Lett. Math. Phys.* 18, 97 (1989).
- Y. Hua, T. K. Sarkar, "Matrix pencil method for estimating parameters of exponentially damped/undamped sinusoids in noise," *IEEE Trans. ASSP* 38, 814 (1990).
- G. E. Volovik, *The Universe in a Helium Droplet* (Oxford, 2003) — the vacuum as a quantum fluid.
