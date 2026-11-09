# The Claude electron: the toroidal photon held by a knotted ether (H44)

Hypothesis H44 (CLAUDE.md §1): the electron is H23's toroidal photon (one
Compton-wavelength photon circulating twice round a torus, Williamson &
van der Mark 1997) held together by H29's nonlinear ether (a Faddeev knot
whose magnetic term is knotted light). A Hestenes zitterbewegung electron
is the standard twin.

```
go test -short ./examples/claude_electron       # start degrees, gradient with winding, bench numbers (< 2 s)
go test ./examples/claude_electron              # + a coarse relaxation: steps 1, 2, 4 (~10 s)
go run ./examples/claude_electron               # h = 0.1 ℓ (~2 min); -h 0.15 for a quick run
```

**Result: bounded.** The ether holds a single-loop knot, but not H23's double
loop. The double loop has Hopf charge 2, so it could not be a fermion even if
it were held. The single loop can carry spin ½, but no size for it satisfies
the measured bounds: spin and Coulomb energy need ℓ ≳ 444 fm, while
scattering, excited-electron searches and hydrogen need ℓ ≲ 0.05 fm.
Charge is still not derived.

## Why this combination

The electron test bench (book chapter 9) left H23 ahead. Its gaps, and the
other models that address them:

| Gap in H23 | What another model contributes | Where |
|---|---|---|
| Charge 0.910 e drives the anomaly and hydrogen misses | With q = e the anomaly is α/2π and the atom gives Bohr's levels | `light_hydrogen` (H32) |
| Not self-confined under linear Maxwell | A unit-vector ether holds a knot (E = 547.5 √(κ2κ4)) whose magnetic term is knotted light | `ether_knot` (H29) |
| Spin ½ of a soliton | Hopf solitons of odd charge can be quantized as fermions (Krusch & Speight 2006) | here, step 2 |
| Compton size against the < 10⁻¹⁸ m scattering bound | Hestenes: a point charge circulating at c on a Compton helix | here, step 0 |

## How it is computed

`axi.go` is `ether_knot`'s axisymmetric relaxer, generalized to Battye &
Sutcliffe's A_{M,K} family: n = R_z(Mφ) m(ρ, z), with the field turning M
times round the axis and K times round the core. The Hopf charge is Q = M·K.
The preimage curves (the field lines of B_n) lie on tori and close after
K turns the long way round and M turns through the hole. So:
- Williamson & van der Mark's double loop (twice round the ring, once
  round the tube) is **A_{1,2}, Q = 2**;
- Rañada's knot is A_{1,1}.

The start is the rational map W = Z1^M / Z0^K on the 3-sphere. The charge is
followed as the degree (1/4π)∫∫ m·(m_ρ × m_z) dρ dz. Relaxation uses arrested
Newton flow, as in `ether_knot`.

`spin.go` gives Faddeev's kinetic energy T = κ2∫|ṅ|² + κ4∫Σ|ṅ × ∂_i n|² for
rigid rotations: I∥ about the axis, where rotation equals isorotation, and I⊥
averaged over φ. Isorotations that move the vacuum have infinite inertia and
are frozen.

Finkelstein–Rubinstein quantization allows a fermion only for odd Q (Krusch &
Speight). The lowest state then has J = ½, |J3| = ½, with
E_rot = ħ²/(4I⊥) + ħ²/(8I∥).

Units: energy E0 = √(κ2κ4), length ℓ = √(κ4/κ2).

## Results (h = 0.1 ℓ, `run_output.txt`)

**Step 0: standard twin.** Hestenes' zitter radius ħ/2mc = 193.08 fm is
exactly H23's loop radius λ/4π. A charge e moving at c on that circle has
μ = ecr/2 = μ_B/2, so with S = ħ/2 its **g is 1.000**, not 2. Classical
circulation gives g = 2 only on radius ħ/mc, where an orbiting mass m would
carry ħ rather than ħ/2. Dirac's g = 2 comes from the spin current (Gordon
decomposition), not from circulation. H23's eq. 12, g = 2(1 + α′/2π), is
asserted in the paper, not obtained from its own loop.

**Step 1: confinement.**

| Knot | Q | Degree start → end | E / E0 | E2/E4 | Core / ℓ |
|---|---|---|---|---|---|
| A11 single loop | 1 | 0.985 → 0.988 | 541.08 | 0.933 | 1.150 |
| A12 double loop (H23) | 2 | 1.934 → 0.143 after 8000 steps (still falling) | 1144 → 685 | — | falls to 0.15 (the grid) |
| A21 (windings swapped) | 2 | 0.974 → 0.989 (×M = 2) | 884.64 | 0.957 | 1.750 |

- **The double loop is not held.** Its core falls onto the axis and its
  charge leaves through the grid, both at h = 0.15 (back to Q = 1 within
  1000 steps) and at h = 0.1. Within the axial ansatz with M = 1 it has no
  minimum. The held Q = 2 knot is A21, which turns twice round the axis
  and once round the core: the opposite winding to the double loop. In
  full 3-D the double loop is free to turn into A21; that run was not
  made here.
- A11 equals `ether_knot` at the same grid (541.08) to the printed digits.
  E(A21)/E(A11) = 1.635, near the Q^¾ scaling (1.682).
- Control: the double loop as knotted light alone (κ2 = 0) spreads, rms
  radius 1.02 → 4.23 ℓ.

**Step 2: spin.**

| Knot | Fermion (FR) | I∥ | I⊥ (E0ℓ²/c²) | least ℓ for spin ½ within m_e c² |
|---|---|---|---|---|
| A11 | yes | 473.3 | 721.1 | **444 fm = 1.15 λ̄** |
| A21, A12 | no (Q = 2) | 2132 | 1866 | — |

m_e c² = E·E0 + c_r(ħc/ℓ)²/E0 has a real E0 only for ℓ ≥ 2√(E·c_r)·λ̄.
At that size half the mass is rotation and the rim turns at 0.65 c. A
spin-½ knot of the electron's mass cannot be smaller than its Compton
length.

**Step 3: charge.**

| Option | q / e | Anomaly (meas. 0.00115965) | Lyman-α (meas. 121.567 nm) | Derived |
|---|---|---|---|---|
| H23 eq. 5, √(3ε0ħc)/2π | 0.9103 | 0.000962 | 146.70 nm | yes |
| Rañada's topological unit √(ε0ħc) | 3.3023 | 0.0127 | 11.15 nm | yes |
| declared e | 1 | 0.0011614 | 121.57 nm | no |

The Faddeev ether carries no charge. In the gauged version (Babaev, Faddeev &
Niemi) the pair charge is an input. None of the three options derives e.

**Step 4: κ2, κ4 fixed by m_e c² and a size (A11).**

| Reading | ℓ | Spin energy / m_e c² | Excited electron | Coulomb / m_e c² | 1S–2S shift | g (rigid) |
|---|---|---|---|---|---|---|
| A: ring = λ̄/2 | 167.9 fm | 1.75 (no real E0), rim 1.72 c | 1.60 MeV | 0.007 | 127 GHz | 2.81 |
| B: ℓ = 10⁻¹⁸ m | 0.001 fm | 4.9×10¹⁰, rim 2.9×10⁵ c | 268 GeV | 1230 | 4.5 Hz | 2.81 |
| ℓ_min (spin) | 444 fm | 0.25 (½ the mass at the boundary), rim 0.65 c | 0.60 MeV | 0.003 | 887 GHz | 2.81 |
| Measured | — | — | none below 5.6 TeV | — | QED + proton to ~kHz | 2.00232 |

The rigid g assumes the charge is spread like the knot's energy and turns
with it: g = m⟨ρ²⟩/I∥, independent of ℓ. It is 2.81, so the ESR line would
be at 39.35 GHz/T against the measured 28.025. A mechanical ring with charge
and mass distributed alike would give 1.

**Windows on ℓ:**

| Constraint | Allows ℓ |
|---|---|
| scattering, rms < 10⁻¹⁸ m | ≤ 5.5×10⁻⁴ fm |
| no excited electron below 5.6 TeV (CMS, JHEP 05 (2020) 052) | ≤ 4.8×10⁻⁵ fm |
| hydrogen 1S–2S finite size < 10 kHz (declared) | ≤ 0.047 fm |
| Coulomb energy of e at the ring < m_e c² | ≥ 1.23 fm |
| spin ½ within m_e c² | ≥ 444 fm |

The overlap is empty: the lower bounds exceed the upper bounds by 9×10⁶.
Reading B, the small knot on Hestenes' helix, fails on the knot's own spin.
A knot that small turns faster than light, unless its spin is moved off the
knot and into the orbit, which then carries ħ, not ħ/2 (step 0).

**Step 5: register.** H44 is *bounded*. The bounding numbers are: Q = 2 for
the double loop (bosonic, and not held in the axial ansatz); ℓ ≥ 1.15 λ̄
for a spin-½ knot of the electron's mass, against ℓ ≤ 5×10⁻⁵ fm from
excited-state searches; and charge not derived.

What would move it:
- a confining medium in which the double loop has a minimum and that also
  lets Q = 2 carry spin ½ (Faddeev's does neither);
- or a mechanism that gives the electron spin without rigid rotation;
- and a derivation of e.

## Standard prediction

The Dirac electron with QED: pointlike below 10⁻¹⁸ m, g = 2.002 319 304 36,
no excited states below 5.6 TeV. Charge is an input in the Standard Model
too; its quantization is tied to anomaly cancellation or to grand
unification.

## Tests

| Test | Checks | Tolerance |
|---|---|---|
| `TestStartDegree` | A_{M,K} starts carry degree K on the grid | 8 % (grid) |
| `TestGradientWinding` | energy gradient with winding M = 1, 2 vs finite differences | 1e-5 |
| `TestBenchValues` | g = 1 for a charge at c on λ/4π; WvdM charge 0.910318; Rañada 3.30227 e; proton 1S–2S shift ≈ 1 MHz | 1e-12 / 2e-6 / 1e-4 / 0.03 |
| `TestElectronKnot` | A11 keeps degree 1; I∥ both ways agree; ℓ_min > λ̄; windows empty; reading A has no real E0; double loop unties | as stated |

## Files

| File | Contents |
|---|---|
| `axi.go` | A_{M,K} axisymmetric Faddeev knots (from `ether_knot`, winding M added), rational-map starts, degree, relaxation |
| `spin.go` | inertia tensors, energy-weighted radii, FR spin coefficient |
| `bench.go` | step 0 twin, charge options, readings, windows, hydrogen finite-size shift |
| `main.go` | steps 0–5 |
| `main_test.go` | tests above |
| `run_output.txt` | the h = 0.1 run |

## Sources

- J. G. Williamson, M. B. van der Mark, *Ann. Fond. Louis de Broglie* 22, 133 (1997).
- L. Faddeev, A. J. Niemi, *Nature* 387, 58 (1997); R. A. Battye, P. M. Sutcliffe, *Phys. Rev. Lett.* 81, 4798 (1998) and *Proc. R. Soc. A* 455, 4305 (1999) (A_{M,K} family).
- E. Babaev, L. D. Faddeev, A. J. Niemi, *Phys. Rev. B* 65, 100512 (2002).
- S. Krusch, J. M. Speight, "Fermionic quantization of Hopf solitons", *Commun. Math. Phys.* 264, 391 (2006), arXiv:hep-th/0503067; D. Finkelstein, J. Rubinstein, *J. Math. Phys.* 9, 1762 (1968).
- A. F. Vakulenko, L. V. Kapitanski, *Sov. Phys. Dokl.* 24, 433 (1979): E ≥ c·Q^¾.
- D. Hestenes, *Found. Phys.* 20, 1213 (1990).
- A. F. Rañada, *Lett. Math. Phys.* 18, 97 (1989); M. Arrayás, J. L. Trueba, A. F. Rañada, "Topological electromagnetism: knots and quantization rules", IntechOpen (2012).
- CMS Collaboration, "Search for an excited lepton that decays via a contact interaction to a lepton and two jets", *JHEP* 05 (2020) 052, arXiv:2001.04521.
- X. Fan et al., *Phys. Rev. Lett.* 130, 071801 (2023): measured g.
- G. S. Adkins, C. R. Nappi, E. Witten, *Nucl. Phys. B* 228, 552 (1983): collective-coordinate quantization of a soliton.
