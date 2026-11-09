# claude_electron (H44): issues report

Date: 2026-11-08. Covers `examples/claude_electron` as it stands in the
working tree (not yet committed): `axi.go`, `spin.go`, `bench.go`,
`main.go`, `main_test.go`, `README.md`, `run_output.txt`, and the H44 text
in `CLAUDE.md` and `book/09-the-electron-test-bench.md` §9.8.

State at the time of writing: `go vet` is clean and all four tests pass
(`go test ./examples/claude_electron`, 7.6 s). The example's verdict,
*bounded*, still stands after the checks below. Several of the numbers
quoted for it do not: they are wrong, or they have not converged.

The issues come in three groups:

- **A.** What the hypothesis itself ran into. These are the findings that
  bound H44.
- **B.** Numerical issues: results that depend on the grid, the box or
  the ansatz.
- **C.** Errors and inconsistencies in the code and the write-up, with
  fixes.

---

## A. Issues with the hypothesis (physics findings)

| # | Issue | Evidence | Consequence |
|---|---|---|---|
| A1 | **H23's double loop is not held by the ether.** In the field's own winding numbers, Williamson & van der Mark's double loop (twice round the ring, once round the tube) is A_{1,2}, Hopf charge Q = 2. | Degree 1.934 → 0.143 in 8000 steps at h = 0.1 ℓ. The core falls to 0.15 ℓ (the grid) within 1000 steps. At h = 0.15 it returns to Q = 1 within 1000 steps. | The combination H23 + H29 fails at its first step for the configuration H23 actually describes. |
| A2 | **The charge-2 knot that does hold winds the other way.** A_{2,1} goes twice round the axis and once round the core. | E = 884.64 E0 = 1.635 × A11 (Vakulenko–Kapitanski Q^¾ gives 1.682). | It is not H23's geometry. |
| A3 | **Any Q = 2 knot is a boson.** Finkelstein–Rubinstein quantization allows spin ½ only for odd Q (Krusch & Speight 2006). | Applied by rule in `main.go` (`Fermion: q%2 == 1`). | Even if A1 were fixed, the double loop could not be the electron. |
| A4 | **The single loop (A11) can be a fermion but cannot be small.** Its J = ½ rigid rotation fits inside m_e c² only for ℓ ≥ 2√(E·c_r)·λ̄. | ℓ_min = 444 fm = 1.15 λ̄ (h = 0.1, box 10 ℓ). | — |
| A5 | **Size windows do not overlap.** | Lower bounds: spin ℓ ≥ 444 fm; Coulomb ℓ ≥ 1.23 fm. Upper bounds: excited electron ℓ ≤ 4.8×10⁻⁵ fm, scattering ≤ 5.5×10⁻⁴ fm, hydrogen 1S–2S ≤ 0.047 fm. | Gap ≈ 9×10⁶. No choice of κ2, κ4 survives. |
| A6 | **Charge is not derived.** | H23 eq. 5 gives 0.910 e (anomaly 17 % low, Lyman-α 146.7 nm). Rañada's topological unit gives 3.30 e. The gauged ether (Babaev–Faddeev–Niemi) takes the charge as an input. | Only a declared e reproduces the anomaly and Bohr's levels. |
| A7 | **H23's g = 2 is asserted, not computed (step 0).** A charge e at c on H23's loop radius λ̄/2 (equal to Hestenes' zitter radius) has μ = μ_B/2, so with S = ħ/2, g = 1.000. | `GCirculation`, `TestBenchValues` (1e-12). | Classical circulation gives g = 2 only at radius λ̄, where an orbiting mass would carry ħ rather than ħ/2. |
| A8 | **The rigid-rotor g is wrong too.** Charge spread like the energy and turning with the knot gives g = E⟨ρ²⟩/I∥. | 2.81 at h = 0.1, box 10, which puts ESR at 39.35 GHz/T against the measured 28.025. Not converged (B4). | — |
| A9 | **Small-knot reading B (ℓ = 10⁻¹⁸ m) is not physical.** | Spin energy 4.9×10¹⁰ m_e c², rim 2.9×10⁵ c, Coulomb 1230 m_e c². | Spin would have to move from the knot into an orbit, and the orbit then carries ħ (A7). |

None of the examples tested 3-D relaxation, a different confining medium,
or a mechanism for spin without rigid rotation. These are the open
directions recorded in the README.

---

## B. Numerical issues

### B1. The double loop's "untying" is a lattice process, and no h → 0 study was made

In the continuum the Hopf charge cannot change. On the grid it changes
only when the core shrinks to the grid spacing. The run shows exactly
this: the core sits at 0.15 ℓ at every printout. The physical statement
behind it is that **A_{1,2} has no minimum in the axial ansatz: it
collapses**. Only two grids were tried (h = 0.15 and 0.10), and the run
stopped while E was still falling (685 at 8000 steps). The README wording
"has no minimum" is supported, but only qualitatively. No study shows the
collapse rate is independent of h.

### B2. Only the axial ansatz was tested for A12

In full 3-D the double loop may relax into A21, or into a twisted
unknot or link. Battye & Sutcliffe's Q = 2 minimum is axial A21, which
agrees with this run. The README says this 3-D run "was not made here".
Until it is made, A1 is a statement about the ansatz, not about the
Faddeev model.

### B3. Energy and degree do not fall monotonically during A12's collapse

E 892.0 → 895.1 (steps 1000 → 2000), and degree 0.179 → 0.181 (steps
6000 → 7000). Arrested Newton flow keeps the step on which the energy
rises and zeroes the velocity only afterwards, and every `Relax` call
starts again from rest. This is consistent with a lattice untying, not a
smooth descent. It is harmless, but the README table shows only the end
points.

### B4. The inertias, and everything derived from them, are not converged in box or grid

The README reports only one grid (h = 0.1) and one box (10 ℓ). The
coarse test grid (h = 0.15, box 8) gives I∥ = 385 against 473, a 19 %
difference. A scratch convergence run made for this report (same code,
A11, tol 1e-7) gives:

| h / ℓ | box / ℓ | E / E0 | degree | I∥ | I⊥ | ℓ_min (fm) | g (rigid) |
|---|---|---|---|---|---|---|---|
| 0.150 | 10 | 530.45 | 0.967 | 415.9 | 624.9 | 470.8 | 2.697 |
| 0.100 | 10 | 541.08 | 0.988 | 473.3 | 721.1 | 444.0 | 2.811 |
| 0.075 | 10 | 544.39 | 0.993 | 489.7 | 749.3 | 437.3 | 2.843 |
| 0.100 | 14 | 540.38 | 0.988 | 515.5 | 746.8 | 431.3 | 2.689 |
| 0.150 | 20 | 529.69 | 0.968 | 475.0 | 659.7 | 450.4 | 2.520 |

- **The box matters more than the grid.** Going from box 10 to 14 raises
  I∥ by 9 %. Going from box 10 to 20 at h = 0.15 raises it by 14 %. The
  ether's field falls off as a power law, so the ρ²-weighted inertia
  integrand has a slow tail, and the vacuum held at 0.9 × box cuts it off.
- **ℓ_min is robust.** It goes as I^(−½) and moves only 431–471 fm, so
  the 9×10⁶ gap in A5 is unaffected.
- **g (rigid) is not robust.** It drifts 2.52–2.84 with box and grid. Its
  sign of disagreement with 2.0023 is clear, but "2.81, independent of ℓ"
  should be quoted with a box-extrapolated value or a range.
- **E at h = 0.1 is 541.08, against ether_knot's extrapolated 547.5.**
  The README quotes 541.08 without saying it is a finite-h value.

### B5. The knot's tone is borrowed, not recomputed

`toneW = 1.36` (squashing mode, ω = 1.36 c/ℓ, Q ≈ 6–7) comes from
`examples/ether_knot`. It sets the tightest upper bound,
ℓ ≤ 4.8×10⁻⁵ fm. Two further assumptions sit inside that bound:

- reading the mode as an excited electron, m_e* = ħω/c², which a broad
  resonance (Q ≈ 6) makes approximate;
- taking CMS's contact-interaction limit with Λ = m_e*, which is
  model-dependent.

Even if this window is dropped, scattering (5.5×10⁻⁴ fm) and hydrogen
(0.047 fm) still leave a gap of ≥ 10⁴.

### B6. Two discretizations of the energy

`Axi.Energy` uses link differences, while `Inertia` and `Moments.E` use
central differences at nodes (as ether_knot's `Size`). The energy weights
of ⟨ρ²⟩ and ⟨r²⟩ therefore come from a slightly different energy than
the relaxed total. The difference is O(h²) and inside B4's spread, but it
is not stated.

### B7. The degree is not an integer on the grid

A11 relaxes to 0.988 (0.993 at h = 0.075). The A12 start is 1.934.
`TestStartDegree` allows 8 %. Q is taken as M·round(degree), which is
fine as long as the rounding is never close. It would not be if a
partially untied state were relaxed.

---

## C. Errors and inconsistencies in the code and write-up

### C1. The "rim at 0.65 c" at ℓ_min is computed with the wrong E0 (error)

`Knot.At` computes the spin energy, rim speed and tone with the
**static** E0 = m_e c²/E, ignoring the rotation. At ℓ_min the
self-consistent root of m_e c² = E·E0 + c_r(ħc/ℓ)²/E0 is E0 = m_e c²/(2E),
half the static value. The code's own doc comment on `MinSizeSpin` says
"half the mass is rotation", but the printed "spin energy 0.25 m_e c²"
contradicts it. With the self-consistent E0:

| Quantity at ℓ_min | Printed (static E0) | Self-consistent |
|---|---|---|
| spin energy / m_e c² | 0.25 | 0.50 |
| rim speed, ω⊥ from √(J(J+1)) = √¾ | 0.65 c | **1.30 c** |
| rim speed, ω⊥ from √(J(J+1) − J3²) = √½ (C2) | 0.53 c | **1.06 c** |
| rim speed, axial, ω∥ = ħJ3/I∥ | 0.57 c | **1.14 c** |

So at ℓ_min the knot's rim moves faster than light. Requiring a rim
speed below c raises ℓ_min only slightly, to ≈ 460 fm with the √¾
definition (≈ 445 fm with √½). **The verdict holds and gets slightly
stronger, but the statement "rim 0.65 c" is wrong.** It appears in three
places:

- `README.md` (step 2 paragraph and step 4 table);
- `CLAUDE.md` H44 row ("rim 0.65 c");
- `book/09-the-electron-test-bench.md` line 325–326 ("with the rim at
  0.65 c").

**Fix:** in `Knot.At`, use `E0At(ℓ)` when it exists and the static E0
only when it does not, and say which one is used. Then re-run and update
the three texts. The ℓ_min row should read "½ the mass is rotation, rim
≈ 1.1–1.3 c".

### C2. Ambiguous rim-speed definition

`RimSpeedC` uses ω⊥ = ħ√(J(J+1))/I⊥, which takes the total J = √¾ ħ as
if it all lay perpendicular. For J = ½, |J3| = ½, the perpendicular
angular momentum is √(J(J+1) − J3²) ħ = ħ/√2. Both are crude
(semiclassical) for a spin-½ state. The report should state the
definition, or give both the axial and the perpendicular rim speeds.

### C3. Run-time claims disagree

- The `main.go` header says "h = 0.1 (about five minutes)".
- The README says "h = 0.1 ℓ (~2 min)".
- `run_output.txt` shows 20 s (A11) + 41 s (A21) plus the A12 segments
  and control.

**Fix:** pick one measured figure.

### C4. "Independent of ℓ" for the rigid g

This is true analytically: g = E⟨ρ²⟩/I∥ contains no ℓ. But the value
depends on box and grid (B4). Rigid rotation is also unphysical in both
readings A (rim 1.72 c) and B, where g is quoted. The README should
qualify this.

### C5. Declared, uncited tolerance

`h1s2sTolMHz = 0.010` ("theory–experiment agreement taken as 10 kHz") is
marked as declared. It should cite the 1S–2S theory budget (the proton
radius term dominates the uncertainty), or state that 10 kHz is a
deliberately loose choice. The window it sets (0.047 fm) is not the
binding one.

### C6. The Coulomb window uses a thin-shell model at the core radius

The window uses αħc/(2·ring). A charge spread like the energy (as in the
g calculation) would give a different factor of order 1. This is
harmless (1.23 fm against 444 fm) but inconsistent within the example.

### C7. The test asserts a lattice artifact

`TestElectronKnot` requires the double loop's degree to fall below 1.5
within 2000 steps at h = 0.15. That checks the collapse only through its
lattice signature. It would fail if a later relaxer held the core off the
grid longer, for example with a smaller dt or a finer grid. Consider
asserting instead that the core radius falls monotonically below a few h.

### C8. Production and test grids differ by more than the tolerances quoted elsewhere

The test (h = 0.15, box 8) gives I∥ 385.2, ℓ_min 484 fm, g 2.788. The
run (h = 0.1, box 10) gives 473.3, 444 fm, 2.811. The test checks only
inequalities, so it passes, but the README's tests table does not say
the test runs on a different grid from the reported numbers.

---

## D. Bookkeeping

| Item | State |
|---|---|
| Git | the whole directory is untracked (`?? examples/claude_electron/`); not committed |
| `examples/README.md` row | present (line 61) |
| Book citation | present, `book/09` §9.8; carries the C1 error |
| `CLAUDE.md` | H44 register row and the layout entry present; register row carries the C1 error |
| SPDX headers | present on all `.go` files |
| Figures / PDF report | none, unlike sibling examples (`fdamador_neutrino`, `rodin_resonance`); not required |
| `internal/export` | no circuits in this example; nothing to export |
| Benchmarks (§8 rule 4) | none; this is an example, not a device or analysis, so the rule does not apply strictly |

---

## E. Recommended actions, in order

1. **Fix C1.** Use the self-consistent E0 in `Knot.At`, re-run, and
   correct "rim 0.65 c" in the README, `CLAUDE.md` H44 and `book/09`
   §9.8. Optionally fix C2 in the same change.
2. **Quote converged numbers (B4).** Add a box study (10, 14, 20 ℓ) and
   an h study (0.15, 0.1, 0.075) for A11. Quote E extrapolated to h → 0
   (against 547.5), and I∥, I⊥, ℓ_min and g with their spread.
3. **Make the A12 collapse less grid-dependent (B1).** Run A12 at
   h = 0.075 and 0.05 and report the core radius against steps.
4. **Run the double loop in 3-D (B2).** Use `ether_knot`'s 3-D rebuild
   as the starting point, to see whether it turns into A21.
5. **Smaller items.** Fix the run-time statements (C3), the g caveat
   (C4), the 1S–2S citation (C5) and the test assertion (C7). Then
   commit the example.
