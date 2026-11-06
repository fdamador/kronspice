# Hydrogen poles: Russell's "four pairs of poles" read as hydrogen (H34)

Hypothesis H34 (CLAUDE.md §1, user hypothesis 2026-11-05) reads Walter Russell's painting "Analysis of the Octave Bar Magnet" as hydrogen. Beneath it Russell wrote *"a bar magnet has four pairs of poles — not one pair"*. H34 takes the painting's two touching spheres as two hydrogen atoms. The bar reading of the same words is H33 ([octave_magnet](../octave_magnet)).

```
go run ./examples/hydrogen_poles                # about a minute
go run ./examples/hydrogen_poles -svg viz       # also the figures and animations
go test ./examples/hydrogen_poles
```

Reports: [hydrogen_atom_report.html](hydrogen_atom_report.html) (one atom), [h2_report.html](h2_report.html) (two atoms) and [hydrogen_spectrum_report.html](hydrogen_spectrum_report.html) (spectral resonances). All Russell studies: [doc/russell-octave-reports.html](../../doc/russell-octave-reports.html).

## Why hydrogen, and what "four pairs" means for an atom

- **Poles come in powers of two.** A dipole has one pair; an octupole (2³ poles) has four. So "four pairs, not one pair" names the octupole against the dipole.
- **No two-pair state.** A state of definite parity has no magnetic quadrupole. An atom's magnetic poles therefore come as one pair or four pairs, never two.
- **Four pairs need higher states.** A magnetic 2^k-pole needs 2j ≥ k. The ground state (1s, j = ½) can have only one pair; four pairs first become possible in 2p3⁄2.
- **On a sphere round the atom.** An axial octupole shows as four alternating zones of B_r (N S N S from pole to pole); a dipole shows two.

## One atom (run Parts 1–2)

**The model** (`atom.go`, `series.go`):
- **The electron.** Pauli spinors |n l j m⟩ built with Clebsch–Gordan coefficients from the hydrogen radial functions. Its field comes from two sources: the orbital current J = −2μB Σ m_l|ψ|²/ρ, and the spin magnetization M = −(g_s/2)μB ψ†σψ with g_s = 2.002319.
- **The field.** Both enter as one azimuthal current, J_orb + ∇×M, whose field is computed exactly from the Legendre-series Green's function (Jackson §5.5).
- **The proton.** Its spin dipole (1.52×10⁻³ μB) is added at the centre.

| Level | Dipole (μB) | Landé −g_j m_j | Octupole (μB·a0²) | Pole pattern round the atom |
|---|---|---|---|---|
| 1s½, m = ½ | −1.001160 | −1.001160 | 0 | 2 zones (one pair) at every radius |
| 2p½, m = ½ | −0.332947 | −0.332947 | 0 | 2 zones at every radius |
| 2p3⁄2, m = 3⁄2 | −2.001160 | −2.001160 | +18.0209 | 2 zones at every radius (r_c = 4.25 a0 falls inside the cloud) |
| 2p3⁄2, m = ½ | −0.667053 | −0.667053 | −54.0626 | **4 zones out to 15.6 a0**, then 2 |

**Checks:**
- **Normalization:** ∫\|ψ\|² = 1 to 10⁻⁸.
- **Dipole:** it equals the Landé value to 10⁻⁸.
- **Selection rules:** even-l moments are zero to 10⁻¹⁰ (parity), and the j = ½ octupoles are zero to 10⁻¹¹ (2j ≥ k).
- **Wigner–Eckart:** the j = 3⁄2 octupoles stand in the ratio −3.000000, as the Clebsch–Gordan coefficient m(5m² − 3j(j+1) + 1) predicts.
- **Field twin:** the series field equals a direct sum of 1.28 million current loops to 2.4×10⁻⁴.
- **Where the octupole comes from:** entirely from the spin magnetization. The l = 1 orbital current carries none.

**The bench readings:**
- **Spheres.** N and S zones of B_r on spheres round the atom (Figs. H2, H3). For an exterior field the octupole leads inside r_c = √(3m3/m1), or √(−2m3/m1) when the ratio is negative.
- **Line beside the axis.** B_z along a line beside the spin axis, the octave bar's own test (Fig. H4). At 5 a0 from the axis, 2p3⁄2 m = ½ gives 4 sign changes; the others give 2 or none.

## Two atoms, H₂ (run Part 3)

- **The four magnets.** The painting's two spheres, read as a molecule, hold four small magnets: two electron spins and two proton spins.
- **The model** (`h2.go`). The Heitler–London wavefunction (Heitler & London 1927) with the exact one- and two-centre integrals (Sugiura 1927 for the exchange integral).
  - **Check:** the singlet binds at R_e = 1.6425 a0 with D_e = 3.156 eV (0.1160 Eh), the classic Heitler–London result (measured: 1.4011 a0, 4.748 eV).
  - **The triplet** lies above two separate atoms at every R, so it never binds.

| Spin arrangement | Spin density | Moment (μB) | Bound? | Zones round the molecule |
|---|---|---|---|---|
| singlet electrons + para protons | 0 | 0 | yes | none — no field at all |
| singlet electrons + ortho protons | 0 | +0.00304 | yes | 4 zones from 0.3 to 1.6 a0 (two proton dipoles end to end), 2 beyond |
| triplet electrons + para protons | 2.0001 | −2.00232 | no | 2 zones at every radius |
| triplet electrons + ortho protons (all four aligned) | 2.0001 | −1.99928 | no | 2 zones, except 4 in a thin shell round the protons (0.6–0.8 a0) |

- **The electron field.** It comes from the triplet spin magnetization M_z = −(g_s/2)μB s_z, summed as current loops. Inside the cloud it is evaluated on cell corners, where the neighbouring loops' near fields cancel, and interpolated.
  - **Check:** h = 0.05 against 0.025 agrees to 3.4×10⁻⁴.
- **The proton spins.** Point dipoles at the nuclei, spins quantized along the bond.

## Spectral resonances found by kronspice (run Part 4)

The model (`spectrum.go`):
- **Atomic levels.** They are resonances of Kron's 1945 radial network (*Phys. Rev.* 67, 39), solved by kronspice's AC analysis in its dual form, as in [light_hydrogen](../light_hydrogen). Each level is where the driving-point reactance changes sign. Grids of h = 0.1 and 0.05 a0 are Richardson-extrapolated, and the reduced mass is applied.
- **Fine structure.** The p networks carry the spin–orbit term (g_s − 1)(α²/2)⟨L·S⟩/r³ (Bethe & Salpeter §12), so the j = ½ and j = 3⁄2 networks resonate apart.
- **H₂ vibration.** An LC circuit: the masses are capacitors and the Heitler–London spring is an inductor 1/k. Its single transfer pole is the vibration.
- **Twin:** every kronspice resonance equals the network matrix's Sturm eigenvalue to 5×10⁻¹⁴ Eh.

| Line | Found | Measured | ppm | How |
|---|---|---|---|---|
| Lyman α (2p → 1s) | 121.5689 nm | 121.5670 nm | −16 | kronspice |
| Balmer α (3d → 2p) | 656.4697 nm | 656.4610 nm | −13 | kronspice |
| Paschen α (4f → 3d) | 1875.627 nm | 1875.61 nm | −9 | kronspice |
| Brackett α (5g → 4f) | 4052.28 nm | 4052.27 nm | −3 | kronspice |
| fine structure 2p3⁄2 − 2p½ | 10 968.9 MHz | ≈ 10 969.1 MHz | −18 | kronspice |
| hyperfine 21 cm | 1420.486 MHz | 1420.406 MHz | +56 | model (ψ(0) of Part 1) |
| ESR of H, per tesla | 28 024.95 MHz | 28 024.45 MHz | +18 | model (dipole of Part 1) |
| H₂ / HD / D₂ vibration ω_e | 3805 / 3296 / 2691 cm⁻¹ | 4401 / 3813 / 3116 | −13.5 % | kronspice LC |

The remainders, and what is missing for each:
- **Optical lines:** the Dirac level shifts and the Lamb shift.
- **21 cm:** the Breit, recoil and proton-structure terms.
- **ESR:** the bound g-factor correction, −α²/3.
- **H₂:** the soft Heitler–London spring. Its isotope ratios, set by the masses alone, agree to 0.07 %: HD/H₂ 0.8661 against 0.8664, D₂/H₂ 0.7073 against 0.7079. The rotational constant at the measured bond length is 60.89 cm⁻¹ against 60.853.

## Result

**H34: bounded.**
- **One atom.** The standard calculation gives Russell's "four pairs, not one" a precise meaning. Ground-state hydrogen has one pair, and so does every j = ½ level. Four pairs (an octupole) first appear in 2p3⁄2, which with m = ±½ shows four N S N S zones round the atom out to 15.6 a0. The octupole comes from the electron's spin magnetization.
- **Two atoms.** The bond forms only with the electron spins opposite, so the electron magnets cancel exactly. That leaves the proton magnets, 658 times weaker. Ortho-H₂ shows a four-zone pattern near the protons. Para-H₂, the commonest form at low temperature, has no field at all.
- **What would move it.** An observation of an excited atom's magnetic octupole, or of H₂'s nuclear four-zone field. Neither has been measured, to our knowledge. Hydrogen's spin-½ proton has no octupole for the electron's to couple to.

## Tests

| Test | Checks | Tolerance |
|---|---|---|
| `TestNormAndDipole` | ∫\|ψ\|² = 1; dipole = −g_j m_j μB, all four levels | 1e-6 |
| `TestOctupoleSelection` | even-l moments zero; j = ½ octupole zero; j = 3⁄2 ratio −3 | 1e-8; 1e-6 |
| `TestSeriesTwin` | series field = direct loop sum | 1e-3 |
| `TestZones` | 1s 2 zones at 1–40 a0; 2p3⁄2 m = ½ 4 zones to 14 a0, 2 beyond 18 a0; r_c = 15.61 a0 | exact; 0.05 a0 |
| `TestE1` | E1(x) against Abramowitz & Stegun Table 5.1 | 1e-12 |
| `TestHeitlerLondon` | R_e = 1.6425 a0, D_e = 0.1160 Eh; triplet unbound at every R | 0.002 a0, 5e-4 Eh |
| `TestSpinDensity` | ∫ spin density 0 (singlet), 2 (triplet) | exact; 1e-3 |
| `TestSpectrum` | kronspice resonance = matrix; Lyman α to 30 ppm; fine structure, 21 cm to 100 ppm; ESR to 30 ppm; LC pole = √(k/μ); D₂/H₂ ratio | 1e-12 Eh; as listed; 1e-9; 2e-3 |
| `TestMoleculeField` | triplet far field → dipole of −g_s μB (corrections ∝ 1/r²); ortho protons 4 zones at 1 a0; para moment 0 | 10/r²; exact |

## Files

| File | Contents |
|---|---|
| `atom.go` | levels, spinors, orbital current, spin magnetization, Landé factor |
| `series.go` | Legendre-series field of an azimuthal current, moments, Gauss–Legendre |
| `loop.go` | current loops (exact elliptic field), point dipoles, the loop twin |
| `poles.go` | atom field with proton, zone and line scans, octupole radius |
| `h2.go` | Heitler–London energies, densities, spin field of H₂, corner grid |
| `spectrum.go`, `specfigs.go` | Part 4: Kron's radial network and the LC molecule in kronspice, the lines, Figs. S1–S3 |
| `main.go` | Parts 1–4 of the run |
| `figures.go`, `svg.go` | Figs. H1–H4, M1–M4 |
| `anim.go` | rasterizer and the two animations |
| `main_test.go` | the tests above |
| `hydrogen_atom_report.html`, `h2_report.html`, `hydrogen_spectrum_report.html` | the reports |
| `viz/` | `figH1_atom_fieldlines.svg`, `figH2_sphere_scans.svg`, `figH3_zones.svg`, `figH4_line_scan.svg`, `figM1_energies.svg`, `figM2_spin_density.svg`, `figM3_h2_fieldlines.svg`, `figM4_h2_zones.svg`, `figS1_fine_structure.svg`, `figS2_lines.svg`, `figS3_vibration.svg`, `atom_spheres.gif`, `h2_bond.gif` |
| `run_output.txt` | the recorded run |

## Sources

- W. Russell, "Analysis of the Octave Bar Magnet" (painting); see [doc/russell-octave-report.md](../../doc/russell-octave-report.md). The hypothesis source.
- J. D. Jackson, *Classical Electrodynamics*, 3rd ed., §5.5 (Green's-function expansion), §5.6.
- J. J. Sakurai, *Modern Quantum Mechanics*, §3.8 (coupling l and ½), §3.11 (Wigner–Eckart, projection theorem).
- W. Heitler, F. London, *Z. Phys.* 44, 455 (1927); Y. Sugiura, *Z. Phys.* 45, 484 (1927); J. C. Slater, *Quantum Theory of Molecules and Solids* I (1963), App. 6.
- G. Kron, *Phys. Rev.* 67, 39 (1945); H. A. Bethe, E. E. Salpeter, *Quantum Mechanics of One- and Two-Electron Atoms* (1957); NIST Atomic Spectra Database; H. Hellwig et al., *IEEE Trans. Instrum. Meas.* 19, 200 (1970); J. S. Tiedemann, H. G. Robinson, *Phys. Rev. Lett.* 39, 602 (1977); E. W. Hagley, F. M. Pipkin, *Phys. Rev. Lett.* 72, 1172 (1994); K. P. Huber, G. Herzberg, *Constants of Diatomic Molecules* (1979).
- M. Abramowitz, I. A. Stegun, *Handbook of Mathematical Functions*, Table 5.1.
- CODATA 2018 (μB, a0, g_s, μp); W. Kołos, L. Wolniewicz, *J. Chem. Phys.* 49, 404 (1968), the H₂ bond to spectroscopic accuracy.
