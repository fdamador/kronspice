# Light hydrogen: the hydrogen atom from the light models (H32)

Hypothesis H32 (CLAUDE.md §1, user hypothesis 2026-11-05): the hydrogen
atom is built from the repository's light models. The electron is the
toroidal photon of Williamson & van der Mark ([photon_electron](../photon_electron),
H23). The proton is a circuit of light ([nucleon_circuits](../nucleon_circuits),
H31): one light ring, three quark rings, or Skyrme's knot. The knot of
[knotted_light](../knotted_light) and [ether_knot](../ether_knot) (H29) is
what would hold the electron's light together. Every prediction is set
beside measured hydrogen and beside the standard account: the Dirac
electron with QED, and the measured proton.

```
go run ./examples/light_hydrogen       # about 40 seconds
go run ./examples/light_hydrogen -svg examples/light_hydrogen/viz   # also the report figures and the 3-D scene
go test ./examples/light_hydrogen
KRONSPICE_EXPORT=1 go test ./examples/light_hydrogen -run TestExport   # export/: netlists, KiCad schematics, results sheets
```

## The constituents, as the three examples give them

| Part | Model | Numbers used |
|---|---|---|
| Electron (H23) | one Compton wavelength of light on a double loop | ring radius λ_C/4π = 193.1 fm; charge q = √(3ε0ħc)/2π = 0.910 e (eq. 5); g = 2(1 + α′/2π), α′ = (q/e)²α (eqs. 8, 12) |
| Electron, charge set to e | the same, with the paper's eq. 13 comparison | q = e, g = 2(1 + α/2π) |
| Proton, light ring (H31) | +e carried at c round radius ħ/(m_p c) | 0.2103 fm, 1 μN |
| Proton, quark rings (H31) | uud rings fitted to the proton's moment | 2.7928 μN, rms charge radius 0.678 fm |
| Proton, Skyrme knot (H31) | B = 1 hedgehog, computed in `nucleon_circuits/skyrme` | 1.880 μN; no charge radius computed |
| Confinement (H29) | the ether knot | energy 547.5√(κ2κ4), core radius 1.2ℓ, squashing tone ω = 1.416 c/ℓ |

## Results

**1. Scales.** The electron ring is 193 fm, the proton ring 0.21 fm, and
the Bohr radius 52,918 fm. Both rings are far smaller than the atom
(0.0036 and 4×10⁻⁶ of a0). At the atom's scale they are point charges and
point moments, so the atom is the Coulomb problem with charges q_e and
q_p.

**2. The bound state as Kron's network.** Kron's 1945 Schrödinger network
(*Phys. Rev.* 67, 39) for s states is a radial ladder:
- series inductors carry the kinetic energy;
- shunt inductors carry the potential −Z/r;
- unit capacitors carry the energy;
- the wave function u = rψ is the node voltage, with short circuits at r = 0 and at the edge.

| Z = q_e q_p / e² | E_1s (network, h → 0) | Exact −Z²/2 | E_2s | Exact −Z²/8 |
|---|---|---|---|---|
| 1 | −0.499 999 997 Ha | −0.5 | −0.125 000 000 | −0.125 |
| 0.910 318 (H23) | −0.414 339 490 Ha | −0.414 339 492 | −0.103 584 873 | −0.103 584 873 |

The same network is also built in kronspice. Like `kron_schrodinger`, it uses the dual form (capacitors and inductors swapped), because AC analysis needs a DC operating point. With 299 nodes, its highest resonance, found by bisection on the sign of the driving-point reactance, equals the network matrix's ground state to 8×10⁻¹⁵ Ha. The light models change only Z.

**3. Spectrum, neutrality and proton size.** The levels are Bohr's with
reduced mass, scaled by Z². The last column is the proton-size shift of
1S–2S, (2/3n³)(Zα)⁴m_r c²(r/ƛ_r)², relative to the measured radius.

| Electron + proton | Net charge | 1S binding | Lyman-α | 1S–2S size shift |
|---|---|---|---|---|
| toroidal photon + light ring | +0.0897 e | 11.269 eV | 146.70 nm | −0.926 MHz |
| toroidal photon + quark rings | +0.0897 e | 11.269 eV | 146.70 nm | −0.535 MHz |
| q set to e + light ring | 0 | 13.598 eV | 121.568 nm | −0.907 MHz |
| q set to e + quark rings | 0 | 13.598 eV | 121.568 nm | −0.338 MHz |
| measured electron + proton | 0 | 13.598 eV | 121.568 nm | 0 |
| **Measured hydrogen** | \|q\| < 10⁻²¹ e | 13.5984 eV | 121.5668 nm | 1S–2S known to 10 Hz |

The measured columns come from these sources:
- **Neutrality:** Bressi et al., *Phys. Rev. A* 83, 052101 (2011).
- **Ionization energy and Lyman-α:** NIST Atomic Spectra Database.
- **1S–2S:** Parthey et al., *Phys. Rev. Lett.* 107, 203001 (2011).

The remaining 23.0 GHz between the reduced-mass Bohr 1S–2S and the measured value is the Dirac and QED (Lamb) part, which no row models.

**4. The 21 cm line as a circuit.** The electron's spin, spread over the
1s cloud, is a magnetization M(r) = μ_e|ψ|² ẑ. A magnetization is a
current, J = ∇×M = (2/a)M sinθ φ̂, so the cloud is 32 coaxial current
loops round the proton ring. Each loop and the ring is a kronspice
inductor driven to its current, coupled to the ring by Maxwell's
coaxial-circle mutual inductance.

The interaction ΣM·I_p·I_k, taken from the solved currents, equals
Fermi's contact energy μ_p·(2/3)μ0 M(0) to 1.3×10⁻⁵, the remainder being
the ring's finite size. The ledger's stored energy equals ½ΣLI² + ΣMII to
4×10⁻¹⁶.

The circuit gives the aligned and anti-aligned (zz) energy ±U. Spin ½
(S·S rather than S_zS_z) makes the splitting 4U.

| Electron + proton | μ_e | μ_p | 21 cm line | / measured |
|---|---|---|---|---|
| toroidal photon + light ring | 0.9112 μB | 1 μN | 349.2 MHz | 0.246 |
| toroidal photon + quark rings | 0.9112 μB | 2.7928 μN | 975.3 MHz | 0.687 |
| toroidal photon + Skyrme knot | 0.9112 μB | 1.880 μN | 656.5 MHz | 0.462 |
| q set to e + light ring | 1.00116 μB | 1 μN | 508.6 MHz | 0.358 |
| q set to e + quark rings | 1.00116 μB | 2.7928 μN | 1420.51 MHz | 1.000 07 |
| q set to e + Skyrme knot | 1.00116 μB | 1.880 μN | 956.2 MHz | 0.673 |
| measured electron + proton | 1.00116 μB | 2.7928 μN | 1420.50 MHz | 1.000 07 |
| **Measured** | | | 1420.405 751 768 MHz | 1 |

The measured line is from Hellwig et al., *IEEE Trans. Instrum. Meas.* 19, 200 (1970). The 0.007 % left over is the Breit, recoil and proton-structure terms beyond Fermi's formula.

**5. What would hold the electron's light.** The ether knot of H29 is sized to the electron:
- **Scale:** its core radius 1.2ℓ is set to the ring's 193.1 fm, giving ℓ = 160.9 fm.
- **Energy:** 547.5√(κ2κ4) is set to m_e c², giving κ2 = 9.29×10⁻⁴ J/m and κ4 = 2.41×10⁻²⁹ J·m.
- **Prediction:** its squashing tone, ω = 1.416 c/ℓ, is then ħω = 1.737 MeV. That is an excited electron at 2.248 MeV.

Against measurement:
- No excited electron has been seen. Collider searches exclude excited electrons up to the TeV scale (PDG 2024, "Searches for excited leptons").
- g − 2 limits the electron's size to about 10⁻¹⁸ m, against the knot's 10⁻¹³ m.
- The knot carries no charge and no spin ½. Those still come from H23's assumption.

**6. The spectrum from the frequency response.** An AC sweep of the whole network shows each level as a pole of the driving-point impedance.
- **Pinning:** each pole is pinned by bisection, then two network spacings are extrapolated to zero.
- **Levels:** the poles equal the network matrix's eigenvalues to all ten printed digits, and reach −Z²/2n² to 3×10⁻⁶.
- **Lines, q = e:** their differences give the Lyman and Balmer lines within the 0.01 nm of the measured values.
- **Lines, q = 0.910 e:** every line comes out 1/Z² = 1.207 times longer.
- **1S–2S:** 13 ppm low. About 9 ppm of that is the fine structure and Lamb shift the Bohr levels leave out.

[REPORT.md](REPORT.md) has the full answer, with seven figures, to "does the spectral analysis match the frequency response?". The interactive 3-D view is [viz/atom3d.html](viz/atom3d.html). Open it in a browser after cloning; it loads three.js from cdnjs, so it needs a network connection, and GitHub shows its source rather than running it. A hosted copy is at https://claude.ai/artifact/Hb8rksK7WtBHxBkM31ESD4.

**7. Hydrogen under pressure.** The network's outer short circuit is a hard wall, and moving it in compresses the atom (Michels, de Boer & Bijl, *Physica* 4, 981, 1937):
- **Lines:** every level rises and every line shortens. Lyman α goes from 121.6 nm free to 26.8 nm at a 2 a0 wall, which takes 375 GPa, about the pressure at Earth's centre.
- **Degeneracy:** the free atom's 2s = 2p splits.
- **Exact checks:** at 2 a0 the ground state is exactly −1/8 Ha; at 1.835 a0 it reaches zero and the atom is ionized (Fig. 8).

This is Walter Russell's "high pressure, short waves" read as the confined atom, compared in [doc/russell-octave-report.md](../../doc/russell-octave-report.md).

**8. The two-lobed 2p_z state.** A centrifugal term l(l+1)/2r² in the shunts gives the p states:
- **Level:** E_2p = −1/8 Ha to 1e-9; its kronspice twin agrees with the matrix to 2×10⁻¹⁴.
- **Shape:** two lobes on the axis and a nodal plane between them (Fig. 9, and the interactive 3-D view's middle tab).

## What it means

**With the toroidal-photon charge, the light atom fails three measurements.**
- It is not neutral: it carries 0.090 e, where the measured limit is 10⁻²¹ e.
- Every line moves by Z² = 0.829: Lyman-α comes out at 146.7 nm.
- The 21 cm line falls to 0.69 of the proton factor.

**With the charge set to e, the levels are Bohr's**, and what is left
tests the proton:
- **Light-ring proton:** misses 21 cm by its moment (0.36 of the measured line) and 1S–2S by its size (0.91 MHz).
- **Quark-ring proton:** fitted to the moment, it gets 21 cm. It misses 1S–2S by 0.34 MHz because it is smaller than the measured proton (0.678 fm against 0.841). The measurement is good to 10 Hz.

**The electron's light is never held.** The knot of H29 that would hold it predicts an excited electron at 2.25 MeV, which has not been seen.

**Standard prediction:** the Dirac electron with QED and the measured proton. It gives all of these columns, including the 23 GHz and the 0.007 % that the Bohr and Fermi formulas leave out.

## Tests

| Test | Checks | Tolerance |
|---|---|---|
| `TestNetworkCoulomb` | network 1s and 2s against −Z²/2 and −Z²/8 for Z = 1 and 0.910 | 1e-7 |
| `TestNetworkTwin` | kronspice resonance against the network matrix | 1e-9 Ha |
| `TestToroidalElectron` | charge 0.910318 e, g from eq. 12, ring radius 193.08 fm | 1e-6 |
| `TestStandardHydrogen` | ionization, Lyman-α and the 21 cm line of the standard atom against measurement | 1e-4 |
| `TestHyperfineCircuit` | loop interaction = Fermi's contact energy; ledger = ½ΣLI² + ΣMII; 4U/h = Fermi's formula | 1e-4, 1e-12 |
| `TestLightAtom` | the toroidal-photon atom's net charge (0.09 e ≫ 10⁻²¹ e); levels scale as Z² | 1e-12 |
| `TestKnotExcitation` | knot energy = m_e c²; squashing tone 1.737 MeV | 0.002 MeV |
| `TestSpectrum` | network poles = matrix (1e-9); levels → −Z²/2n² (5e-6); lines within 0.01 nm of measured; light-atom lines 1/Z² longer (1e-5) | as stated |
| `TestConfined` | E_1s = −1/8 at a 2 a0 wall; zero crossing in (1.83, 1.84) a0; 2p below 2s under compression; free 2s = 2p = −1/8 | 1e-6 |
| `TestTwin2p` | the l = 1 network in kronspice = its matrix | 1e-9 |
| `TestExport` | the radial network and the hyperfine loops export, and their results sheets solve | — |

## Files

| File | Contents |
|---|---|
| `main.go` | constituents, atoms, Bohr levels, size shift, Fermi formula, the five parts |
| `network.go` | Kron's radial Schrödinger network: Sturm-count eigenvalues, Richardson extrapolation, the kronspice dual-form twin |
| `hyperfine.go` | the electron's spin cloud as coaxial loops round the proton ring; the knot sized to the electron |
| `spectrum.go` | the frequency response of the whole network: AC sweep, poles, extrapolation, Lyman and Balmer lines |
| `pressure.go` | the confined atom (levels against the wall radius, pressure on the wall) and the 2p_z cloud sampler |
| `svg.go`, `figures.go`, `figures2.go`, `atom3d.go` | the report's SVG figures (standard library only), the 3-D figure and the 3-D scene data |
| `REPORT.md`, `viz/` | the report and its figures; `viz/atom3d.html` is the interactive 3-D view, written with `viz/scene3d.json` by the `-svg` run |
| `atom3d.html.tmpl` | the 3-D page's template, embedded in the program |
| `main_test.go`, `export_test.go` | the tests above |
| `export/` | `kron_radial_network` and `hyperfine_loops`: `.cir`, `.kicad_sch`, results sheets |
| `run_output.txt` | the recorded run |

## Sources

- J. G. Williamson, M. B. van der Mark, "Is the electron a photon with
  toroidal topology?", *Ann. Fond. Louis de Broglie* 22, 133 (1997): the
  electron model (H23).
- `examples/nucleon_circuits` (H31): the light-ring and quark-ring
  protons; G. S. Adkins, C. R. Nappi, E. Witten, *Nucl. Phys. B* 228, 552
  (1983), the Skyrme proton.
- `examples/knotted_light`, `examples/ether_knot` (H29): A. F. Rañada,
  *Lett. Math. Phys.* 18, 97 (1989); L. Faddeev, A. J. Niemi, *Nature* 387,
  58 (1997); R. A. Battye, P. M. Sutcliffe, *Phys. Rev. Lett.* 81, 4798
  (1998).
- G. Kron, "Electric circuit models of the Schrödinger equation," *Phys.
  Rev.* 67, 39–43 (1945).
- E. Fermi, *Z. Phys.* 60, 320 (1930); H. A. Bethe, E. E. Salpeter,
  *Quantum Mechanics of One- and Two-Electron Atoms*, §22; M. I. Eides,
  H. Grotch, V. A. Shelyuto, *Phys. Rep.* 342, 63 (2001), finite-size and
  hyperfine corrections.
- CODATA 2022 constants; PDG 2024; NIST Atomic Spectra Database.
- H. Hellwig et al., *IEEE Trans. Instrum. Meas.* 19, 200 (1970);
  C. G. Parthey et al., *Phys. Rev. Lett.* 107, 203001 (2011); G. Bressi
  et al., *Phys. Rev. A* 83, 052101 (2011).
