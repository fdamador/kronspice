# RNA bases and base pairs as circuits

The four RNA bases, the two Watson–Crick pairs and the G·U wobble pair as
electric circuits, with the infrared light each one absorbs and releases
found by kronspice's AC analysis. The RNA twin of
[dna_vibration](../dna_vibration): every experiment there is repeated here
on the RNA force field, with uracil in place of thymine and RNA's own
non-Watson–Crick pair, G·U, added
([doc/chemistry-plan.md](../../doc/chemistry-plan.md)).

```
go run ./examples/rna_vibration [-out examples/rna_vibration] [-stacks=false] [-only-stacks] [-stack-keys AU2,GC2]
go test ./examples/rna_vibration
```

## The construction

**Molecules.** 9-methyladenine, 1-methyluracil, 9-methylguanine and
1-methylcytosine (the methyl stands in for the ribose, as in gas-phase
spectroscopy), the Watson–Crick pairs A·U and G·C, and the G·U wobble pair
(G N1–H···O2 U and U N3–H···O6 G). **Starting structure:** the solution
NMR structure of the RNA duplex r(GAGUGCUC)₂, PDB
[1EKA](https://www.rcsb.org/structure/1EKA), which has all hydrogens and
is near A-form (X. Chen, J. A. McDowell, R. Kierzek, T. R. Krugh,
D. H. Turner, *Biochemistry* 39, 8970, 2000). A·U is A2·U7, G·C is G1·C8
and G·U is G5·U4 of that duplex; the single bases are taken from the A·U
and G·C pairs. Everything is then minimized in the AMBER RNA force field
OL3 (Cornell et al. base terms with the bsc0 and OL3 χ corrections), as
converted by OpenMM. `prep/prep.py` writes `molecules.json` from OpenMM's
`amber14/RNA.OL3.xml` and the 1EKA file (URLs, the OpenMM commit and the
download date are in its docstring): atoms with masses, charges and
coordinates, and every force-field term the Go code evaluates. The methyl
caps are made exactly as in the DNA example (C1′ becomes the methyl
carbon, O4′ and C2′ become hydrogens), this example's choice, stated in
the output.

**Circuit** (`circuit.go`), **outputs** (`analysis.go`, `spectra.go`),
force field (`forcefield.go`) and modes (`modes.go`) are those of
dna_vibration: atoms are capacitors (mass); every AMBER bond and angle is
a spring branch (VCVS chain for Wilson's B row, an inductor for the spring
force, CCCS back to the atoms); everything else is the residual stiffness
R = K − Σ k·BᵀB through displacement integrators and one VCCS per non-zero
entry; the forces between the two bases of a pair are the group G_HB and,
in a stack, the forces between bases of different pairs the group G_ST.
Light enters as a current source at every atom. Absorption is the power
the field puts in during the AC sweep (0–3700 cm⁻¹, 2 cm⁻¹ steps,
averaged over x, y and z); release per mode comes by the Einstein A
coefficient and by Maxwell's radiated power of the vibrating dipole
(`em.LarmorPower`), which must agree; release spectra are given for one
quantum per mode and for thermal equilibrium at 300 K.

## Results

| Molecule | Atoms | Vibrations | Range (cm⁻¹) | d/λ | Thermal release at 300 K |
|---|---|---|---|---|---|
| A (9-methyladenine) | 18 | 48 | 48–3367 | 2.5e-4 | 2.52e-20 W |
| U (1-methyluracil) | 15 | 39 | 52–3313 | 2.0e-4 | 1.24e-20 W |
| G (9-methylguanine) | 19 | 51 | 31–3363 | 2.5e-4 | 2.97e-20 W |
| C (1-methylcytosine) | 16 | 42 | 45–3366 | 2.2e-4 | 2.52e-20 W |
| A·U pair | 33 | 93 | 12–3365 | 4.3e-4 | 3.60e-20 W |
| G·C pair | 35 | 99 | 18–3370 | 4.2e-4 | 5.16e-20 W |
| G·U wobble pair | 34 | 96 | 12–3363 | 4.3e-4 | 4.13e-20 W |

Every molecule is at a true minimum (RMS gradient ≤ 7.0e-10 kcal/mol/Å,
all 3N−6 vibrations real) and lumped by the λ/10 rule. The circuit's AC
solution and a direct solve of the same admittance give the same absorbed
power on and off resonance, and the two release routes agree for every
mode (`rna_test.go`).

**Against DNA** (both examples' run outputs): adenine, guanine and
cytosine are the same molecules with the same bonded base terms, while
the partial charges are each force field's own. Their N–H lines agree to
within half a wavenumber (A NH₂ antisymmetric stretch 3367.0 here, 3367.1
in DNA; G N1–H 3308.3 and 3308.7; C NH₂ antisymmetric 3365.6 and 3365.9);
the soft methyl and out-of-plane modes differ more (C lowest mode 45 cm⁻¹
here, 13 in DNA), from the charges and the different starting geometries.
Uracil (15 atoms, 39 vibrations) releases 1.24e-20 W at 300 K against
thymine's 1.28e-20 W; A·U releases 3.60e-20 W against A·T's 3.68e-20 W.

**Model against measurement** (`summary/measured.csv`, written by
`prep/measured.py`; harmonic model against gas-phase, matrix or
helium-droplet data). The adenine, guanine and cytosine values and sources
are reused from dna_vibration, because the molecules are the same; only
the model column is this run's. The uracil values are new.

| Band | Measured | Model | Difference |
|---|---|---|---|
| A NH₂ antisymmetric stretch | 3569 | 3367.0 | −5.7 % |
| A NH₂ symmetric stretch | 3452 | 3246.7 | −5.9 % |
| U N3–H stretch (1-methyluracil, Ar matrix) | 3428.6 | 3312.7 | −3.4 % |
| U N3–H stretch (uracil, gas) | 3436 | 3312.7 | −3.6 % |
| U C2=O stretch | 1737.5 | 1789.3 | +3.0 % |
| U C4=O stretch | 1700 | 1458.0 | −14.2 % |
| G NH₂ antisymmetric stretch | 3543 | 3363.3 | −5.1 % |
| G N1–H stretch | 3440 | 3308.3 | −3.8 % |
| G NH₂ symmetric stretch | 3430 | 3244.0 | −5.4 % |
| G C6=O stretch | 1752 | 1544.2 | −11.9 % |
| C NH₂ antisymmetric stretch | 3572.7 | 3365.6 | −5.8 % |
| C NH₂ symmetric stretch | 3451.7 | 3244.6 | −6.0 % |
| C C2=O stretch | 1740 | 1685.2 | −3.1 % |

The N–H stretches sit 3–6 % below the data, the AMBER harmonic value
against anharmonic measurements, as for DNA. The C=O stretches are spread
over several ring modes in this force field: the mode with the largest
share of the C4=O bond's energy holds only 19 % of it and sits at
1458.0 cm⁻¹, 14.2 % below the measured band, so the model has no clean
C4=O line. No gas-phase or droplet infrared spectrum of an A·U or G·U pair
was found (the closest is the gas-phase A·T pair of the DNA example, a
non-Watson–Crick isomer), and no measured low-frequency mode of an A-form
RNA duplex, so neither is compared. As for DNA, the model's H-bonded
stretches move up on pairing (U N3–H3 +26.3 in A·U, G N1–H1 +21.8 in G·C
and +30.4 in G·U, G N2–H21 +125.9 cm⁻¹ in G·C), where H-bonded N–H
stretches are measured to move down: a fixed-charge model cannot give that
red shift. The UV absorption maxima (UMP 262 nm, Cavaluzzi & Borer 2004;
A, G, C reused from the DNA example) and fluorescence maxima (A, G, C
reused; none found for UMP) are listed as measured only; they are
electronic and not in the vibration circuit.

**Pairing.** A·U has two hydrogen bonds (N6–H···O4 1.87 Å, N3–H···N1
1.91 Å H···A) and G·C three (1.92, 1.81, 1.83 Å), donor–acceptor 2.83–2.92 Å.
Modes with ≥ 25 % of their energy in the inter-base forces: A·U at 28, 45,
53, 56, 72 and 134 cm⁻¹; G·C at 61, 62, 115 and 152 cm⁻¹.

**The wobble pair.** G·U keeps exactly the two wobble hydrogen bonds after
minimization (`TestWobble`): G N1–H1···O2 U, H···A 1.80 Å, D···A 2.80 Å,
178°, and U N3–H3···O6 G, 1.88 Å, 2.89 Å, 170°. Its modes living in the
pairing are at 21, 43, 52, 55, 68 and 131 cm⁻¹, the strongest H-bond
mode (79 % in the H-bonds) at 131 cm⁻¹, against 134 cm⁻¹ for A·U (84 %) and
152 cm⁻¹ for G·C (80 %): two hydrogen bonds like A·U, and a pairing about
as stiff, well below G·C's three. Its thermal release, 4.13e-20 W, sits
between A·U and G·C. On pairing, G N1–H1 moves +30.4 cm⁻¹ (bonded to
U O2 instead of C N3) and U N3–H3 +3.1 cm⁻¹ (bonded to G O6 instead of
A N1).

**Stacks** (`stacks/summary*.csv`, `stacks/run_output.txt`): 2–4 stacked
A·U pairs (strand I reads AA, AAU, AAUU), G·C pairs (CG, CGC, CGCG),
poly(A)·poly(U) runs (AA, AAA, AAAA), the ApU step and the tandem wobble
pairs U·G / G·U (GU2). PDB 1EKA has no AA/UU, GG/CC or CG step, so every
built stack repeats one real Watson–Crick step of the duplex, G1·C8 → A2·U7
(5′-GA-3′/5′-UC-3′: 3.77 Å along the pair normal, 29.4° about it, in a
simplified Olson et al. 2001 frame); U·A and C·G are A·U and G·C turned
180° about the dyad. AU2 and polyA2 are therefore the same molecule (both
AA/UU) and give the same numbers. GU2 alone is cut from 1EKA as it is.

| Stack | Pairs | Atoms | Vibrations | Lowest (cm⁻¹) | Modes below 200 cm⁻¹ | of which rigid-pair |
|---|---|---|---|---|---|---|
| AU | 1 | 33 | 93 | 12.5 | 11 | 0 |
| AU2 | 2 | 66 | 192 | 8.5 | 28 | 4 |
| AU3 | 3 | 99 | 291 | 8.6 | 43 | 9 |
| AU4 | 4 | 132 | 390 | 6.4 | 59 | 12 |
| GC | 1 | 35 | 99 | 17.9 | 11 | 0 |
| GC2 | 2 | 70 | 204 | 16.7 | 28 | 6 |
| GC3 | 3 | 105 | 309 | 8.2 | 45 | 7 |
| GC4 | 4 | 140 | 414 | 10.0 | 60 | 13 |
| polyA3 | 3 | 99 | 291 | 5.9 | 44 | 9 |
| polyA4 | 4 | 132 | 390 | 4.8 | 60 | 14 |
| ApU | 2 | 66 | 192 | 13.0 | 28 | 5 |
| GU | 1 | 34 | 96 | 12.4 | 12 | 0 |
| GU2 | 2 | 68 | 198 | 10.6 | 28 | 3 |

As for DNA, stacking turns each single-pair line into a band of nearby
lines (`summary/stack_band_splitting.png`) and adds rigid-pair modes below
200 cm⁻¹; the lowest mode falls from 12.5 cm⁻¹ (one A·U) to 6.4 cm⁻¹ (four,
AAUU) and 4.8 cm⁻¹ (AAAA), against 12.7 → 3.4 cm⁻¹ for A·T in the DNA
example. With the backbone absent, the free stacks do not keep the A-form
step on minimization: AU2 rise 4.10 → 5.83 Å, GC2 twist 35.00° → 77.18°,
GU2 twist 50.72° → 93.09° (rise and twist here are the simplified
centroid and C1′–C1′ measures of `StackGeometry`), which the summaries
record. The 2-pair stacks' modal spectra agree with the direct solve and
with kronspice's AC analysis of the circuit to ≤ 1.2e-5, the size of the
series spring resistance the modal sum leaves out. No measured
low-frequency band of an A-form RNA duplex was found to compare with.

## Files

| File or folder | What it holds |
|---|---|
| `main.go` | the run: bases and pairs, netlists, spectra, pairing report, summary tables |
| `circuit.go`, `forcefield.go`, `modes.go`, `analysis.go`, `spectra.go`, `stacks.go` | see above (adapted from dna_vibration) |
| `rna_test.go`, `ff_test.go` | minima, stiffness split, circuit twin, release twin, passivity, hydrogen bonds (the wobble pair's included), modal twin, stack parts; analytic gradient against differences |
| `molecules.json` | the AMBER input (atoms and force-field terms), embedded in the binary |
| `netlists/` | each base and pair as a SPICE netlist (`a.cir` … `gu.cir`) |
| `spectra/` | absorption and release spectra per molecule (CSV) |
| `modes/` | every mode with its intensity, release rate, lifetime, pairing share and assignment |
| `summary/` | measured values with sources, bands, circuits, hydrogen bonds; atoms, bonds and mode shapes for the renderers; pictures |
| `stacks/` | stacked-pair netlists, modes, spectra, low modes, summaries and the stack run's output |
| `sound/` | the spectra as audio, 1 cm⁻¹ → 3 Hz, loudness by √(IR intensity) |
| `prep/` | `prep.py` (molecules.json from OpenMM's RNA.OL3.xml and PDB 1EKA), `measured.py` (summary/measured.csv), `render_bases.py`, `render_pairs.py`, `render_stacks.py` (pictures and sound files) |
| `viz/make_report.py`, `rna_vibration_report.pdf` | the PDF report |
| `run_output.txt` | the run above (bases and pairs) |

## Sources

- AMBER RNA force field OL3: W. D. Cornell et al., *JACS* 117, 5179 (1995); A. Pérez et al., *Biophys. J.* 92, 3817 (2007); M. Zgarbová, M. Otyepka, J. Šponer, A. Mládek, P. Banáš, T. E. Cheatham, P. Jurečka, *J. Chem. Theory Comput.* 7, 2886 (2011); as converted by OpenMM (`amber14/RNA.OL3.xml`).
- X. Chen, J. A. McDowell, R. Kierzek, T. R. Krugh, D. H. Turner, *Biochemistry* 39, 8970 (2000), doi:10.1021/bi992938e — PDB 1EKA, the starting structure.
- A. Yu. Ivanov, S. A. Krasnokutski, G. G. Sheina, Yu. P. Blagoi, *Biophysical Bulletin* (Visnyk Kharkiv Univ.) 2(11), 14 (2002) — 1-methyluracil in an argon matrix (read from the page images of the Russian original).
- P. Colarusso, K. Zhang, B. Guo, P. F. Bernath, *Chem. Phys. Lett.* 269, 39 (1997) — uracil in the gas phase (value read from Colarusso's thesis, Univ. Waterloo 1997, Table 4.1).
- Reused from dna_vibration: C. Plützer, M. Nir, M. S. de Vries, K. Kleinermanns, *PCCP* 3, 5466 (2001); M. J. Nowak et al., *PCCP* 21, 22857 (2019); A. Min et al., *Bull. Korean Chem. Soc.* 30, 3039 (2009); M. Nir et al., *PCCP* 4, 732 (2002); D. Onidas et al., *J. Phys. Chem. B* 106, 11367 (2002).
- M. J. Cavaluzzi, P. N. Borer, *Nucleic Acids Res.* 32, e13 (2004) — UV absorption maxima (UMP 262 nm).
- G. Varani, W. H. McClain, *EMBO Rep.* 1, 18 (2000) — the G·U wobble pair.
- W. Saenger, *Principles of Nucleic Acid Structure* (1984), Ch. 6 — hydrogen-bond geometry.
- W. K. Olson et al., *J. Mol. Biol.* 313, 229 (2001) — the base-pair reference frame, simplified here for placing stacks.
- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations* (1955) — the B matrix.
