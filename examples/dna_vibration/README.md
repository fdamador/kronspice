# DNA bases and base pairs as circuits

The four DNA bases and the two Watson–Crick pairs as electric circuits,
with the infrared light each one absorbs and releases found by kronspice's
AC analysis. The chemistry branch's step after
[water_vibration](../water_vibration)
([doc/chemistry-plan.md](../../doc/chemistry-plan.md)): the same mobility
analogy, now on an all-atom force field, with the base pairing and the
stacking between pairs as named coupling groups of the circuit.

```
go run ./examples/dna_vibration [-out examples/dna_vibration] [-stacks=false] [-only-stacks] [-stack-keys AT2,GC2]
go test ./examples/dna_vibration
```

## The construction

**Molecules.** 9-methyladenine, 1-methylthymine, 9-methylguanine and
1-methylcytosine (the methyl stands in for the sugar, as in gas-phase
spectroscopy) and the pairs A·T and G·C, cut from the all-atom B-DNA
dodecamer d(CGCGAATTCGCG) and minimized in the AMBER DNA force field (OL15;
parm99/bsc0 base terms). `prep/prep.py` writes `molecules.json` from
OpenMM's `DNA.OL15.xml` and its `bdna.xyz` structure: atoms with masses,
charges and coordinates, and every force-field term the Go code evaluates
(bonds, angles, torsions, impropers, nonbonded pairs with 1-4 scaling).
The methyl caps are this example's choice, stated in the output.

**Circuit** (`circuit.go`). Atoms are capacitors (mass); every AMBER bond
and angle is a spring branch: a chain of VCVS forms the rate of change of
its internal coordinate (Wilson's B row), an inductor carries the spring
force and CCCS hand it back to each atom, Kron's connection tensor between
the bond frame and the atom frame. Everything else in the force field
(torsions, out-of-plane, electrostatics and van der Waals, the pre-stress
at the minimum) is the residual stiffness R = K − Σ k·BᵀB, applied through
displacement integrators and one VCCS per non-zero entry. In a pair, every
force between the two bases, the hydrogen bonds and the rest of their
electrostatics, is its own group G_HB; in a stack the forces between bases
of different pairs are a third group G_ST. The circuit's stiffness is K
exactly up to entries below 1e-9 of the largest. Light enters as a current
source at every atom (the field's force on the atom's charge).

**Outputs** (`analysis.go`, `spectra.go`). The absorption spectrum is the
power the field puts in during the AC sweep (0–3700 cm⁻¹, 2 cm⁻¹ steps,
averaged over x, y and z). Release per mode comes by two routes that must
agree: the Einstein A coefficient from the IR intensity, and Maxwell's
radiated power of the vibrating dipole (`em.LarmorPower`, itself checked
against the Poynting flux of the full dipole fields). Release spectra are
given for one quantum per mode and for thermal equilibrium at 300 K
(Kirchhoff). `modes.go` holds the normal-mode solve and unit conversions;
`forcefield.go` the AMBER energy and its analytic gradient and Hessian;
`stacks.go` the 2-, 3- and 4-pair stacks.

## Results

| Molecule | Atoms | Vibrations | Range (cm⁻¹) | d/λ | Thermal release at 300 K |
|---|---|---|---|---|---|
| A (9-methyladenine) | 18 | 48 | 49–3367 | 2.5e-4 | 2.58e-20 W |
| T (1-methylthymine) | 18 | 48 | 32–3314 | 2.1e-4 | 1.28e-20 W |
| G (9-methylguanine) | 19 | 51 | 30–3363 | 2.5e-4 | 2.90e-20 W |
| C (1-methylcytosine) | 16 | 42 | 13–3366 | 2.2e-4 | 2.65e-20 W |
| A·T pair | 36 | 102 | 13–3365 | 4.4e-4 | 3.68e-20 W |
| G·C pair | 35 | 99 | 18–3370 | 4.2e-4 | 5.20e-20 W |

Every molecule is at a true minimum (RMS gradient ≤ 4.3e-10 kcal/mol/Å,
all 3N−6 vibrations real) and lumped by the λ/10 rule. The circuit's AC
solution and a direct solve of the same admittance give the same absorbed
power on and off resonance, and the two release routes agree for every
mode (`dna_test.go`).

**Model against measurement** (`summary/measured.csv`; harmonic model
against gas-phase, matrix or helium-droplet data):

| Band | Measured | Model | Difference |
|---|---|---|---|
| A NH₂ antisymmetric stretch | 3569 | 3367.1 | −5.7 % |
| A NH₂ symmetric stretch | 3452 | 3246.8 | −5.9 % |
| T N3–H stretch | 3435 | 3314.1 | −3.5 % |
| T C4=O stretch | 1720 | 1745.7 | +1.5 % |
| G NH₂ antisymmetric stretch | 3543 | 3363.2 | −5.1 % |
| G N1–H stretch | 3440 | 3308.7 | −3.8 % |
| G C6=O stretch | 1752 | 1544.0 | −11.9 % |
| C NH₂ antisymmetric stretch | 3572.7 | 3365.9 | −5.8 % |
| C C2=O stretch | 1740 | 1685.0 | −3.2 % |
| A·T inter-base H-bond stretch | 94 (duplex in water) | 123.5 (one pair in vacuum) | +31 % |

The N–H stretches sit 3–6 % below the data, the AMBER harmonic value
against anharmonic measurements. The measured red shifts of the N–H
stretches on pairing (−39 to −186 cm⁻¹ in the gas-phase A·T isomer) are not
reproduced: with fixed charges the model's H-bonded stretches move up
(T N3–H +27.6, G N2–H21 +126.3 cm⁻¹). UV absorption and fluorescence maxima
are listed as measured only; they are electronic and not in the vibration
circuit.

**Pairing.** A·T has two hydrogen bonds (N6–H···O4 1.87 Å, N3–H···N1
1.90 Å H···A) and G·C three (1.92, 1.81, 1.82 Å), donor–acceptor 2.82–2.92 Å
near the crystal values. Modes with ≥ 25 % of their energy in the
inter-base forces: A·T at 24, 45, 54, 72 and 123 cm⁻¹; G·C at 62, 63, 116
and 154 cm⁻¹.

**Stacks** (`stacks/summary*.csv`): 2–4 stacked A·T and G·C pairs, poly(dA)
runs and ApT. Stacking turns each single-pair line into a band of nearby
lines and adds rigid-pair modes below 200 cm⁻¹ (5 of 30 for AT2, 14 of 66
for AT4, 15 of 60 for GC4); the lowest mode falls from 12.7 cm⁻¹ (one
A·T) to 3.4 cm⁻¹ (four). With the backbone absent, the rise and twist of
the free stacks move away from the duplex's on minimization (for example
AT2 twist 45.7° → 12.7°), which the summaries record.

## Files

| File or folder | What it holds |
|---|---|
| `main.go` | the run: bases and pairs, netlists, spectra, pairing report, summary tables |
| `circuit.go`, `forcefield.go`, `modes.go`, `analysis.go`, `spectra.go`, `stacks.go` | see above |
| `dna_test.go`, `ff_test.go` | minima, stiffness split, circuit twin, release twin, passivity, hydrogen bonds, modal twin, stack parts; analytic gradient against differences |
| `molecules.json` | the AMBER input (atoms and force-field terms), embedded in the binary |
| `netlists/` | each base and pair as a SPICE netlist (`a.cir` … `gc.cir`) |
| `spectra/` | absorption and release spectra per molecule (CSV) |
| `modes/` | every mode with its intensity, release rate, lifetime, pairing share and assignment |
| `summary/` | measured values with sources, bands, circuits, hydrogen bonds; atoms, bonds and mode shapes for the renderers |
| `stacks/` | stacked-pair netlists, modes, spectra, low modes and summaries |
| `sound/` | the spectra as audio, 1 cm⁻¹ → 3 Hz, loudness by √(IR intensity) |
| `prep/` | `prep.py` (molecules.json from OpenMM's files), `render_pairs.py` and `render_stacks.py` (pictures and the stack sound files) |
| `run_output.txt` | the run above (bases and pairs) |

## Sources

- AMBER DNA force field OL15 (W. D. Cornell et al., *JACS* 117, 5179, 1995; A. Pérez et al., *Biophys. J.* 92, 3817, 2007; M. Zgarbová et al., *JCTC* 11, 5723, 2015), as converted by OpenMM.
- C. Plützer, M. Nir, M. S. de Vries, K. Kleinermanns, *PCCP* 3, 5466 (2001); C. Plützer, I. Hünig, K. Kleinermanns, M. Nir, M. S. de Vries, *ChemPhysChem* 4, 838 (2003) — adenine, thymine and the gas-phase A·T pair.
- M. J. Nowak, I. Reva, A. J. Lopes Jesus, L. Lapinski, R. Fausto, *PCCP* 21, 22857 (2019) — 9-methylguanine in an argon matrix.
- A. Min, S. Lee, M. Y. Choi, R. E. Miller, *Bull. Korean Chem. Soc.* 30, 3039 (2009); M. Nir et al., *PCCP* 4, 732 (2002) — cytosine.
- M. González-Jiménez et al., *Nat. Commun.* 7, 11799 (2016) — inter-strand mode of DNA in water.
- M. J. Cavaluzzi, P. N. Borer, *Nucleic Acids Res.* 32, e13 (2004); D. Onidas et al., *J. Phys. Chem. B* 106, 11367 (2002) — UV absorption and fluorescence (measured only).
- W. Saenger, *Principles of Nucleic Acid Structure* (1984), Ch. 6 — hydrogen-bond geometry.
- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations* (1955) — the B matrix.
