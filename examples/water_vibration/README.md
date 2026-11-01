# Water's bonds as a circuit

The vibration spectrum of a water molecule from kronspice's AC analysis:
atoms are mass capacitors, bonds are spring inductors, and Wilson's B matrix
is an E/F transformer between the bond frame and the atom frame. The first
molecule of the chemistry plan ([doc/chemistry-plan.md](../../doc/chemistry-plan.md));
the electronic side of the same molecule is in [kron_schrodinger](../kron_schrodinger)
and the DNA bases built the same way are in [dna_vibration](../dna_vibration).

```
go run ./examples/water_vibration [-netlists dir]
go test ./examples/water_vibration
```

## The construction

The mobility (force–current) analogy of Firestone 1933, the one Kron and
Branin use for network analogies (CLAUDE.md §3a.1):

| Mechanics | Circuit |
|---|---|
| velocity of a coordinate | node voltage |
| force | current |
| mass m | capacitor C = m to ground |
| spring of stiffness k | inductor L = 1/k |
| damper | resistor to ground |

The six in-plane Cartesian coordinates of O, H1 and H2 are nodes with the
atomic mass as a capacitor. The three bond coordinates (two O–H stretches
and the H–O–H bend) are springs. Each spring's rate of change is a linear
combination of the atomic velocities, s′ = Σ B_kj v_j (Wilson's B matrix),
built from a chain of voltage-controlled voltage sources; the spring force
flows in an inductor; and that force acts back on every atom through
current-controlled current sources. The E/F pair is an ideal multiport
transformer, Kron's connection tensor between the two frames. The force
constants form the inductance matrix L = F⁻¹, so the stretch–stretch
interaction f_rr is a mutual inductance (K element).

Circuit units: capacitance in daltons, inductance in Å/mdyn, so the
circuit's ω = 1/√(LC) maps to the molecule's through
√(100 N/m / 1 u) and is reported in cm⁻¹.

Geometry r_e = 0.958 Å, ∠HOH = 104.4776° (Hoy & Bunker 1979). The three
force constants (f_r, f_rr, f_a; no stretch–bend term, so three constants
for three H₂O lines) are fitted by Newton's method to H₂O's three measured
fundamentals. The same constants then predict D₂O, HDO and H₂¹⁸O with no
refit: isotopes change masses, not bonds (Born–Oppenheimer).

Force field fitted to the H₂O fundamentals: f_r 7.6712, f_rr −0.0942,
f_a 0.7030 mdyn/Å (to the harmonic frequencies: 8.4388, −0.1191, 0.7514).

## Results

The circuit's AC peaks (coarse sweep, fine sweep, parabola through the top
three points) equal the GF eigenvalues to 1e-6: the circuit *is* the
molecule's equations, not an approximation of them.

| Molecule | Mode | Circuit AC peak (cm⁻¹) | Measured (cm⁻¹) | Error |
|---|---|---|---|---|
| H₂O | bend | 1595.00 | 1595.0 | 0 (fitted) |
| H₂O | stretch | 3657.00 | 3657.0 | 0 (fitted) |
| H₂O | stretch | 3756.00 | 3756.0 | 0 (fitted) |
| D₂O | bend | 1162.51 | 1178.0 | −1.32 % |
| D₂O | stretch | 2647.26 | 2671.0 | −0.89 % |
| D₂O | stretch | 2752.05 | 2788.0 | −1.29 % |
| HDO | bend | 1395.34 | 1402.0 | −0.48 % |
| HDO | stretch | 2697.36 | 2727.0 | −1.09 % |
| HDO | stretch | 3708.98 | 3707.0 | +0.05 % |
| H₂¹⁸O | bend | 1589.15 | — | prediction only |
| H₂¹⁸O | stretch | 3647.51 | — | prediction only |
| H₂¹⁸O | stretch | 3740.70 | — | prediction only |

Measured values are the gas-phase fundamentals in NIST CCCBDB (Shimanouchi;
Huber & Herzberg). H₂¹⁸O is not in CCCBDB and its band origins have not
been checked against a primary source, so it is shown as a prediction.

The remaining error is anharmonicity. A measured fundamental sits below
the harmonic frequency by a shift that shrinks with heavier atoms
(ω·x ∝ 1/μ against ω ∝ 1/√μ). A harmonic force field fitted to H₂O's
fundamentals absorbs hydrogen's larger shift, so it comes out soft and
under-predicts the deuterium modes by about 1 %; the O–H stretch of HDO,
closest to H₂O, is within 0.1 %. The next step in the plan is Morse
(nonlinear) inductors.

The tests also check the three zero modes (two translations, one in-plane
rotation) and that the force field is physically sensible: stretch stiffer
than bend, a weak negative stretch–stretch interaction.

## Files

| File | What it holds |
|---|---|
| `main.go` | masses, geometry, Wilson's B, GF eigenproblem, the fit, the netlist writer, the AC sweep and peak finder |
| `main_test.go` | fit reproduces H₂O; AC peaks equal eigenvalues; isotopes within 1.5 % |
| `netlists/h2o.cir`, `d2o.cir`, `hdo.cir`, `h218o.cir` | each isotopologue as a SPICE netlist (E, F, K elements), runnable in kronspice or any SPICE |
| `export/` | the same netlists with KiCad schematics and results sheets (`TestExport`) |
| `run_output.txt` | the run above |

## Sources

- F. A. Firestone, "A new analogy between mechanical and electrical
  systems," *J. Acoust. Soc. Am.* 4, 249 (1933) — the mobility analogy.
- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations*,
  McGraw-Hill (1955), Ch. 4 — the B matrix and the GF method.
- A. R. Hoy, P. R. Bunker, *J. Mol. Spectrosc.* 74, 1 (1979) — equilibrium
  geometry of water (via NIST CCCBDB).
- NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB):
  measured fundamentals (T. Shimanouchi, *Tables of Molecular Vibrational
  Frequencies*; K. P. Huber, G. Herzberg, 1979) and harmonic frequencies.
- NIST "Atomic Weights and Isotopic Compositions" — atomic masses.
