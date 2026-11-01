# Chemistry in kronspice — plan

User request, 2026-10-27: simulate chemical bonds as circuits, starting with
water, building a molecule library up to a DNA strand, with the vibration
(infrared/Raman), electromagnetic and magnetic (spin, NMR) spectra coming out
of each model and checked against measured data. Kron's diakoptics tears
large molecules into pieces (atoms, groups, base pairs).

## Two circuit models, both done for the first molecules

| Model | What the circuit is | Spectrum it gives | Status |
|---|---|---|---|
| Bond vibrations (`examples/water_vibration`) | Mobility analogy (Firestone 1933; Branin): atom masses are capacitors, bond stiffnesses are inductors, Wilson's B matrix is a transformer (E/F pairs) between bond and atom frames | Infrared/Raman fundamentals as AC resonances | H₂O: force field fitted to the three measured fundamentals. D₂O and HDO predicted with the same bonds, within 1.3 % (anharmonicity). Netlists in `netlists/`. |
| Electronic structure (`examples/kron_schrodinger`) | Kron 1945: the Schrödinger equation as an LC network; ψ is node voltage, energy levels are resonances | Bond length, bond energy, electronic levels | H atom and H₂⁺ (one electron): R_e 1.998 a₀ against exact 1.997; D_e 2.786 eV against 2.79 eV. |
| DNA bases and pairs (`examples/dna_vibration`) | AMBER OL15 force field as a circuit: masses = capacitors, bond/angle springs = inductors via Wilson's B, the rest as a coupling network; AC sweep = IR absorption; release from Einstein A and em.LarmorPower | IR lines, release rates, H-bond (pairing) modes | N–H stretches 3–6 % below gas/matrix data; H-bond red shifts on pairing not reproduced (fixed charges); A·T H-bond stretch 123.5 cm⁻¹ against 94 cm⁻¹ in DNA in water; UV and fluorescence measured only (summary/measured.csv). |

## Next steps

1. **Water's electrons.** Ten interacting electrons need a self-consistent
   field (Hartree–Fock): one Kron network per orbital, a Poisson network for
   electron repulsion, iterated, on a 3-D grid (water is not axially
   symmetric). Validation targets:
   - geometry (0.958 Å, 104.48°);
   - dipole moment (1.85 D);
   - Hartree–Fock limit energy.
   Diakoptics tears the 3-D grid into blocks.
2. **Anharmonic bonds.** Replace the harmonic springs with nonlinear
   inductors (Morse potential). This gives fundamentals and overtones
   directly, and removes the ~1 % isotope error.
3. **Peak intensities.** Each mode's infrared strength comes from the dipole
   derivative: the charge carried by each atom's motion. This is a
   controlled-source output in the circuit, checked against measured
   band intensities.
4. **Molecule library.** Add the following in this order, each with its
   measured IR/Raman spectrum as the validation case:
   - adenine, thymine, guanine and cytosine;
   - the sugar and phosphate backbone pieces;
   - the A–T pair (two hydrogen bonds) and the G–C pair (three). G–C must
     come out as the stiffer link.
5. **DNA strand.** Base-pair blocks joined by backbone and stacking
   springs. Kron's tearing works per base pair, so a mutation (a base
   change) refactors one block only. Low-frequency collective modes of the
   strand are the check against measured terahertz/Raman spectra.
6. **Magnetic spectra.**
   - Nuclear spins as coupled resonators: Larmor frequency × (1 − shielding),
     with J-couplings as mutual coupling between resonators. This gives
     NMR-style spectra, checked against measured chemical shifts.
   - Electron spin (paired or unpaired) for radicals.
7. **Fields.** The electric field and potential of each molecule come from
   the electronic model (the potentials on the Kron grid). Any scalar field
   beyond the electric potential is not supported by measurement. It enters
   only as a hypothesis in the register (§1 of `CLAUDE.md`), with the
   standard prediction computed beside it.

## Rules

- Every molecule is validated against measured spectra (NIST CCCBDB, the
  Shimanouchi tables, NMR databases) before it enters the library.
- Numbers that cannot be traced to a primary source are left out, not
  guessed. For example, the H₂¹⁸O band origins are shown as a prediction
  only until a verified measurement is added.
- Keep the Go model and the `.cir` netlist of each molecule together, so
  every result can be checked in kronspice and in other SPICE simulators.
