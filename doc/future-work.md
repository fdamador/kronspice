# Future work requested

Items the user has asked to come back to, in the order they were raised.
Each is recorded here so it is not lost between milestones.

## The water molecule as a network (requested during M9b)

A circuit-analysis review of the water molecule (H₂O, and D₂O alongside
it): its vibrations, torsion and translation mechanics, treated as a
network so that harmonic resonances can be read off and compared.

What is already in place for it:

- **Energy domains** (`doc/domains.md`): mechanical networks in the
  mobility analogy (velocity across, force through), with masses as
  capacitors to ground and springs as inductors, side by side with
  electrical and acoustic networks and on the same ledger.
- **Kron's frames** (`kron/`): normal modes are a change of reference
  frame of the primitive network (Kron 1939 treats coupled oscillators and
  rotating machines the same way), and the mesh, node and torn forms all
  apply to a linear mechanical network.
- **H₂O and D₂O media** (`fluid/media.go`), with cited properties.

What it needs:

- A small-oscillation model of the molecule: three point masses (O, H, H or
  O, D, D) with bond-stretch and bend force constants from a cited force
  field; translation (3 modes) and rotation (3 modes, torsion about the
  principal axes) as zero-stiffness modes; the three vibrational modes
  (symmetric stretch, bend, antisymmetric stretch) as resonances.
- Validation against the measured fundamentals (gas phase, about
  3657, 1595 and 3756 cm⁻¹ for H₂O and 2671, 1178 and 2788 cm⁻¹ for D₂O;
  values to be taken from a cited spectroscopy table, e.g. NIST, before
  use), and the isotope shift H₂O → D₂O as a check of the mass terms.
- Rotational mechanics as a network: the inertia tensor in the rotational
  mobility analogy (angular velocity across, torque through), and how
  vibration–rotation coupling enters as mutual terms (chords).
- A harmonic-resonance review: the mode frequencies, their ratios and
  overtones, and where the liquid's acoustic modes (the flask of
  `examples/acoustic_drive`) sit relative to them — about 25 kHz against
  about 10¹⁴ Hz, some ten orders of magnitude apart, which the review
  should state plainly.

## Magnetic monopole evidence (pending, from the user)

The user will supply monopole search results later. Context recorded on
2026-10-25: the 2026 Nobel Prize in Physics went to Francis Halzen for
high-energy neutrino astronomy with the IceCube Neutrino Observatory
(nobelprize.org/prizes/physics/2026), and IceCube also searches for
relativistic magnetic monopoles. Until published evidence of a monopole
exists, `q_m = 0` in normal mode and the dual pair of §3f stays in the
hypothetical mode; a positive result would enter as a cited model with its
own known-value tests (CLAUDE.md §1 rule 3 applies in spirit), and a null
result as a recorded bound on H2/H3.
