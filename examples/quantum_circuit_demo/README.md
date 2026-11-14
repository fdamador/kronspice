# Quantum circuit demo: hydrogen as a circuit, and the quark test (H22)

This folder is one self-contained Python script with the two simulations from the
2026-11-12 discussion. It needs only `numpy` and `matplotlib` (no Go, no scipy). It
re-derives in a few hundred lines what `examples/kron_schrodinger`,
`examples/light_hydrogen` and `examples/quantumtricity` do in Go, so the method can
be read in one place.

```
python3 examples/quantum_circuit_demo/quantum_circuit_demo.py    # about 15 s
```

Evaluation report: [report.html](report.html).

Outputs: `run_output.txt`, `hydrogen_ladder.png`, `quark_test.png` and
`hydrogen_ladder.cir` (a plain SPICE netlist of the coarse ladder).

## Part A: hydrogen as a resonant ladder circuit (Kron 1945)

The radial Schrödinger equation, −½u″ + V_l(r)u = Eu, is built as an LC ladder in
its dual form, the same construction as `kron_schrodinger`:

| Circuit element | Quantum meaning |
|---|---|
| series capacitor 1/(2h²) between neighbouring nodes | kinetic energy |
| shunt capacitor V_l(r) + V0 at each node | potential energy (−1/r, plus l(l+1)/2r² for p states) |
| shunt inductor 1 H at each node | the energy |
| short circuit at r = 0 and at the outer edge | u = 0 boundary conditions |
| node voltage | wavefunction u(r) = rψ(r) |
| resonance at ω² = 1/(E + V0) | energy level E |

A 1 A AC source drives the node at r = 1 a0. The script sweeps frequency, finds the
poles of the driving-point reactance, and bisects each one to full precision.
A Sturm-sequence eigen-solve of the same network matrix is the independent check.
Three grids (h = 0.04, 0.02, 0.01 a0) are extrapolated to h → 0.

| State | Circuit (h → 0) | Exact | AC vs eigen check |
|---|---|---|---|
| 1s | −0.500000010 Ha (−13.60569 eV) | −0.5 | ≤ 8e-13 |
| 2s | −0.125000000 | −0.125 | ≤ 8e-13 |
| 3s | −0.055555556 | −1/18 | ≤ 8e-13 |
| 2p | −0.125000000 | −0.125 | ≤ 6e-13 |
| 3p | −0.055555556 | −1/18 | ≤ 6e-13 |

- **Lyman α** (2p → 1s, reduced mass applied): 121.5684 nm (measured 121.5670 nm).
- **Balmer α** (3p → 2s): 656.470 nm (measured 656.46 nm).
- **What's missing:** the remaining ~12 ppm is fine structure and the Lamb shift,
  which this network does not contain. `examples/hydrogen_poles` adds spin–orbit.
- **Wavefunctions:** the node voltages near each resonance lie on the exact 1s, 2s
  and 2p wavefunctions (`hydrogen_ladder.png`, right panel).

`hydrogen_ladder.cir` is the s-wave ladder at h = 0.25 a0 (79 nodes) in normalized
units (farads, henries). Its 1s resonance on that coarse grid is E = −0.492423 Ha,
at f = 0.083794 Hz. The `.ac` card brackets that frequency. It has not yet been run
through ngspice or kronspice's netlist parser.

## Part B: Quantumtricity (H22), circuit electricity against the quark model

**Step 1 (U(1), standard).** The circuit is a two-mesh RLC at 1 kHz: V1 1 V,
R1 100 Ω, L1 10 mH, C1 1 µF shared, L2 20 mH, R2 50 Ω.

| Gauge-theory identity | Circuit law | Residual |
|---|---|---|
| Curl∘Grad = 0 | A·C = 0 (incidence × loop matrix) | 0 |
| zero holonomy of a pure-gauge field | KVL, Cᵀv = 0 | 1e-16 |
| field = gradient of a potential | v = Aᵀφ | 8e-16 |
| Gauss's law | KCL, A·i = 0 | 0 |
| gauge invariance | φ → φ + const leaves v | 1e-15 |
| energy conservation | source power = resistor power | 2e-19 W |

**Step 2 (SU(3), the hypothesis).** A random 3×3 colour matrix goes on every branch
of the same graph. Both loops are based at the node they share, and the script
compares W₁W₂ with W₂W₁:

| Link variables | Composition defect ‖W₁W₂ − W₂W₁‖/3 |
|---|---|
| random SU(3), 2000 fields | mean 0.758 (min 0.223, max 0.969) |
| diagonal (Cartan) SU(3), the abelian part | 6e-17 |
| voltage, current, charge, flux, complex power as phases exp(ix) | 0 |

Gauge transformations leave the Wilson-loop traces unchanged (≤ 3e-15), which
checks the SU(3) construction itself.

**Verdict, the same as `examples/quantumtricity`.** Circuits are the abelian U(1)
limit of a lattice gauge theory, and the quark model is the non-abelian SU(3) case
on the same kind of graph. H22 stays **bounded**. What would move it is a circuit
quantity that carries an index along its branches and whose loop products do not
commute. The next construction, listed as open in `quantumtricity`, is the matter
sector: colour-triplet charges on the nodes, with hopping currents along the links.

## Sources

- G. Kron, "Electric circuit models of the Schrödinger equation," *Phys. Rev.* 67, 39 (1945).
- K. G. Wilson, *Phys. Rev. D* 10, 2445 (1974).
- J. Kogut, L. Susskind, *Phys. Rev. D* 11, 395 (1975).
- CODATA 2018 (hartree, proton/electron mass ratio).
