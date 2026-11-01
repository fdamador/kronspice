# Kron's Schrödinger network: the electron of a chemical bond

The Schrödinger equation of a chemical bond solved as an electric network,
after Kron's *Phys. Rev.* 67, 39 (1945): series inductors carry the kinetic
energy, shunt inductors to ground the potential energy and shunt capacitors
the energy eigenvalue; the wave function ψ is the junction voltage to ground
and an eigenstate is a resonance of the network (CLAUDE.md §3a.1). The
molecule is H₂⁺, the simplest bond (two protons, one electron), with the H
atom as its check; both have exact solutions. The electronic side of the
chemistry plan ([doc/chemistry-plan.md](../../doc/chemistry-plan.md)); the
vibrational side is [water_vibration](../water_vibration).

```
go run ./examples/kron_schrodinger
go test ./examples/kron_schrodinger
```

## The construction

The ground state is axially symmetric (m = 0), so the network is laid out
on the (ρ, z) half-plane with the nuclei on the axis. In atomic units, the
equation multiplied by ρ so the network is reciprocal is

```
−½ [∂ρ(ρ ∂ρψ) + ρ ∂z²ψ] + ρ V ψ = E ρ ψ,    V = −Σ Z_a / |r − R_a|
```

On a grid of spacing h with junctions at ρ_i = (i + ½)h (no junction on
the axis, so it needs no boundary condition) this is A ψ = E B ψ with
B = diag(ρ_i). In Kron's network: series branches of admittance
ρ_face/(2h²) between neighbours, a shunt branch ρ_i·V_i at each junction
and a shunt capacitor ρ_i. The outer boundary is a short circuit (ψ = 0)
far from the nuclei, Kron's termination "by a short circuit at a point of
zero voltage".

Two ways of finding the eigenvalue:

- **Kron's driven network.** The network is driven at a trial energy σ by
  sources B·x and its response is the next x (inverse iteration with the
  factors computed once); the response converges to the resonance nearest
  σ and the Rayleigh quotient gives E. Energies at h, h/2 and h/4 are
  extrapolated to h → 0 (Richardson, with the order estimated from the
  three runs).
- **A kronspice LC circuit.** The same network written as a SPICE netlist
  in its dual form (capacitors and inductors exchanged, because Kron's
  all-inductor network has loops of shorts at DC): series capacitors
  between junctions, shunt capacitors ρ(V + V0) and shunt inductors 1/ρ,
  with V0 shifting the potential so every element is positive. Then every
  eigenstate is a resonance at ω² = 1/(E + V0), the ground state the
  highest, and kronspice's AC analysis finds it.

## Results

Finite-difference energies lie below the exact ones (a difference scheme
is not variational; the grid samples the Coulomb wells too deeply) and
converge with an order set by the cusps, about 1.45 on these grids.

**H atom**

| h (a₀) | Junctions | E (Ha) |
|---|---|---|
| 0.200 | 9,800 | −0.51087970 |
| 0.100 | 39,200 | −0.50382944 |
| 0.050 | 156,800 | −0.50123643 |
| h → 0 (order 1.44) | | −0.499728 = −13.5983 eV (exact −0.5 Ha = −13.6057 eV) |

**H₂⁺ at R = 2 a₀, electronic energy**

| h (a₀) | Junctions | E (Ha) |
|---|---|---|
| 0.200 | 7,800 | −1.11745397 |
| 0.100 | 31,200 | −1.10780517 |
| 0.050 | 124,800 | −1.10429632 |
| h → 0 (order 1.46) | | −1.102291 (exact −1.1026342) |

**Bond curve** E(R) + 1/R, each R extrapolated from h = 0.2 and 0.1, with a
least-squares cubic through the minimum:

| | Model | Exact |
|---|---|---|
| Bond length R_e | 1.998 a₀ = 1.057 Å | 1.997 a₀ = 1.057 Å |
| Bond energy D_e | 0.1024 Ha = 2.786 eV | 0.1026 Ha = 2.79 eV |

**Circuit twin** (coarse grid, h = 0.4 a₀, 1,950 junctions): the kronspice
AC resonance is −1.13848455 Ha and the driven-network eigenvalue on the
same grid is −1.13848455 Ha; the test requires agreement to 1e-7. The
exported results sheet also solves this 1,950-node circuit by automatic
tearing into 4 pieces joined by 179 links, equal to MNA.

Water has ten interacting electrons: the next step is a self-consistent
field (Hartree–Fock) on a 3-D network, with a Poisson network for the
repulsion, which is not done here.

## Files

| File | What it holds |
|---|---|
| `main.go` | the (ρ, z) grid, Kron's network, the driven-network eigensolver, Richardson extrapolation, the LC netlist and its AC search, the bond curve |
| `main_test.go` | convergence to the exact energies; circuit resonance equals the eigenvalue; cubic minimum |
| `netlists/h2plus_coarse.cir` | the coarse H₂⁺ network as an LC netlist |
| `export/` | the same netlist with a KiCad schematic and results sheet (MNA and torn, `TestExport`) |
| `run_output.txt` | the run above |

## Sources

- G. Kron, "Electric circuit models of the Schrödinger equation," *Phys.
  Rev.* 67, 39–43 (1945) (`reference/`).
- D. R. Bates, K. Ledsham, A. L. Stewart, *Phil. Trans. R. Soc. A* 246,
  215 (1953); J. M. Peek, *J. Chem. Phys.* 43, 3004 (1965) — exact H₂⁺
  energies.
- CODATA 2018 — the hartree and the bohr.
