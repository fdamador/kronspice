# Nucleon circuits: the proton and the neutron as circuits of light (H31)

Hypothesis H31 (CLAUDE.md §1, user hypothesis 2026-11-05): the proton and
the neutron are circuits of light, as the electron (H23) and the neutrino
(H28) are pictured in this repository. Each picture is checked against
what any model of the nucleons has to reproduce, with the standard
accounts, both as quark rings and as Skyrme's knot, beside it.

```
go run ./examples/nucleon_circuits      # under a second
go test ./examples/nucleon_circuits
KRONSPICE_EXPORT=1 go test ./examples/nucleon_circuits   # export/: netlists, KiCad schematics, results sheets
```

| | Proton | Neutron |
|---|---|---|
| Charge | +e | 0 |
| Spin | ½ | ½ |
| Mass | 938.272 MeV | 939.565 MeV (1.293 MeV heavier) |
| Magnetic moment | +2.7928 μN | −1.9130 μN |
| Size | charge radius 0.84075 fm | ⟨r²⟩ = −0.1155 fm² |
| Lifetime | stable | 878.4 s (n → p e ν̄, 0.782 MeV released) |

**The light ring.** Light of energy mc² going once round a ring of
circumference h/(mc) sits at radius R = ħ/(mc). A charge q carried round it
at c is a current qc/(2πR) on an area πR², so its moment is qcR/2 = qħ/(2m):
the Dirac moment, exactly (`TestLightRingDirac`). Each ring is a kronspice
circuit: an inductor (a thin ring of core radius R/2,
L = μ0R(ln 16 − 7/4)) carrying the ring's current from a current source,
coupled to its neighbours by Maxwell's coaxial-circle mutual inductance.
The ledger books every ½LI² and every M·I1·I2, matching the analytic sum to
1e-6 (`TestCircuitTwin`).

## 1. The proton as one light ring

- Radius ħ/(m_p c) = 0.210 fm, against the measured 0.841 fm (4.0×).
- Moment 1.000 μN, against the measured 2.793 μN (2.79×).
- A ring carrying +e at c with the measured moment needs radius 0.587 fm,
  and light on a ring that size carries 336 MeV, a third of the proton.

**One ring cannot be the whole proton.** It reproduces the moment of a
point particle; the measured moment is 2.79 times that.

## 2. The neutron as the light rings of its decay products

The neutron decays to a proton, an electron and an antineutrino. In the
light model these are three rings: the proton ring, the electron's photon
ring (H23) and the inverted neutrino ring (H28).

| Requirement | Result |
|---|---|
| Charge 0 | +e − e + 0 = 0: passes |
| Spin ½ | three spin-½ rings can combine to ½: passes. (Rutherford's 1920 neutron, a bare p + e pair, could not.) |
| Unstable | m_n − (m_p + m_e) = +0.782 MeV: the neutron holds more energy than its parts, so it can decay: passes |
| Moment −1.913 μN | each ring at its own size (electron ring 386 fm): −1835.5 μN, 959× too large: **fails** |
| Moment, refitted | an electron ring of radius 0.990 fm gives exactly −1.913 μN (kronspice: currents +13 kA and −7.7 kA, magnetic energy 0.130 MeV) |
| Energy of that ring | light on a 0.990 fm ring carries 199.4 MeV, against 1.293 MeV available: 154× too much: **fails** |
| Binding | the two rings' Coulomb attraction is −1.6 MeV: far too weak to hold it |

The antineutrino ring is neutral, with no current and no moment, so it
changes none of this. **The decay-products picture gets the charge, the
spin and the instability right, and fails on the moment and the energy by
two to three orders of magnitude.** This is why Rutherford's 1920 neutron
was abandoned after Chadwick's 1932 discovery, and the light rings do not
escape it.

## 3. Standard counterpart: three quark light rings

The constituent-quark picture, written in the same terms: each quark is a
light ring (u: +2/3 e, d: −1/3 e), and the nucleon's spin-½ state is a
weighted mix of spin arrangements (for the proton, u↑u↑d↓ with weight 2/3,
u↑u↓d↑ and u↓u↑d↑ with 1/6 each; the SU(6) wavefunction). Each arrangement
is a three-ring kronspice circuit, and its moment is the I·πR² of the
solved currents.

| | Quark rings | Measured |
|---|---|---|
| Quark ring energy (the one free number, from μ_p) | 335.96 MeV, radius 0.587 fm | — |
| Proton moment | +2.7928 μN (fitted) | +2.7928 |
| **Neutron moment** | **−1.8619 μN (predicted)** | −1.9130 (−2.7 %) |
| μn/μp | −2/3 | −0.685 |
| Three rings' energy | 1007.9 MeV | 938.3 / 939.6 (7.4 % high; binding not modeled) |
| m_n − m_p from electromagnetism | −1.05 MeV (rings one radius apart; −0.69 to −1.29 MeV from half to four radii apart) | lattice QED part −1.00 MeV |

The ring that fits the proton's moment in Part 1 is exactly one quark ring.
With that one number the neutron's moment follows to 2.7 %. Electromagnetism
alone makes the neutron about 1 MeV *lighter*, in line with the lattice's
QED part; the measured +1.293 MeV needs the d quark heavier than the u
(lattice QCD part +2.52 MeV; Borsanyi et al., *Science* 347, 1452, 2015).

## 4. The knot picture: Skyrme's model

In Skyrme's model the proton and the neutron are the same knot (winding
number 1), in two states of internal rotation. Adkins, Nappi & Witten
computed its properties, with two constants fitted to the nucleon and Δ
masses (*Nucl. Phys. B* 228, 552, 1983):

| | Knot (ANW) | Measured |
|---|---|---|
| Proton moment (μN) | +1.87 | +2.79 |
| Neutron moment (μN) | −1.31 | −1.91 |
| \|μp/μn\| | 1.43 | 1.46 |
| Isoscalar radius (fm) | 0.59 | 0.72 |
| g_A | 0.61 | 1.23 |

These are quoted, not computed here. The knot gives the neutron a negative
moment and the moment ratio to 2 % with no charged parts at all, but the
absolute moments and sizes come out about 30 % low. Our ether knot (H29,
`examples/ether_knot`) is Faddeev's two-component cousin of this model;
the nucleon needs Skyrme's four-component field. Computing it is a natural
next step.

## 5. Neutron decay

- **Energy budget:** 939.5654 − 938.2721 − 0.5110 = 0.7823 MeV, shared by
  the electron and the antineutrino.
- **Lifetime:** the light rings carry no weak coupling, so the light model
  gives no lifetime. The standard model gives
  τ = 4908.6 s / (V_ud²(1 + 3λ²)) = 880.5 s with PDG's V_ud = 0.97373 and
  λ = 1.2754 (Czarnecki, Marciano & Sirlin, *PRL* 120, 202002, 2018),
  against the measured 878.4 s.

## What it means

- The light ring captures one thing exactly: a ring of light carrying a
  charge has the Dirac moment qħ/2m.
- The proton is not one ring, and the neutron is not its decay products
  held together: that picture fails on moment and energy by factors of 150
  to 1,000.
- Three rings work: the quark-ring circuit predicts the neutron's moment
  to 2.7 % from one number, and gets the electromagnetic part of the mass
  difference in sign and size.
- The knot (Skyrme) gets the same pattern with no charged parts.

What neither the light model nor the circuits supply: the binding that
holds three rings together (in the standard account, gluons), and the weak
interaction that turns d into u in the decay.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestLightRingDirac` | ring of radius ħ/(m_p c) with +e at c: moment 1 μN from the circuit's current; radius 0.2103 fm | 1e-9; 1e-6 |
| `TestCircuitTwin` | three coupled rings, both current signs: ledger stored energy = ½ΣLI² + ΣMI1I2; currents as set; ledger | 1e-6; 1e-9 |
| `TestQuarkRings` | μp fitted, μn/μp = −2/3; electromagnetic m_n − m_p near −1 MeV | 1e-9; ±0.5 MeV |
| `TestNeutronFromDecayRings` | electron ring for the measured moment 0.990 fm; its energy > 100× the 1.293 MeV available | 1e-3 |
| `TestLifetime` | standard-model τ within 0.5 % of 878.4 s | 0.5 % |
| `TestExport` | the four ring circuits as netlists, KiCad schematics and results sheets | — |

## Sources

- Measured values: CODATA 2022 (masses, moments, r_p = 0.84075 fm); PDG 2024 (⟨r²⟩_n, τ_n = 878.4 s, V_ud, λ).
- J. G. Williamson, M. B. van der Mark, *Ann. Fond. Louis de Broglie* 22, 133 (1997) — the light ring (H23).
- E. Rutherford, *Proc. R. Soc. A* 97, 374 (1920); J. Chadwick, *Nature* 129, 312 (1932).
- G. S. Adkins, C. R. Nappi, E. Witten, "Static properties of nucleons in the Skyrme model," *Nucl. Phys. B* 228, 552 (1983); table as reproduced in [arXiv:0912.5142](https://arxiv.org/abs/0912.5142).
- Sz. Borsanyi et al., "Ab initio calculation of the neutron-proton mass difference," *Science* 347, 1452 (2015).
- A. Czarnecki, W. J. Marciano, A. Sirlin, *PRL* 120, 202002 (2018).
