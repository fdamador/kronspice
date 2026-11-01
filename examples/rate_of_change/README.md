# Higher rates of change: in the material and in the field equations

Where do second, third and higher time derivatives enter
electromagnetism? The question came from a premise recorded in
[doc/vassilatos-references.md](../../doc/vassilatos-references.md): that
a fast pulse acts on space the way a fast slap acts on water, so that
above some rate the medium answers like a solid. In kronspice the answer
lives in two places, kept apart:

1. **In the material** (standard physics, `em/dispersion.go`). Maxwell's
   constitutive relation 𝔇 = (1/4π)K𝔈 (*Treatise* Art. 619, CLAUDE.md
   §3c) with a permittivity that depends on how fast the field changes.
   Water's molecules need about 8 ps to turn with a field, so water
   answers a slow field with ε ≈ 78.6 and a fast one with ε ≈ 5.4:
   electrically it stiffens when hit fast. The slap-on-water picture is
   real, and it is Debye relaxation.
2. **In the field equations** (hypothesis H19, `em/podolsky.go`).
   Bopp–Podolsky electrodynamics, (1 + a²□)□A = μ₀J, fourth order in
   time and space, with one new length a. Every result is printed beside
   Maxwell (a = 0), as the register requires.

```
go run ./examples/rate_of_change
go test ./examples/rate_of_change ./em
KRONSPICE_EXPORT=1 go test ./examples/rate_of_change -run TestExport   # the water cell's netlist and results sheet
```

## The construction

**Polarization equations.** A Debye term carries the first rate of change
of the polarization, τ·dP/dt + P = ε₀Δε·E, and a Lorentz term the second,
d²P/dt² + γ·dP/dt + ω₀²P = ε₀Δε·ω₀²E, with D = ε₀ε∞E + ΣP_k (Debye 1929;
Lorentz 1909; the auxiliary differential equations of Taflove & Hagness,
ch. 9). In the frequency domain they are Δε/(1 + jωτ) and
Δε·ω₀²/(ω₀² − ω² + jγω). Water at 25 °C is one Debye term: ε(0) = 78.6,
ε∞ = 5.4, τ = 8.3 ps, valid below about 100 GHz (the faster relaxation
near 1 ps and the infrared bands are left out). The time-domain solve
uses an exponential integrator that is exact for a field varying
linearly between grid points.

**The circuit twin.** Each term is an exact circuit through the
Kron/Branin network analogy (§3a.1): a plate cell of area A and gap d
filled with the material is C∞ = ε₀ε∞A/d in parallel with one branch per
term, a series R–C for Debye (C_k = ε₀ΔεA/d, R_k = τ/C_k) and a series
R–L–C for Lorentz. For the 1 cm², 1 mm cell used here: C∞ = 4.78 pF,
C_D = 64.8 pF, R_D = 0.128 Ω. Its admittance is exactly jω(ε₀A/d)·ε(ω),
so the netlist solved by kronspice AC is the third form of the same
model.

**Podolsky.** The Lagrangian adds (a²/2)∂_μF^μσ∂^νF_νσ to Maxwell's,
giving the fourth-order field equation above in Lorenz gauge. The code
holds the static point-charge potential and field, the finite
self-energy, the two free-wave modes and the Coulomb-test bound
(`PodolskyPotential`, `PodolskyField`, `PodolskySelfEnergy`,
`PodolskyModes`, `PodolskyBound`, `PodolskyLength`).

## 1. Water at 25 °C (single Debye: ε(0) 78.6, ε∞ 5.4, τ 8.3 ps)

| Frequency | ε′ | ε″ | Loss tangent |
|---|---|---|---|
| 1 MHz | 78.60 | 0.00 | 0.000 |
| 100 MHz | 78.60 | 0.38 | 0.005 |
| 1 GHz | 78.40 | 3.81 | 0.049 |
| 10 GHz | 62.95 | 30.01 | 0.477 |
| 19.2 GHz (1/2πτ) | 41.95 | 36.60 | 0.872 |
| 50 GHz | 14.79 | 24.47 | 1.655 |
| 100 GHz | 8.00 | 13.54 | 1.693 |

"Slapping" the water with a field ramped from 0 to 1 over a rise time T:

| Ramp time T | ε seen at the end of the ramp | After it settles (20 τ) |
|---|---|---|
| 1 µs | 78.6 | 78.6 |
| 1 ns | 78.0 | 78.6 |
| 100 ps | 72.5 | 78.6 |
| 10 ps | 36.1 | 78.6 |
| 1 ps | 9.6 | 78.6 |

Slow, the molecules follow; fast, they cannot turn in time and the water
answers like a stiff solid, ε → 5.4. The analytic value at the end of the
ramp is ε∞ + Δε·(1 − (τ/T)(1 − e^(−T/τ))), which the integrator meets to
10⁻⁶. Three forms of one model agree:

| Comparison | Result |
|---|---|
| Circuit twin (kronspice AC) against ε(ω) at 1 kHz to 100 GHz | worst difference 7.9×10⁻¹³ |
| Time-domain response against ε(ω) (`em/dispersion_test.go`) | twin test passes |
| Causality, Kramers–Kronig: ε′(0) − ε∞ = (2/π)∫ε″/ω dω | 73.2 recovered |

The exported results sheet tears the cell automatically into two pieces
joined by the Debye branch; the torn and full solves agree to 2.2×10⁻¹⁶
and the power is the same in both (Roth 1959).

## 2. Bopp–Podolsky (hypothetical, H19)

Point charge: φ = q/(4πε₀r)·(1 − e^(−r/a)), finite at the centre
(q/(4πε₀a)), below Coulomb by the fraction e^(−r/a). For an electron with
a = 1 fm (illustration):

| r / a | Maxwell φ | Podolsky φ | Fraction below Coulomb |
|---|---|---|---|
| 0 | infinite | 1.440 MV | 1 |
| 0.1 | 14.40 MV | 1.370 MV | 0.905 |
| 1 | 1.440 MV | 0.910 MV | 0.368 |
| 3 | 480.0 kV | 456.1 kV | 0.050 |
| 10 | 144.0 kV | 144.0 kV | 4.5×10⁻⁵ |
| 30 | 48.0 kV | 48.0 kV | 9×10⁻¹⁴ |

- The field energy of a point charge is finite, q²/(8πε₀a), and equals
  m_e c² when a = 1.409×10⁻¹⁵ m (half the classical electron radius),
  Podolsky's motivation.
- The free field splits into the ordinary Maxwell wave, ω = ck, and a
  massive mode ω = c√(k² + 1/a²) starting at c/a: 3×10²³ rad/s
  (4.8×10²² Hz) for a = 1 fm. The Podolsky mass is ħ/(ac), and
  `PodolskyLength` converts a mass in eV back to a.
- A Coulomb-law test good to 10⁻⁶ at 1 cm bounds only a < r/ln(1/ε) =
  0.72 mm. Measured Bhabha (e⁺e⁻) scattering bounds the Podolsky mass
  far more tightly (Bufalo, Pimentel & Soto 2014); that value is to be
  entered from the paper's results, not assumed (CLAUDE.md H19).

Standard prediction beside every line: a = 0, Coulomb's potential,
infinite self-energy, one wave mode. The effect is gone a few a away, so
only very short-range experiments can see it. Status: untested.

## Tests

| Test | Checks |
|---|---|
| `TestCircuitTwin` | the water cell and a Lorentz test material solved by kronspice AC equal jω(ε₀A/d)ε(ω) to 10⁻⁹ (§3a rule 4) |
| `TestSlap` | ramp responses equal the analytic end value to 10⁻⁶ and settle to ε(0) |
| `em/dispersion_test.go` | water's ε(ω), Kramers–Kronig, time–frequency twin; `TestPodolsky`: potential, field, self-energy and modes |
| `TestRuns` | the whole report runs |

## Files

| File | What it holds |
|---|---|
| `main.go` | water table, ramp test, circuit twin, Podolsky table |
| `main_test.go` | circuit twin and ramp tests |
| `export_test.go` | the water cell's netlist, KiCad schematic and results sheet in `export/` |
| `run_output.txt` | the run above |

## Sources

- J. C. Maxwell, *A Treatise on Electricity and Magnetism* (1873), Vol. II Art. 619: 𝔇 = (1/4π)K𝔈.
- P. Debye, *Polar Molecules* (Chemical Catalog Co., 1929); H. A. Lorentz, *The Theory of Electrons* (Teubner, 1909): the polarization equations.
- A. Taflove, S. C. Hagness, *Computational Electrodynamics: The Finite-Difference Time-Domain Method*, 3rd ed. (2005), ch. 9: auxiliary differential equations for dispersive media.
- D. C. Elton, "The origin of the Debye relaxation in liquid water and fitting the high frequency excess response," arXiv:1704.01667, §II.B and Table I: ε(0), ε∞, τ.
- J. Barthel et al., *Chem. Phys. Lett.* 165, 369 (1990), τ = 8.3 ps; R. Buchner, J. Barthel, J. Stauber, *Chem. Phys. Lett.* 306, 57 (1999), τ = 8.4 ps.
- F. Bopp, *Ann. Phys.* 430, 345 (1940); B. Podolsky, "A generalized electrodynamics," *Phys. Rev.* 62, 68 (1942).
- R. Bufalo, B. M. Pimentel, D. E. Soto, *Phys. Rev. D* 90, 085012 (2014), arXiv:1407.1476: the free modes (Sec. 2.2) and the Bhabha-scattering bound.
- CODATA 2018: e, m_e c², ħc.
- G. Vassilatos, *Secrets of Cold War Technology*, p. 27: the slap-on-water analogy (history only, not a physics source; see `doc/vassilatos-references.md`).
