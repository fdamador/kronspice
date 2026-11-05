# The fdamador neutrino: the toroidal field inverted (H28)

Hypothesis H28 (CLAUDE.md §1, user hypothesis 2026-11-04): take the
toroidal field of the photon electron ([photon_electron](../photon_electron),
H23) and turn it inside out, so that the **electric field is the ring** and
the **magnetic field circles round it**. The proposal is that this inverted
structure is the neutrino. This example analyzes the inversion with
Maxwell's equations, builds it as kronspice circuits, and puts it on a
neutrino test bench next to the measured values and the standard
prediction.

```
go run ./examples/fdamador_neutrino                                     # about 1 minute
go run ./examples/fdamador_neutrino -viz examples/fdamador_neutrino/viz # also the chart and sound data (about 2 minutes)
python3 examples/fdamador_neutrino/viz/render.py                       # resonance.png and the .wav files
go test ./examples/fdamador_neutrino                                   # -short skips the ring-line sweeps
KRONSPICE_EXPORT=1 go test ./examples/fdamador_neutrino -run TestExport # netlists, KiCad schematics, results sheets
```

## What the inversion is, from Maxwell alone

In the charged ring and the photon electron the ring carries a current and
B wraps round the tube. Here E runs round the ring (toroidal) and B wraps
round it (poloidal). Before any model:

1. **It has no charge.** A closed loop of E has ∇·E = 0, so by Gauss's law
   it encloses no charge. It is neutral, as a neutrino is, and so it can
   never be the electron.
2. **It cannot be static.** ∮E·dl round the ring is not zero, which
   Faraday's law allows only while the flux through the hole changes. The
   ring sustains itself as an oscillation: E round the ring is a ring of
   displacement current (Ampère–Maxwell), its B threads the hole, and the
   changing flux drives E. That is a resonance.
3. **Its magnetic moment oscillates.** The displacement-current ring is a
   magnetic dipole m(t) = I(t)·πa² whose time average is zero, and a
   time-varying magnetic dipole radiates.

Standard physics already has this exact pattern: the TE01δ mode of a
dielectric ring resonator, E round the ring and H looping through the hole
(Kajfez & Guillon, *Dielectric Resonators*, 1986), used in microwave
filters. There the dielectric is what holds the field in place.

## 1. The E-ring resonator in kronspice

The ring as a lumped series loop of 24 sections: the displacement path
along the ring as series capacitors (C = ε₀εr·πr²/2πa), the ring's
self-inductance (L = μ₀a(ln(8a/r) − 7/4), Maxwell, *Treatise* Art. 706), the
radiation resistance of a small magnetic-dipole loop
(R_rad = (η₀π/6)(ka)⁴, Balanis eq. 5-24) and the dielectric loss
(tan δ·ω₀L). With radiation alone the Q has a closed form:

Q_rad = 6Λ/(π(ka)³), with ka = √2·(a/r)/√(εrΛ) and Λ = ln(8a/r) − 7/4,

so **Q_rad ∝ εr^(3/2): the permittivity is the confinement.**

Bench ring: a = 10 mm, r = 2.5 mm, εr = 100, tan δ = 1×10⁻⁴.

| Quantity | Formula | kronspice |
|---|---|---|
| Resonance | 2.0606 GHz | 2.0606 GHz (AC) |
| Q | 40.52 (radiation alone 40.68) | 40.52 (AC bandwidth), 40.53 (ring-down) |
| Ring-down ledger | | kick 3.151×10⁻¹³ J = radiated + heat 3.144×10⁻¹³ J + still stored 6.7×10⁻¹⁶ J; time-stepping error 2.5×10⁻¹⁸ J |
| Net charge | 0 | 0 by construction; 9×10⁻¹⁵ of a section's charge to ground through the DC bleeders |
| Magnetic moment | oscillates | peak 1.65×10⁻⁶ A·m² per volt of kick; mean over the last 10 periods −0.007 of their peak |
| Electrical size | | diameter 0.137 λ (em.CheckQuasiStatic) |

| εr | f₀ | ka | Q_rad | Q (with loss) | kronspice AC |
|---|---|---|---|---|---|
| 30 | 3.762 GHz | 0.788 | 6.69 | 6.68 | 6.68 |
| 100 | 2.061 GHz | 0.432 | 40.7 | 40.5 | 40.5 |
| 300 | 1.190 GHz | 0.249 | 211 | 207 | 207 |
| 1000 | 0.652 GHz | 0.137 | 1,286 | 1,140 | 1,140 |

At εr = 30 the ring is no longer small against the wavelength (ka ≈ 0.8),
so that row is at the edge of the lumped model.

## 2. The neutrino test bench

In vacuum (εr = 1, a/r = 4) the ring is electrically large at its own
resonance, ka = 4.3, so the lumped model fails and gives Q ≈ 0.04: the
field is gone in under one period. That agrees with the field solve of
[photon_electron](../photon_electron) Part 3, where light on a closed loop
leaves within half a period. Nothing in vacuum holds a ring of field.

So the bench asks what holding it would take. Neutrinos from SN1987A
(about 10 MeV) crossed 168,000 light-years. With a declared, illustrative
mass of 0.1 eV (KATRIN bounds m < 0.45 eV) the Lorentz factor is 10⁸ and
the proper time 5.3×10⁴ s, so the ring must keep its energy for
ω₀τ₀ = 8×10¹⁸ radians.

| Quantity | E-ring model | Measured / standard |
|---|---|---|
| Charge | 0 (∇·E = 0) | consistent with 0 |
| Spin ½ | the crossed inverted ring keeps half-integer modes (Part 3) | ½ |
| Static magnetic moment | 0 (it oscillates) | < 2.8×10⁻¹¹ μB (Borexino 2017); standard ≈ 3×10⁻²⁰ μB·(m/0.1 eV) |
| Oscillating moment, one quantum ħω₀ stored | 0.073 μB at ω₀ = 1.5×10¹⁴ rad/s | no radiative decay observed |
| Lifetime in vacuum | ≈ 3×10⁻¹⁶ s | ≥ 5.3×10⁴ s proper time (SN1987A) |
| Confinement needed to survive | effective εr = 3.4×10¹³ at a = 1.5×10⁻¹² m | none known |
| Mass | not given by the model (declared 0.1 eV) | < 0.45 eV (KATRIN 2025) |
| Weak interaction | not in the model | the only way neutrinos are detected |

The inversion gets the two properties a field picture can get: no charge
and no static moment. It fails on stability by about 20 orders of
magnitude in lifetime, and it says nothing about mass or the weak
interaction. Status: *bounded*.

## 3. The inverted travelling wave as a circuit

The ring line of [photon_electron](../photon_electron) with E and B
exchanged is its circuit dual (Kron, Branin): series capacitors carry E
along the path and shunt inductors carry B round it. That is a high-pass,
backward-wave ladder: the longest wavelength has the highest frequency.
Differential modes follow ω = 1/(2√(C·L_g)·sin(πk/N)). Built in kronspice
with 100 sections, C = 1 nF, L_g = 1 µH (0.1 Ω loss), driven at one point:

| Mode | Plain ring, AC | Theory | Wavelengths | Crossed ring, AC | Theory | Wavelengths |
|---|---|---|---|---|---|---|
| 1 | 80.114 MHz | 80.115 MHz | 1 | 160.21 MHz | 160.21 MHz | ½ |
| 2 | 40.078 MHz | 40.077 MHz | 2 | 53.420 MHz | 53.421 MHz | 1½ |
| 3 | 26.740 MHz | 26.740 MHz | 3 | 32.073 MHz | 32.074 MHz | 2½ |
| 4 | 20.079 MHz | 20.078 MHz | 4 | 22.931 MHz | 22.932 MHz | 3½ |

Agreement is within 5×10⁻⁵ and both AC ledgers balance. Crossing the
conductors still gives half-integer modes, so the 720° (spin ½) property of
the photon electron survives the inversion.

## 4. Inflow: the zero-point field balances the outflow

The ring loses energy through its radiation resistance. By the
fluctuation–dissipation theorem (Callen & Welton 1951) that resistance is
also a noise source: the field the ring radiates into drives it back. At
the neutrino ring's 0.1 eV the thermal part is negligible even in the
2.7 K cosmic background (occupation e⁻⁴³⁰), so the inflow is the
zero-point part, S_V = ħωR_rad/π. It is booked on the magnetic-field
channel, and R_rad is a resistor whose power leaves on that same channel,
so inflow and outflow sit in one ledger column.

Scaled loop (L = 1 H, C = 1 F, R_rad = 1/Q, ħ → 1 J·s; the circuit is
linear, so energies scale with ħ), 4 realizations at the bench ring's
radiation Q of 40.7:

| Quantity | kronspice | Expected |
|---|---|---|
| Stored energy | 0.4941 ± 0.0067 ħω₀ | 0.4916 (½ħω₀ for the noise modes present) |
| Magnetic channel | zero-point in 391.77 J = radiated out 389.36 J + stored 2.15 J | in = out + stored; time-stepping 6.7×10⁻⁴ of the outflow |
| Energy fluctuation σ/⟨E⟩ | 0.97 | ≈ 1 for a random field |
| Mean magnetic moment | 1.9×10⁻⁴ of its rms | 0 |

**Does the Casimir effect add inflow?** No, but it changes the outflow. The
Casimir force comes from boundaries changing which zero-point modes fit
between them. That is a fixed energy shift, a one-shot release at most
(see [forward_casimir](../forward_casimir), H25), not a continuous supply.
The same boundaries also change how fast the ring radiates (the Purcell
factor F, R_rad → F·R_rad), and the zero-point noise changes with it, so
the held energy should not depend on F:

| Purcell factor F | Q | Stored (kronspice) | Expected | Lifetime |
|---|---|---|---|---|
| 4 (enhanced) | 10.2 | 0.4794 ± 0.0023 ħω₀ | 0.4817 | ×¼ |
| 1 (free space) | 40.7 | 0.4941 ± 0.0067 ħω₀ | 0.4916 | ×1 |
| 0.2 (suppressed) | 203 | 0.4843 ± 0.0318 ħω₀ | 0.4942 | ×5 |

So the zero-point field does balance the outflow: the ring never decays
to zero and the books close. But what it holds is half a quantum of
random noise, ½ħω₀ = 0.05 eV for the 0.1 eV ring, with fluctuations as
large as the mean. Every mode of every oscillator gets the same, so it
singles out no particle. A coherent ring lasting SN1987A's proper time
would need F ≈ 5×10⁻²¹; the classic measurement of inhibited emission
found about 20× (Hulet, Hilfer & Kleppner 1985), and a free neutrino has
no walls.

One lesson from building it: the first version read 8 % high at low Q.
The noise started at a non-zero voltage, so the t = 0 operating point
charged the series capacitors through the 1 GΩ DC bleeders, an offset that
relaxes over 4×10⁹ s and sat in the stored energy as a bias proportional
to R_rad. Ramping the noise on from zero during the discarded settling time
removed it (both points then within 0.5σ). The expected value computed
beside the result is what exposed it.

## Chart and sound

`viz/resonance.png`: (1) the E-ring's resonance at εr 30–1000, narrowing
as Q rises; (2) Q against εr, formula and kronspice, on the εr^1.5 line;
(3, 4) the plain and crossed inverted ring lines' spectra with the modes
labelled in wavelengths round the ring; (5) one zero-point run, the
stored energy fluctuating about ½ħω₀.

`viz/fdamador_neutrino.wav` (15 s, 44.1 kHz) and `viz/sound_*.wav`: every
real frequency is far above hearing (2.06 GHz, 20–160 MHz, 24 THz), so
each sound is slowed by one factor, keeping every frequency ratio, decay
and Q:

| Part | Source | Audible |
|---|---|---|
| 1. Ceramic ring, εr 100 (Q 40.5) | kronspice ring-down | 440 Hz, a 0.14 s tick |
| 2. Ceramic ring, εr 1000 (Q 1140) | kronspice ring-down | 440 Hz, a 3.4 s bell |
| 3. Ring held by the zero-point field | kronspice noise run | a hiss centred on 440 Hz, 2.8 s |
| 4. Plain inverted ring line | four highest AC modes, synthesized | 880, 440, 294, 221 Hz: the undertone series 1, ½, ⅓, ¼ |
| 5. Crossed inverted ring line | the same | 880, 293, 176, 126 Hz: the odd undertones 1, ⅓, ⅕, ⅐ |

Parts 1–3 are the simulated waveforms themselves; parts 4 and 5 are tones
at the measured mode frequencies (equal loudness, 1.2 s decay). The
crossing removes every even member of the series, which is the spin-½
property heard rather than seen.

`viz/harmonics.wav` (14 s) adds those modes one at a time, 1.6 s each: the
plain ring 880 → +440 → +294 → +221 Hz, then the crossed ring 880 → +293 →
+176 → +126 Hz. Parts 1–3 have no harmonics of their own: each is a single
linear resonance, one tone at 440 Hz.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestERingAC` | AC resonance and Q against 1/√(LC) and ω₀L/R; ledger | 2e-4, 2e-3; 1e-9 |
| `TestRingDown` | ring-down Q; kick = radiated + heat + stored; no net charge; mean moment ≈ 0 | 1 %; 1e-3; 1e-9; 0.02 |
| `TestQScaling` | Q_rad ∝ εr^1.5; closed form = ω₀L/R_rad; AC at εr = 300 | 1e-12; 2e-3 |
| `TestBench` | proper time; the required ring resonates at mc²/ħ with exactly the required Q; vacuum ring electrically large | 1e-9 |
| `TestDualRingModes` | plain and crossed CL-ring modes against the ladder dispersion; ledgers | 1e-4; 1e-9 |
| `TestZeroPointInflow` | stored energy = expected for the noise modes; magnetic channel in = out + stored; σ/⟨E⟩ ≈ 1; mean moment ≈ 0 | 4σ or 0.03; 2e-3; 0.7–1.3; 0.01 |
| `TestExport` | E-ring and both ring lines exported (torn-solution tolerance 1e-6 for the E-ring, from 1 GΩ bleeders beside 0.29 Ω) | |

## Files

| File | What it holds |
|---|---|
| `main.go` | `ERing` (L, C, R_rad, Q, the kronspice loop, AC and ring-down), `Bench` (the neutrino test), `DualRing` (inverted ring lines) and the report |
| `main_test.go` | the tests above |
| `zpf.go` | Part 4: `ZPFRing`, the scaled loop in its radiation resistance's zero-point noise, on the magnetic channel |
| `viz.go`, `viz/render.py` | `-viz` data; the chart `viz/resonance.png` and the sound files |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt` | the run above |

## Sources

- J. G. Williamson, M. B. van der Mark, "Is the electron a photon with toroidal topology?", *Ann. Fond. Louis de Broglie* 22, 133 (1997) — the field that is inverted here (H23).
- J. C. Maxwell, *Treatise on Electricity and Magnetism*, Vol. II, Art. 706 — self-inductance of a ring.
- C. A. Balanis, *Antenna Theory*, 3rd ed., Wiley 2005, eq. 5-24 — radiation resistance of a small loop.
- D. Kajfez, P. Guillon (eds.), *Dielectric Resonators*, Artech House 1986 — the TE01δ mode.
- M. Agostini et al. (Borexino), *Phys. Rev. D* 96, 091103 (2017) — neutrino magnetic moment < 2.8×10⁻¹¹ μB.
- M. Aker et al. (KATRIN), *Science* 388, 180 (2025) — m_ν < 0.45 eV.
- SN1987A neutrino burst: Kamiokande-II (Hirata et al., *PRL* 58, 1490, 1987) and IMB (Bionta et al., *PRL* 58, 1494, 1987); LMC distance about 51 kpc.
- H. B. Callen, T. A. Welton, "Irreversibility and generalized noise," *Phys. Rev.* 83, 34 (1951) — the zero-point noise of a resistance.
- E. M. Purcell, *Phys. Rev.* 69, 681 (1946) — boundaries change the emission rate; R. G. Hulet, E. S. Hilfer, D. Kleppner, "Inhibited spontaneous emission by a Rydberg atom," *PRL* 55, 2137 (1985).
- H. E. Puthoff, *Phys. Rev. D* 35, 3266 (1987) — the zero-point balance of [zpf_hydrogen](../zpf_hydrogen) (H24), whose noise construction Part 4 follows.
- B. W. Lee, R. E. Shrock, *Phys. Rev. D* 16, 1444 (1977); K. Fujikawa, R. E. Shrock, *PRL* 45, 963 (1980) — the standard Dirac-neutrino moment 3eG_F m_ν/(8√2π²).
