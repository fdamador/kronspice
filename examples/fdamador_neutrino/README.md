# The fdamador neutrino: the toroidal field inverted (H28)

Hypothesis H28 (CLAUDE.md §1, user hypothesis 2026-11-04): take the
toroidal field of the photon electron ([photon_electron](../photon_electron),
H23) and turn it inside out, so that the **electric field is the ring** and
the **magnetic field circles round it**. The proposal is that this inverted
structure is the neutrino. This example analyzes the inversion with
Maxwell's equations, builds it as kronspice circuits, and puts it on a
neutrino test bench next to the measured values and the standard
prediction. Parts 5 and 6 implode and explode it (H30) and test the
user's reading of the three neutrino flavours as three wave types (H39).

```
go run ./examples/fdamador_neutrino                                     # Parts 1–4, 5a, 6: about 4 minutes
go run ./examples/fdamador_neutrino -field                              # also Part 5b, the field solver on grids r/2 and r/3: about 13 minutes (r/2 2 min, r/3 8 min)
go run ./examples/fdamador_neutrino -field -viz examples/fdamador_neutrino/viz # and the chart and sound data
python3 examples/fdamador_neutrino/viz/render.py                       # resonance.png and the .wav files
python3 examples/fdamador_neutrino/viz/render_part56.py                # implosion.png and flavours.png
go test ./examples/fdamador_neutrino                                   # about 3 minutes; -short skips the long runs
KRONSPICE_EXPORT=1 go test ./examples/fdamador_neutrino -run TestExport # netlists, KiCad schematics, results sheets
```

`run_output.txt` is the output of the `-field` run.

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

The ramp was added after the mismatch was seen, so it was then tested for
whether it steers the result. On identical noise (12 realizations each),
ramps of 2Q, 5Q and 8Q give 0.4815, 0.4815 and 0.4814 ħω₀ at Q = 10.2
(expected 0.4817) and 0.4900, 0.4900, 0.4899 at Q = 40.7 (expected
0.4916), while an abrupt start brings back exactly the biased 0.5230 and
0.5004. Raw rows in `ramp_check.txt`; the full account, with the equations
and the record of events, is in the reference document
[The 8 % zero-point bias](https://claude.ai/artifact/RvCYRSsSTSUjKchKxgeToa).

## 5. Implosion and explosion

H30 (implosion and explosion, user hypothesis 2026-11-05) applied to the
E-ring: a particle as a cycle of contraction and expansion whose books
balance only when the inward half is counted. In standard physics the
inward half of a radiating resonator is its ring-down run backwards: an
incoming wave that grows as e^{+αt} at the ring's frequency. Every joule of
it is declared and booked; nothing here supplies it.

### 5a. Circuit

The bench ring of Part 1 (a = 10 mm, r = 2.5 mm, εr = 100, tan δ = 1e-4;
2.0606 GHz, Q = 40.52) with its radiation resistance R_rad treated as a
port to the outside: a line of impedance R_rad on which an incoming wave
v_inc is a source 2·v_inc behind R_rad. The outgoing wave on that line is
v_ref = v_inc − R_rad·i: the reflection while the ring is driven, the
explosion after. For a one-port resonator a drive e^{st} is reflected with
amplitude 1 − 2R_rad/Z(s); it vanishes at the growth rate
(R_rad − R_d)/2L, so that drive is absorbed with no reflection at all
(coherent perfect absorption, the time-reversed laser: Chong, Ge, Cao &
Stone, *PRL* 105, 053901, 2010). The literal reversal of the lossy
ring-down, rate (R_rad + R_d)/2L, is predicted to reflect
(R_d/(R_rad + R_d))² = 1.6×10⁻⁵ of its power.

All drives carry the same incident energy (2.298×10⁻¹⁰ J) over 8 τ_E
(τ_E = Q/ω₀ = 3.129 ns) and are then switched off; the explosion is followed
as long again. "Absorbed" is the share of the incident energy not sent back
while the drive lasts; "peak" the most the ring held:

| Drive | Absorbed | Peak stored | Explosion | Ledger | Wave books |
|---|---|---|---|---|---|
| time-reversed ring-down, CPA rate (R_rad − R_d)/2L | 0.999656 | 0.99536 | 0.99123 | 5.1×10⁻⁶ | 1.0×10⁻⁵ |
| time-reversed ring-down, exact rate (R_rad + R_d)/2L | 0.999640 | 0.99538 | 0.99124 | 5.1×10⁻⁶ | 1.0×10⁻⁵ |
| growing at twice the rate | 0.889452 | 0.88742 | 0.88373 | 6.7×10⁻⁶ | 8.9×10⁻⁶ |
| growing at half the rate | 0.901874 | 0.89483 | 0.89112 | 3.6×10⁻⁶ | 9.2×10⁻⁶ |
| the ring-down itself (decaying, not reversed) | 0.028366 | 0.53939 | 0.02031 | 1.1×10⁻⁷ | 2.1×10⁻⁷ |
| constant burst at f₀, 2 τ_E | 0.799747 | 0.79681 | 0.79350 | 5.5×10⁻⁶ | 8.1×10⁻⁶ |
| constant burst at f₀, 8 τ_E | 0.486970 | 0.47658 | 0.47459 | 1.6×10⁻⁶ | 4.9×10⁻⁶ |

- The time-reversed drive's whole reflection, 3.4×10⁻⁴, is its switch-on
  (it starts at e⁻⁸ of its final power: e⁻⁸ = 3.4×10⁻⁴). The difference
  between the CPA rate and the exact reversal is 1.6×10⁻⁵, the predicted
  value.
- CPA run, kronspice ledger: source work 4.58343×10⁻¹⁰ J = R_rad
  4.56406×10⁻¹⁰ J + R_d 1.857×10⁻¹² J + stored at the end 7.77×10⁻¹⁴ J,
  time-stepping 5.1×10⁻⁶. Wave books: incident 2.29791×10⁻¹⁰ J = reflected
  7.914×10⁻¹⁴ J + explosion 2.27775×10⁻¹⁰ J + heat 1.857×10⁻¹² J + stored
  7.77×10⁻¹⁴ J.
- The decaying drive (the ring-down not reversed) fills the ring to 0.54
  early on, but the ring gives almost all of it back while the drive is
  still arriving. A steady drive at resonance is mostly reflected (R_rad ≫
  R_d: the ring is strongly over-coupled); the best constant burst holds
  0.80.

**Pulsing ring** (as knotted_light Part 2): an implosion of 3 τ_E of CPA
growth every 40 periods, six cycles, drives locked to the free ring left
from the cycle before, or in antiphase to it:

| Cycle (locked) | In | Out | Heat (R_d) | Stored start → end | Peak |
|---|---|---|---|---|---|
| 1 | 2.175×10⁻¹⁰ J | 2.082×10⁻¹⁰ J | 1.44×10⁻¹² J | 0 → 7.91×10⁻¹² J | 2.052×10⁻¹⁰ J |
| 2 | 2.175×10⁻¹⁰ | 2.151×10⁻¹⁰ | 1.70×10⁻¹² | 7.91 → 8.63×10⁻¹² | 2.241×10⁻¹⁰ |
| 3–5 | 2.175×10⁻¹⁰ | 2.158×10⁻¹⁰ | 1.72×10⁻¹² | 8.67 → 8.67×10⁻¹² | 2.250×10⁻¹⁰ |

In the steady cycle the outside supplies 2.175×10⁻¹⁰ J and takes back
2.158×10⁻¹⁰ J; the difference, 1.718×10⁻¹² J (0.79 %), is the dielectric
heat (1.716×10⁻¹² J). In antiphase the second cycle peaks lower (1.871×10⁻¹⁰
J, the residual ring cancels part of the drive) and the train settles to
the same steady cycle by cycle 4. Ledgers 6.6×10⁻⁶ and 6.8×10⁻⁶. The pulsing
ring is paid for, every cycle, by the declared incoming wave.

### 5b. Field

The same torus in the field solver (`em/fdtd.go`): εr = 100, lossless,
staircased on grids of r/2 and r/3 (dx = 1.25 and 0.833 mm), a closed
Huygens surface at ±2a and the absorbing layer from ±3a. Method:

1. A first run starts from an azimuthal displacement field in the tube
   (D_φ = D₀·cos²(πs/2r), B = 0: E round the ring, no charge). That start
   excites the TE01δ-like mode near 2.05 GHz and, weakly, n = 0 modes of
   the tube's cross-section near 7 GHz. The Fourier transform of every
   field component at the lumped f₀ over 15 periods is the low mode alone.
2. The ring-down proper starts from that filtered field: f₀ from the zero
   crossings of E_φ in the tube, Q from the decay of the energy inside the
   surface, and the ring's outgoing field on the surface as a phasor at
   f₀.
3. The implosion: the surface injects the time reverse of that outgoing
   field (Ê → conj Ê, B̂ → −conj B̂), growing as e^{+γt} (γ = ω₀/2Q) for
   6 τ_E, and is switched off at an E maximum in the ring. The explosion
   then crosses the surface into the absorber. Energy outside the surface
   or absorbed before the peak is what the ring sent back (or the declared
   wave failed to match).
4. For comparison a constant burst of the same incoming pattern and
   incident energy, 2 τ_E long with one-period raised-cosine edges.

| Quantity | Grid r/2 (389,017 nodes, 128 s) | Grid r/3 (912,673 nodes, 453 s) | Lumped model (Part 1, lossless) / 5a |
|---|---|---|---|
| f₀ (ring-down) | 2.101 GHz | 2.0571 GHz | 2.0606 GHz |
| Q (energy decay; field peaks) | 41.69; 41.56 | 38.9; 39.04 | Q_rad 40.68 |
| delivered by the surface / declared incident (2γ·W₀·∫g²dt) | 1.0253 | 0.9966 | — |
| held at the peak / delivered | 0.9976 | 0.9974 | 0.99536 of the incident (with R_d) |
| outside the surface + absorbed at the peak / delivered (reflection + mismatch) | 2.40×10⁻³ | 2.61×10⁻³ | 3.4×10⁻⁴ (switch-on) |
| explosion afterwards / delivered (followed 3 τ_E) | 0.9473 | 0.9466 | 0.99123 (followed 8 τ_E) |
| still inside at the end / delivered | 5.03×10⁻² | 5.08×10⁻² | — |
| E_φ share of the electric energy in the tube, at the peak | 0.9966 | 0.9971 | 1 by construction |
| B_φ share of the magnetic energy in the tube, T/4 later | 2.26×10⁻⁴ | 7.36×10⁻⁵ | 0 |
| B_z at the centre; at (2a, 0, 0), T/4 later | −2.38×10⁻⁹ T; 1.74×10⁻¹⁰ T | −2.45×10⁻⁹ T; 1.74×10⁻¹⁰ T | opposite signs: the flux returns outside |
| correlation of E at the peak with the free ring-down's E | 0.9991 | 0.9996 | — |
| field ledger residual / incident (implosion run; ring-down) | 3.3×10⁻¹³; 5.1×10⁻¹⁴ | 4.8×10⁻¹³; 1.3×10⁻¹² | 5.1×10⁻⁶ (time-stepping) |
| constant burst, same pattern and incident energy: peak inside; out when it ends | 0.7923; 0.2379 | 0.7649; 0.2389 | 0.79681 (2 τ_E burst) |

Between the grids f₀ moves by 2 % and Q by 7 %; the lumped model of Part 1
lies between them, so the field solver confirms it to within the
discretization at these grids (the tube is only 2 and 3 cells in radius; a
finer grid r/4, estimated at about 40 minutes, was not run). The absorbed
fraction does not depend on the grid: the time-reversed outgoing mode is
taken in, with 2.4–2.6×10⁻³ of the delivered energy outside the surface at
the peak, against 0.24 for the constant burst, which reaches 0.76–0.79 —
the same contrast as the circuit (0.99536 against 0.79681). The rebuilt
field is the TE01δ pattern: E round the ring (0.997 of the electric
energy), B threading the hole and closing round the tube a quarter period
later, the same field as the free ring-down (correlation 0.9996). Its
energy then leaves as the ordinary ring-down. The ±2.5 % gap between the
declared incident energy and what the surface delivered is the error of
W₀ and γ from the ring-down fit, used to normalise the declared wave; the
ratios above are taken against the delivered energy, which the ledger
measures directly. Approximations: staircased torus, lossless dielectric,
graded absorber at 3a (k·r = 1.3, where the near field is weak but not zero),
single-mode phasor (the 7 GHz cross-section modes filtered out).

Standard prediction, beside both parts: the implosion is the explosion
reversed; it needs a declared incoming wave carrying all its energy, is
absorbed without reflection only when it is the time reverse of the
ring-down, and is followed by the ordinary ring-down. That is what both
parts find. Nothing in the ring supplies the incoming wave: for a neutrino
(Part 2) it would have to arrive from outside, coherent, every period.

## 6. Flavours as wave types (H39)

The user's hypothesis (2026-11-07): "for the flavors of the neutrino, it
seems to me to be three different types of standing wave oscillations on
the wave: constructive, destructive and standing wave. Check for that and
see what the implosion/explosion would do to this coherent unit."

**Where the types live.** The uniform E-ring (TE01δ, azimuthal order 0) is
one state, with no travelling or standing form. The types live on the
ring's azimuthal harmonics n ≥ 1, i.e. on the inverted ring lines of Part
3 (here with Part 3's element values, C = 1 nF and L_g = 1 µH, on 24
sections, 1 mΩ shunt loss and 1 pF stray capacitance to ground, which lets
the state be imposed exactly at t = 0; it lowers ω by 1.5 % and is in every
formula). Harmonic n is a degenerate pair, z_j = u_j + iωL_g·g_j =
A·e^{−iκj} + B·e^{iκj}: A runs one way (CW), B the other (CCW).

**One fixed definition**, relative to a reference point c (where the state
is made or detected):

| Type | Definition | Meaning |
|---|---|---|
| C, constructive | travelling, B = 0 | every point reaches the full amplitude; nothing cancels anywhere |
| S, standing | \|A\| = \|B\|, antinode of u at c | fixed nodes and antinodes |
| D, destructive | \|A\| = \|B\|, node of u at c | the two waves cancel at c |

A state that cancels everywhere does not exist: Σ_j|z_j|² = N(|A|² + |B|²)
exactly (the two waves are orthogonal; the smallest singular value is √N =
4.8990), so only zero cancels everywhere and "destructive" can only mean
cancelling at a chosen point. On a perfect ring S and D are the same state
turned by a quarter wavelength, and C = (S + iD)/√2: **the three types
span two dimensions, not three.**

### Test 1: persistence

Each type imposed on the perfect ring (harmonic 1, c = ½) and followed for
100 periods in a kronspice transient:

| Ring | Type | Own share, min | Other share, max | Harmonic purity, min | Turns per period | Energy at end | Ledger |
|---|---|---|---|---|---|---|---|
| plain | C | 1.0000000000 | 1.73×10⁻¹¹ | 1.0000000000 | 0.99909 | 0.99469 | 9.5×10⁻⁵ |
| plain | S | 1.0000000000 | 7.0×10⁻²⁵ | 1.0000000000 | 0 | 0.99469 | 9.5×10⁻⁵ |
| plain | D | 1.0000000000 | 2.5×10⁻²⁵ | 1.0000000000 | 0 | 0.99469 | 9.5×10⁻⁵ |
| crossed | C | 1.0000000000 | 4.5×10⁻¹² | 1.0000000000 | 0.99909 | 0.99723 | 1.0×10⁻⁴ |
| crossed | S | 1.0000000000 | 1.2×10⁻²³ | 1.0000000000 | 0 | 0.99723 | 1.0×10⁻⁴ |
| crossed | D | 1.0000000000 | 1.8×10⁻²³ | 1.0000000000 | 0 | 0.99723 | 1.0×10⁻⁴ |

All three are stationary: C turns once a period (0.99909: trapezoidal
frequency warping at 60 steps a period), S and D keep their nodes, and the
energy falls only by the shunt loss. Standard prediction: every
superposition of a degenerate pair is stationary. None decays or
converts.

### Test 2: flavour oscillation

One section's capacitor made C(1 + δ) on both wires lifts the degeneracy:
the eigenmodes become the standing waves with node and antinode at the
defect, split by Δf ≈ δ·f/(N(1 + δ))·4C sin²(κ/2)/(4C sin²(κ/2) + C_g)
(f = 19.139 MHz; the last factor is the stray capacitance's share). Started as CW, the
wave converts to CCW and back, P = sin²2θ·sin²(πΔf·t):

| Ring | δ | Δf exact (eigenmodes) | Δf first order | Δf kronspice | sin²2θ | θ | Largest residual | Ledger |
|---|---|---|---|---|---|---|---|---|
| plain | 0.05 | 37,428 Hz | 37,426 Hz | 37,327 Hz | 1.0000 | 45.00° | 5.7×10⁻⁴ | 2.0×10⁻⁵ |
| plain | 0.10 | 71,454 | 71,449 | 71,266 | 1.0000 | 45.00° | 1.1×10⁻³ | 6.2×10⁻⁵ |
| plain | 0.20 | 131,000 | 130,990 | 130,660 | 1.0001 | 45.00° | 2.0×10⁻³ | 8.3×10⁻⁵ |
| plain | 0.40 | 224,530 | 224,560 | 223,990 | 1.0001 | 45.00° | 3.6×10⁻³ | 9.3×10⁻⁵ |
| crossed | 0.20 | 245,510 | 245,370 | 244,860 | 1.0000 | 45.00° | 1.8×10⁻³ | 9.5×10⁻⁵ |

(kronspice's Δf is 0.27 % below the exact value: the same trapezoidal
warping.) In the CW/CCW basis θ is 45° for every reciprocal defect: the
eigenmodes of a reciprocal, lossless ring are real standing waves, equal
mixtures of CW and CCW. A different angle needs a basis tied to a point:
started as the standing wave with its antinode s sections from the
defect's centre, θ = κ·s:

| s | θ predicted | sin²2θ predicted | kronspice sin²2θ | θ |
|---|---|---|---|---|
| 1 | 15.0° | 0.2500 | 0.2507 | 15.02° |
| 2 | 30.0° | 0.7500 | 0.7512 | 30.04° |
| 3 | 45.0° | 1.0000 | 1.0000 | 44.94° |

**Three flavours.** One harmonic gives two states, so three flavours need
a third state near the pair. The plain ring has one at its band edge: the
π-mode (κ = π, neighbouring nodes in antiphase: a standing wave that
cannot travel) lies 0.86 % below the pair at κ = 11π/12, because the CL
dispersion is flat there. With ring flavours f1, f2 = the standing cos and
sin at c = ½ on κ = 11π/12 and f3 = the π-mode, two capacitor defects mix
the three. Their sizes (≤ 0.5) and sections and the assignment of f1–f3 to
e, μ, τ were searched for the closest match to the measured mixing
(NuFIT 6.0, normal ordering with Super-Kamiokande atmospheric data:
Esteban et al., *JHEP* 12 (2024) 216, arXiv:2410.05380, Table 1):

| Quantity | Ring, tuned | Measured |
|---|---|---|
| defects | δ = 0.2064 on section 2, 0.1628 on section 23 | — |
| e, μ, τ | f1, f2, f3 | — |
| θ₁₂ | 33.66° | 33.68° ± 0.72 |
| θ₂₃ | 42.28° | 43.3° ± 0.9 |
| θ₁₃ | 8.56° | 8.56° ± 0.11 |
| splitting ratio | Δω₃₁/Δω₂₁ = 33.56 | Δm²₃₁/Δm²₂₁ = 2.513×10⁻³/7.49×10⁻⁵ = 33.55 |
| CP phase | 0 or 180° only (U real) | δ_CP = 212° (+26/−41), consistent with CP conservation at about 1σ |
| cost | 1.28 (sum of squared deviations in measured σ) | |

The mixing matrix (rows f1, f2, f3; columns the modes by ω, 1.54798,
1.58623, 1.58737×10⁷ rad/s) is

```
-0.1423  +0.7687  +0.5205
-0.6357  +0.3093  -0.6599
+0.6991  +0.4889  -0.5004
```

and up to 8.7 % of a mode lies outside the three ring states (the defects
these angles need mix in other harmonics). A kronspice run started in e
(f1) for one slow beat (1,399 periods) follows the ring's full mode
expansion (all 24 modes, same readout) to 1.1×10⁻⁴; against the 3×3
formula P(α→β) = |Σᵢ U_αi U_βi e^{iωᵢt}|² it deviates by up to 0.277, the
fast beats of that 8.7 %. Ledger 2.6×10⁻⁴; with U real, P(e→μ)(t) and
P(μ→e)(t) agree to 1.4×10⁻² (two separate runs).

**What matches and what does not.** The beat form matches: two-state
conversion sin²2θ·sin²(πΔf·t) with Δf set by the splitting, and three-state
mixing through a unitary (here real orthogonal) matrix. With two defects
chosen for the purpose the ring reproduces the measured angles and
splitting ratio. That is tuning, not prediction: four numbers fitted with
four knobs. What does not match: the neutrino's flavour is defined by the
weak interaction (the charged lepton produced or detected), the splittings
are mass differences entering as Δm²L/4E, matter changes them (MSW), and a
CP phase may be present; the ring's flavours are positions on a ring, its
splittings are defect sizes in hertz, it has no weak production or
detection, no mass, no matter effect, and a reciprocal lossless ring
cannot carry a CP phase. And the user's three types are not three states:
C, S and D span two dimensions, so they cannot be three orthogonal
flavours. Standard prediction: the oscillation formula is the generic
mathematics of coupled modes.

### Test 3: implosion and explosion of the coherent unit

A port on the ring line (R = 811.7 Ω, differential, at node 0, giving a
standing wave with its antinode there Q_ext = 40.5, like the E-ring), a
second port a quarter wavelength on at node 6 where needed. Reference c =
node 0. Free explosion through port 0, 80 periods:

| Type | Energy left | End state | Ledger |
|---|---|---|---|
| S (antinode at the port) | 0.0000 | — | 2.3×10⁻³ |
| D (node at the port) | 0.9957 | D 1.0000 | 9.8×10⁻⁵ |
| C (travelling) | 0.4979 | D 1.0000 (S 0.0000) | 1.2×10⁻³ |

S leaves completely, D is dark to the port and stays (the 0.4 % is the
shunt loss), and C loses its S half and is left as D: **the explosion
converts a travelling wave into a standing wave with its node at the
port.**

Rebuild: the explosion of the time-reversed target (CW ↔ CCW) is recorded
through passive ports and played back reversed in time as each port's
incident wave (PWL, 20 points a period, 400 steps a period):

| Target, ports | Absorbed | Ring at the peak / incident | Fidelity | Explosion | Ledger | Time-stepping |
|---|---|---|---|---|---|---|
| S, port 0 | 0.99867 | 0.9980 | 1.00000 | 0.9956 | 1.6×10⁻⁴ | 3.1×10⁻⁴ |
| C, ports 0 and 6 | 0.99864 | 0.9979 | 1.00000 | 0.9954 | 1.6×10⁻⁴ | 3.1×10⁻⁴ |
| C, port 0 only | 0.99865 | 0.9979 | 0.50000 (it is S) | 0.9956 | 1.6×10⁻⁴ | 3.1×10⁻⁴ |

D cannot be built through port 0 at all (its explosion there is zero). A
travelling wave needs two ports in quadrature; one port builds only its
standing half.

**A flavour superposition.** Defect δ = 0.2 on section 0 and the port at
node 2: the state made at the port, F (antinode at the port), is at θ =
22.5° to the defect's modes. Its free explosion (120 periods) leaves 0.0221
of the energy. F itself leaves at once; what remains is almost all F⊥
(its share reaches 0.9904), which is dark to the port, but the defect's
beat keeps turning it back into F, so it too leaves, more slowly: the
beat goes on through the explosion. Rebuilt
by the reversed explosion: absorbed 0.99082, ring 0.9901 of the incident
energy, fidelity 0.98440 (F 0.9844, F⊥ 0.0156), limited by the 0.0221 the
recording left behind. Rebuilt by a single-tone growth at the port:
absorbed 0.99841, ring 0.9980, fidelity 0.96903. To rebuild a flavour state
the incoming wave must carry both mode frequencies with the right
amplitudes and phases (the time-reversed beat); after the peak the beat
resumes as in the free explosion.

**Status (H39): bounded.** The wave types are stationary states of a ring
and none converts by itself; a perturbation makes the generic two- and
three-mode beat, which reproduces the form of neutrino oscillation and,
when tuned, its angles and splitting ratio; the three types are not three
independent states; implosion and explosion return a type unchanged only
if the port couples to it, and turn a travelling wave into a standing one
otherwise. Nothing in the ring supplies flavour, mass, weak interaction or
CP violation.

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

`viz/implosion.png` (Part 5, from `render_part56.py`): 5a — stored energy
for the time-reversed drive against bursts and the un-reversed ring-down;
incident and outgoing power of the CPA run (the outgoing wave is ~0 while
driven, the explosion after); the outgoing power of drives that are
reflected; the pulsing ring's books per cycle. 5b — energy inside the
Huygens surface and outside it on both grids; the burst; E_φ on the y = 0
plane at the peak and the B lines a quarter period later.

`viz/flavours.png` (Part 6): the three types' patterns on the ring at five
instants; the CW → CCW beats for four defects (time in beat periods) and
the standing starts at θ = 15°, 30°, 45°; Δf against δ (exact, first
order, kronspice); the three-state run against the full mode expansion and
the 3×3 formula; the free explosions of S, D, C; the rebuilds; the flavour
superposition's conversion through the explosion and the rebuild.

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestERingAC` | AC resonance and Q against 1/√(LC) and ω₀L/R; ledger | 2e-4, 2e-3; 1e-9 |
| `TestRingDown` | ring-down Q; kick = radiated + heat + stored; no net charge; mean moment ≈ 0 | 1 %; 1e-3; 1e-9; 0.02 |
| `TestQScaling` | Q_rad ∝ εr^1.5; closed form = ω₀L/R_rad; AC at εr = 300 | 1e-12; 2e-3 |
| `TestBench` | proper time; the required ring resonates at mc²/ħ with exactly the required Q; vacuum ring electrically large | 1e-9 |
| `TestDualRingModes` | plain and crossed CL-ring modes against the ladder dispersion; ledgers | 1e-4; 1e-9 |
| `TestZeroPointInflow` | stored energy = expected for the noise modes; magnetic channel in = out + stored; σ/⟨E⟩ ≈ 1; mean moment ≈ 0 | 4σ or 0.03; 2e-3; 0.7–1.3; 0.01 |
| `TestRampIndependence` | ramps 2Q and 8Q agree on identical noise, at the expected value; an abrupt start shows the start-up excess | 0.002; 4σ or 0.01; > 0.02 |
| `TestImplosionCPA` | 5a: the time-reversed drive reflects only its switch-on (≤ 2e⁻⁸); exact reversal adds (R_d/(R_rad+R_d))²; peak ≥ 0.99; other drives reflected, the decaying drive absorbs < 0.1; ledger and wave books | 1e-5; 1e-4 |
| `TestTrainBooks` | 5a pulsing ring: steady cycle in − out = heat + Δstored; ledger | 1e-4 |
| `TestFieldImplosionCoarse` | 5b on grid r/2, short windows: f₀ and Q near the lumped values; held at the peak / delivered ≥ 1 − 1.5e⁻⁴ (the shortened 4 τ_E window's switch-on); E_φ ≥ 0.99 and B_φ ≤ 0.01 of their energies; field ledgers; burst holds less | 5 %, 10 %; 1e-10 |
| `TestNoFullCancellation` | the CW and CCW waves of a harmonic are orthogonal (nothing non-zero cancels everywhere) | 1e-12 |
| `TestTypesPersist` | C, S, D stationary on the plain and crossed rings; C turns once a period; ledger | 1e-9; 2e-3; 1e-3 |
| `TestFlavourBeat` | Δf from kronspice = exact eigen-splitting; = first order; CW → CCW fully; standing start s = 1 gives sin²30° | 1 %; 2 %; 0.01; 0.01 |
| `TestThreeFlavour` | the fit reaches the measured angles (1.5σ) and ratio (5 %); kronspice = full mode expansion; P(e→μ) = P(μ→e) | 1e-3; 0.03 |
| `TestImplodeTypes` | through one port S leaves, D stays, C is left as D with half its energy; reversed explosion rebuilds S (one port) and C (two ports), one port gives half of C; ledgers | 1e-3, 0.99, 0.01; 0.9999; 1e-3 |
| `TestExport` | E-ring and both ring lines exported (torn-solution tolerance 1e-6 for the E-ring, from 1 GΩ bleeders beside 0.29 Ω) | |

## Files

| File | What it holds |
|---|---|
| `main.go` | `ERing` (L, C, R_rad, Q, the kronspice loop, AC and ring-down), `Bench` (the neutrino test), `DualRing` (inverted ring lines) and the report |
| `main_test.go` | the Part 1–4 tests above |
| `zpf.go` | Part 4: `ZPFRing`, the scaled loop in its radiation resistance's zero-point noise, on the magnetic channel |
| `viz.go`, `viz/render.py` | `-viz` data; the chart `viz/resonance.png` and the sound files |
| `implode.go` | Part 5a: `ERing.BuildPort` (R_rad as a port with SIN-segment incident-wave sources), `RunPort` (wave and ledger books), the drives (`Growth`, `Burst`), `Implosions`, the pulsing `Train` |
| `field.go` | Part 5b: `FieldRing` — the torus on the Yee grid, the Huygens surface stencils, filtering run, ring-down, implosion with the time-reversed mode, burst |
| `flavour.go` | Part 6: `FRing` (CL ring line with defects and ports, imposed start), wave types and projections, `Persistence`, exact modes (Jacobi), `FlavourBeat`, the band-edge three-state cluster (`ClusterMix`, `FitPMNS`, `ThreeFlavour`), explosions and rebuilds |
| `part56.go` | the report for Parts 5 and 6 |
| `viz56.go`, `viz/render_part56.py` | `-viz` data for Parts 5 and 6; `viz/implosion.png`, `viz/flavours.png` |
| `part56_test.go` | the Part 5 and 6 tests above |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |
| `ramp_check.txt` | the ramp-length independence check, raw rows |
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
- Y. D. Chong, L. Ge, H. Cao, A. D. Stone, "Coherent perfect absorbers: time-reversed lasers," *PRL* 105, 053901 (2010) — the no-reflection condition of Part 5a.
- A. Taflove, S. C. Hagness, *Computational Electrodynamics*, 3rd ed., Artech 2005, ch. 5 — the total-field/scattered-field (Huygens) surface of Part 5b.
- I. Esteban, M. C. Gonzalez-Garcia, M. Maltoni, I. Martinez-Soler, J. P. Pinheiro, T. Schwetz, "NuFit-6.0: updated global analysis of three-flavor neutrino oscillations," *JHEP* 12 (2024) 216, arXiv:2410.05380 — the measured angles, splittings and δ_CP of Part 6 (Table 1, normal ordering with SK-atm).
- L. Wolfenstein, *Phys. Rev. D* 17, 2369 (1978); S. P. Mikheyev, A. Yu. Smirnov, *Sov. J. Nucl. Phys.* 42, 913 (1985) — matter (MSW) effects, which the ring has no counterpart of.
- B. W. Lee, R. E. Shrock, *Phys. Rev. D* 16, 1444 (1977); K. Fujikawa, R. E. Shrock, *PRL* 45, 963 (1980) — the standard Dirac-neutrino moment 3eG_F m_ν/(8√2π²).
