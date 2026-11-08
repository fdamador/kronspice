# The fdamador neutrino: the toroidal field inverted (H28)

Hypothesis H28 (CLAUDE.md §1, user hypothesis 2026-11-04): take the
toroidal field of the photon electron ([photon_electron](../photon_electron),
H23) and turn it inside out, so that the **electric field is the ring** and
the **magnetic field circles round it**. The proposal is that this inverted
structure is the neutrino. This example analyzes the inversion with
Maxwell's equations, builds it as kronspice circuits, and puts it on a
neutrino test bench next to the measured values and the standard
prediction. Parts 5 and 6 implode and explode it (H30) and test the
user's reading of the three neutrino flavours as three wave types (H39);
Part 7 spins it (H40).

```
go run ./examples/fdamador_neutrino                                     # Parts 1–4, 5a, 6, 7a: about 6 minutes
go run ./examples/fdamador_neutrino -field                              # also Part 5b, the field solver on grids r/2 and r/3: about 13 minutes (r/2 2 min, r/3 8 min)
go run ./examples/fdamador_neutrino -field -viz examples/fdamador_neutrino/viz # and the chart and sound data
python3 examples/fdamador_neutrino/viz/render.py                       # resonance.png and the .wav files
python3 examples/fdamador_neutrino/viz/render_part56.py                # implosion.png and flavours.png
python3 examples/fdamador_neutrino/viz/render_part7.py                 # spin.png (Part 7)
go test ./examples/fdamador_neutrino                                   # about 9 minutes; -short skips the long runs (15 s)
KRONSPICE_EXPORT=1 go test ./examples/fdamador_neutrino -run TestExport # netlists, KiCad schematics, results sheets
```

`run_output.txt` is the output of the `-field` run; Part 7 was appended from
`-only 7 -field -spincells 2,3`.

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

## 7. A spinning ring (H40)

The user's question (2026-11-07): "If the torus spins about its axis and
spins around the E axis (the ring), would this also add to the missing
information of neutrino data?" After Parts 1–6 the missing pieces are mass,
the weak interaction, a CP phase (the reciprocal ring gives only a real
mixing matrix, δ_CP = 0 or 180°), left-handed helicity, three generations
and confinement.

"Spinning" has two readings, kept apart:

- **7a, a material ring that turns.** The moving medium makes the ring
  non-reciprocal: the co- and counter-rotating waves split (Sagnac).
- **7b, a field that spins.** In vacuum a field "spins" only by carrying
  circulating energy flow, E×H. Spin about the axis needs B round the ring
  beside the E-ring (H23's B-ring with H28's E-ring) or a pattern that turns;
  spin "around the ring" is poloidal circulation, the twist of H23's (2, 1)
  path.

```
go run ./examples/fdamador_neutrino -only 7                  # Part 7a alone: about 1.5 minutes
go run ./examples/fdamador_neutrino -only 7 -field -spincells 2,3 -viz examples/fdamador_neutrino/viz   # with 7b on grids r/2 and r/3: about 27 minutes (r/2 6, r/3 20); the chart uses the finest grid
python3 examples/fdamador_neutrino/viz/render_part7.py      # viz/spin.png
```

### 7a. Rotation on the ring line

A ring turning at Ω delays the co-rotating wave by Δt = 4AΩ/c² per turn.
For a co-moving medium, Fresnel drag makes this delay the same in any medium
(Sagnac, *C. R. Acad. Sci.* 157, 708, 1913; Post, *Rev. Mod. Phys.* 39, 475,
1967). A resonance therefore splits by

Δω/ω = Δt/T_round trip = 2aΩ·v_g/c²,

which is `spacetime.SagnacSplitting` × v_g/c (the ring-laser value for
v_g = c). In the ring's own frame this is a gyroscopic (Coriolis) term,
antisymmetric and workless: Landau & Lifshitz, *Mechanics* §39, the same
structure as Kron's speed tensor G. In kronspice it is a gyrator across every
section of the CL ring line of Part 6:

Cmat·du/dt = −g + 2·Gy·Q·u, Qᵀ = −Q.

On a uniform ring the CW and CCW waves then satisfy c_κω² ± 4·Gy·sin κ·ω −
1/Lg = 0. They split by exactly 4·Gy·sin κ/c_κ. The gyrator transconductance
Gy is set from the non-reciprocity ν = Δω/ω wanted on a pair. The complex
modes are solved exactly, from −A² for the real antisymmetric generator ẏ =
A·y.

**Test 1, the splitting** (harmonic 1, ν = 0.01, started standing, 120
periods, 60 steps a period; Δω = ω_CW − ω_CCW in rad/s):

| Ring | Gy | Δω exact (quadratic) | complex modes | 4Gy sin κ/c_κ | exact at the trapezoidal frequencies | kronspice | S → D amplitude | ledger |
|---|---|---|---|---|---|---|---|---|
| plain | 8.032×10⁻⁵ S | 1.20257×10⁶ | 1.20257×10⁶ | 1.20257×10⁶ | 1.19928×10⁶ | 1.19927×10⁶ | 0.99999 | 9.3×10⁻⁵ |
| crossed | 8.151×10⁻⁵ S | 2.34984×10⁶ | 2.34984×10⁶ | 2.34984×10⁶ | 2.34341×10⁶ | 2.3434×10⁶ | 0.99999 | 1.0×10⁻⁴ |

kronspice measures the splitting at its own trapezoidal frequencies, ω →
(2/h)·atan(ωh/2), to within 3.8×10⁻⁶. Rotation alone turns the standing
wave S into D and back with amplitude 1: the pattern turns. Each gyrator passes
up to 3.2×10⁻³ ω·W (plain ring) from section to section and absorbs 5×10⁻¹⁹
ω·W. The net work of all of them over the run is 6.4×10⁻¹⁹ of the stored
energy. So at constant Ω the rotation does no work on the waves, and its
mechanical ledger column is zero. The rim-speed conversion matches
`spacetime.SagnacSplitting`·v_g/c to 3.5×10⁻¹⁶.

**Test 2, defect and rotation together.** A defect δ = 0.1 on section 0 and
a rotation that splits the pair by the same amount (71,454.2 Hz each). The
combined splitting is 100,862 Hz from the exact modes, against
√(Δd² + Δr²) = 101,052 Hz. Started CW, kronspice gives sin²2θ = 0.5030
(θ = 22.59°) against the predicted Δd²/(Δd² + Δr²) = 0.5000 (22.50°), with
ledger 7.6×10⁻⁵. Part 6 found θ = 45° in the CW/CCW basis for every
reciprocal defect. Rotation gives any θ, because the modes become elliptical:
complex mixtures of CW and CCW.

**Test 3, the CP phase.** This is Part 6's three-state cluster: defects 0.2064
on section 2 and 0.1628 on section 23; e, μ, τ = f1, f2, f3. Rotation is
added with ν on the uniform pair at κ = 11π/12. J = Im(U_e1 U_μ2 U*_e2 U*_μ1)
is the Jarlskog invariant (Jarlskog, *PRL* 55, 1039, 1985), taken from the
unitary matrix closest to the projected U (Löwdin). Reversing the rotation
gives 360° − δ_CP.

| ν | δ_CP | J | J/J_max | θ₁₂ | θ₂₃ | θ₁₃ | Δ31/Δ21 | cost |
|---|---|---|---|---|---|---|---|---|
| 0 | 180.0° | 0 | 0 | 33.66° | 42.28° | 8.56° | 33.570 | 1.285 |
| −1×10⁻⁵ | 180.6° | −0.00036 | −0.0098 | 33.66° | 42.28° | 8.56° | 33.568 | 1.285 |
| −1×10⁻⁴ | 185.6° | −0.00359 | −0.0973 | 33.72° | 42.28° | 8.56° | 33.402 | 1.314 |
| −3×10⁻⁴ | 196.2° | −0.01037 | −0.2792 | 34.17° | 42.28° | 8.57° | 32.144 | 4.341 |
| −1×10⁻³ | 221.5° | −0.02555 | −0.6631 | 37.17° | 42.29° | 8.60° | 23.645 | 198.2 |
| −3×10⁻³ | 236.5° | −0.03426 | −0.8342 | 41.95° | 42.40° | 8.92° | 10.356 | 2,098 |
| +1×10⁻³ | 138.5° | +0.02555 | +0.6631 | 37.17° | 42.29° | 8.60° | 23.645 | 198.2 |

(Beyond ν ≈ 5×10⁻³ the cluster's modes re-order and the mass labels change;
at −0.01 the run reads 40.5°.) The rotation that makes the phase also
splits the close pair, which spoils θ₁₂ and the splitting ratio. Refitting
the two defect sizes and ν on Part 6's sections gives:

| Fit | δ (sections 2, 23) | ν | δ_CP | J | θ₁₂ | θ₂₃ | θ₁₃ | Δ31/Δ21 | cost |
|---|---|---|---|---|---|---|---|---|---|
| δ_CP free within NuFIT's band | 0.2063, 0.1629 | −1.404×10⁻⁴ | 187.9° | −0.00503 | 33.67° | 42.28° | 8.54° | 33.41 | 1.33 |
| δ_CP pinned at 212° | 0.2040, 0.1651 | −5.567×10⁻⁴ | 210.9° | −0.01815 | 33.75° | 42.29° | 8.21° | 31.37 | 17.7 |
| measured (NuFIT 6.0, NO) | | | 212° (+26/−41) | | 33.68 ± 0.72 | 43.3 ± 0.9 | 8.56 ± 0.11 | 33.55 | |

So a phase inside the measured 1σ band (188°) costs nothing in the angles.
The central 212° cannot be reached with θ₁₃ (3.2σ low) and the ratio
(−6.5 %) held at once. A search over every section pair (`-cpall`, 5.4
minutes, run once) found only the mirror pair (1, 22) at the same cost. At
the pinned point:

```
|U|∠arg U  (rows e, μ, τ; columns ν1, ν2, ν3)
0.8168∠+158.4°  0.5553∠+121.2°  0.1563∠ +90.5°
0.3384∠+116.5°  0.6568∠ −90.3°  0.6739∠ +87.2°
0.4672∠+130.2°  0.5102∠ −96.6°  0.7221∠ −92.7°
```

**T violation in kronspice.** At the pinned point, separate runs were started
in e and in μ, 1,313 periods each (one slow beat). They follow the exact
complex mode expansion to 2.1×10⁻⁵ and 2.3×10⁻⁵ (ledgers 2.5×10⁻⁴). The
largest |P(e→μ)(t) − P(μ→e)(t)| is 0.1730 (exact expansion 0.1731). With the
rotation off and the same defects it is 0.0130 (0.0130 exact): that
difference is the readout's, since the flavour projections are not exact
energy projections, and the 3×3 formula is off by 0.397 for the same 8.6 %
leak as Part 6. **So the rotating ring does violate T, and so CP.**

**What the rotation needs** (bench ring a = 10 mm; aΩ = ν·c²/2v_g). The
burst limit is the T1000G carbon-fibre hoop, 1,881 m/s
(`spacetime.BurstSpeed`). It is the strongest material in the package; a
ceramic ring bursts far lower.

| Case | Guide (v_g) | ν | rim speed aΩ | (v/c)² | × burst speed | ½v² per kg |
|---|---|---|---|---|---|---|
| δ_CP 188° (fit) | c (most favourable) | 1.40×10⁻⁴ | 2.11×10⁴ m/s | 4.9×10⁻⁹ | 11.2 | 2.2×10⁸ J |
| δ_CP 188° (fit) | c/10 (εr = 100 ring) | 1.40×10⁻⁴ | 2.11×10⁵ m/s | 4.9×10⁻⁷ | 112 | 2.2×10¹⁰ J |
| δ_CP 188° (fit) | the circuit line's own wave (2,750 m/s) | 1.40×10⁻⁴ | 2.3×10⁹ m/s | ≥ c: impossible | — | — |
| δ_CP 211° (pinned) | c | 5.57×10⁻⁴ | 8.34×10⁴ m/s | 7.8×10⁻⁸ | 44.4 | 3.5×10⁹ J |
| δ_CP 211° (pinned) | c/10 | 5.57×10⁻⁴ | 8.34×10⁵ m/s | 7.8×10⁻⁶ | 444 | 3.5×10¹¹ J |
| δ_CP 211° (pinned) | circuit line | 5.57×10⁻⁴ | 9.1×10⁹ m/s | ≥ c: impossible | — | — |

At the burst speed, with the most favourable guide, the reachable ν is
1.25×10⁻⁵: 44× below the 212° point. δ_CP there is 0.80° from 180°. **The
rotation the CP phase needs is unphysical for any real ring:** 11 to 444
times the burst speed of the strongest fibre, and faster than light for the
slow wave the circuit line carries. Energy: the gyrators do no net work at
constant Ω. The ring's kinetic energy (½v² per kilogram) and the extra ΩJ_z
of a co-rotating wave seen from the lab are both paid by whatever turns the
ring. Nothing enters from the rotation itself.

**Standard prediction.** Rotation-induced non-reciprocity is ordinary Sagnac
physics (ring-laser gyroscopes). In the Standard Model, CP violation in
neutrino mixing is a complex phase of the weak-interaction (PMNS) matrix,
not rotation.

### 7b. The spinning field

This is the torus of Part 5b in `em/fdtd.go` (εr = 100, lossless, staircased,
absorber from 3a, diagnostics inside the ±2a box). Each mode is made as in
5b: a start field excites it, a probe run finds its frequency (the spectral
peak), and a second run's Fourier transform at that frequency is the mode's
field. The modes are:

- **E-ring**: TE01δ, E round the ring, m = 0 (H28).
- **B-ring**: B round the ring, E poloidal, m = 0 (the dual pattern; H23's
  ring-current picture in the same torus).
- **m = 1**: the lowest odd mode the start field D_φ ∝ cos φ excites. Its
  sin φ partner is the same field turned 90° on the grid, which is exact
  because the grid and the staircased torus share that symmetry. "Travelling"
  is cos φ ± i·sin φ: a pattern that turns.

Quantities:

- Angular momentum J_z = ε0∫ρ(E×B)_φ dV in the Abraham (vacuum) form asked
  for, and ∫ρ(D×B)_φ dV in the Minkowski form beside it. Inside a medium the
  two differ; the Abraham–Minkowski question is not settled here.
- A pseudoscalar χ = c·ε0·∫E·B dV / U.
- Helicity h = (1/2c)∫(A·B/μ0 − C·D/ε0) dV, with A and C carried beside the
  solver (Cameron, Barnett & Yao, *New J. Phys.* 14, 053050, 2012). It is
  normalized so that h·ω/U = ±1 for circular light in vacuum (tested); inside
  εr = 100 the same definition gives up to (1 + n²)/2n = 5.05.

For modes of azimuthal order m_k, J_z = Σ m_k U_k/ω_k. So J_z·ω/U is an
energy-weighted mean of whole numbers, and **½ per quantum is never a
property of one coherent classical mode**, only of a mixture.

Grid r/2 (dx = 1.25 mm, 389,017 nodes, 364 s for nine ring-downs of 20
periods). Modes: E-ring 2.1008 GHz, B-ring 4.8245 GHz, m = 1 4.4396 GHz
(97.1 % of its energy in the torus; overlap with itself turned 90° 0.0000,
turned 180° −1.0000).

| Field | J_z·ω/U Abraham: start [range] | Minkowski: start [range] | χ: start [range] | Q | field ledger |
|---|---|---|---|---|---|
| E-ring alone | 0 | 0 | 0 | 41.7 | 3.2×10⁻¹³ |
| B-ring alone | 0 | 0 | — | 112 | 5.8×10⁻¹⁴ |
| E + B, quadrature (+90°) | 0 | 0 | −0.137 [−0.195, +0.198] | 63.3 | 2.8×10⁻¹⁴ |
| E + B, quadrature (−90°) | 0 | 0 | +0.137 [−0.198, +0.195] | 63.3 | 2.8×10⁻¹⁴ |
| E + B, in phase | 0 | 0 | −0.025 [−0.193, +0.199] | 63.3 | 1.0×10⁻¹³ |
| m = 1 standing | 0 | 0 | — | 1,380 | 2.6×10⁻¹³ |
| m = 1 travelling (+) | +0.0112 [+0.0089, +0.0123] | +1.890 [+1.712, +1.987] | 0 | 1,380 | 1.7×10⁻¹³ |
| m = 1 travelling (−) | −0.0112 [−0.0123, −0.0089] | −1.890 [−1.987, −1.712] | 0 | 1,380 | 1.7×10⁻¹³ |
| ½ E-ring + ½ m = 1 travelling | +0.0037 [+0.0037, +0.0072] | +0.622 [+0.622, +1.184] | — | 239 | 2.7×10⁻¹³ |

("0" means below 10⁻⁴; ranges while the energy is above 10⁻³ of its start.)
The helicity h of E + B starts at ±2.10 and swings through [−4.3, +4.3]. But
the potentials, built from the filtered phasors (A = iÊ/ω), miss curl A = B
by 0.69 of max |B| at the start (0.13 for m = 1). So h is read only for its
sign and its swings, and χ is the cleaner number.

What the field solver shows:

- **Spin about the axis as E-ring + B-ring (m = 0).** The azimuthal E×B
  circulation is there: ε0∫ρ|(E×B)_φ|·ω/U = 0.022. It integrates to J_z = 0
  in both forms, because a field unchanged by a rotation about the axis
  carries no angular momentum about it. Its handedness (χ, h) starts with the
  sign of the relative phase. Then, as the two modes (2.10 and 4.82 GHz)
  beat, it swings through both signs within a period. It is not conserved,
  as expected where the dielectric breaks electric–magnetic duality
  symmetry (Fernandez-Corbaton et al., *PRL* 111, 060401, 2013).
- **A pattern that turns (m = 1).** This field has net J_z. Its sign is set
  by the phase (the ± runs are exact mirror images) and held while it
  radiates (Minkowski 1.71–1.99, never changing sign). Its size is not a
  clean 1. The Abraham form is about 0.01 because the field's momentum sits
  in the dielectric. The Minkowski form is 1.9 because the cubic grid has only
  fourfold symmetry, which couples m to m ± 4, and the lattice "m = 1" mode is
  only 58 % m = 1 on the ring circle (27 % m = −3, 7 % m = 5, 8 % m = −7). On
  the staircased torus J_z is conserved only modulo 4.
- **Decay.** The travelling and standing m = 1 patterns ring down with the
  same Q (1,380 on this grid). Spinning neither lengthens nor shortens the
  decay.
- **Handedness against measurement.** Neutrinos are left-handed (helicity
  −½) and antineutrinos right-handed. Maxwell's equations are parity-even,
  so the sign of J_z or χ is set only by how the field is started. Nothing
  in the model selects the left-handed one, and a ring at rest has no
  momentum to measure helicity against.

The "twist" reading (spin around the ring, poloidal circulation) is H23's
(2, 1) double loop. It was run in the field solver in
[photon_electron](../photon_electron), where it leaves within half a period
in vacuum. In the dielectric torus a poloidal travelling wave would be one of
the cross-section modes seen near 7 GHz in 5b. The torus has no poloidal
rotation symmetry (its inner and outer sides differ), so the angular momentum
of that circulation is not conserved either.

**Standard prediction.** A field's helicity and angular momentum are
ordinary electromagnetism. Neutrino handedness comes from the V−A weak
interaction.

**Grid r/3** (dx = 0.833 mm, 912,673 nodes, 1,219 s). Modes: E-ring
2.0571 GHz, B-ring 4.7021 GHz, m = 1 4.4354 GHz (93.1 % inside; turned 90°
0.0000, turned 180° −1.0000). The m = 1 lattice mode is now 95.6 % m = 1 on
the ring circle (3.1 % m = −3, 1.3 % m = 5). The potentials miss curl A = B
by 0.07 (E + B) and 0.03 (m = 1) of max |B| at the start.

| Field | J_z·ω/U Abraham: start [range] | Minkowski: start [range] | χ: start [range] | Q | field ledger |
|---|---|---|---|---|---|
| E-ring alone | 0 | 0 | 0 | 38.9 | 9.5×10⁻¹³ |
| B-ring alone | 0 | 0 | — | 103 | 9.5×10⁻¹³ |
| E + B, quadrature (+90°) | 0 | 0 | −0.149 [−0.217, +0.214] | 58.8 | 6.3×10⁻¹³ |
| E + B, quadrature (−90°) | 0 | 0 | +0.149 [−0.214, +0.217] | 58.8 | 6.3×10⁻¹³ |
| E + B, in phase | 0 | 0 | −0.019 [−0.217, +0.215] | 58.7 | 3.1×10⁻¹³ |
| m = 1 standing | 0 | 0 | — | 201 | 5.7×10⁻¹³ |
| m = 1 travelling (cos φ + i sin φ) | −0.0230 [−0.0238, −0.0225] | −0.793 [−0.837, −0.790] | 0 | 201 | 3.5×10⁻¹³ |
| m = 1 travelling (cos φ − i sin φ) | +0.0230 [+0.0225, +0.0238] | +0.793 [+0.790, +0.837] | 0 | 201 | 3.5×10⁻¹³ |
| ½ E-ring + ½ m = 1 travelling | −0.0074 [−0.0129, −0.0074] | −0.256 [−0.437, −0.256] | — | 92.7 | 4.3×10⁻¹³ |

On the finer grid the lattice mode is nearly pure m = 1, and its Minkowski
J_z·ω/U settles at −0.79. Its sign is the one that pattern's sense of turning
calls for: with time as e^{iωt}, cos φ + i·sin φ = e^{iφ} turns toward −φ.
On r/2 the same phase read +1.89, because there 27 % of the pattern was
m = −3 and 8 % m = −7, which turn the other way and carry 3 and 7 times as
much per quantum. So **the turning pattern's angular momentum converges
toward one unit per quantum (m = 1) as the grid resolves the mode.** It
does not converge toward ½. The Abraham form stays at about 0.02, the share
of the momentum outside the dielectric.

The m = 1 Q moves from 1,383 to 201 between the grids (that mode is poorly
resolved at 2 cells per tube radius). On each grid, though, travelling and
standing decay alike, which is the point being tested. The E + B handedness
behaves the same on both grids: χ starts at ∓0.14–0.15 with the relative
phase and swings through both signs (±0.20–0.22).

### 7c. Mass and confinement

Rotation adds energy, and so mass E/c², but no binding. A material ring at
the burst speed carries ½v² = 1.77×10⁶ J/kg, ΔM/M = v²/2c² = 2.0×10⁻¹¹. In
the lab, a wave of order m in a ring turning at Ω has energy U(1 + mΩ/ω).
At the burst Ω of the bench ring (1.88×10⁵ rad/s, against ω = 1.295×10¹⁰
rad/s) that is mΩ/ω = 1.45×10⁻⁵ for m = 1. In the field solver the turning
pattern decays as the standing one: Q = 1,380 for both on r/2 (Q/ω = 49 ns)
and 201 for both on r/3 (7.2 ns). A centrifugal pull outward needs a wall to push against. The
dielectric is that wall here; a neutrino has none.

### Status (H40): bounded

What rotation supplies, and what it cannot:

| Missing piece | What rotation gives | Measured / standard |
|---|---|---|
| CP phase | Yes, in form: a turning ring is non-reciprocal, U becomes complex, J ≠ 0 and P(e→μ) ≠ P(μ→e) in kronspice (0.173 against 0.013 without). δ_CP = 188° is reached with the angles kept. 212° is reached only with θ₁₃ at −3.2σ and the ratio −6.5 %. Either needs a rim speed 11–444× the burst speed of the strongest fibre (≥ c for the circuit's own slow wave). At the burst speed δ_CP moves 0.8° from 180°. | δ_CP = 212° (+26/−41), a phase of the weak-interaction (PMNS) matrix |
| Left-handed helicity | A sign, chosen by how the field is started; not selected, and not conserved for E + B in the dielectric | helicity −½ (V−A) |
| Spin ½ | No: J_z·ω/U is a mean of whole m; the turning pattern tends to one unit (−0.79 on r/3, 95.6 % m = 1); ½ only for a mixture | ½ |
| Mass | ½v² of the ring and ΩJ_z of the wave, paid by the drive (ΔM/M ≤ 2×10⁻¹¹ at the burst speed) | m < 0.45 eV |
| Confinement | None: the turning pattern decays as the standing one | lifetime ≥ 5.3×10⁴ s (SN1987A) |
| Weak interaction, three generations | Nothing | the content of the Standard Model |

Limits: the circuit's gyrator is the first-order Sagnac term set on one pair
(its ω² dependence across the cluster's 2.5 % spread is not followed). The
field solves are staircased (tube radius 2–3 cells, fourfold lattice
symmetry). The helicity's potentials come from filtered phasors. J_z is
given in both Abraham and Minkowski forms without settling between them.

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

`viz/spin.png` (Part 7, from `render_part7.py`): δ_CP of the rotating
three-state ring against the non-reciprocity ν, with the rim speed it needs
and the burst limit; the mixing angles' drift in measured σ; P(e→μ) and
P(μ→e) in kronspice at δ_CP ≈ 211° and their difference with and without
rotation; 7b — J_z·ω/U of the spinning fields (Minkowski, with Abraham
dotted) and the helicity's swings as the E- and B-rings beat.

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
| `TestSagnacSplit` | 7a: gyrators split harmonic 1 by 4·Gy·sin κ/c_κ (quadratic = complex modes = first order); kronspice = exact at its trapezoidal frequencies; S → D amplitude 1; gyrators absorb nothing, pass energy; ledger (plain and crossed) | 1e-9, 1e-12; 1e-3; 0.01; 1e-12; 1e-3 |
| `TestRimTwin` | rim speed aΩ = νc²/2v_g equals `spacetime.SagnacSplitting`·v_g/c; speeds ≥ c flagged impossible | 1e-12 |
| `TestJarlskog` | no rotation: J = 0, δ_CP = 180°, Part 6's angles; with rotation J ≠ 0; reversing it reverses J (δ → 360° − δ); mode residual | 1e-12, 1e-6, 0.01; 1e-9; 1e-10 |
| `TestTwoStateRot` | defect + equal rotation: CW → CCW sin²2θ = ½; Δf kronspice = exact; exact ≈ √(Δd² + Δr²) | 0.01; 1 %; 2 % |
| `TestTViolation` | pinned CP point: kronspice = exact complex modes; P(e→μ) − P(μ→e) matches the prediction and is > 3× its no-rotation (readout) value; gyrator work | 2e-3; 0.01; 1e-9 |
| `TestHelicityNorm` | h·ω/u = ±1 for circular light, 0 for linear, at every instant | 1e-12 |
| `TestFieldIntegrals` | grid J_z of E ∝ (x, y, 0), B ∝ ẑ against the closed form −ε0E0B0·16b⁵/3L; Minkowski path; ∫E·B = 0 | 0.3(dx/b)²; 1e-12 |
| `TestSpinFieldCoarse` | 7b on r/2, short runs: m = 1 orthogonal to itself turned 90°, −1 at 180°; travelling ± are mirror images with net Minkowski J_z > 0.5, Abraham of the same sign; E + B (m = 0) J_z ≈ 0, helicity sign flips with the phase; field ledgers | 0.01; 1e-9; 0.02; 1e-10 |
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
| `spin.go` | Part 7a: gyrators on `FRing` (the ring's rotation), exact complex modes (`RotModes`), the rotating three-state cluster with J and δ_CP (`RotClusterMix`, Löwdin unitarization), `FitCP`, the Sagnac, two-state and T-violation runs, rim speed against the burst limit (`Rim`) |
| `spinfield.go` | Part 7b: `SpinField` — modes of the torus found and filtered (E-ring, B-ring, m = 1), turned 90° on the grid, combined; ring-downs with A and C carried beside the solver; J_z (Abraham, Minkowski), ∫E·B, helicity, Q; azimuthal content |
| `part7.go` | the report for Part 7 |
| `viz7.go`, `viz/render_part7.py` | `-viz` data for Part 7; `viz/spin.png` |
| `part7_test.go` | the Part 7 tests above |
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
