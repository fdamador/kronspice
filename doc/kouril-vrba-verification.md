# Kouřil & Vrba — transcription and checks

F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits: Principles,
Theory and Applications*, Ellis Horwood / John Wiley / SNTL, 1988, 413 pp.
(ISBN 0-85312-606-2). Page photographs taken from the user's copy are in
`reference/kouril-vrba/` (`kv_pNNN-MMM.jpg`, one per two-page spread); page
numbers below are the printed ones. Equations are transcribed from the
photographs; symbols follow the book.

## Batch 1: §4.4–4.5, pp. 128–142

### §4.4.1 Circuits with controlled linear elements (pp. 128–130)

- **(4.63)** `y(t) = P(t)·x(t)` — the output of a controlled (parametric)
  linear element is the product of the time-varying parameter P(t) and
  the input x(t).
- **(4.64), (4.65)** `P(t) = Σ a_k1·e^(jk1ω1t)`,
  `x(t) = Σ b_k2·e^(jk2ω2t)`; **(4.66)**
  `y(t) = Σ Σ c_k1,k2·e^(j(k1ω1 + k2ω2)t)`, `c = a_k1·b_k2`: the output
  contains the combination frequencies `k1ω1 + k2ω2`.
- **p. 129, energy (quoted in substance):** in a circuit with time-varying
  elements, the energy delivered to the loads and the energy turned into
  heat come not only from the voltage or current source in the circuit
  but also from the source that changes the element's parameter; how the
  two shares divide depends on the conditions in the circuit.
  → *Engine:* this is why the N element's parameter variation is driven
  through a mechanical port and booked there (`device/machine.go`), and
  why the other conversation's time-varying C and L book `ṡ·∫Q` on a
  mechanical port.
- **(4.67)** if y is the charge q (controlled capacitor) or the flux
  linkage ψ (controlled inductor), the current or voltage is
  `z = dy/dt = Σ Σ j(k1ω1 + k2ω2)·c_k1,k2·e^(j(k1ω1+k2ω2)t)`: amplitudes
  proportional to the combination frequency, phases shifted by π/2.

### §4.4.2 Autoparametric non-linear circuits (pp. 130–133)

- **(4.68)** `y + Δy = f(x + Δx) ≈ f(x) + P_d(x)·Δx`, `P_d = dy/dx` the
  differential parameter.
- **(4.69)** with `x = X0 + X1·cos α1`, `Δy = P_d(t)·Δx`: a small signal
  sees a linear element with a time-varying parameter — superposition
  holds for small signals while the spectrum is transformed (p. 131,
  Fig. 77: `Δi = G_d(V0 + V1·cos α1)·Δv`).
- **(4.70)–(4.77)** `P_d(α1) = P0 + P1·cos α1 + P2·cos 2α1 + …` with
  `P0 = (1/2π)∫P_d dα1`, `Pk = (1/π)∫P_d·cos kα1 dα1`; a small harmonic
  `Δx = X2·cos(ω2t + φ2)` produces the side frequencies `ω2 ± kω1`
  (4.76), and for storage elements the derived current or voltage
  (4.77).
  → *Engine:* this is the linear periodically time-varying small-signal
  model ("periodic AC"). Not implemented; it is the natural next analysis
  for parametric amplifiers and mixers.

### §4.5 General relations for powers (pp. 133–142)

The book cites Manley and Rowe [43] and Pantell [47] (reference list,
p. 396): J. M. Manley, H. E. Rowe, "Some General Properties of Nonlinear
Elements — Part I," *Proc. IRE* 44 (1956), No. 7, pp. 904–913; R. H.
Pantell, "General Power Relationship for Positive and Negative Nonlinear
Resistive Elements," *Proc. IRE* 46 (1958), No. 12, pp. 1910–1913.

- **(4.78)–(4.81)** `i, v, q, ψ` as two-dimensional Fourier series with
  coefficients `J, U, Q, S` (two-sided, `J_−k1,−k2 = J*_k1,k2` etc.);
  `J = j(k1ω1 + k2ω2)·Q`, `U = j(k1ω1 + k2ω2)·S`.
- **(4.82)** active power at `k1ω1 + k2ω2`: `P = U·J* + U*·J`;
  **(4.83)** reactive power `W = −j(U·J* − U*·J)`. With peak phasors
  `V = 2U`, `I = 2J` this is `P = ½·Re(V·I*)`, the engine's AC convention.
- **(4.84)** `α1 = ω1t`, `α2 = ω2t`.

**Non-linear resistor (Pantell), pp. 135–137**

- **(4.85)–(4.88)** `(∂i/∂αn)(∂v/∂αn) = G_d(v)·(∂v/∂αn)²`, averaged over
  both phases: `½ Σ Σ k_n²·P_k1,k2 = g_n` (two-sided),
  `g_n = (1/(2π)²)∫∫ G_d(v)·(∂v/∂αn)² dα1 dα2`.
- **(4.89), (4.90)** on the one-sided spectrum `k1ω1 + k2ω2 > 0`:
  `Σ Σ k_n²·P_k1,k2 = g_n` (n = 1, 2); `g_n ≥ 0` unless the differential
  conductance is negative (otherwise the circuit is potentially
  unstable).
- **(4.91)** frequency multiplier: `P1 − c²·Pc ≥ 0`, efficiency
  `η = Pc/P1 ≤ 1/c²` — at most 25 % to the 2nd harmonic, 11 % to the 3rd
  (Fig. 78: 100 mW at 1 MHz → 25 mW at 2 MHz, 11 mW at 3 MHz).
- **(4.92)** mixer: `η = Po/(P1 + P2) ≤ 1/(c1² + c2²)`, 50 % for simple
  mixing (Fig. 79).
- → **`analysis.TestPantellResistiveMultiplier`**: a piecewise-linear
  resistor (B element, 1 mS / 21 mS knee) driven through 50 Ω at 1 kHz;
  from the engine's transient waveforms, `Σ k²·P_k = 1.1877 mW` against
  `⟨G_d·(dv/dα)²⟩ = 1.1853 mW` (0.2 %, 200 harmonics); the 2nd harmonic
  takes 4.9 % ≤ 25 % and the 3rd 0.3 % ≤ 11 %.

**Non-linear storage element (Manley–Rowe), pp. 138–141**

- **(4.93)** `(1/(2π)²)∫∫ q·∂v/∂αn dα1 dα2 = 0` (q and v periodic in each
  phase).
- **(4.94), (4.95)** `Σ Σ k_n/(k1ω1 + k2ω2)·P_k1,k2 = 0` (n = 1, 2), on
  the two-sided or the one-sided spectrum. The weighted powers depend on
  the combination frequencies and coefficients only — not on the drive
  level nor on the shape of the characteristic. By duality the same holds
  for a non-linear inductive element.
- **(4.96)** storage frequency multiplier: `P1 = Pc`, `η = 1` (Fig. 80).
- **(4.97), (4.98)** mixer: `P1 : P2 : Po = c1ω1 : c2ω2 : (c1ω1 + c2ω2)`,
  `P1 + P2 = Po`, `η = 1`.
- **(4.99), (4.100)** up-converter (`c1 = c2 = 1`):
  `P1 : P2 : Po = ω1 : ω2 : (ω1 + ω2)`, power gain `K_P = 1 + ω2/ω1`.
- **(4.101)** `c1 = −1, c2 = 1`: `−P1 : P2 : Po = ω1 : ω2 : (ω2 − ω1)`. The
  signal source at ω1 *receives* power; only the pump at ω2 delivers; the
  element appears as a negative resistance at ω1 and the circuit is
  potentially unstable (p. 141) — self-excited oscillation at ω1 and
  ω2 − ω1 (subharmonic when rational), or an amplifier when an auxiliary
  selective load absorbs ω2 − ω1 (the idler).
- → **`analysis.TestManleyRowePumpedInductor`**: the N element as a
  pumped inductance `L(t) = L0·(1 + 0.2·cos ωp·t)` (pump 11 kHz through
  the mechanical port), signal source at 3 kHz and idler load at 8 kHz
  behind series-tuned branches. For an inductance pumped mechanically
  rather than electrically, (4.93) still gives the n = 1 relation from the
  terminal powers alone (`∫∫ L(α2)·i·∂i/∂α1 = 0`), and energy conservation
  turns the n = 2 relation into `Σ k2·P/(k1ωs + k2ωp) + Pp/ωp = 0`, with Pp
  the mechanical pump power. From the engine's transient: both sums
  vanish to 3.3e−4 and 1.9e−4 of the signal term over every combination
  frequency up to 60 kHz; the signal source receives power
  (`Ps = −42.1 µW`) while the pump delivers `Pp = 158.8 µW` and the idler
  takes 113.8 µW; per rad/s 2.232e−9 : 2.297e−9 : 2.264e−9 (eq. 4.101's
  equality within 3 %, the rest carried by ωp ± ωs and higher products,
  which the sums include).

**Reactive powers, pp. 141–142**

- **(4.102)** non-linear resistor: `Σ Σ k_n·W_k1,k2 = 0` — detuning one
  resonant circuit introduces reactances in the others.
- **(4.103)–(4.106)** storage element:
  `Σ Σ k_n²/(k1ω1 + k2ω2)·W = h_n`, `h_n = ⟨C_d(v)·(∂v/∂αn)²⟩`, or for an
  inductor `⟨L_d(i)·(∂i/∂αn)²⟩`. Not yet tested.

### §4.6 (start, pp. 142–143)

Polyphase symmetrized systems (push-pull, ring, bridge, balanced
circuits); Fig. 82: the basic configuration — an active linear circuit
feeding a non-linear or parametric element controlled by p(t). To be read
with the rest of §4.6.

## Batch 1.3: §6.4.4–6.4.5, pp. 222–228

### §6.4.4 Parametric amplifier with one resonant circuit (pp. 222–225)

- **(6.23)** `C(t) = C0/(1 + m·sin 2ω1t)`, pumped at `ωs = 2ω1`, circuit
  tuned to `ω1 = 1/√(LC0)` (Fig. 130a).
- **(6.24)** the pumped capacitor adds `R_e = m·cos 2φ/(2ω1C0)` (from eq.
  10.81), negative for `φ ∈ (π/4, 3π/4)`, most negative at `2φ = π`.
- **(6.25)–(6.27)** `A_v = V2/V1 = Q/(1 + (m/2)·Q·cos 2φ)`;
  `m_c = −2/(Q·cos 2φ)`; at `2φ = π`, `A_v = Q/(1 − m/m_c)`. The amplifier
  is *synchronous*: signal and pump phases must be locked.
- **(6.28), (6.29)** asynchronous pump `ωs = 2(ω1 + Ω)`: `A_v(t)` swings
  between `Q/(1 − mQ/2)` and `Q/(1 + mQ/2)` with period `π/Ω` (Fig. 130c:
  `mQ/2 = 0.8` gives `5Q` and `5Q/9`); Fig. 130b: varactor realization.
- → **`analysis.TestSynchronousParametricAmplifier`**: series R, L, C(t),
  Q = 50, m = m_c/2, signal phase 90° and 0° to the pump: voltage gain
  across L 99.95 and 33.35 against eq. 6.27's 99.99 and 33.33.

### §6.4.5 Parametric amplifier with two resonant circuits (pp. 225–228)

- Tanks `C1L1` at `ω1` and `C2L2` at `ω2 = ωs − ω1`, coupled by C(t) only
  (Fig. 131a). By Manley–Rowe (§4.5.3), for `ωs ≫ ω1` the idler output at
  `ω2` carries much more power than the signal input.
- **(6.30)–(6.36)** `C(t) = C0 − C1·cos ωs·t` (the book reuses C1 for the
  modulation amplitude; called ΔC below); the ω2 current develops
  `v2 = ½·V1·ω2·ΔC·R2·sin(ω2t − φ1)` (6.34) across the idler tank, which
  returns `½·V2·ω1·ΔC·cos(ω1t + φ1)` at ω1 (6.36).
- **(6.37)–(6.39)** `Y_e(ω1) = jω1C0 − ¼·ΔC²·ω1·ω2·R2`: capacitance C0 and
  negative conductance `G_e = −¼·ω1·ω2·ΔC²·R2`.
- **(6.40)** stable while `G1 + G_e > 0`, i.e. `¼·ω1·ω2·ΔC²·R2 < 1/R1`.
- → **`analysis.TestTwoCircuitParametricAmplifier`**: tanks at 100 and
  230 kHz, Q = 200, C0 = 3 nF of 100 nF per tank, ΔC set for
  `¼·ω1ω2ΔC²R1R2 = 0.5` (eq. 6.39 gain 2). Engine gain 1.785 and idler
  3.81 V; an exact harmonic balance over the sidebands `ω1 + n·ωs`
  (|n| ≤ 6, in the test) gives 1.759 and 3.739 V, and an independent
  fine-step RK4 integration 1.760. The book's 2 assumes each tank is a
  short circuit at the other's frequency; here an off-tune tank is ≈ 5 Ω
  against the coupling capacitor's 145–330 Ω, enough to lower the gain by
  12 % because the gain `1/(1 − k)` amplifies any change in k. The engine's
  remaining 1.5 % is trapezoidal integration error at the default
  tolerance on two high-Q tanks, amplified the same way (a fixed-step
  trapezoidal integration in Python gives 1.80 at 100 steps per pump
  period).

## Batch 1.2: §10.1 end, §10.2, pp. 370–393

### §10.1 end (pp. 370–371)

- **(10.24)–(10.29)** inverting non-linear amplifier with active loads
  (Fig. 227): transfer characteristics, reference voltage `V_r`, and the
  curved knee when `v2X < v3X`. Not modelled.

### §10.2.1 Resonance in non-linear resonant circuits (pp. 371–377)

- Definition: a second-order circuit with at least one non-linear storage
  element; analysed by equivalent linearization (§5.2.1) or slowly
  changing amplitudes, assuming the response is almost harmonic (high Q).
- **(10.30)–(10.34)** capacitive diode
  `C_d(v) = C_d(0)/√(1 − v/Vφ)`; about bias V0,
  `C_d(v1) = C0(1 + √a·v1)^(−1/2)`, `C0 = C_d(0)·√(Vφ/(Vφ − V0))`,
  `√a = 1/(V0 − Vφ)`; `C_d ≈ C0 + a1·v1 + a2·v1²`, `a1 = C0√a/2` (printed
  without the sign the expansion gives: `−C0√a/2`), `a2 = 3C0a/8`.
- **(10.35)–(10.37)** with `v1 = V1·sin ωt` the book takes the first
  harmonic of `C_d(v1)·dv1/dt` as `(C0 + a2V1²/2)·ωV1·cos ωt`, so
  `C_e = C0(1 + 3aV1²/16)`.
- **(10.38)** `ω_r(V1) = ω0/√(1 + 3aV1²/16)`; **(10.39)–(10.46)** driven
  parallel tank: amplitude (10.40)–(10.44) and phase (10.45), (10.46)
  characteristics; Figs. 229–230: maxima tilt to lower frequencies; jump
  phenomena and hysteresis between `α_L` and `α_H` (`V_max` above ≈ 1 V).
- **Discrepancy found.** The fundamental of eq. 10.35 is
  `(C0 + a2V1²/2 − a2V1²/4)·ωV1·cos ωt`: the term `−(a2/2)V1²·cos 2ωt`
  times `cos ωt` contributes `−a2V1²/4`, which the book drops, so the
  describing-function result is `C_e = C0(1 + 3aV1²/32)`. Further, the
  method assumes harmonic voltages are negligible; for the circulating
  current of a free or high-Q tank, the a1 term makes a second harmonic
  of ≈ √a·V1²/6, not small, which mixes back to the fundamental. For this
  curve the charge is exactly `q = (2C0/√a)(√(1 + √a·v) − 1)`, so
  `v(q) = q/C0 + √a·q²/(4C0²)` is quadratic, and the Lindstedt–Poincaré
  expansion of `L·q̈ + v(q) = 0` gives
  `ω/ω0 = 1 − (5/192)·a·V1² + O(V1⁴)`, against the book's
  `1 − (18/192)·aV1²` (and `9/192` with the dropped term restored). The
  qualitative conclusions (lower frequency with amplitude, tilted
  curves, jumps) are unchanged.
- → **`analysis.TestNonlinearTankBackbone`**: free oscillation of 1 mH
  with this varactor curve (V0 = −4 V, Vφ = 0.7 V) at V1 = 1 V
  (aV1² = 0.045): engine `ω/ω0 = 0.998818`, i.e. 5.03/192·aV1²;
  Lindstedt 0.998825; an independent RK4 of the ODE (Python) 5.033/192;
  eq. 10.38 0.995798. At V1 = 3 V the book's formula is off by 2.4 % in
  frequency.

### §10.2.2 Resonance phenomena in parametric resonant circuits (pp. 377–393)

- Pp. 377–378: a parametric storage element draws energy at one frequency
  and supplies it at another; the energy comes from the source that
  changes the parameter, which may be mechanical, light or thermal — the
  circuit then converts that energy to electrical (the same point as
  p. 129; the engine books it on the element's scale port).
- **(10.47)–(10.50)** step change of C at constant charge Q1:
  `W2 = (C1/C2)·W1`, `V2 = (C1/C2)·V1`; `C1/C2 = 5` gives `W2 = 5W1`, paid
  by the work against the field.
- **(10.51)–(10.54)** square-wave pumping (Fig. 232), capacitance reduced
  at voltage peaks, restored at zeros: oscillation grows if
  `ΔC/C0 > ω0·C0·R·π`, i.e. `ΔC/C0 > π/Q`.
- **(10.55)–(10.65)** `C(t) = C0/(1 + m·cos 2ωt)`:
  `q̈ + 2αq̇ + ω0²(1 + m·cos 2ωt)·q = 0`, `α = R/2L`, `ω0² = 1/(LC0)`;
  with `τ = ωt`, `c = ω0²/ω²`, `2b = m·ω0²/ω²`, `2d = 2α/ω` and
  `q = y·e^(−dτ)`: Mathieu's equation `ÿ + (a + 2b·cos 2τ)·y = 0`,
  `a = c − d²`.
- **(10.66)–(10.70)** `y = A1·e^(μτ)F(τ) + A2·e^(−μτ)F(−τ)`; lossless
  instability when μ is real (10.67); with losses `|μ| > d` (10.70).
  Fig. 233: Ince–Strutt chart, tongues at `a = n²` (`ω1 = 2ω0/n`), most
  favourable at `a = 1` (`ω1 = 2ω0`, energy delivered twice per period).
  The same holds for time-varying inductance (p. 385).
- **(10.71)–(10.92)** amplitude limiting: elastance
  `D(q, t) = D0(1 + m·sin 2ω1t) + a·D0·q²` (Fig. 234), equivalent
  impedance `Z_e(ω1) = R_e + 1/(jω1C_e)`,
  `R_e = m·cos 2φ/(2ω1C0)` (10.81), `C_e` (10.82), maximum power
  delivered `P_Cmax = −¼·m·ω1·C0·V1²` (10.87); with `x = ω1/ω0`,
  `d = ω0C0R`, `A² = a(I1/ω1)²`:
  `2dx + m·cos 2φ = 0`, `x² − (1 + ¾A²) − (m/2)·sin 2φ = 0` (10.90), so
  `x² = (1 − d²/2 + ¾A²) ± ½√(m² − 4d²(1 + ¾A²) + d⁴)` (10.92)
  (checked: the algebra from 10.90 is exact). Figs. 235–236: curves tilt
  right for a > 0, jumps at D and F.
- **(10.93)–(10.103)** loci of maxima (10.94) and of vertical tangents
  (10.96, `A² = (4/3)(x² − 1)`); bandwidth `2Δ ≈ ½√(m² − 4d²)` (10.100);
  critical `m_c = 2d = 2/Q` (10.101); detuned
  `m_c = 2√(d² + (2Δω/ω0)²)` (10.103); Fig. 237.
- → **`analysis.TestParametricSelfExcitation`**: series R, L, C(t),
  Q = 50, pumped at 2ω1 with ω1 = ω0 and ω1 = 1.02·ω0. Envelope growth
  `σ = −α + √((mω0/4)² − Δω²)` (the averaged Mathieu exponent, zero at
  eq. 10.103's `m_c`): decay at 0.9·m_c and growth at 1.1·m_c in both
  cases; tuned at 2·m_c σ/ω0 = 0.01000 (predicted 0.01000), detuned at
  1.5·m_c 0.01663 (0.01693).
- → **`analysis.TestParametricOscillatorAmplitude`**: the book's Fig. 234
  circuit (pumped linear capacitor in series with a cubic-elastance
  capacitor), m = 0.1, d = 0.02, tuned: engine steady `A² = 0.060899`,
  eq. 10.92 at x = 1 gives 0.06092 (0.03 %).

## Batch 1.4: references, pp. 394–397

62 entries. Those cited from the pages transcribed so far: [4] Angot
(Ince–Strutt chart), [8] Blackwell–Kotzebue, *Semiconductor-Diode
Parametric Amplifiers* (1961), [9] Blaquière, [19] Cunningham, [33]
Kharkevich, *Non-linear and Parametric Phenomena in Radio-Electronics*
(1956), [39] Louisell, *Coupled Mode and Parametric Electronics* (1960),
[43] Manley–Rowe and [47] Pantell (full citations under §4.5 above), [45]
Minorsky, *Nonlinear Oscillations* (1962). Also relevant to this project:
[14] Chua, "Memristor — The Missing Circuit Element" (1971), [16] Chua,
"The Rotator" (1967), [58] Tucker, *Circuits with Periodically-Varying
Parameters* (1964), [60] Vrba, "Polyphase Systems of Non-linear and
Parametric Circuits" (1970).

## Next batches

The second list, less what is now done: §1.4, §2.1.4–2.1.6, §3.7, §4.6,
§5.2 (equivalent linearization, used throughout §10.2) and §6.4.2.
