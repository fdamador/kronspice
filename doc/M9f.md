# M9f — quaternion polyphase networks and GA power theory

CLAUDE.md §3b items 2 and 3. Files: `kron/polyphase.go`,
`analysis/power.go`, `examples/power_filter/`.

## Quaternion polyphase (`kron/polyphase.go`)

Three phase values are the pure quaternion x_a·i + x_b·j + x_c·k.

- **Clarke and Park as rotors.** The power-invariant Clarke matrix is
  orthogonal with determinant +1, so it is a rotation; Park adds a rotation
  by −θ about the zero-sequence axis. `ClarkeRotor()` and `ParkRotor(θ)`
  are the unit quaternions, and `Rotate3` applies r·x·r*. The matrices
  (`ClarkeMatrix`, `ParkMatrix`) are kept for the twin tests.
- **Quaternion power.** `InstantPower(v, i)` = v·i*: the scalar part is
  p = Σ v·i, the vector part −v×i (Akagi's instantaneous imaginary power,
  as a vector). p is invariant under every frame change and the vector
  turns with the frame.
- **Kron's rotating frame, three-phase.** `ParkFrame(n, triples, p)` is
  the moving frame C(θ) = Pᵀ(pθ) built from the Park rotor, with its exact
  derivative, for `RotatingFrame` (M9c.1). It closes the M9c.1 gap
  (`RotorFrame` handled two-phase pairs only).

## GA power (`analysis/power.go`)

- `HarmonicsOf(t, v, i, f0, periods, hmax)` resamples the last whole
  periods uniformly and takes RMS cosine and sine components by DFT.
- `GeometricPower(h)` forms M = V·I in Cl(2H) (`hyper/ga`): scalar part P;
  same-harmonic bivectors Q_k = Im(V_k·I_k*) (positive inductive);
  cross-harmonic bivectors grouped by pair of harmonics, largest first,
  which say which pairs of harmonics cause the distortion power; and the
  IEEE 1459 groups S₁ = |V₁I₁|, D_I = |V₁I_H|, D_V = |V_HI₁|, S_H = |V_HI_H|.
- `IEEE1459Of(t, v, i, f0, periods)` computes IEEE Std 1459-2010
  single-phase quantities conventionally: time-domain RMS and mean power,
  fundamental phasors, harmonic parts by subtraction, THD and power
  factors.

## Checks

| Check | Result |
|---|---|
| Clarke and Park rotors = matrices; r·x·r* = matrix product | 500 random angles and vectors, 1e−14 |
| Balanced set → constant d, q (√3·V·e^(jφ)); zero sequence → 0 axis only | 1e−12 |
| s = v·i*: scalar Σv·i, vector −v×i; p invariant and vector co-rotating under Clarke and Park; balanced p = 3VI cos φ, \|q\| = 3VI sin φ | 1e−13, 1e−10 |
| Salient-pole stator in Kron's Park frame: L′ = diag(Ls + Ms + 1.5Lm, Ls + Ms − 1.5Lm, Ls − 2Ms); speed voltages e_d = −ωLq·i_q, e_q = ωLd·i_d against the abc flux differentiated numerically; D = dC/dθ | 1e−15; 1e−6 |
| **GA vs IEEE 1459** on random band-limited waveforms (9 harmonics, DC offset): P, S, S₁, D_I, D_V, S_H, P₁, Q₁; S² = P² + Q² + D²; S² = S₁² + D_I² + D_V² + S_H²; Σ pair² = D² | 50 cases, 1e−9 |
| R–L from the engine's transient: P = I²R, Q₁ = I²ωL, no distortion | 1e−4 |
| Single-phase bridge rectifier from the engine: GA and IEEE 1459 agree (S to 1 %: harmonics above the 32nd); largest distortion pair V₁ × I₃ | PF 0.56, THD_I 140 % |

## Example: `examples/power_filter`

A 60 Hz, 120 V six-pulse diode rectifier (0.5 mH line inductance,
2200 µF and 10 Ω on the DC side), two series L–C branches per phase at the
PCC, four E12 genes, fitness 1 − PF from the GA power of the three phases.

| | PF | THD_I | largest distortion pairs (phase a) |
|---|---|---|---|
| no filter | 0.903 | 40 % | 1–5 961 VA, 1–7 333 VA, 1–11 192 VA |
| evolved (30 × 20, 562 evaluations) | 0.957 | 28 % | 1–5 611 VA, 1–7 402 VA, 1–13 151 VA |

The evolved branches sit at 284 Hz and 673 Hz rather than exactly at the
5th and 7th harmonics: maximizing PF trades harmonic absorption against
the leading fundamental reactive power the capacitors add (Q₁ from
+515 var to −226 var per phase). The quaternion instantaneous power of the
three phases averages to ΣP (7389 W and 7736 W, equal to the GA scalar
parts). The test (`go test ./examples/power_filter`) runs a short GA and
requires PF to improve by at least 0.01.

## Not done here

- IEEE 1459 three-phase effective quantities (S_e, the positive-sequence
  P₁⁺ and unbalance power) and GA power for three-phase sets.
- An FFT: the DFT over the needed harmonics costs O(N·H) and is not a
  bottleneck at these sizes.
