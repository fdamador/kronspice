# Tewari's Space Power Generator test through the energy ledger

P. Tewari's Space Power Generator (SPG) is a homopolar machine: an
electromagnet on a rotating iron core, a static yoke, output brushes on
slip rings, with the Faraday-disk EMF V = ½ωBR² (UPSM eq. 6.16). His book
reports one loaded test and reads it as an efficiency above one. This
example runs that test, number by number, through the energy ledger
(hypothesis H13, CLAUDE.md §1; `doc/tewari-references.md`). Whether the
shaft bearings could supply the missing torque is asked in
[bearing_motor](../bearing_motor); the electron model of the same
hypothesis is on the bench in [space_vortex](../space_vortex).

```
go run ./examples/space_power_generator
go test ./examples/space_power_generator
```

## The published test (UPSM pp. 159–160)

A dc drive motor (DM) runs the SPG at 2700 r/min.

| Condition | Reported |
|---|---|
| No load | DM input W0 = 2943 W; excitation 650 W; EMF E = 3.41 V |
| Loaded, same speed | I = 3438 A; DM input W_L = 8620 W |

From these the book computes W_r = W_L − W0 = 5677 W (6.20),
W_E = E·I = 11 724 W (6.21), an incremental efficiency of 2.066 (6.22) and
2.49 for the SPG alone after dividing by the DM's efficiency 0.83 (6.23),
and calls it a violation of energy conservation (p. 161).

## The ledger

1. **Arithmetic.** 3.41 V × 3438 A = 11 724 W; eq. 6.22 prints 11 727 W,
   a typo. The ratios 2.065 and 2.49 follow.
2. **All inputs.** DM 8620 W + excitation 650 W = 9270 W in, against the
   claimed 11 724 W out: ratio 1.26. The excitation is left out of the
   book's ratio.
3. **The assumption.** W_E multiplies the load current by the EMF measured
   at no load. Standard electrodynamics requires the mechanical power the
   rotor takes from the shaft to equal E_load·I (Poynting's theorem; the
   back torque is k·I when E = k·ω, wherever the reaction lands). The
   book's own numbers give the extra shaft power under load:
   0.83 × 5677 W = 4712 W, a torque increase of 16.7 N·m at 282.7 rad/s.

| | Loaded EMF | Back torque at 3438 A |
|---|---|---|
| Standard prediction | ≤ 1.37 V (≤ 1.65 V even with a lossless DM) | 16.7 N·m, as the drive motor delivered |
| The claim | 3.41 V, unchanged from no load | needs k·I = 41.5 N·m |

So the standard prediction is that the EMF falls by at least 60 % at
3438 A (armature reaction, saturation, speed sag at the rotor; the load
loop carries 3438 ampere-turns against the excitation). The implied loop
resistance E0/I is 0.992 mΩ. The book reports no loaded-EMF measurement;
it says the drops "can also be" measured (p. 160).

**Decisive measurement, not in the book:** the sum of the voltage drops
around the loaded loop (brushes, rings, conductor, shunt, load) at
3438 A, or the shaft torque from a transducer rather than inferred from
the DM's input and an assumed efficiency. A closed-loop, self-running test
would settle it outright. Status: *untested*; the claim rests on an
unmeasured loaded EMF.

## Files

| File | What it holds |
|---|---|
| `main.go` | `Published` test values, `Analyse` (eqs. 6.20–6.23, all inputs, the loaded-EMF bound) |
| `main_test.go` | the book's numbers and the bound (1.3706 V, 1.6513 V lossless), torque consistency to 1e-12 |

## Sources

- P. Tewari, *Universal Principles of Space and Matter*, Grain and Chaff, Toronto, 2017 (1st ed. 2009), §6.7, pp. 149–161 (`reference/tewari_upsm_2017.pdf`). Hypothesis source.
- P. Tewari, "Structural relation between the vacuum space and the electron," *Physics Essays* 31, 108 (2018).
- `doc/tewari-references.md` for the reading notes and the status of H13.
