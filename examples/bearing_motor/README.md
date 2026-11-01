# Ball-bearing (Huber) motor: the measurement registry and the SPG

A current through a rotating ball bearing drives the shaft with a torque
that scales as I², aids the rotation and is zero at rest. kronspice carries
it as the P element (`device/bearing.go`, `doc/bearing-motor.md`). This
example holds the public record of what has been measured about such
motors, the fits that would set the element's coefficients, and one
question: could bearing torque, aiding the drive motor, account for the
Space Power Generator's published figures
([space_power_generator](../space_power_generator), H13)? Marinov's claim
of heat beyond the input is register item H14 (CLAUDE.md §1).

```
go run ./examples/bearing_motor
go test ./examples/bearing_motor
```

## The registry

Every number is a record with its source and a status: *measured*
(reported from an experiment), *estimate* (no instrument named) or *claim*
(asserted, no data shown). Only measured values enter a calculation.

| Finding | Source | Status |
|---|---|---|
| Steady currents 43.5–70.15 A; self-rotation at high current, small efficiency; a drive motor needed at low current | Moyssides & Hatzikonstantinou, *IEEE Trans. Magn.* 26, 1274 (1990), abstract | measured |
| Back-EMF that depends on shaft speed; no generator action | same | measured |
| Self-rotation 70–140 A, aided 18–140 A (two discs) | Moyssides 1992, via Lauterbach, Soong & Abbott (2001) | measured, secondary |
| Torque ∝ I² | Choo, Soong & Abbott, *Proc. SPIE* 5649 (2005), per McDonald (2020) | measured |
| Mechanical output / electrical input = 0.10 | Marinov, *Electronics & Wireless World*, April 1989 (calorimeter) | measured, data unpublished |
| Heat 8 % higher with rotation than without | Marinov 1989 | claim (H14), not used |
| 1000 r/min (Milroy 1959); "at least 1500 rev/min" (Ogden 1989) | Milroy, *J. Appl. Mech.* 34, 525 (1967); Ogden 1989 | measured / estimate |
| Needs a hand spin to start; thermal effects start it, the Lorentz force cannot | Harrison (electricstuff); McDonald (2020) fn. 6 | measured |
| Bearing temperature 80 °C at 30 A | Univ. Adelaide 2022 | measured |

## Coefficients: none invented

The element's torque law is T = κ·I²·tanh(ω/ω0). `FitTorque` fits κ and
ω0 to a measured (I, ω, T) table and `FitBackEMF` to a measured back-EMF
table; in EM mode the two must agree, since E·I = ω·T. No such table is
open access, so `Datasets` is empty, the program lists the nine papers
whose tables would fix the coefficients (Gruenberg 1978, Moyssides 1990 and
1997, Watson 1993 and 1999, Choo 2005, Shen 1999, Polivanov 1973, Huber
1959), and the element refuses to run without fitted KAPPA and W0.
`TestFitRecovers` checks the fitting code on synthetic points
(κ and ω0 recovered to 1e-6), not a bearing.

## Could bearing torque explain the SPG?

All inputs are measured: the SPG test (Tewari, UPSM pp. 159–160: 2700
r/min, 3438 A, no-load EMF 3.41 V, drive motor 2943 W → 8620 W at
η = 0.83) and the bearing efficiency above.

| | Value |
|---|---|
| Shaft torque the claimed E0·I needs | 41.5 N·m |
| Shaft torque the drive motor supplied | 16.7 N·m |
| Gap a bearing torque would have to fill | 24.8 N·m, 7.01 kW |
| κ·ψ(ω) needed if all 3438 A cross the bearings | ≥ 2.1×10⁻⁶ N·m/A² |

- **EM mode** (the measured mode: a back-EMF exists). The bearings then
  drop 2.04 V of the 3.41 V, and the rest of the loop gets 4.71 kW, not
  11.72 kW. The 1.37 V left is exactly the loaded-EMF bound of the SPG
  ledger: a bearing motor moves shaft work around and adds none.
- **Thermal mode** at the measured η = 0.10: the bearings would need 70 kW
  of electrical input, 20.4 V across them, more than the 3.41 V the
  machine makes. Excluded.

The engine run (`SimulateSPGLoop`) puts the no-load EMF, the P element at
the SPG's speed with the derived κψ and the rest of the loop sized for
3438 A into one circuit, and sweeps the unknown contact resistance:

| Rb | Load | Bearing shaft work | E0·I | Ledger residual |
|---|---|---|---|---|
| 1 µΩ | 4.70 kW | 7.01 kW | 11.72 kW | 0 |
| 10 µΩ | 4.59 kW | 7.01 kW | 11.72 kW | 0 |
| 100 µΩ | 3.53 kW | 7.01 kW | 11.72 kW | 0 |

The load never gets more than (E0 − back-EMF)·I, whatever Rb.

**Decisive bench tests:** insulate the SPG bearings (ceramic balls or
insulated housings) and repeat the loaded run; measure the voltage across
each bearing under load; put a torque transducer between the drive motor
and the rotor.

## Files

| File | What it holds |
|---|---|
| `main.go` | registry, `Needed` papers, `FitTorque` / `FitBackEMF`, `CheckSPG`, `SimulateSPGLoop` |
| `main_test.go` | fit recovery, registry sourcing, the SPG bound (24.8 N·m, 1.37 V) and the swept engine run |
| `export_test.go` | netlist, KiCad schematic and results sheet in `export/` |

## Sources

- J. Huber, *Elektrotechnik und Maschinenbau* 76(8), 169–174 (1959).
- R. A. Milroy, *J. Appl. Mech.* 34, 525 (1967).
- H. Gruenberg, "The ball bearing as a motor," *Am. J. Phys.* 46, 1213 (1978).
- P. G. Moyssides, P. Hatzikonstantinou, *IEEE Trans. Magn.* 26, 1274 (1990) and 33, 4566 (1997).
- D. B. Watson et al., *IEE Proc. A* 140 (1993); *IEE Proc. Sci. Meas. Technol.* 146, 83 (1999).
- J. L. Choo, W. L. Soong, D. Abbott, *Proc. SPIE* 5649, 700 (2005).
- S. Marinov, "The intriguing ball-bearing motor," and F. Ogden, "Great balls of fire!," *Electronics & Wireless World*, April 1989.
- K. T. McDonald, "Ball-Bearing Motor," Princeton (2020).
- P. Tewari, *Universal Principles of Space and Matter*, 2017, §6.7 (`reference/tewari_upsm_2017.pdf`).
