# Tesla's homopolar with the load current through its ball bearings

Tesla's two-disk homopolar generator (U.S. Patent 406,968,
[homopolar](../homopolar)) with its output current taken through the shaft
ball bearings instead of brushes, so that each bearing is also a Huber
ball-bearing motor (the P element, [bearing_motor](../bearing_motor);
H13, H14 in CLAUDE.md §1). Circuit simulation only: every part is a lumped
element and the shafts are mechanical nodes (across = ω, through =
torque). The question is the least drive torque per watt delivered to an
external resistor or capacitor.

```
go run ./examples/homopolar_bearing [-p 20] [-gens 25] [-pop 30]
go test ./examples/homopolar_bearing
```

## The trade-off and the measured inputs

The bearing bore sets the axle radius r_i, so bigger bearings remove
magnetic area from each disk, EMF = ½Bω(r_o² − r_i²); but bigger, more
heavily loaded bearings carry more current before they are damaged, and
have more friction. No coefficient is assumed; unmeasured ones are swept.

- **Bearing voltage under DC:** almost constant, 0.73 V average at 4 A
  (0.5–1.5 V reported), from repeated breakdown of the oxide and lubricant
  film (Safdarzadeh, Farahi, Binder et al., *Lubricants* 12, 148, 2024;
  greased 51208 bearings, 40 mm bore). The P element carries it as VF.
- **Damage limit:** 0.1–0.3 A/mm² of Hertzian contact area, A_Hz = 1.0 mm²
  at 400 N and 2.5 mm² at 1600 N for that bearing (same paper); other
  sizes scale by Hertz theory, A ∝ (F·R)^(2/3). `TestHertzArea` recovers
  the paper's second point from the first within 2 %.
- **Friction torque** M = ½·f·F·d with f = 0.0015 (SKF simplified model).
- **Rim-speed limit** of a solid copper disk: hoop stress (3+ν)/8·ρv² at
  the yield of annealed copper, 69 MPa, with a safety factor of 1.5
  (Timoshenko; ASM Handbook): 111 m/s.

Swept: the Huber coefficient κψ in EM and thermal modes, the belt's
rolling-contact resistance, and the shaft collector's resistance for the
insulated-bearing comparison.

## Results (15 cm copper disks, 0.8 T, from `run_output.txt`)

**Reference unit**, 40 mm bearings × 2 per shaft at 400 N, 600 rad/s, 2 Ω
load: EMF 10.608 V, 4.571 A, 41.80 W into the load from 77.29 W at the
shaft (efficiency 0.541, 0.00308 N·m/W). The bearing current is 7.62× the
0.3 A/mm² damage limit.

**Bore against magnetic area**, load matched to the bearing bank's current
limit:

| Bore | Area kept | EMF | I_max | P_load | N·m/W |
|---|---|---|---|---|---|
| 10 mm | 99.9 % | 10.788 V | 0.238 A | 2.22 W | 0.0073 |
| 40 mm | 98.2 % | 10.608 V | 0.600 A | 5.49 W | 0.0107 |
| 120 mm | 84.0 % | 9.072 V | 1.248 A | 9.50 W | 0.0171 |

At the damage limit a bearing bank passes 0.2–1.3 A, so a through-bearing
unit is a watts-to-tens-of-watts machine.

**Huber torque, κψ swept 0 to 0.1 N·m/A²:** in EM mode it lowers the drive
torque only by lowering the output through its back-EMF (24.7 W at 1e-3,
2.6 W at 1e-2), so it never improves the torque per watt; in thermal mode
it recovers at most η·H, about 10 % of the bearing heat, leaving 41.8 W at
0.0031 N·m/W. `TestNoFreeTorque` checks that the shaft power is never
below the load plus the measured losses.

**Bearing voltage range and the insulated alternative:** through the
bearings at 0.5, 0.73 and 1.5 V each: 46.1, 41.8 and 28.9 W (efficiency
0.578, 0.541, 0.418); with the bearings insulated and a shaft collector of
10 µΩ to 1 mΩ: 56.2 W at efficiency 0.661, independent of the collector.

**Capacitor load**, 1 F through 2 Ω: electrical conversion 42.8 % with a
stepped field, 65.6 % with a 20 s ramp, 73.3 % with 100 s, against the
75.8 % ceiling set by the bearing drops (2·V_b·Q whatever the ramp), while
friction runs throughout.

**Evolved** for at least 20 W with the least drive torque per watt (genes:
bore, preload, thickness, speed up to the rim limit, load, bearings per
shaft; 702 evaluations each):

| | Through the bearings | Bearings insulated |
|---|---|---|
| Design | 108 mm × 6 per shaft, 739 rad/s, 4.84 Ω | 10 mm × 1 per shaft, 739 rad/s, 0.256 Ω |
| Load | 21.1 W at 2.09 A (1.00× the current limit) | 689 W at 51.9 A |
| Efficiency | 0.134 | 0.998 |
| Drive torque per watt | 0.0101 N·m/W | 0.00135 N·m/W (7.4× less) |

The insulated unit's torque per watt is 1/(η·ω); in every case the drive
torque per load watt is at least 1/ω. Ledger residual 0 at every operating
point; transient mismatch below 1e-4 of the shaft energy.

**Bench data still needed:** the Huber coefficient κψ, the belt's
rolling-contact resistance, and the real collector's contact drop (a
carbon brush drops about 1 V, which would narrow the insulated design's
lead; a liquid-metal collector would not).

## Files

| File | What it holds |
|---|---|
| `main.go` | `Design`, measured constants, bearing bank, trade-off sweeps, capacitor charge, GA problem |
| `main_test.go` | Hertz scaling, no free torque, bore monotonicity, capacitor conversion ceiling |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt` | the run above |

## Sources

- O. Safdarzadeh, A. Farahi, A. Binder et al., *Lubricants* 12(5), 148 (2024), DOI 10.3390/lubricants12050148.
- K. L. Johnson, *Contact Mechanics*, Cambridge 1985, §4.2; T. A. Harris, *Rolling Bearing Analysis*.
- SKF rolling-bearing catalogue, simplified friction model.
- S. Timoshenko, *Strength of Materials* II, §44; ASM Handbook vol. 2 (C10100 annealed).
- N. Tesla, U.S. Patent 406,968 (1889); `doc/bearing-motor.md`.
