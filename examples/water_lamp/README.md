# The flashing neon lamp, with water instead of neon

The classic flashing-lamp (relaxation) oscillator: a supply charges a
capacitor through a resistor until the lamp strikes and empties it, then it
repeats. Neon lamps strike reliably because a trace of radioactive gas
(krypton-85, a beta emitter) supplies the first free electron; without it,
in the dark, the lamp waits for a cosmic ray (the "dark effect"). This
example puts water in place of the neon and the radioactive starter, three
ways, beside the neon lamp itself, then lets each one drive the Tesla
bifilar coil of [pulsed_bifilar](../pulsed_bifilar), tunes the flash
circuit to the coil as a spark-gap Tesla coil is tuned, and tries plates in
parallel. The misted gap is the §3j front end, and its faster recovery is
register hypothesis H5, carried here as an untested parameter.

```
go run ./examples/water_lamp -viz examples/water_lamp/viz   # about 3 minutes
python3 examples/water_lamp/viz/render.py
go test ./examples/water_lamp
KRONSPICE_EXPORT=1 go test ./examples/water_lamp -run TestExport   # the coil network's netlist and results sheet
```

## The construction

**The circuit.** A supply Vs charges C through R; the lamp sits across C.
When it strikes it discharges C down to its maintaining voltage through
its on-resistance and goes out when its current falls below the
extinction current; then C charges again. With no statistical wait the
period is the textbook R·C·ln((Vs − Vm)/(Vs − Vb)) plus the short on-time.

**The starting electron.** A gas lamp strikes only when a free electron is
there to start the avalanche. The wait is a random process: free
electrons appear at a declared rate S, and each one starts a discharge
with probability 1 − Vb/V once V is past breakdown (declared form), so
the mean wait at the supply voltage is 1/(S·(1 − Vb/Vs)). All randomness
comes from one seeded generator (CLAUDE.md §4), so runs repeat exactly.
Cosmic rays give about 2 ion pairs per cm³ per second in air at 1 atm
(declared), scaled with gas density.

**The lamps** (`lamps.go`; cited or declared values):

| Lamp | Supply, R, C | Strike | Maintain | Starting electrons |
|---|---|---|---|---|
| Neon NE-2, radioactive starter | 120 V, 10 MΩ, 10 nF | 76 V (TekWiki) | 58 V | 10⁴/s (declared) |
| Neon NE-2, no starter, dark | same | 76 V | 58 V | cosmic rays in ~0.03 cm³ at ~60 Torr: 0.0047/s |
| Water-vapour tube, dark | 600 V, 22 MΩ, 10 nF | 480 V, Paschen minimum of water vapour at 0.6 Torr·cm (Škoro 2011) | 0.76 × strike (NE-2 ratio, declared) | cosmic rays in 2 cm³ at 1.2 Torr: 0.0063/s |
| Water-vapour tube, UV LED starter | same | 480 V | 365 V | 10³/s photoelectrons (declared) |
| Spark gap in pure liquid water | 3,000 V, 1 MΩ, 10 nF | E·τ^(1/3) = 0.135 MV/cm·µs^(1/3) ± 12 % over 1 mm (Stygar 2006) | 50 V arc | none needed |
| Air spark gap, dry | 10,000 V, 1 MΩ, 10 nF | Paschen's law, 1 mm at 1 atm: 5,035 V (Raizer) | 30 V arc | 10³/s room-air ions; 1 ms recovery (declared) |
| Air spark gap with water mist | same | dry × 29.4/21.7 = 6,822 V (Liu 2018) | 30 V | as dry; recovery as dry |

Paschen's law for air is V = B·pd/ln(A·pd/ln(1 + 1/γ)) with A = 15 /(cm·Torr),
B = 365 V/(cm·Torr) (Raizer, Table 4.1 and eq. 7.3) and γ = 0.01 declared;
for water vapour the curve V = B·pd/ln(pd/x₀) is pinned to the measured
minimum, 480 V at 0.6 Torr·cm. The liquid gap's breakdown relation was
measured for pulses of 0.011 to 0.6 µs; the millisecond charging here is
far outside that range, so its strike voltage is an extrapolation and is
labelled so. Pure water conducts (0.055 µS/cm), so a leak current flows
all the time through 1 mm × 0.2 cm² and splits water: the ledger books
1.229 V × charge (Faraday's law) on the chemical channel and the rest as
heat. The fog factor is Liu's 5 cm rod–plane result applied to 1 mm, also
an extrapolation, labelled.

**The engine.** One RK4 integrator for Parts 1 and 2 on the state
[V, I_primary, I_coil, V_coil, supply energy, R heat, leak, lamp energy,
coil heat], with the energies integrated as states so every ledger closes
to RK4 accuracy. Part 2 drives the coil's first mode: Ls 4.853 mH, f1
41.61 kHz, Q 83 from `pulsed_bifilar/run_output.txt`, its 5-turn primary
Lp 5.97 µH computed from the same coaxial-circle formulas, and the
measured primary-to-coil mutual 28.1 µH (k = 0.165). A one-mode model;
`Twin` checks it against a kronspice transient of the same network.

## 1. Seven lamps in the same circuit

| Lamp | Supply | Strikes at | Flashes every | Energy per flash |
|---|---|---|---|---|
| Neon NE-2 with radioactive starter | 120 V | 78 ± 1 V | 0.039 s (±6 %) | 13.5 µJ |
| Neon, no starter, in the dark | 120 V | 120 V | 540 s ± 100 % (expected wait 576 s) | 55 µJ |
| Water-vapour tube, in the dark | 600 V | 600 V | 662 s ± 108 % (expected wait 792 s) | 1.13 mJ |
| Water-vapour tube, ultraviolet LED starter | 600 V | 492 ± 4 V | 0.174 s (±5 %) | 0.55 mJ |
| Spark gap in pure liquid water | 3,000 V | 1,076 ± 91 V | 4.4 ms (±11 %) | 5.8 mJ |
| Air spark gap, dry | 10,000 V | 6,798 ± 894 V | 11.8 ms (±27 %) | 0.235 J |
| Air spark gap with water mist | 10,000 V | 8,224 ± 643 V | 18.1 ms (±24 %) | 0.34 J |

- **Water-vapour tube**: low-pressure water vapour at its measured Paschen
  minimum. Without a starter it has the same dark-effect problem as neon:
  in 10 hours it flashes 54 times at the full supply voltage, about every
  11 minutes, at random. An ultraviolet LED freeing electrons from the
  cathode does the radioactive gas's job (10 flashes in 2 s, ±5 %).
- **Liquid water**: strikes near 1.08 kV with the 12 % scatter of the
  breakdown relation; the leak splits 0.17 mJ worth of water into
  hydrogen and oxygen in 2 s (7×10⁻¹⁰ mol H₂).
- **Water mist**: the fog raises the strike voltage by 1.36× and lengthens
  the interval. Faster recovery with mist (H5) would shorten it; it is a
  parameter here, untested.

Every ledger: supply = heat in R + lamp (light and heat) + leakage +
stored, residual ≤ 10⁻⁵ relative (≤ 10⁻¹⁰ J where no cosmic-ray run is
involved). Picture: `viz/lamps.png`.

## 2. Driving the bifilar coil

The reduced coil model agrees with the kronspice transient to 8×10⁻⁴ of
the peak.

| Lamp | Flashes | Coil peak | Energy to the coil | Share of flash energy |
|---|---|---|---|---|
| Neon, glow | 11 | 0.03 V | 1.3×10⁻¹¹ J | ≈0 |
| Water-vapour glow, UV starter | 7 | 0.05 V | 2.3×10⁻¹¹ J | ≈0 |
| Liquid-water spark | 41 | 44 V | 96 µJ | 0.2 % |
| Mist spark | 14 | 313 V | 1.9 mJ | 0.6 % |

Glow discharges are too slow and gentle; sparks ring the coil, but the
10 nF flash capacitor rings with the 6 µH primary at about 650 kHz,
16 times the coil's 41.6 kHz, so almost nothing reaches it. Picture:
`viz/coil.png`.

## 3. Tuning the flash circuit to the coil (`tuned.go`)

Lp·Cp = Ls·Cs, as in a spark-gap Tesla coil: the primary rings at the
coil's own frequency and the energy swings from primary to coil and back
at the beat rate, the first complete transfer after half a beat,
1/(2k·f1) = 72.8 µs with the 5-turn primary (Cp = 2.45 µF). An ordinary
arc keeps conducting through the current zeros until the current has
stayed small for its de-ionization time (20 µs declared), so the energy
swings back and burns in the arc; a quenched gap opens at the first
current zero after the transfer and leaves the energy in the coil. The
charging resistor is kept large enough that the supply cannot hold the
arc by itself (Vs/R ≤ I_ext/2), or the arc never goes out.

| Gap and circuit | Flashes | Spark current | To the coil per flash | Of the flash | Of the supply | Coil peak |
|---|---|---|---|---|---|---|
| Liquid water, untuned (10 nF, Part 2) | 41 | 48 A | 2.3 µJ | 0.2 % | 0.00 % | 44 V |
| Liquid water, tuned, ordinary arc | 47 | 173 A | 6.2 mJ | 4.4 % | 0.23 % | 2,401 V |
| Liquid water, tuned, quenched gap | 47 | 175 A | 6.2 mJ | 4.4 % | 0.23 % | 2,421 V |
| Mist, untuned (10 nF, Part 2) | 14 | 340 A | 0.14 mJ | 0.6 % | 0.01 % | 313 V |
| Mist, tuned (2.45 µF, 5 turns), ordinary arc | 5 | 3,490 A | 1.9 J | 3.3 % | 1.08 % | 68.5 kV |
| Mist, tuned, quenched gap | 5 | 3,486 A | 5.4 J | 9.5 % | 3.16 % | 68.4 kV |
| Mist, evolved: 0.249 µF, 20 turns, quench at 0.65 of the transfer time | | 438 A | 4.35 J | **74.7 %** | 25.5 % | 51.6 kV |

The genetic algorithm (`evo.Run`, 24 designs × 12 generations, seed 1)
chooses the flash capacitor (1 nF to 20 µF), the primary turns (1 to 20;
Lp from the turns' self and mutual inductances, M scaled with the turns,
declared) and the quench time (0.2 to 3 beat halves) for the mist gap,
minimizing 1 − share with at least three flashes in six charging time
constants; its ledger closes to 2.9×10⁻⁶. 20 turns is the top of the
allowed range. Charging a capacitor through a resistor from a fixed supply
loses half the supply energy in the resistor, so "of the supply" is at
most 50 %.

**Caution:** the evolved design rings the coil to about 52 kV, so
neighbouring bifilar turns carry about 26 kV, far beyond enamelled wire
(a few kV, declared). It would need much thicker insulation or a gentler
drive (§3e rule).

## 4. Plates in parallel (`parallel.go`)

| Several mist gaps on one capacitor | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| Strike voltage | 6,914 ± 29 V | 5,922 ± 2,394 V | 4,254 ± 2,726 V | 3,602 ± 2,452 V |
| Flash interval | 0.532 s ± 1 % | 0.468 s ± 60 % | 0.24 s ± 140 % | 0.15 s ± 192 % |
| Spark current | 438 A | 438 A | 438 A | 436 A |
| Energy to the coil per flash | 4.35 J (74.7 %) | 3.23 J (71.5 %) | 1.34 J (59.7 %) | 0.68 J (45.3 %) |

Only the first gap to break down carries the flash; more gaps give more
starting electrons, so a gap also re-strikes before the others have
recovered: smaller, more frequent flashes, and the spark current does not
add.

| Capacitor plates in parallel (primary retuned) | 1 | 2 | 4 | 8 |
|---|---|---|---|---|
| Capacitance; primary turns | 0.249 µF; 20 | 0.498 µF; 13 | 0.997 µF; 8 | 1.99 µF; 6 |
| Detuning | −0.5 % | −1.0 % | 4.7 % | −5.4 % |
| Spark current | 438 A | 846 A | 1,695 A | 2,818 A |
| Energy per flash | 5.8 J | 11.5 J | 23.6 J | 46.2 J |
| To the coil per flash | 4.35 J (75 %) | 6.66 J (58 %) | 7.85 J (33 %) | 8.3 J (18 %) |

More plates raise the energy per flash (×n) and the spark current nearly
as much; with fewer primary turns the coupling falls, so the share
reaching the coil drops, and a bigger bank needs the coil and primary
redesigned together.

Every joule comes from the supply; the water adds no energy.

## Tests

| Test | Checks |
|---|---|
| `TestPaschen` | water vapour's 480 V minimum at 0.6 Torr·cm; air at 1 atm, 1 mm near 4.5 to 5 kV |
| `TestRelaxationPeriod` | a primed lamp follows R·C·ln((Vs − Vm)/(Vs − Vb)); every ledger closes |
| `TestDarkWait` | the mean dark wait equals 1/(S·(1 − Vb/Vs)) |
| `TestCoil` | the reduced coil matches kronspice; every lamp's coil energy is on the ledger |
| `TestTuning` | tuning moves far more energy into the coil, a quenched gap more than an arc; ledgers close |
| `TestParallel` | parallel gaps do not add spark current; parallel plates do |

## Files

| File | What it holds |
|---|---|
| `lamps.go` | the lamps, Paschen curves, the flashing-lamp cases and the RK4 engine |
| `coil.go` | the bifilar coil's first mode, its primary and its kronspice twin |
| `tuned.go` | tuning, arc and quenched gaps, the genetic-algorithm search |
| `parallel.go` | parallel gaps and parallel capacitor plates |
| `main.go`, `report34.go` | the report and the `-viz` traces |
| `main_test.go` | the tests above |
| `export_test.go` | the coil network's netlist, KiCad schematic and results sheet in `export/` |
| `run_output.txt` | the run above |
| `viz/render.py`, `viz/lamps.png`, `viz/coil.png` | the pictures |

## Sources

- N. Škoro, D. Marić, G. Malović, W. G. Graham, Z. Lj. Petrović, "Electrical breakdown in water vapor," *Phys. Rev. E* 84, 055401(R) (2011): the Paschen minimum of water vapour.
- TekWiki, "Neon Lamps": NE-2 strike and maintaining voltages, the krypton-85 starter, the dark effect.
- W. A. Stygar et al., "Water-dielectric-breakdown relation for the design of large-area multimegavolt pulsed-power systems," *Phys. Rev. ST Accel. Beams* 9, 070401 (2006): E·τ^(1/3) = 0.135 ± 12 %.
- Yu. P. Raizer, *Gas Discharge Physics* (Springer, 1991), Table 4.1 and eq. 7.3: Paschen's law for air.
- Y. Liu, S. Jiang, Z. Zhong, J. Geng, F. Lv, "Alternating current discharge characteristics and simulation analysis of rod-plane short air gaps under salt fog conditions," *Energies* 11, 2251 (2018): fog raises the discharge voltage of 5 cm gaps from 21.7 to 29.4 kV.
- C. S. R. Bhoomireddy, *Recovery times of a spark gap with high pressure gases as insulation suitable for repetitive applications*, PhD thesis, HBNI/BARC (2016): spark-gap recovery.
- F. E. Terman, *Radio Engineers' Handbook* (1943), §3, coupled circuits: energy transfer between two tuned circuits in 1/(2k·f).
- CODATA 2018: the Faraday constant; 1.229 V reversible water-splitting voltage at 25 °C.
- The quenched gap (Tesla's quenched gap; Wien's quenched spark) is named in `tuned.go` without a citation.
