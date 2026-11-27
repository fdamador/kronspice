# Pulse-driven resonance in ordinary and Tesla bifilar coils

One-sign voltage pulses drive seven coils of the same 0.5 mm copper wire
(about 126 m each) through a 5-turn primary, and the example finds the
pulse rates at which the energy builds up into a large standing wave
(CLAUDE.md §3e step 3, "for pulse drive, the equivalent condition is
repeating the pulse at the round-trip time"). Three of the seven are the
bifilar solenoid wired other ways: its strands in series with opposing
currents, and in parallel with currents the same way and opposing
([Opposing currents, series and parallel](#opposing-currents-series-and-parallel)).
An AC sweep gives volts and amps against frequency for the five
solenoid windings. A secondary pickup then
takes power off at low voltage and high current. The bifilar solenoid's
first mode is the load that [water_lamp](../water_lamp) drives with its
spark gaps, and the same coil formulas serve [coiled_tether](../coiled_tether)
and [bial_bifilar](../bial_bifilar).

```
cd examples/pulsed_bifilar
go run .                                  # about 15 minutes, writes viz/*.csv (relative to the working directory)
go run . -only-pickup                     # the pickup study alone
go run . -only-vi                         # the V–I sweep alone (viz/vi_sweep.csv)
python3 viz/render.py                     # pictures and animations
python3 viz/make_report.py                # pulsed_bifilar_report.pdf
python3 viz/make_current_report.py        # pulsed_bifilar_current_report.pdf (the four connections)
go test ./examples/pulsed_bifilar
KRONSPICE_EXPORT=1 go test ./examples/pulsed_bifilar -run TestExport   # AC netlists, schematics, results sheets
```

## The construction

**Wire and geometry** (declared in `main.go`): 0.5 mm copper, 0.7 mm
between neighbouring wire centres (0.2 mm of enamel, ε = 3 declared), and
windings cut so that each holds about the same length of wire as the
400-turn solenoid:

| Coil | Turns | Radius | Height | Wire | Surroundings |
|---|---|---|---|---|---|
| Ordinary cylindrical solenoid | 400 | 5.0 cm | 27.9 cm | 125.7 m | grounded can of twice its radius |
| Tesla bifilar solenoid (US 512,340) | 400 | 5.0 cm | 27.9 cm | 125.7 m | same can |
| Tesla bifilar cone, 45° opening upward | 248 | 2.0 to 14.2 cm | 12.2 cm | 126.4 m | grounded plane 5 cm below the base |
| Tesla bifilar flat (Archimedean) spiral | 214 | 2.0 to 16.9 cm | 0 | 127.1 m | grounded plane 5 cm below |
| Series bifilar solenoid, opposite currents | 400 | 5.0 cm | 27.9 cm | 125.7 m | same can |
| Parallel bifilar solenoid, same direction | 400 | 5.0 cm | 27.9 cm | 2 × 62.8 m | same can |
| Parallel bifilar solenoid, opposite currents | 400 | 5.0 cm | 27.9 cm | 2 × 62.8 m | same can |

In a Tesla bifilar coil the wire runs up the even positions (strand 1);
its end is joined to the start of strand 2, which runs up the odd
positions between them in the same sense of rotation, so both strands
carry the current the same way round, neighbouring turns belong to the
two strands and sit half a wire length apart in potential (the patent's
V/2). The other three solenoid windings use the same positions:

- **Series, opposite currents** (`NewCounterCoil`): the wire folded in
  half and the doubled wire wound side by side, so strand 2 comes back
  down the odd positions with its current circulating the other way (the
  non-inductive winding of wire-wound resistors). Neighbouring turns are
  n − 1 − 2k wire turns apart: almost the whole coil voltage at the base,
  none at the fold.
- **Parallel, same direction** and **parallel, opposite currents**
  (`NewParallelCoil`): both strands start at the base terminal and end at
  the far terminal, strand 2 wound the same way or the other way round.

Each turn carries a sign, +1 or −1 for its sense of circulation, and the
turn-to-turn mutual inductances take the product of the signs. kronspice's
K element takes only positive coupling, so `Build` gives each section on a
returning strand the reversed reference direction (the dot convention);
the wire current is unchanged. The parallel windings' terminal inductance
is (L₁L₂ − M²)/(L₁ + L₂ − 2M) from the strands' block sums, and
`PrimaryM` gives their coupling to the primary at the terminals and round
the internal loop up one strand and down the other.

**The distributed network.** Each coil is 20 sections, not a lumped L
and C. Every turn is a circle: the inductance matrix of the sections is
summed turn by turn from Maxwell's formula for coaxial circles
(`em.MutualCircles`, *Treatise* art. 701, complete elliptic integrals)
with a single turn's high-frequency self-inductance μ₀R(ln(8R/r_w) − 2)
(Grover 1946), and the sections are coupled in kronspice with K
elements. Capacitance between physically adjacent turns is the two-wire
formula π·ε₀ε_r·πD/acosh(p/d) with each pair's mean diameter (12.1 nF
summed over all neighbouring turns); capacitance to the surroundings
treats the winding as a sheet facing the nearest grounded surface
(2πε₀/ln(R_can/R) per metre for the can, ε₀·area/height per ring for the
plane), split half to each end node of a section. Wire resistance is
Terman's skin-effect value at f1. A 1 GΩ bleeder gives the floating coil
a DC reference.

**The drive.** The pulses reach the coil only through magnetic coupling,
as in a Tesla coil: the source (100 V, 50 Ω) drives a 5-turn primary of
1 mm wire at 2 mm pitch just below the base (radius 6 cm around the
solenoid, 1 cm outside the base of the cone and the spiral), coupled to
every section by the same coaxial-circle mutuals. Pulses are PULSE
sources with edges 2 % and width 10 % of 1/f1, repeated at rates r/f1 =
1/3, 1/2, 0.9, 0.97, 1, 1.03, 1.1, 1.5 and 2, for min(600, max(60,
5Q/π + 20)) cycles of f1 (about five ring-up time constants Q/(πf1)) at
50 steps per period. The first resonance f1 and its Q come from an AC
sweep of the terminal voltage per volt of drive, from 1 kHz up to the
frequency at which one section's wire is λ/10 (4.77 MHz for the
solenoids): the first peak with gain above 2, or, for the
opposite-current windings, which have none, the strongest peak, used as
the pulse rate and marked "(peak)".

**Approximations declared:** sections are equipotential inside; the
bifilar junction lead and radiation are left out; proximity effect is not
in the wire resistance; the absolute tolerances are loosened to the
circuit's scale (1 µA, 1 mV, 1 pC) so the step control does not chase
picoamp ringing inside the sectioned network.

## The coils

| Coil | f1 | Wire / λ | AC gain at top | Q | Inductance | Pulse width | Cycles per run |
|---|---|---|---|---|---|---|---|
| Ordinary solenoid | 1.167 MHz | 0.489 | ×338.6 | 210 | 4.853 mH (Nagaoka sheet 4.873 mH) | 85.7 ns | 354 |
| Bifilar solenoid | 41.61 kHz | 0.017 | ×12.1 | 83 | 4.853 mH | 2.40 µs | 151 |
| Bifilar cone | 35.81 kHz | 0.015 | ×6.6 | 96 | 6.542 mH | 2.79 µs | 172 |
| Bifilar flat spiral | 34.93 kHz | 0.015 | ×10.8 | 94 | 6.842 mH | 2.86 µs | 169 |
| Series bifilar, opposite currents | 593.0 kHz (peak) | 0.249 | ×0.559 | 71 | 14.70 µH | 169 ns | 132 |
| Parallel bifilar, same direction | 1.995 MHz | 0.836 | ×233.3 | 238 | 1.213 mH | 50.1 ns | 398 |
| Parallel bifilar, opposite currents | 528.1 kHz (peak) | 0.221 | ×0.105 | 148 | 3.674 µH | 189 ns | 255 |

The two solenoids have the same inductance to 10⁻¹² (the strands of a
series bifilar coil carry the same current the same way), within 0.4 %
of Nagaoka's current-sheet formula; the flat spiral is within 1.7 % of
Wheeler's flat-spiral formula. The ordinary solenoid's wire is 0.489
wavelengths at f1: it rings as a free half-wave line, the first resonance
Knight reports for a free single-layer coil (§3e), with the voltage
maxima at its ends and the current maximum at its middle. The bifilar
windings, with 12 nF between their strands, resonate 28 to 33 times lower.
The other three solenoid windings are described below.

## Build-up only in step with resonance

A one-sign pulse train at rate r contains only the frequencies r, 2r,
3r, …, so it builds the coil up only when one of these lands on f1: at
r = f1, f1/2, f1/3. Coil-end voltage after the run (100 V pulses from
50 Ω):

| Coil | r = f1/3 | f1/2 | 0.97 f1 | f1 | 1.03 f1 | 1.5 f1 | 2 f1 |
|---|---|---|---|---|---|---|---|
| Ordinary solenoid | 2,632 V | 3,868 V | 893 V | **7,426 V (×74)** | 657 V | 480 V | 356 V |
| Bifilar solenoid | 97 V | 143 V | 55 V | **282 V** | 59 V | 7 V | 5 V |
| Bifilar cone | 53 V | 79 V | 27 V | **155 V** | 28 V | 4 V | 3 V |
| Bifilar flat spiral | 87 V | 128 V | 44 V | **252 V** | 46 V | 6 V | 5 V |
| Series bifilar, opposite currents | 5.1 V | 7.4 V | 4.9 V | **13.9 V** | 5.4 V | 2.0 V | 2.0 V |
| Parallel bifilar, same direction | 1,862 V | 2,665 V | 444 V | **4,918 V** | 398 V | 78 V | 54 V |
| Parallel bifilar, opposite currents | 1.6 V | 2.2 V | 1.3 V | **2.6 V** | 1.0 V | 1.8 V | 1.3 V |

Between neighbouring turns at r = f1: 38 V in the ordinary solenoid, 141 V
in the bifilar solenoid, 78 V in the cone and 127 V in the spiral, half
the coil voltage for every bifilar winding, as the patent says. The
ledger shows every joule in the coil came from the pulse source:

| Coil, at f1 | In | Heat | Stored | Time-stepping error |
|---|---|---|---|---|
| Ordinary solenoid | 2.352 mJ | 2.204 mJ | 0.132 mJ | 1.6×10⁻⁵ J |
| Bifilar solenoid | 79.82 mJ | 79.71 mJ | 0.115 mJ | 1.2×10⁻⁶ J |
| Bifilar cone | 108.3 mJ | 108.3 mJ | 0.035 mJ | 1.4×10⁻⁷ J |
| Bifilar flat spiral | 108.5 mJ | 108.4 mJ | 0.093 mJ | 9.9×10⁻⁸ J |
| Series bifilar, opposite currents | 2.526 mJ | 2.445 mJ | 0.075 mJ | 6.1×10⁻⁷ J |
| Parallel bifilar, same direction | 1.699 mJ | 1.622 mJ | 0.073 mJ | 4.0×10⁻⁶ J |
| Parallel bifilar, opposite currents | 6.204 mJ | 6.079 mJ | 0.125 mJ | −2.6×10⁻⁶ J |

The bifilar coils take far more energy from the source (their wide pulses
at low f1 spend most of it as heat in the primary circuit) and store less:
their low-frequency, low-Q mode builds up less than the ordinary coil's
half-wave resonance at 1.17 MHz.

## Opposing currents, series and parallel

The same 400-turn bifilar solenoid wired four ways, with pulses at each
winding's own f1 (peaks in the last 10 % of the run, mean power over the
run; the pulse width is 0.1/f1, so low-f1 windings take wider pulses):

| Winding | L | M to primary | f1 | Coil ends | Neighbouring turns | Wire current | Primary current | Power in | Into heat |
|---|---|---|---|---|---|---|---|---|---|
| Series, same direction (Tesla) | 4.853 mH | 28.14 µH | 41.61 kHz | 282 V | 141 V | 0.428 A | 1.97 A | 22.0 W | 22.0 W |
| Series, opposite currents | 14.70 µH | 0.269 µH | 593.0 kHz (peak) | 13.9 V | 15.0 V | 0.514 A | 1.39 A | 11.3 W | 11.0 W |
| Parallel, same direction | 1.213 mH | 14.07 µH | 1.995 MHz | 4,918 V | 26.8 V | 0.405 A | 1.03 A | 8.51 W | 8.13 W |
| Parallel, opposite currents | 3.674 µH | 0.135 µH (loop 28.14 µH) | 528.1 kHz (peak) | 2.6 V | 209 V | 0.497 A | 1.55 A | 12.8 W | 12.6 W |
| Ordinary (reference) | 4.853 mH | 28.14 µH | 1.167 MHz | 7,426 V | 37.9 V | 0.453 A | 0.903 A | 7.75 W | 7.26 W |

- **Series, opposite currents.** The fields cancel: L falls 330×, the
  coupling to the primary 100×, and there is no resonance with gain above
  2. The current still flows (0.51 A), but the coil holds almost no
  field and builds no voltage.
- **Parallel, same direction.** L is exactly a quarter of the series
  winding's (for identical strands (L₁² − M²)/(2(L₁ − M)) = (L₁ + M)/2
  against 2(L₁ + M) in series). The strands meet at the same point of
  the wire, so the strand capacitance holds no voltage: f1 rises to
  1.995 MHz, 1.7× the ordinary coil and 48× the series bifilar, with the
  highest Q (238). This is the way to raise a bifilar coil's resonance.
- **Parallel, opposite currents.** L at the terminals is 3.7 µH and the
  terminals show 2.6 V, but current can run up one strand and down the
  other: a loop whose turns all circulate the same way, a shorted coil
  coupled to the primary by the full 28.14 µH. Its EMF appears between
  neighbouring turns, 209 V, the largest turn-to-turn voltage of any
  winding, with nothing at the terminals to show it.
- **No connection adds energy**: each ledger closes, and each winding
  only changes how L and C share the stored energy and where the current
  flows.

**The two ends joined** (`Coil.Short`, `short_test.go`): with the
Tesla bifilar solenoid's ends joined through 1 mΩ and pulses at its
open-ends f1 (41.61 kHz), the coil becomes a shorted transformer
secondary. Its current is (M/L)·I_primary = 0.0058 × I_primary, set by
the flux linkage M·I_p + L·I = 0 and not by any resonance (the 41.6 kHz
resonance disappears). Over 151 periods:

| Ends | Coil ends | Neighbouring turns | Wire current | Current through the join | Stored at the end | Power in |
|---|---|---|---|---|---|---|
| Open | 282 V | 141 V | 0.43 A | — | 115 µJ | 22.0 W |
| Joined (1 mΩ) | 29 µV | 1.5 V | 0.036 A | 0.029 A | 4.6 nJ | 22.3 W |
| Joined (1 Ω) | 28 mV | 1.5 V | 0.035 A | 0.028 A | 4.6 nJ | 22.3 W |

Joining the ends lowers the voltage **and** the current (12× less wire
current, 25,000× less stored energy). The source still delivers about
22 W, but almost all of it goes to heat in its own 50 Ω and the primary:
the shorted coil reflects the drive instead of storing it. Volts are
traded for amps only at equal power or energy (P = V·I, ½LI² = ½CV²), by
an impedance transformation such as the tuned pickup below (3.8 A at
0.77 V from a coil ringing at 3.6 kV), not by shorting. The ordinary
solenoid joined the same way rings at its closed-loop resonance
(2.007 MHz) with 0.32 A in the wire and 0.071 A through the join.

Why: an open coil at resonance is a high-impedance source. Its Thevenin
impedance at 41.61 kHz is V_open / I_joined = 12.05 V / 0.115 mA =
105 kΩ per volt of drive, which is the tank's Q·ωL (83 × 1,271 Ω). Joining
the ends draws V_open / 105 kΩ. The 18 mA that circulates inside the open
coil at resonance is reactive current sloshing between its L and its strand
C; it is not available at the terminals.

**Checked against ngspice** (`twin_ngspice_test.go`, `viz/ngspice_twin.py`;
ngspice-39 in the `kronspice-tools` image): the same netlists, ends open
and joined. AC at 10 kHz, 41.61 kHz and 100 kHz agrees to ≤ 1×10⁻⁵ (the
precision ngspice prints); the pulse transient agrees to 0.4–0.5 % open and
3.6–5 % joined (µV across the 1 mΩ join, different time steps), with the
same conclusion in both: joining the ends lowers voltage and current.

```
KRONSPICE_NGSPICE_TWIN=/tmp/twin go test ./examples/pulsed_bifilar -run TestNgspiceTwinShort
docker run --rm -v /tmp/twin:/w -w /w kronspice-tools sh -c 'for f in open_ac joined_ac open_tran joined_tran; do ngspice -b $f.cir > $f.log; done'
python3 examples/pulsed_bifilar/viz/ngspice_twin.py /tmp/twin
```

**Volts and amps against frequency** (`visweep.go`, `viz/vi_sweep.csv`,
`viz/vi_frequency.png`): an AC sweep per volt of drive on the primary,
1 kHz to 4.77 MHz, for the five solenoid windings.

| Winding | Peak | Coil ends | Wire current | V/I |
|---|---|---|---|---|
| Ordinary | 1.168 MHz | 297 V | 12.2 mA | 24.3 kΩ |
| Series, same direction (Tesla) | 41.69 kHz | 6.42 V | 9.64 mA | 666 Ω |
| Series, opposite currents | 592.2 kHz | 0.422 V | 18.8 mA | 22.4 Ω |
| Parallel, same direction | 1.995 MHz | 241 V | 14.7 mA | 16.4 kΩ |
| Parallel, opposite currents | 527.8 kHz | 0.0892 V | 17.7 mA | 5.04 Ω |

Lowering f1 by adding capacitance at fixed L leaves the current at
resonance about the same and lowers the voltage: at equal stored energy
½LI² = ½CV², I is unchanged and V = I·√(L/C) falls by the square root of
the capacitance ratio (Tesla bifilar against ordinary: about 800× the C,
28× lower f1, 46× lower V per volt of drive, 0.8× the current). Compared
at equal voltage instead, the bifilar carries √(C ratio) more current and
stores the C ratio more energy. Below 1 MHz both opposite-current
windings show a comb of wire-current peaks (largest 24.2 mA series,
20.9 mA parallel, per volt of drive, against 9.6 mA for the Tesla
bifilar at its resonance) with almost no voltage at the terminals: modes
of the strand pair whose fields cancel outside. Above about 1.5 MHz every
winding has modes of its own, up to the λ/10 limit. The full set of volts,
amps and power is in `pulsed_bifilar_current_report.pdf`.

## Low-voltage, high-current pickup (`secondary.go`)

A 4-turn pickup of 2 mm copper (radius 6.5 cm, inside the can; 4.082 µH,
0.072 Ω at f1) wound round the ordinary coil's middle, where the current
is largest, tuned with a series capacitor to the coil's resonance,
C = 1/((2πf1)²L) = 4.56 nF, so that only the load and the wire resistance
remain in the pickup loop. Its coupling to every section and to the
primary comes from the same coaxial-circle formulas. The load (and, where
stated, the capacitor) is chosen from a grid for the largest load power
per V² of drive, and the pulse rate is the coupled system's best
frequency, searched from 0.7 to 1.4 f1.

| Pickup | C | Load | Rate / f1 | Load current | Load voltage | Load power | Past the source's 50 Ω | Share | Coil |
|---|---|---|---|---|---|---|---|---|---|
| Middle, tuned 1/(ω1²L) | 4.56 nF | 0.2 Ω | 1.2975 | 3.83 A | 0.77 V | 1.40 W | 2.11 W | 66.5 % | 3,651 V |
| Middle, capacitor optimized | 1.14 nF | 10 Ω | 0.9717 | 0.67 A | 6.67 V | 1.20 W | 1.31 W | 91.3 % | 3,442 V |
| Middle, untuned | none | 1 Ω | 1.1073 | 1.70 A | 1.70 V | 1.31 W | 1.64 W | 79.7 % | 3,502 V |
| Top (high-voltage end), tuned | 4.56 nF | 0.5 Ω | 1.0487 | 2.04 A | 1.02 V | 1.01 W | 1.16 W | 86.8 % | 2,230 V |

The pickup converts voltage into current; its load never gets more than
the power that reaches the primary past the source's own resistance
(ledgers close to 10⁻⁵ J). A tightly coupled tuned pickup splits the
coil's resonance in two (about 0.83 and 1.30 f1), so the best pulse rate
moves with it.

## Pictures

`WriteViz` writes each coil's geometry (with each turn's place along the
wire and sense of circulation), voltage envelope, peak and ring-up
"glow" data, and its currents (`*_ienv.csv`: wire and primary current per
period; `*_current.csv`: each section's current at the moment of largest
current) to `viz/*.csv` from the full-resolution run at f1;
`viz/render.py` turns them into:

| Picture | Shows |
|---|---|
| `coil_voltage.png` | the ordinary and Tesla bifilar solenoids, every turn coloured by its voltage at the moment of largest voltage, with the voltage along the wire and between neighbouring turns |
| `coil_voltage_shapes.png` | the same for the cone and the spiral |
| `coil_voltage_counter.png`, `coil_voltage_parallel.png` | the same for the series and parallel bifilars, currents the same way and opposing |
| `coil_current.png`, `coil_current_parallel.png` | every turn coloured by the current circulating in it (red one way round, blue the other), with each section's current along the wire and in the sense of its turns |
| `coil_amps.png` | wire and primary current building up pulse by pulse at f1, all seven coils |
| `coil_vip.png` | volts, amps and power at f1, all seven coils side by side |
| `vi_frequency.png` | V, I and V/I against frequency for the five solenoid windings (AC, per volt of drive) |
| `pulse_sweep.png` | coil-end voltage against rate / f1 |
| `ring_up.png` | the coil-end voltage building up pulse by pulse |
| `ring_up.gif`, `ring_up_shapes.gif` | the coils glowing as the pulses build them up |

`python3 viz/render.py viz run_output.txt vi_frequency coil_vip` redraws
only the named pictures.

## Tests

| Test | Checks |
|---|---|
| `TestCoilNetwork` | every winding uses each position once; L within 1 % of Nagaoka and equal for both solenoids; symmetric capacitance matrix; the bifilar coil holds (N − 1) strand-to-strand gaps |
| `TestFlatSpiralInductance` | the spiral against Wheeler's a²N²/(8a + 11c), held to 10 % |
| `TestPulseBuildUp` | pulses at f1 build more than three times higher than at 1.5 f1; ledger closes to 1 % |
| `TestPickup` | the tuning capacitor resonates the pickup at f1 to 10⁻¹²; ledger to 1 %; load power below the power past the source |
| `TestCounterCoil` | series opposite-current winding: L below 1 % and primary coupling below 2 % of the Tesla winding's; no section spans the fold, so every K is positive in `Build`'s reference directions; (N − 1) strand gaps; no resonance above gain 2; ledger to 1 % |
| `TestParallelCoils` | parallel windings: same direction L = ¼ of the series winding's to 10⁻³ and f1 above the ordinary coil's; opposite currents L below 1 % and no resonance above gain 2; ledgers to 1 % |
| `TestShortedEnds` | ends joined through 1 mΩ: loop current = (M/L)·I_primary at 10 kHz to 1 % (ordinary and Tesla bifilar); pulsed at the bifilar's f1, wire current below a fifth and stored energy below a thousandth of the open coil's; ledgers to 1 % |
| `TestNgspiceTwinShort` | with `KRONSPICE_NGSPICE_TWIN=<dir>`: writes the open and joined Tesla bifilar as ngspice netlists (AC and pulse transient) and kronspice's values beside them, for `viz/ngspice_twin.py` |
| `TestExport` | the four coils of the first study and the series opposite-current coil, swept around resonance, and the tuned pickup, as netlists and results sheets |

## Files

| File | What it holds |
|---|---|
| `main.go` | geometries, the seven windings (series and parallel, same and opposite currents), the sectioned network, resonance sweep, pulse runs with volts, amps and power, ledger, picture data |
| `visweep.go` | the AC V–I sweep of the five solenoid windings (`-only-vi`) |
| `secondary.go` | the tuned pickup, its AC search and pulse runs |
| `main_test.go` | inductance against Nagaoka and Wheeler, the opposite-current and parallel windings, build-up and pickup tests |
| `export_test.go` | AC netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt`, `pickup_output.txt` | the runs above |
| `twin_ngspice_test.go`, `viz/ngspice_twin.py`, `ngspice_twin_output.txt` | the ends-open / ends-joined netlists for ngspice and the kronspice–ngspice comparison |
| `viz/render.py`, `viz/*.csv` | the pictures and animations and their data |
| `viz/make_report.py`, `pulsed_bifilar_report.pdf` | report of the first four coils and the pickup |
| `viz/make_current_report.py`, `pulsed_bifilar_current_report.pdf` | report of the four connections: volts, amps and power, V and I against frequency |

## Sources

- N. Tesla, "Coil for Electro-Magnets," US 512,340 (1894): the bifilar winding, V/2 between neighbouring turns; "Electrical Transformer," US 593,138 (1897): the flat spiral and the quarter-wave secondary (CLAUDE.md §3e).
- J. C. Maxwell, *A Treatise on Electricity and Magnetism* (1873), Vol. II art. 701: mutual inductance of coaxial circles.
- F. W. Grover, *Inductance Calculations* (Van Nostrand, 1946): single-turn self-inductance.
- H. Nagaoka, "The inductance coefficients of solenoids," *J. Coll. Sci. Tokyo* 27, 1 (1909): current-sheet solenoid check.
- H. A. Wheeler, "Simple inductance formulas for radio coils," *Proc. IRE* 16, 1398 (1928): flat-spiral check.
- S. Ramo, J. R. Whinnery, T. Van Duzer, *Fields and Waves in Communication Electronics*, ch. 4: two-wire capacitance between turns.
- F. E. Terman, *Radio Engineers' Handbook* (1943), sec. 2: skin-effect resistance.
- D. W. Knight, "The self-resonance and self-capacitance of solenoid coils," v1.00 (2016), §13: a free coil's first resonance at half a wavelength of wire.
