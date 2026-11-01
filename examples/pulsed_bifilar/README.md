# Pulse-driven resonance in ordinary and Tesla bifilar coils

One-sign voltage pulses drive four coils of the same 0.5 mm copper wire
(about 126 m each) through a 5-turn primary, and the example finds the
pulse rates at which the energy builds up into a large standing wave
(CLAUDE.md §3e step 3, "for pulse drive, the equivalent condition is
repeating the pulse at the round-trip time"). A secondary pickup then
takes power off at low voltage and high current. The bifilar solenoid's
first mode is the load that [water_lamp](../water_lamp) drives with its
spark gaps, and the same coil formulas serve [coiled_tether](../coiled_tether)
and [bial_bifilar](../bial_bifilar).

```
go run ./examples/pulsed_bifilar                    # about 10 minutes, writes viz/*.csv
go run ./examples/pulsed_bifilar -only-pickup       # the pickup study alone
python3 examples/pulsed_bifilar/viz/render.py       # pictures and animations
go test ./examples/pulsed_bifilar
KRONSPICE_EXPORT=1 go test ./examples/pulsed_bifilar -run TestExport   # AC netlists, schematics, results sheets
```

## The construction

**Wire and geometry** (declared in `main.go`): 0.5 mm copper, 0.7 mm
between neighbouring wire centres (0.2 mm of enamel, ε = 3 declared), and
four windings cut so that each holds about the same length of wire as the
400-turn solenoid:

| Coil | Turns | Radius | Height | Wire | Surroundings |
|---|---|---|---|---|---|
| Ordinary cylindrical solenoid | 400 | 5.0 cm | 27.9 cm | 125.7 m | grounded can of twice its radius |
| Tesla bifilar solenoid (US 512,340) | 400 | 5.0 cm | 27.9 cm | 125.7 m | same can |
| Tesla bifilar cone, 45° opening upward | 248 | 2.0 to 14.2 cm | 12.2 cm | 126.4 m | grounded plane 5 cm below the base |
| Tesla bifilar flat (Archimedean) spiral | 214 | 2.0 to 16.9 cm | 0 | 127.1 m | grounded plane 5 cm below |

In a bifilar coil the wire runs out along the even positions (strand 1)
and back along the odd ones (strand 2), the end of the first joined to the
start of the second, so neighbouring turns belong to the two strands and
sit half a wire length apart in potential (the patent's V/2).

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
sweep of the terminal voltage per volt of drive, 1 kHz to 100 MHz.

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

The two solenoids have the same inductance to 10⁻¹² (the strands of a
series bifilar coil carry the same current the same way), within 0.4 %
of Nagaoka's current-sheet formula; the flat spiral is within 1.7 % of
Wheeler's flat-spiral formula. The ordinary solenoid's wire is 0.489
wavelengths at f1: it rings as a free half-wave line, the first resonance
Knight reports for a free single-layer coil (§3e), with the voltage
maxima at its ends and the current maximum at its middle. The bifilar
windings, with 12 nF between their strands, resonate 28 to 33 times lower.

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

The bifilar coils take far more energy from the source (their wide pulses
at low f1 spend most of it as heat in the primary circuit) and store less:
their low-frequency, low-Q mode builds up less than the ordinary coil's
half-wave resonance at 1.17 MHz.

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

`WriteViz` writes each coil's geometry, voltage envelope, peak and
ring-up "glow" data to `viz/*.csv` from the full-resolution run at f1;
`viz/render.py` turns them into `coil_voltage.png` (the two solenoids,
every turn coloured by its voltage at the moment of largest voltage, with
the voltage along the wire and between neighbouring turns),
`coil_voltage_shapes.png` (the cone and the spiral), `pulse_sweep.png`
(coil-end voltage against rate / f1), `ring_up.png` (the coil-end voltage
building up pulse by pulse) and the animations `ring_up.gif` and
`ring_up_shapes.gif`.

## Tests

| Test | Checks |
|---|---|
| `TestCoilNetwork` | every winding uses each position once; L within 1 % of Nagaoka and equal for both solenoids; symmetric capacitance matrix; the bifilar coil holds (N − 1) strand-to-strand gaps |
| `TestFlatSpiralInductance` | the spiral against Wheeler's a²N²/(8a + 11c), held to 10 % |
| `TestPulseBuildUp` | pulses at f1 build more than three times higher than at 1.5 f1; ledger closes to 1 % |
| `TestPickup` | the tuning capacitor resonates the pickup at f1 to 10⁻¹²; ledger to 1 %; load power below the power past the source |
| `TestExport` | the four coils swept around resonance, and the tuned pickup, as netlists and results sheets |

## Files

| File | What it holds |
|---|---|
| `main.go` | geometries, the sectioned network, resonance sweep, pulse runs, ledger, picture data |
| `secondary.go` | the tuned pickup, its AC search and pulse runs |
| `main_test.go` | inductance against Nagaoka and Wheeler, build-up and pickup tests |
| `export_test.go` | AC netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt`, `pickup_output.txt` | the runs above |
| `viz/render.py`, `viz/*.csv` | the pictures and animations and their data |

## Sources

- N. Tesla, "Coil for Electro-Magnets," US 512,340 (1894): the bifilar winding, V/2 between neighbouring turns; "Electrical Transformer," US 593,138 (1897): the flat spiral and the quarter-wave secondary (CLAUDE.md §3e).
- J. C. Maxwell, *A Treatise on Electricity and Magnetism* (1873), Vol. II art. 701: mutual inductance of coaxial circles.
- F. W. Grover, *Inductance Calculations* (Van Nostrand, 1946): single-turn self-inductance.
- H. Nagaoka, "The inductance coefficients of solenoids," *J. Coll. Sci. Tokyo* 27, 1 (1909): current-sheet solenoid check.
- H. A. Wheeler, "Simple inductance formulas for radio coils," *Proc. IRE* 16, 1398 (1928): flat-spiral check.
- S. Ramo, J. R. Whinnery, T. Van Duzer, *Fields and Waves in Communication Electronics*, ch. 4: two-wire capacitance between turns.
- F. E. Terman, *Radio Engineers' Handbook* (1943), sec. 2: skin-effect resistance.
- D. W. Knight, "The self-resonance and self-capacitance of solenoid coils," v1.00 (2016), §13: a free coil's first resonance at half a wavelength of wire.
