# 6. Coils, pulses and gaps

Nikola Tesla's patents set three of this chapter's questions. A coil
wound with two conductors side by side, the end of one joined to the
start of the other, stores far more energy between its turns than an
ordinary winding (US 512,340, 1894). A coil driven by sharp one-sided
pulses from a disruptive discharge builds up a large standing wave when
the pulses arrive in step with it (US 462,418, US 593,138). And a spark
gap, the element that makes those pulses, is a nonlinear device in its
own right, whose behaviour depends on what fills the gap. kronspice
treats all three with the ledger on, so that every volt of build-up is
traced to the source that paid for it (CLAUDE.md §3e, §3g, §3j). The
examples are [bial_bifilar](../examples/bial_bifilar), the project's
first experiment, [pulsed_bifilar](../examples/pulsed_bifilar),
[water_lamp](../examples/water_lamp), and the two spark-gap diode
studies [sparkgap_diode](../examples/sparkgap_diode) and
[sparkgap_diodebridge](../examples/sparkgap_diodebridge).

## 6.1 Tesla's comparison, and what a coil really does

Tesla's example is a 1,000-turn coil with 100 V across it. Wound
normally, adjacent turns differ by about 0.1 V; wound bifilar, the turn
next to any given turn belongs to the other conductor, half a coil
away, so adjacent turns differ by about 50 V. Since the energy in a
capacitance goes with the square of the voltage, he states that the
energy stored in the coil becomes (50/0.1)² = 250,000 times as great
(US 512,340, p. 2, ll. 3–7). His purpose was to build into the coil the
capacitance needed to cancel its self-inductance, instead of adding
condensers.

The engine reproduces the idealized comparison exactly. `em/coil.go`
computes the turn-to-turn capacitance with V/N per gap for an ordinary
winding and V/2 for a bifilar one, and `bial_bifilar`'s
`TestTeslaComparison` checks that the bifilar coil stores (N/2)² =
250,000 times the ordinary coil's turn-to-turn energy at 1,000 turns.
CLAUDE.md §3e then draws a line that this chapter keeps: the idealized
ratio and a real coil's energy are separate things, and the first is
never presented as the second.

The caveat comes from Knight (2016, §6 and §13). Turn-to-turn
capacitance models have no predictive power for a single-layer
solenoid's self-capacitance or self-resonance; Medhurst's 1947
measurements already showed them underestimating it, and a coil's
end-to-end capacitance does not change with the number of turns. The
V/N model is therefore kept only for the comparison above. Predictions
of an ordinary coil's resonance come from `em/helix.go`, Knight's
helical transmission-line model, and the coupling between bifilar
strands, a real capacitance between two conductors at different
potentials, from the coupled-line model of §5.2. The sectioned coils of
§6.3 are built that way, and the first-mode result there is what the
later examples use.

## 6.2 A Bi–Al element charging a bifilar coil

The project's first experiment put a bismuth–aluminium capacitor element
in parallel with a bifilar coil, self-driven, with the element as the
only energy source, and let evolution trade the tank's Q at resonance
against the energy stored per volume (milestone M5a). The element comes
in three variants, each with a declared port:

| Variant | Element | Energy enters through |
|---|---|---|
| dry | Bi and Al plates with a solid dielectric; the Bi–Al leads are a thermocouple, EMF 75.5 µV/K × ΔT, in series with the plate capacitance | thermal port |
| thermopile | the dry element charged by n Bi–Al pairs in series, with the legs' resistance, volume and conductive heat leak | thermal port |
| wet | Al\|electrolyte\|Bi galvanic cell (at most 1.970 V), electrolyte resistance, Randles electrodes, capacity by Faraday's law | chemical port |

The genes are the coil (turns, diameter, wire, pitch, insulation) and
the plates (side, gap, dielectric), with the thermopile adding pairs and
leg dimensions; a wire longer than λ/10 at resonance is rejected by the
check of §5.2. The run kept in
[doc/experiments/bial_bifilar_run2.md](../doc/experiments/bial_bifilar_run2.md),
5,710 evaluations per variant, gives:

| Variant | Highest Q | Highest energy per volume |
|---|---|---|
| dry | 39 turns, 150 mm, 311.7 kHz, Q 1,050, 5.4 pJ/m³ | 15 turns, 5 mm, 10 µm gap, 78.9 kHz, Q 0.75, 1.254 µJ/m³ |
| wet | 23 turns, 141 mm, 11.16 kHz, Q 6.4, 4.7 mJ/m³ | 10 turns, 5.0 mm, 15.63 kHz, Q 2.4, 89.3 J/m³ |
| thermopile | 3 pairs, 11 turns, 21.15 MHz, Q 1,250, 86 pJ/m³ | 5,000 pairs, 19 turns, 27.75 MHz, Q 43.3, 31.25 J/m³ |

What the fronts say is sobering and checkable. In the dry element the
only store is the field in the plate gap, ½ε(EMF/gap)² diluted by the
plates: 1.255 µJ/m³ by that formula against 1.254 µJ/m³ simulated. With
a millivolt Seebeck EMF the stored energy is necessarily microscopic. The
wet element stores mostly ½LI² from its steady current, 4.079 A in the
best design, a current a real aluminium electrode would not sustain, so
those are upper bounds until a cell is measured. The thermopile's 5,000
pairs raise the EMF to 37.75 V and the energy by EMF², but holding
100 K across the legs conducts 47.83 W of heat to store 1.4 mJ. At DC no
voltage stands across the turn-to-turn capacitance, so the bifilar
winding stores the same energy as an ordinary one of the same geometry
and raises Q on none of the 59 dry designs. The ledger over the six
switch-on transients has a worst per-point residual of 1.6×10⁻¹⁴ of peak
power, and model-supplied energy is exactly 0 J: every joule came from
heat or chemistry through the declared port. The Q and resonance values
rest on the superseded coil model of §6.1 and are marked in the README as
needing a re-run.

The element is also the bench for register hypothesis H4, Jack's temporal
field, which predicts an EMF from the time variation of temperature under
rapid heating with an imposed anisotropy. The standard prediction, the
Seebeck EMF from the instantaneous temperature difference alone, is what
this example computes; the difference is the measurement target.

## 6.3 Pulses on coils: build-up only in step with resonance

`pulsed_bifilar` leaves lumped coils behind. Four windings of the same
0.5 mm copper wire, about 126 m each, are each cut into twenty sections:
an ordinary 400-turn solenoid, the same solenoid wound bifilar, a
bifilar cone opening upward at 45°, and a bifilar flat spiral. Every
turn is a circle, so the inductance matrix of the sections is summed
turn by turn from Maxwell's formula for coaxial circles (*Treatise*
Art. 701, `em.MutualCircles`) and the sections are coupled in the
circuit with K elements; capacitance between physically adjacent turns
uses the two-wire formula (12.1 nF in all), capacitance to the grounded
can or plane is added per section, and the wire resistance is Terman's
skin-effect value at the first resonance. The drive reaches the coil
only by magnetic coupling, as in a Tesla coil: 100 V one-sign pulses
from 50 Ω into a 5-turn primary below the base.

An AC sweep first finds each coil's first resonance:

| Coil | f1 | Wire / λ | Q | Inductance |
|---|---|---|---|---|
| Ordinary solenoid | 1.167 MHz | 0.489 | 210 | 4.853 mH (Nagaoka 4.873 mH) |
| Bifilar solenoid | 41.61 kHz | 0.017 | 83 | 4.853 mH |
| Bifilar cone | 35.81 kHz | 0.015 | 96 | 6.542 mH |
| Bifilar flat spiral | 34.93 kHz | 0.015 | 94 | 6.842 mH |

The two solenoids have the same inductance to 10⁻¹², because the strands
of a series bifilar coil carry the same current the same way, and the
flat spiral is within 1.7 % of Wheeler's formula. The ordinary solenoid's
wire is 0.489 wavelengths at f1: it rings as a free half-wave line, the
first resonance Knight reports for a free single-layer coil, with voltage
maxima at its ends. The bifilar windings, with 12 nF between their
strands, resonate 28 to 33 times lower.

The pulse runs then sweep the repetition rate. A one-sign pulse train at
rate r contains only the frequencies r, 2r, 3r, and so on, so it can
build the coil up only when one of those lands on f1. The coil-end
voltage after each run, from `run_output.txt` (the full run takes about
ten minutes, so the numbers here are the recorded ones):

| Coil | r = f1/3 | f1/2 | 0.97 f1 | f1 | 1.03 f1 | 2 f1 |
|---|---|---|---|---|---|---|
| Ordinary solenoid | 2,632 V | 3,868 V | 893 V | **7,426 V** | 657 V | 356 V |
| Bifilar solenoid | 97 V | 143 V | 55 V | **282 V** | 59 V | 5 V |
| Bifilar cone | 53 V | 79 V | 27 V | **155 V** | 28 V | 3 V |
| Bifilar flat spiral | 87 V | 128 V | 44 V | **252 V** | 46 V | 5 V |

Build-up happens at f1, f1/2 and f1/3 and nowhere else, which is the
pulse-drive condition of CLAUDE.md §3e step 3 seen directly. Between
neighbouring turns at f1 the bifilar windings carry half the coil
voltage, 141 V in the solenoid against 38 V for the ordinary winding, as
the patent says. And the ledger shows where the energy went:

| Coil, at f1 | In | Heat | Stored | Time-stepping error |
|---|---|---|---|---|
| Ordinary solenoid | 2.352 mJ | 2.204 mJ | 0.132 mJ | 1.6×10⁻⁵ J |
| Bifilar solenoid | 79.82 mJ | 79.71 mJ | 0.115 mJ | 1.2×10⁻⁶ J |
| Bifilar cone | 108.3 mJ | 108.3 mJ | 0.035 mJ | 1.4×10⁻⁷ J |
| Bifilar flat spiral | 108.5 mJ | 108.4 mJ | 0.093 mJ | 9.9×10⁻⁸ J |

The bifilar coils take far more from the source, because their wide
pulses at low f1 spend most of it as heat in the primary circuit, and
store less: their low-frequency, low-Q mode builds up less than the
ordinary coil's half-wave resonance at 1.17 MHz. Every joule in any coil
came from the pulse source.

The bifilar solenoid's first mode is what the later examples drive, and
`water_lamp` exports it as a netlist a reader can run: the mode as its
inductance, the capacitance that resonates it and a resistance for its Q,
coupled by K to the 5-turn primary, with one 1 A pulse of 2 µs:

```
bifilar coil, first mode, pulsed primary
IP 0 p PULSE(0 1 0 2e-07 2e-07 2e-06 1)
LP p 0 5.974707646431556e-06
LS s 0 0.004853
CS s 0 3.014633635100319e-09
RS s 0 105309.11616312526
K1 LP LS 0.16502237933600086
.tran 0.1u 240u
.end
```

From L and C, 1/(2π√LC) = 41.61 kHz, the f1 above. The report says the
coil end is ringing at 2.64 V after 240 µs, and its ledger is worth
reading closely:

```
energy in            3.02933 µJ
energy out           2.79479 µJ (returned to sources and ports)
energy dissipated    13.7236 nJ
stored: start → end  0 J → 12.3915 nJ
integration mismatch 208.421 nJ (time-stepping error, not a physical channel)
```

Most of the pulse's energy goes in and comes straight back to the source
through the primary; only 26 nJ stays as heat and stored energy. Against
that small net, a mismatch of 208 nJ is large, and the label on it is
doing its job. Rerun with a 10 ns maximum step (`.tran 0.1u 240u 0 10n`)
and the mismatch falls to 29 nJ while heat and stored energy are
unchanged to four digits, 13.72 nJ and 12.39 nJ: the error was in the
time stepping across the pulse edges, not in the physics. The example's
own test holds its one-mode model against the full sectioned coil to
8×10⁻⁴ of the peak.

A pickup closes the study: four turns of 2 mm copper round the ordinary
coil's middle, where its current is largest, tuned by a series capacitor
C = 1/(ω1²L) = 4.56 nF so only the load and the wire remain in the loop.
With a 0.2 Ω load it delivers 3.83 A at 0.77 V, 1.40 W, which is 66.5 % of
the 2.11 W that gets past the source's own 50 Ω; a 10 Ω load with the
capacitor re-optimized takes 91.3 %. The pickup converts voltage into
current. Its load never gets more than the power that reached the
primary, and the ledgers close to 10⁻⁵ J.

## 6.4 The flashing lamp, with water instead of neon

The neon relaxation oscillator is the simplest pulse source: a supply
charges a capacitor through a resistor until the lamp strikes and empties
it. Neon lamps strike reliably because a trace of krypton-85 supplies the
first free electron; without it, in the dark, a lamp waits for a cosmic
ray. `water_lamp` puts water in place of the neon and the starter three
ways, beside the neon lamp itself, in the same circuit. The engine is a
seeded RK4 integrator with every energy carried as a state, so each
ledger closes to the integrator's accuracy; the starting electron is a
random process with a declared rate, so the mean wait once past
breakdown is 1/(S(1 − Vb/Vs)), and the test checks that identity.

| Lamp | Supply | Strikes at | Flashes every | Energy per flash |
|---|---|---|---|---|
| Neon NE-2 with radioactive starter | 120 V | 78 ± 1 V | 0.039 s (±6 %) | 13.5 µJ |
| Neon, no starter, in the dark | 120 V | 120 V | 540 s ± 100 % (expected 576 s) | 55 µJ |
| Water-vapour tube, in the dark | 600 V | 600 V | 662 s ± 108 % (expected 792 s) | 1.13 mJ |
| Water-vapour tube, ultraviolet LED starter | 600 V | 492 ± 4 V | 0.174 s (±5 %) | 0.55 mJ |
| Spark gap in pure liquid water | 3,000 V | 1,076 ± 91 V | 4.4 ms (±11 %) | 5.8 mJ |
| Air spark gap, dry | 10,000 V | 6,798 ± 894 V | 11.8 ms (±27 %) | 0.235 J |
| Air spark gap with water mist | 10,000 V | 8,224 ± 643 V | 18.1 ms (±24 %) | 0.34 J |

The water-vapour tube runs at the measured Paschen minimum of water
vapour, 480 V at 0.6 Torr·cm (Škoro et al. 2011). Without a starter it has
the neon lamp's dark-effect problem, 54 flashes in ten hours at the full
supply voltage; an ultraviolet LED freeing electrons from the cathode does
the radioactive gas's job. The liquid gap's strike voltage is Stygar's
water-breakdown relation extrapolated from sub-microsecond pulses to
millisecond charging, and is labelled so; its leak current splits a
little water, booked on the chemical channel by Faraday's law. The mist
raises the strike voltage by Liu's fog factor, 29.4/21.7, also an
extrapolation from 5 cm gaps to 1 mm. Faster recovery with mist is
register hypothesis H5, carried here as an untested parameter: the run
uses the dry recovery.

Each lamp then drives the bifilar coil's first mode through the 5-turn
primary. Glow discharges are too slow and gentle, and even the sparks
hand over under 1 % of their energy, because the 10 nF flash capacitor
rings with the 6 µH primary at about 650 kHz, sixteen times the coil's
41.6 kHz. Tuning fixes that, Lp·Cp = Ls·Cs as in a spark-gap Tesla coil,
and a quenched gap that opens at the first current zero after the
transfer keeps the energy from swinging back into the arc:

| Mist gap | To the coil per flash | Of the flash | Coil peak |
|---|---|---|---|
| Untuned (10 nF) | 0.14 mJ | 0.6 % | 313 V |
| Tuned, ordinary arc | 1.9 J | 3.3 % | 68.5 kV |
| Tuned, quenched gap | 5.4 J | 9.5 % | 68.4 kV |
| Evolved: 0.249 µF, 20 turns, quench at 0.65 of the transfer | 4.35 J | 74.7 % | 51.6 kV |

The evolved design is flagged in the same breath: neighbouring bifilar
turns then carry about 26 kV, far beyond enamelled wire, so it needs
thicker insulation or a gentler drive, which is the §3e rule applied.
Parallel gaps on one capacitor do not add spark current, because only
the first to break down carries the flash, while parallel capacitor
plates raise the energy per flash by n and need the primary redesigned to
stay tuned. Through all of it the sentence that ends the program's
output holds: every joule comes from the supply, and the water adds no
energy.

## 6.5 Water as a diode

Can a gap act as a diode because of its electrode materials? It was done
over a century ago. A valve metal, aluminium, grows an insulating oxide
in water when it is the positive electrode and conducts freely when it is
the negative one; paired with lead or iron it is the electrolytic
rectifier of Pollak (1896) and Nodon (1899), about 140 V per cell and
75 to 80 % efficient (*Encyclopaedia Britannica*, 1911). `sparkgap_diode`
models the cell as a piecewise-linear element with the oxide film's
capacitance in parallel, and the export folder holds a full bridge of
four cells on a 60 Hz, 100 V mains with a 2,200 µF filter and a 50 Ω
load. Each cell is one B element:

```
BDA0 ac1 p I=pwl(V(ac1,p), -2800, -1266.77, -140, -0.1, -2, -0.0014, 0, 0, 2, 0, 140, 1380, 2800, 27980)
CDA0 ac1 p 4.378444522813186e-06
```

Run the netlist for its half second of transient. The report ends with
the output between p and n at 43.83 V and −43.92 V, which is the 88 V DC
the example quotes, and a ledger that closes:

```
energy in            105.828 J
energy out           659.572 mJ (returned to sources and ports)
energy dissipated    96.6812 J
stored: start → end  0 J → 8.48738 J
integration mismatch 3.49244 µJ (time-stepping error, not a physical channel)
```

The electrode pairs behave as the history says. Lead, iron or silver
against aluminium gives 980 A forward and 0.07 A reverse through 10 cm
plates 1 cm apart, a ratio of 13,700; iron against copper is symmetric;
aluminium against aluminium blocks both ways. Past 140 V the oxide sparks
over, like a Zener diode, which is why a filtered half-wave rectifier,
which sees twice the peak in reverse, needs two cells in series. In the
bridge, 85 % efficiency sits at the top of the 1911 account's range; at
20 kHz the oxide film's capacitance passes current both ways and the
efficiency collapses to 17 %, the article's warning that these
rectifiers are unsuitable for high frequencies. Silver is a poor positive
electrode: by Faraday's law a 10 µm plating lasts about 18 minutes at
the bridge's current.

At kilovolts the water itself breaks down, and what makes a spark
one-way is the shape of the electrodes, not the metal: a sharp positive
point breaks down at a lower voltage than a negative one. A spark gap is
a switch with memory, which a static curve cannot hold, so the gap
circuits run in a small nodal solver of the example's own, checked
against kronspice below breakdown. A half-wave gap rectifier with a
resistive load works; with a filter capacitor the reverse gap sees twice
the peak and breaks down backwards; a bridge needs a polarity ratio above
2, and the ratio used, 3, is declared and awaits a bench measurement.
The pulse test says what each device is for: the same 10 J through a
water polarity gap gives 4,250 A rising in 0.11 µs, through a cell
470 A in 9.2 µs, because gaps work at kilovolts and cells are limited by
the water's resistance.

`sparkgap_diodebridge` combines cells and gaps for six jobs, each with
its ledger: a mains DC supply (88 V from 100 V peak, 85 %, and 87.96 V
from the same circuit in kronspice, 0.01 % apart), surge protection
where an air-gap arrester and an aluminium–aluminium clamp hold a 4 kV
surge to 117 V on the bridge cells against 221 V unprotected, a spark
pulser that gets 50 % of the mains power into 120 pulses per second only
when its output gap is fired at the mains zero crossings, a Greinacher
doubler that needs two cells per position, the oxide film as a capacitor
whose value is set by its forming voltage (a 1 mH tank tunes from 643 Hz
formed at 10 V to 2,405 Hz at 140 V, kronspice AC against 1/(2π√LC) to
the hertz), and a frequency doubler. The water and the metals add no
energy: every watt at a load comes from the mains or the surge source.

## 6.6 What to keep

- Tesla's 250,000× is an idealized ratio of turn-to-turn energies,
  reproduced exactly and never reported as a coil's energy gain.
- Turn-to-turn capacitance models do not predict a coil's resonance; the
  helical-line and coupled-line models do, and the sectioned coil is the
  honest form.
- One-sign pulses build a coil up only when a harmonic of the pulse rate
  lands on a mode, and the build-up is paid for by the source, mostly as
  heat in the primary circuit.
- A transient ledger's integration mismatch is a number to check against
  the net energy, and shrinks with the step; heat and storage do not
  move.
- Gaps and cells are nonlinear elements with declared, cited parameters,
  and extrapolations beyond the measured range are labelled as such.
- Every design is checked against insulation at its peak turn voltage,
  and a result past that limit is flagged, not celebrated.

## Run it yourself

```
go test ./examples/bial_bifilar -run TestTeslaComparison
go run ./examples/bial_bifilar -quick
go run ./cmd/kronspice examples/water_lamp/export/bifilar_first_mode.cir
go run ./examples/pulsed_bifilar -only-pickup         # the full run takes about ten minutes
go run ./examples/water_lamp                          # about three minutes
go run ./cmd/kronspice examples/sparkgap_diode/export/water_cell_bridge.cir
go run ./examples/sparkgap_diode
go run ./examples/sparkgap_diodebridge
```

## Sources

- N. Tesla, "Coil for Electro-Magnets," US 512,340, 1894; "Method of and
  Apparatus for Electrical Conversion and Distribution," US 462,418,
  1891; "Electrical Transformer," US 593,138, 1897 (`reference/`).
- D. W. Knight, "The self-resonance and self-capacitance of solenoid
  coils," v1.00, 2016, DOI 10.13140/RG.2.1.1472.0887; R. G. Medhurst,
  *Wireless Engineer*, 1947.
- H. A. Wheeler, *Proc. IRE* 16, 1398, 1928; H. Nagaoka, *J. Coll. Sci.
  Tokyo* 27, 1, 1909; F. W. Grover, *Inductance Calculations*, 1946;
  F. E. Terman, *Radio Engineers' Handbook*, 1943; J. C. Maxwell,
  *Treatise*, Vol. II, Art. 701.
- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed.,
  2008; C. P. Steinmetz, *Elementary Lectures on Electric Discharges,
  Waves and Impulses*, 1911 (§3e, §3g).
- Seebeck coefficients: *Electronics Cooling*, Nov. 2006, as recorded in
  `material/`; CRC *Handbook of Chemistry and Physics*, electrochemical
  series; A. J. Bard, L. R. Faulkner, *Electrochemical Methods*, ch. 13.
- P. M. Jack, arXiv:math-ph/0307038, 2003 (hypothesis source for H4).
- N. Škoro, D. Marić, G. Malović, W. G. Graham, Z. Lj. Petrović, *Phys.
  Rev. E* 84, 055401(R), 2011; W. A. Stygar et al., *Phys. Rev. ST Accel.
  Beams* 9, 070401, 2006; Y. Liu et al., *Energies* 11, 2251, 2018;
  Yu. P. Raizer, *Gas Discharge Physics*, 1991; C. S. R. Bhoomireddy, PhD
  thesis, HBNI/BARC, 2016.
- *Encyclopaedia Britannica*, 11th ed., 1911, vol. 27, "Transformers";
  H. Greinacher, *Z. Phys.* 4, 195, 1921; IEC 61000-4-5.
- [doc/experiments/bial_bifilar_run2.md](../doc/experiments/bial_bifilar_run2.md).
