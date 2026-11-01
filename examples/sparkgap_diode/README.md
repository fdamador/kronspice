# A water gap that acts as a diode: spark-gap diode and diode bridge

Can a water gap act like a diode because of its two electrode materials?
Yes — and it was done over a century ago. A **valve metal** (aluminium) in
water grows an insulating oxide when it is the positive electrode and lets
current through freely when it is the negative one. Paired with an ordinary
metal (lead, iron) it is the aluminium electrolytic rectifier of Pollak
(1896) and Nodon (1899): about 140 V per cell, about 75–80 % efficient
(*Encyclopaedia Britannica* 1911, "Transformers"). Past the oxide's
sparking voltage the blocked direction sparks over, like a Zener diode.

At thousands of volts the water itself breaks down in a spark. There, what
makes a spark one-way is the **shape** of the electrodes, not the metal: a
sharp positive point breaks down at a lower voltage than a sharp negative
point. That is the second device, the polarity spark gap.

```
go run ./examples/sparkgap_diode -viz examples/sparkgap_diode/viz   # a few seconds
python3 examples/sparkgap_diode/viz/render.py
go test ./examples/sparkgap_diode
```

## 1. Electrode pairs in water

10 cm × 10 cm plates, 1 cm apart, in a phosphate solution (water 0.1 Ω).

| Pair | +100 V | −100 V | Forward/reverse |
|---|---|---|---|
| lead \| water \| aluminium | 980 A | −0.07 A | 13,700 — a diode |
| iron \| water \| aluminium | 980 A | −0.07 A | 13,700 — a diode |
| copper with silver tips \| water \| aluminium | 98 A | −0.07 A | 1,370 — a diode (tips: 10 % of the plate) |
| silver-plated iron \| water \| aluminium | 980 A | −0.07 A | 13,700 — a diode |
| iron \| water \| copper | 980 A | −980 A | — symmetric (a weak 0.79 V battery) |
| lead \| water \| iron | 980 A | −980 A | — symmetric |
| aluminium \| water \| aluminium | 0.07 A | −0.07 A | — blocks both ways |

Any ordinary (non-valve) metal works as the partner of aluminium,
including silver: the one-way action comes from the aluminium's oxide.
Past 140 V the aluminium's oxide sparks over (−28.7 A at −200 V for
lead|aluminium). Each curve is checked against a kronspice DC sweep of the
same element. The oxide film is also a 4.4 µF capacitor.

## 2. Rectifiers of aluminium–lead water cells (kronspice .TRAN)

Load 50 Ω with a 2,200 µF filter capacitor; source resistance 0.5 Ω.

| Circuit | Output | Ripple | Efficiency | Reverse peak |
|---|---|---|---|---|
| Half-wave, 1 cell, 100 V peak | 76 V | 27 V | 20 % | 167 V — sparks over |
| Half-wave, 2 cells in series | 83 V | 11 V | 84 % | 183 V |
| **Bridge, 1 cell per arm, 100 V peak** | **88 V** | **5 V** | **85 %** | 93 V |
| Bridge, 1 cell per arm, 170 V peak (120 V rms) | 145 V | 10 V | 40 % | 153 V — sparks over |
| Bridge, 2 cells per arm, 170 V peak | 147 V | 8 V | 86 % | 156 V |
| Bridge, 100 V peak, 1 kHz | 88 V | 0.3 V | 84 % | |
| Bridge, 100 V peak, 20 kHz | 86 V | 0.02 V | 17 % | |
| Bridge, copper with silver tips \| aluminium, 100 V peak | 79 V | 4 V | 78 % | 88 V |
| Bridge, silver-plated iron \| aluminium, 100 V peak | 88 V | 5 V | 85 % | 93 V |

**Silver against aluminium.** Both silver versions rectify. The catch is
that when a cell conducts, the silver is the positive electrode, and a
positive metal electrode in water is eaten away. If all the current
dissolved silver (Faraday's law, the worst case — part of it evolves
oxygen or forms a film instead):

| Silver electrode (declared size) | Silver | Lasts at the bridge's current |
|---|---|---|
| Copper with 1 mm solid silver tips on 10 % of the plate | 10.5 g | about 3.3 hours (0.79 A per cell) |
| Iron with 10 µm silver plating | 1.05 g | about 18 minutes (0.88 A per cell) |

After that the copper or iron underneath carries the current and corrodes.
This is why Pollak used lead: as the positive electrode it forms a stable
lead-dioxide skin instead of dissolving. A thin silver plating is
therefore a poor choice for the positive side; solid silver lasts longer
but is still consumed. In the high-voltage spark gap (Part 3), silver tips
make a good point electrode — silver is a usual arc-contact metal — but
the one-way action there comes from the point's shape, not the silver.

- 85 % sits at the top of the 75–80 % the 1911 account gives; the cell's
  leakage and resistance here are declared values.
- A filtered half-wave rectifier sees twice the peak in reverse, so one
  cell sparks over; two in series hold — as the historical practice of
  "two or more valves in series for higher pressures".
- At 20 kHz the oxide film's capacitance passes current both ways and the
  efficiency collapses — the 1911 article's warning that these rectifiers
  are unsuitable for high frequencies.

Every run's energy ledger closes to about 10⁻⁷ or better.

## 3. Polarity spark gaps as rectifiers

A 3 mm point–plane gap in pure water: breaks down at +2.4 kV with the
point positive, −4.8 kV with it negative. **Both numbers are uncertain:**
the positive one is the Sandia water relation (Stygar et al. 2006)
extrapolated from sub-microsecond pulses to 60 Hz, and the ratio of 2 is
declared (J. C. Martin's constants, to be checked).

| Circuit | Polarity ratio | Result |
|---|---|---|
| Half-wave, resistive load, 4 kV peak | 2 | rectifies (pulsed output, 96 %) |
| Half-wave, filter capacitor, 4 kV peak | 2 | breaks down backwards (reverse sees 2× peak) |
| Bridge, filter capacitor, 4 kV peak | 2 | never strikes (two gaps in series need 2 × 2.4 kV) |
| Bridge, filter capacitor, 6 kV peak | 2 | breaks down backwards |
| Bridge, filter capacitor, 6 kV peak | 3 | rectifies, 68 % |

A spark-gap bridge needs a polarity ratio above 2: two gaps in series must
both break down forward while each blocked gap holds the full output. And
a spark gap is not a low-drop diode: it conducts only once the voltage
passes its breakdown, so the output sags and jumps. A bench measurement of
the breakdown voltages and their ratio decides whether a water spark-gap
bridge can work.

The spark gaps are stepped with a small nodal solver of their own (a spark
gap is a switch with memory, which a static kronspice curve cannot hold);
below breakdown it matches kronspice (TestGapSolverTwin).

Picture: `viz/sparkgap_diode.png`.

## 4. Pulse test: which gives the fastest amps?

The same 10 J from a charged capacitor through each device, a 1 µH loop
and a 0.1 Ω load (declared).

| Device | Charge | Peak current | Rise (10–90 %) | To the load | Lost in the device | Left in the capacitor |
|---|---|---|---|---|---|---|
| lead \| water \| aluminium | 100 V, 2 mF | 470 A | 9.2 µs | 4.8 J | 5.2 J | 0 |
| silver-plated iron \| water \| aluminium | 100 V, 2 mF | 470 A | 9.2 µs | 4.8 J | 5.2 J | 0 |
| copper with silver tips \| water \| aluminium | 100 V, 2 mF | 195 A | 2.2 µs | 0.9 J | 8.8 J | 0.3 J |
| **water polarity gap, 3 mm** | 41 kV, 12 nF | **4,250 A** | **0.11 µs** | 0.3 J | 1.6 J | 8.1 J |
| mist-filled air gap, 1 mm | 8 kV, 0.31 µF | 3,510 A | 0.51 µs | 1.1 J | 5.6 J | 3.4 J |

Spark gaps give the highest and fastest current because they work at
kilovolts. The cells give hundreds of amps in microseconds, limited by the
water's resistance, and they hand over all the energy. A gap's arc goes
out at the first current zero, so energy still ringing in the capacitor
stays there. The water gap's 41 kV is Stygar's relation at 1 µs, inside
its measured range; the mist gap's 6.8 kV uses the fog factor from
[water_lamp](../water_lamp) (extrapolated). Picture: `viz/pulses.png`.

## 5. Frequency sweep: amps per volt vs frequency

Every device in the same bridge (source 0.1 Ω, load 1 Ω, no filter), from
50 Hz to 100 kHz. Amps per volt = mean load current ÷ source peak voltage.
A perfect diode bridge gives 0.579 A/V at every frequency.

| Device (peak) | 50 Hz | 100 Hz | 300 Hz | 1 kHz | 3 kHz | 10 kHz | 30 kHz | 100 kHz |
|---|---|---|---|---|---|---|---|---|
| lead \| aluminium (100 V) | 0.459 | 0.459 | 0.459 | 0.459 | 0.460 | 0.463 | 0.481 | 0.508 |
| … efficiency | 73 % | 73 % | 73 % | 73 % | 73 % | 71 % | 60 % | 34 % |
| silver-plated iron \| aluminium (100 V) | same as lead | | | | | | | |
| copper with silver tips \| aluminium (100 V) | 0.192 | 0.192 | 0.192 | 0.192 | 0.193 | 0.193 | 0.195 | 0.191 |
| … efficiency | 31 % | 31 % | 31 % | 30 % | 30 % | 29 % | 22 % | 8 % |
| water gap, ratio 3, fixed 6 kV | 0.247 | 0.198 | 0 | 0 | 0 | 0 | 0 | 0 |
| … efficiency | 47 % | 47 % | | | | | | |
| water gap, ratio 3, voltage raised with frequency | 0.239 | 0.240 | 0 | 0 | 0 | 0 | 0 | 0 |
| … efficiency | 47 % | 47 % | | | | | | |

- **Water cells** hold steady to about 3 kHz. Above that the oxide film (a
  capacitor) passes current both ways: the mean current holds up but the
  efficiency collapses.
- **Water spark gaps** work only at mains frequencies. At a fixed 6 kV the
  breakdown voltage rises with frequency (E·τ^⅓, extrapolated) until, above
  100 Hz, the gaps never strike. With the voltage raised to keep them
  striking, above about 100 Hz the gaps that just carried an arc have not
  recovered and break down backwards: the bridge shorts. Recovery is
  declared (1 − e^(−t/1 ms)), and a hot gap never re-strikes below the
  Paschen minimum of water vapour, 480 V (Škoro et al. 2011) — its channel
  is vapour.

Picture: `viz/frequency_trend.png`.

Bridges that combine these cells and gaps for particular jobs — surge
protection, a spark pulser, a voltage doubler, the oxide film as a
varactor, a frequency doubler — are in [sparkgap_diodebridge](../sparkgap_diodebridge).

The water and the metals add no energy: every watt at the load comes from
the AC source; the rest is heat and a little electrolysis.

## Files

| File | What it holds |
|---|---|
| `cell.go` | electrode materials, the water cell's I–V, its kronspice element |
| `gap.go` | the polarity spark gap |
| `gapsim.go` | the spark-gap rectifier solver |
| `rectifier.go` | half-wave and bridge rectifiers in kronspice; I–V twin |
| `pulse.go` | the 10 J pulse test |
| `sweep.go` | the frequency sweep |
| `main.go` | the report |
| `main_test.go` | pair behaviour, rectifier efficiencies, gap verdicts, pulses, sweep, twins, ledgers |
| `viz/render.py` | `sparkgap_diode.png`, `pulses.png`, `frequency_trend.png` |
| `export_test.go` | the water-cell bridge's netlist, KiCad schematic and results in `export/` |

## Sources

- *Encyclopaedia Britannica*, 11th ed. (1911), vol. 27, "Transformers" — Pollak and Nodon aluminium rectifiers: about 140 V per cell (50–140 V Nodon), efficiency about 75 % (Pollak) and 75–80 % (Nodon), unsuitable for high frequency.
- Anodic-film breakdown rising with the logarithm of electrolyte resistivity (Ikonopisov), as cited in *Periodica Polytechnica Chem. Eng.* (BME), breakdown of barrier films around 500–600 V in borate.
- W. A. Stygar et al., *Phys. Rev. ST Accel. Beams* 9, 070401 (2006) — water breakdown for a field-enhanced anode.
- M. Škoro et al., *Phys. Rev. E* 84, 055401 (2011) — Paschen curve of water vapour (minimum 480 V at 0.6 Torr·cm).
- L. Liu et al. (2018) — fog raising rod–plane breakdown 21.7 → 29.4 kV (as used in [water_lamp](../water_lamp)).
- CRC Handbook of Chemistry and Physics, electrochemical series — standard electrode potentials.
