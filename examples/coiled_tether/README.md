# A five-mile tether, coiled and bifilar, in orbit

A 5-mile (8,046.72 m) copper tether in low Earth orbit is coiled up: does
"multiplication between windings" make it an energy source? Three windings
of the same 2 mm wire (43 Ω) are compared: straight (the ordinary
electrodynamic tether of [tether](../tether)), an ordinary single-layer
coil (2 m across, 1,281 turns) and the same coil wound as Tesla's bifilar
coil (US 512,340, CLAUDE.md §3e), so that each turn lies next to the turn
half a coil away and carries half the coil voltage. Energy comes only from
the sources named: the Sun, and Earth's magnetic and electric fields. The
multipactor check between the turns is the one built for
[double_hull](../double_hull) (§3i), and the burst-speed limit is the
rotor rule of §3l.

```
go run ./examples/coiled_tether
go test ./examples/coiled_tether
KRONSPICE_EXPORT=1 go test ./examples/coiled_tether -run TestExport   # netlists, schematics, results sheets
```

## Declared inputs

Declared, not looked up, and stated in `main.go`:

| Input | Value | Note |
|---|---|---|
| Wire | 5 international miles of 2 mm copper | 8,049 m wound (1,281 whole turns), 43.0 Ω DC |
| Orbit speed | 7,700 m/s | relative to the co-rotating field |
| Earth's field B | 30 µT | order of magnitude in LEO, as in [tether](../tether); a real orbit takes B from IGRF |
| Solar array | 2 m² at 30 %, TSI 1,360.8 W/m² | Kopp & Lean 2011; 816 W to the drive |
| Coil | 2 m across, pitch 3 mm | 1 mm bare gap between turns, vacuum (ε = 1) |
| Gap limit | 10 MV/m | to be replaced by the chosen insulation's data-sheet value |

## The construction

**Windings** (`em.Coil`, the lumped first approximation of §3e): Wheeler's
inductance, 1.362 H for either winding, since the two strands of a series
bifilar coil carry the same current the same way; the two-wire capacitance
between adjacent turns; Terman's skin-effect resistance at the run
frequency; and the patent's turn voltage, V/N for the ordinary coil and
V/2 for the bifilar one. The energy between the turns, referred to the
terminals, is the self-capacitance (N − 1)·C_tt/N² for the ordinary coil
(0.142 pF) and (N − 1)·C_tt/4 for the bifilar coil (58.1 nF). Their ratio
is Tesla's N²/4, 410,240 for 1,281 turns, which is the patent's 250,000
for 1,000 turns (US 512,340, p. 2, ll. 3–7).

**Earth's field.** The motional EMF ∫(v×B)·dl is integrated along the
actual wire path with `em.MotionalEMF` (`em/lorentz.go`, the Lorentz force
by grade selection on a unit test charge, §3f): a straight line for the
tether, a helix of 36 points per turn for the coils, and for the bifilar
coil two interleaved strands with the return lead that joins the end of
the first to the start of the second. A fourth path folds the straight
tether in half, out along one strand and back along the other 1 cm away.
In the craft's frame the same coupling is the electric field E′ = v×B
(Lorentz transformation, `em/frames.go`), not a second source.

**The Sun.** The resonant loop is the drive in series in the winding
loop, as a primary would couple it: a sine source on the light channel,
the wire resistance, the coil inductance and the capacitance across the
coil (the bifilar coil's own inter-winding capacitance, or a tank
capacitor of the same 58.1 nF for the ordinary coil), with an optional
load resistor in series below it. The drive amplitude is set to
√(2(R_wire + R_load)·P_array) so that at resonance the loop draws the
array's whole 816 W. Each run is a kronspice transient at 200 points per
period for seven ring-up time constants (τ = 2L/R), with the energy
ledger on; peaks and mean powers are taken over the last 20 cycles.

**Checks on the turn gap.** The voltage, field and stored energy between
neighbouring turns; and `em.Multipactor` over every amplitude the gap sees
during the ring-up, with the measured as-received copper yield (Baglin et
al. 2002), keeping the worst growth per impact.

## 1. Earth's field: the EMF along the real wire path

| Winding | EMF |
|---|---|
| Straight 5-mile tether | 1,859 V |
| Coiled (ordinary) | 0.888 V, which is E′ × the coil's own 3.84 m length: the turns cancel |
| Coiled (bifilar) | 0.889 V: both strands run the same way through the field |
| Straight, folded out and back | 10⁻¹³ V: out and back cancel |

E′ = v×B = 0.231 V/m. In a uniform field the EMF between a wire's ends is
(v×B)·(r_end − r_start) whatever the path, so coiling the tether throws
its voltage away. The dipole field's gradient across the orbit adds at
most N·A·v·3B/r ≈ 0.4 mV, and all of it is quasi-DC: the field seen by
the craft changes over an orbit (~90 min), far below any LC resonance. A
coil spun at 1 rev/s in the field gives 0.76 V peak at 1 Hz, paid for by
the spin, which slows.

## 2. The Sun drives the coil at resonance (high-Q multiplication)

816 W from the array drive the coil at 565.7 Hz (L 1.362 H, C 58.1 nF,
wire 43.0 Ω, Q 112.5; the wire is 0.015 λ, so the lumped model holds).

| | Ordinary coil + tank capacitor | Tesla bifilar | Bifilar, matched load |
|---|---|---|---|
| Drive amplitude | 265.1 V | 265.1 V | 374.9 V |
| Coil voltage | 29,784 V (×112) | 29,784 V (×112) | 21,064 V (×56) |
| Coil current | 6.15 A | 6.15 A | 4.35 A |
| Between neighbouring turns | 23.3 V, 23 kV/m | 14,892 V, 14.9 MV/m | 10,532 V, 10.5 MV/m |
| Energy between turns | 6.3×10⁻⁵ J | 25.78 J | 12.89 J |
| Total stored energy | 25.78 J | 25.78 J | 12.90 J |
| Multipactor growth per impact | 0.063 | 0.070 | 0.070 |
| Steady state: Sun in = heat (+ load) | 815.2 W = 814.2 W | 815.2 W = 814.2 W | 815.2 W = 406.9 W + 406.9 W load |
| Whole run: in = heat + stored | 309.92 J = 284.13 + 25.78 J (error 0.006 J) | same | 155.13 J = 142.23 + 12.89 J (error 0.003 J) |

Multiplication is real for voltage (by Q), and the bifilar winding puts
the whole stored energy between its turns, N²/4 times the ordinary coil's
turn-to-turn energy, with the same total. All of it comes from the Sun's
array; the matched load takes half the input and never more. The bifilar
turn voltage exceeds the declared 10 MV/m gap limit: it needs insulation
or wider spacing. No electron avalanche starts at 566 Hz.

## 3. and 4. Limits

- Driving the resonance from Earth's field instead would need the coil to
  spin at 566 rev/s: rim speed 3,554 m/s against the 1,881 m/s burst
  speed of T1000G carbon fibre, the strongest rotor material in
  `spacetime/`.
- Where electrons would multiply in the same 1 mm vacuum gap (worst
  amplitude up to 2 kV):

| Frequency | Worst growth per impact | At |
|---|---|---|
| 566 Hz | 0.071 | 480 V |
| 1 MHz | 0.067 | 2,000 V |
| 100 MHz | 0.079 | 2,000 V |
| 1 GHz | 1.657, avalanche | 120 V |
| 3 GHz | 2.029, avalanche | 480 V |

At gigahertz the electrons take their energy from the drive, and the 8 km
wire is many wavelengths long, so it is not a coil resonance.

Picture: `ring_up.png` (the ring-up at resonance).

## Tests

| Test | Checks |
|---|---|
| `TestMotionalEMF` | straight = vBL to 10⁻¹²; each coil = vB × its axial length; folded = 0 |
| `TestResonance` | step-up = Q within 1 %; Sun in = array power = heat within 1 %; ledger closes to 10⁻³; turn-energy ratio = N²/4 to 10⁻⁹ with equal totals; matched load takes half the input |
| `TestRuns` | the whole report runs |

The export sheets carry only the MNA column: automatic tearing finds no
valid cut of this single loop (every torn piece is singular, Aitchison
1983 Lemma 2), which the sheet reports instead of solving.

## Files

| File | What it holds |
|---|---|
| `main.go` | declared inputs, windings, motional EMFs, resonant runs, spin and multipactor limits |
| `main_test.go` | EMF identities, resonance and ledger tests |
| `export_test.go` | the three resonant loops as netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt` | the run above |
| `ring_up.png` | the coil voltage ringing up at resonance |

## Sources

- G. Kopp, J. L. Lean, "A new, lower value of total solar irradiance," *Geophys. Res. Lett.* 38, L01706 (2011): 1,360.8 W/m².
- N. Tesla, "Coil for Electro-Magnets," US 512,340 (1894): the bifilar winding and the N²/4 (250,000×) turn-to-turn energy comparison, p. 2, ll. 3–7.
- H. A. Wheeler, "Simple inductance formulas for radio coils," *Proc. IRE* 16, 1398 (1928): coil inductance.
- F. E. Terman, *Radio Engineers' Handbook* (1943), sec. 2: skin-effect resistance.
- S. Ramo, J. R. Whinnery, T. Van Duzer, *Fields and Waves in Communication Electronics*, ch. 4: two-wire line capacitance between adjacent turns.
- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16 (reference only, CLAUDE.md §3f): the Lorentz force by grade selection behind `em.MotionalEMF`.
- V. Baglin et al., "The secondary electron yield of technical materials and its variation with surface treatments," LHC Project Report 472 (2002), Table 1: copper yield; A. J. Hatch, H. B. Williams (1954) and J. R. M. Vaughan, *IEEE Trans. Electron Devices* 35, 1172 (1988): the two-surface multipactor model of `em/multipactor.go`.
- Toray Composite Materials America, T1000G data sheet (16 Feb 2024): fibre strength and density for the burst speed.
- P. W. Aitchison, *Quart. Appl. Math.* 41, 265 (1983), Lemma 2: why the single loop cannot be torn.
