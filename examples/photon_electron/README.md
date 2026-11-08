# The electron as a photon on a double loop (H23), and electron pumps

Hypothesis H23 (CLAUDE.md §1): J. G. Williamson and M. B. van der Mark,
"Is the electron a photon with toroidal topology?", *Annales de la
Fondation Louis de Broglie* 22, 133–160 (1997). One photon of the Compton
wavelength goes twice round a torus, and the paper derives the electron's
charge, magnetic moment and spin ½ from that. This example puts the paper on
the electron test bench (next to [kron_electron](../kron_electron), H21, and
[space_vortex](../space_vortex), H13), builds its 720° idea as a circuit,
tests its central assumption with the field solver, and then builds the
electron pumps that circuits can actually make: a single-electron pump, a
multipactor pump (§3i, the gap check shared with
[double_hull](../double_hull)) and a resonant ground charge pump (§3j).

```
go run ./examples/photon_electron -viz examples/photon_electron/viz   # about 10 minutes
python3 examples/photon_electron/viz/render.py                         # pictures
go test ./examples/photon_electron                                    # about 1 minute
KRONSPICE_EXPORT=1 go test ./examples/photon_electron -run TestExport  # ring-line netlists and results sheets
```

Every result comes with its standard prediction and an energy ledger; the
source of every joule is named.

## 1. The paper's numbers against measurement

The model: one photon of wavelength λ = h/(m_e c) confined to a closed
path that goes twice round a torus (the paper assumes the confinement; it
adds no field equation). It states q = (1/2π)·√(3ε₀ħc) from the field in a
sphere of diameter λ/2 (eq. 5), α′ = q²/(4πε₀ħc) = (q/e)²α (eq. 8) and
g = 2(1 + α′/2π) (eq. 12, compared with QED's first-order 2(1 + α/2π) in
eq. 13).

| Quantity | Model | Measured (CODATA 2018; Fan et al. 2023) |
|---|---|---|
| Charge, q = (1/2π)·√(3ε₀ħc) (eq. 5) | 0.910318 e | 1 e |
| Coupling α′ = (q/e)²·α (eq. 8) | 0.00604716 | α = 0.00729735 |
| Magnetic anomaly (g − 2)/2 = α′/2π (eq. 12) | 0.000962436 | 0.00115965218 |
| … with q set to e (QED's first-order α/2π) | 0.00116141 | 0.15 % off |
| Spin | ½ (double loop) | ½ |
| Size | λ_C = 2.426×10⁻¹² m; loop radius λ/4π = 1.93×10⁻¹³ m | point-like below 10⁻¹⁸ m |

The anomaly is 17 % low because the charge is 9 % low. Standard
counterpart: the Dirac electron with QED corrections, which gives the
measured values.

## 2. The double loop as a circuit: plain and crossed ring lines

A two-wire ring transmission line (200 LC sections, an LC ladder closed
on itself, a closed path in Kron's sense, as in
[ring_towers](../ring_towers)), once plain and once with its two
conductors swapped at one point. On the plain ring a wave returns to
itself after one trip, so its differential modes have a whole number of
wavelengths round the ring. On the crossed ring the differential voltage
changes sign at the crossing: a wave comes back inverted after one trip
and only matches itself after two, so the modes have half-integer numbers
of wavelengths and the lowest has a wavelength of two circumferences.
Both are solved in kronspice: an AC sweep finds the modes (twin of the
ladder's analytic dispersion) and a transient sends a 2.93 V pulse round.
One trip takes 9.17 µs.

| Mode | Plain ring, AC peak | Ladder theory | Wavelengths | Crossed ring, AC peak | Ladder theory | Wavelengths |
|---|---|---|---|---|---|---|
| 1 | 109,102 Hz | 109,104 Hz | 1 | 54,548 Hz | 54,554 Hz | ½ |
| 2 | 218,181 Hz | 218,182 Hz | 2 | 163,646 Hz | 163,648 Hz | 1½ |
| 3 | 327,205 Hz | 327,206 Hz | 3 | 272,701 Hz | 272,702 Hz | 2½ |
| 4 | 436,148 Hz | 436,149 Hz | 4 | 381,688 Hz | 381,689 Hz | 3½ |

| Trip | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Plain ring, pulse back at the start | +0.90 | +0.69 | +0.64 | +0.57 of launch |
| Crossed ring | −0.91 | +0.70 | −0.59 | +0.52 of launch |

AC peaks agree with the ladder formula within 3×10⁻⁵ (plain) and 1×10⁻⁴
(crossed). On the crossed ring the field returns to itself only after
720°, the paper's spin-½ property. Ledger for each ring: in 3.44×10⁻⁸ J =
heat 1.1×10⁻⁸ J + stored 2.3×10⁻⁸ J, time-stepping error 1.3×10⁻¹¹ J.
Picture: `viz/round_trip.png`.

## 3. Does light stay in the loop? (field solver, `field.go`)

`em/fdtd.go` (Yee's scheme, validated in `em/fdtd_test.go`: TM₁₁₀ cavity
frequency to 10⁻⁷, Hertzian dipole against Larmor within 0.5 %, the exact
Hellwarth–Nouchi "flying doughnut" pulse within 3.8 %, energy ledger to
10⁻¹³, div B at round-off). One wavelength of light is set on the paper's
(2, 1) torus-knot path, twice round the ring of radius R = λ/4π while
once round a tube of radius R/2, as a travelling wave cos(ks − ωt) in a
Gaussian tube of width R/3 with E and B across the path so the energy
flows along it. E and B are built as discrete curls of sampled
potentials, so the start is exactly divergence-free. Grid: λ = 1 m,
12 cells per R, region ±6R (about λ/2, so the absorber sits outside the
near field), followed for 4 periods; the equations are linear, so lengths
scale out.

| Time (periods) | 0 | ⅛ | ¼ | ½ | 1 | 2 |
|---|---|---|---|---|---|---|
| Energy within 2R | 98.9 % | 65.1 % | 5.8 % | 1.6×10⁻⁶ % | 1×10⁻⁴ % | 6×10⁻⁷ % |
| Still in the region | 100 % | 100 % | 100 % | 90.4 % | 0.35 % | 0.012 % |
| Left through the absorber | 0 | 0 | 4×10⁻¹¹ % | 9.6 % | 99.65 % | 99.99 % |

Ledger residual 1.9×10⁻¹³ of the starting energy; div B 7×10⁻¹³. Under
ordinary (linear) Maxwell equations nothing holds light in a loop; it
leaves at the speed of light. The model needs it to stay for ever, so it
needs a confining mechanism the paper does not give. H23 is recorded as
*bounded*. Picture: `viz/knot_light.png`.

## 4. Three electron pumps that work (`pumps.go`)

The model's own frequency, m_e·c²/h = 1.24×10²⁰ Hz, is about 10⁸ times
beyond any circuit, so these pumps run where circuits do. Electrons come
from metal or ground; energy from a declared source.

**1. Single-electron pump** (Pothier et al. 1992). Two aluminium islands
between three tunnel junctions of 200 aF and 200 kΩ, two 50 aF gates
driven 90° apart, at 20 mK; charging energy e²/2C = 0.178 meV (2.1 K).
Solved with the orthodox theory of single-electron tunnelling (Averin &
Likharev 1986; Grabert & Devoret 1992, ch. 2): a master equation over the
islands' electron numbers with rates Γ = −ΔF/(e²R(1 − e^(ΔF/kT))).

| Drive | Electrons per cycle | Current | Heat per electron (paid by the gates) | Ledger |
|---|---|---|---|---|
| 1 MHz | 1.000000000 | 0.1602 pA | 12.7 µeV | −4×10⁻¹⁰ |
| 10 MHz | 0.999999977 | 1.602 pA | 53.1 µeV | −4×10⁻⁸ |
| 30 MHz | 0.999458 | 4.80 pA | 93.2 µeV | −4×10⁻⁷ |
| 100 MHz | 0.950532 | 15.2 pA | 160 µeV | −4×10⁻⁶ |
| 300 MHz | 0.667355 | 32.1 pA | 204 µeV | −1×10⁻⁵ |

Exactly one electron per cycle at low drive (I = e·f); errors grow when
the drive outruns the tunnelling (RC ≈ 40 ps). Uphill against 20 µV at
10 MHz it still pumps 0.999998 electrons per cycle, the gates supplying
72.8 µeV per cycle = 20.0 µeV lifted into the bias source + 52.8 µeV heat.

**2. Multipactor pump in a 1 GHz resonator** (Farnsworth US 2,091,439;
Kishek & Lau 1998). A 1 mm, 1 cm² copper gap (as-received yield, Baglin
2002) is the capacitor of a resonator (C 0.885 pF, L 28.6 nH, Q 1000)
driven by 0.779 mA for 140 V unloaded; a 10 % grid passes electrons to a
collector. The electron sheet loads the gap through Ramo's induced
current, and space charge caps it. The resonator alone agrees with a
kronspice transient to 0.3 % of V₀. In step (140 V) 1,000 starting
electrons grow to 2.0×10⁸, pull the gap voltage down to 39 V and stop
multiplying; the resonator recovers and a stray electron (one per 20 ns
while the gap is empty, declared) restarts it: bursts, as Kishek & Lau
found. Out of step (300 V) no electrons are left after 40 periods. Over
2,000 periods 2.16×10⁸ electrons are collected, 15.2 µA on average, at
1,426 eV of RF drive per collected electron. Ledger (nJ): drive 49.32 =
resonator heat 22.98 + impacts on the plates 31.43 + collector 1.78 +
stored change −6.88 (the bursts also drain what the resonator held at the
start), residual 4×10⁻⁷.

**3. Ground charge pump with a resonant drive** (§3j; Cockcroft & Walton
1932). 10 V at 100 kHz into a series resonator (Q = 20) → 3-stage
negative ladder → 100 nF bank with a 10 MΩ load, 600 cycles, solved by
kronspice. Step-up 10 V → 199.9 V peak (×20.0); bank −1,194.7 V against
the ideal −2·n·peak = −1,199.5 V. The earth lends 7.46×10¹⁴ electrons to
the bank's negative plate (plus 4.0×10¹² through the load) and gets the
same number back from the bank's earthed plate, so its net charge is
unchanged. Ledger: source 3.559 J = heat 3.425 J + stored 0.133 J (bank
0.071 J) + time-stepping error 2.7×10⁻⁵ J.

Picture: `viz/pumps.png`.

## 5. The multipactor pump on a microwave-oven magnetron, evolved (`magnetron.go`)

LG 2M214 oven magnetron (2,450 MHz, 900 W; efficiency 70 % declared, the
sheet gives it only as a chart; 34.6 W heater; 1,320 W from the mains) →
circulator → coupling port → copper re-entrant cavity (Q₀ = 2000) whose
capacitive gap multipacts → grid → collector. The circulator is standard
practice at 2.45 GHz: the magnetron must see a matched load (the sheet
rates VSWR ≤ 1.1), so what the cavity does not take goes to the
circulator's water load. The cavity's circuit twin is C = ε₀A/d, L
resonating it, wall loss R₀ = Q₀/(ωC), the port R_ext = R₀/β, and a
Norton source of the available power. The genetic algorithm (`evo`,
value genome, seed 2450, 40 designs × 30 generations) picks gap, gap area
(side ≤ λ/10 so the gap stays lumped), coupling β, grid transparency and
surface (as-received or conditioned copper, the two measured yields) for
the fewest mains watts per collected microamp, below Kilpatrick's spark
limit and 10 kV (non-relativistic).

| | Best of generation 0 | Evolved |
|---|---|---|
| Gap; area | 2.761 mm; 0.376 cm² | 1.519 mm; 1.500 cm² |
| Coupling β; grid transparency | 0.000139; 0.03 | 0.00153; 0.43 |
| Gap amplitude, unloaded; field / Kilpatrick | 1,038 V; 0.009 | 1,278 V; 0.020 |
| Collected current | 7.5 µA | 659 µA |
| Mains W per µA; electrons per mains joule | 176.6; 3.5×10¹⁰ | 2.00; 3.1×10¹² |
| Into the cavity: wall heat + electrons | 0.33 W: 0.30 + 0.06 | 2.05 W: 1.79 + 0.82 |
| To the water load | 899.7 W | 897.9 W |
| Ledger residual | 4×10⁻⁶ | 4×10⁻⁶ |

The search ran at ≥ 200 steps per period over 1,500 periods and the
results were re-run at 800 steps over 2,000 periods. One lesson from
building it: the first search found a fake 0.3 A that came from time
steps too coarse for fast electrons. The engine now takes at least 40
steps per gap crossing, and results agree to about ±10 % across step
sizes (burst timing). The gap needs only 2.1 W, so one magnetron could
feed about 438 such gaps (321 µA per watt into the cavity). Picture:
`viz/magnetron.png`.

## 6. Scaling to amps: many gaps on one magnetron

The genetic algorithm also chooses how many identical gaps (1 to
1,000,000) share the magnetron through an ideal, lossless, isolated power
divider (declared), each taking P_out/N.

| | LG 2M214, 2,450 MHz, 1 kW | Stellant CWM-100kW, 915 MHz, 100 kW |
|---|---|---|
| Gaps in parallel | 2,686 | 786,127 |
| Gap; area per gap | 0.56 mm; 1.50 cm² | 1.52 mm; 3.96 cm² |
| Total gap area | 0.403 m² | 311 m² |
| Gap voltage, unloaded | 191 V | 196 V |
| Collected current | **1.95 A** (726 µA per gap) | **216 A** (275 µA per gap) |
| Into the gaps / water load | 0.695 / 0.205 kW | 76.4 / 23.6 kW |
| Mains | 1.32 kW | 119 kW |
| Amps per kW of mains; per m² of gap | 1.48; 4.84 | 1.83; 0.70 |

Current follows the copper area in the field and the microwave power.
Scaling the better design: 1 A needs about 0.55 kW from the mains and
1.4 m² of gap plates, 10 A about 5.5 kW and 14 m², 100 A about 55 kW and
144 m², 300 A about 164 kW and 431 m². Feeding thousands of gaps needs a
real power-distribution structure (a large multi-gap cavity rather than a
divider), which is not modelled. The searches were still improving when
they stopped, so these are lower bounds. The electrons come from the
copper; every watt comes from the magnetron.

## The circuit mirror of the fdamador neutrino

The fdamador neutrino ([../fdamador_neutrino](../fdamador_neutrino), H28)
is this example's field turned inside out. In circuit form it is a series
RLC loop: E round the ring is the capacitor, B through the hole the
inductor, radiation a series resistor. Exchanging E and B is circuit
duality (voltage ↔ current, L ↔ C, series ↔ parallel, R ↔ R0²/R), so the
photon electron's mirror is a **parallel RLC tank**: B round the ring
(toroidal flux) is the inductor, E through the hole the capacitor, and
radiation a shunt resistor. With R0 = √(L/C) = 279.1 Ω of the neutrino's
bench loop the dual keeps L = 21.56 nH, C = 0.2767 pF, f0 = 2.0606 GHz and
Q = 40.51; R_rad = 6.862 Ω becomes R_p = 11.36 kΩ and the dielectric loss
27.9 mΩ becomes 2.79 MΩ (`viz/plot_circuit.py` checks that the two Q agree
to 10⁻¹²).

What the mirror is and is not: it is exact circuit algebra between the two
**models**, not a statement about the measured electron or neutrino. The
real dielectric torus is not dual (εr = 100, μr = 1): its B-ring mode rings
at 4.82 GHz with Q 112 in the field solver (fdamador_neutrino Part 7b,
grid r/2), not at the E-ring's 2.06 GHz. The H23 electron itself is a
travelling wave on a ring line, with E and B in phase (§2), and a lumped
tank stands for one of its modes at a time; under linear Maxwell its light
leaves within a quarter to half a period (§3). Neither circuit has mass or
the weak interaction, and neither derives spin ½.

- `viz/circuit.png`: the neutrino's series loop, the duality, and the
  electron's parallel tank, with values.
- `viz/cycle.gif` (`python3 viz/plot_cycle.py`): one period of the light
  (1.24×10²⁰ Hz, slowed about 4×10²⁰× to a 3 s loop): the wave travelling
  along the (2, 1) path, the field along both trips (the value after 360° is
  always minus the value at the start), and the mirror tank's energy moving
  between C and L. Companion of `../fdamador_neutrino/viz/cycle.gif`.

## Tests

| Test | Checks |
|---|---|
| `TestPaperNumbers` | eqs. 5, 8 and 12 of the paper |
| `TestRingModes`, `TestRoundTrips` | whole- and half-wave modes against the ladder; upright and inverted returns; ledgers close |
| `TestKnotLightLeaves` | on a coarse grid the light leaves, the field stays divergence-free, the ledger closes |
| `TestSingleElectronPump` | one electron per cycle at low drive, uphill too, the gates paying |
| `TestMultipactorPump` | resonator matches kronspice; growth in step, death out of step; ledger |
| `TestGroundPump` | the bank reaches −2·n·peak to 1 %; ledger |
| `TestMagnetronEvolution` | cavity matches kronspice; a short evolution improves on its start inside the limits |
| `TestGapArray` | N gaps are N copies of one gap fed P_out/N; the 100 kW tube's declared numbers |

## Files

| File | What it holds |
|---|---|
| `main.go` | Parts 1–2 (paper's numbers, ring lines) and the run order |
| `field.go` | Part 3, light on the double loop in `em/fdtd.go` |
| `pumps.go`, `pumps_report.go` | Part 4 pumps; reports for Parts 4–6 |
| `magnetron.go` | Parts 5–6: magnetrons, cavity, gap arrays, evolution |
| `main_test.go` | known-value, twin and ledger tests for every part |
| `export_test.go` | ring-line netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt` | the run above |
| `viz/render.py`, `viz/*.csv` | `round_trip.png`, `knot_light.png`, `pumps.png`, `magnetron.png` and their data |

## Sources

- J. G. Williamson, M. B. van der Mark, "Is the electron a photon with toroidal topology?", *Ann. Fond. Louis de Broglie* 22, 133 (1997).
- X. Fan et al., "Measurement of the electron magnetic moment", *Phys. Rev. Lett.* 130, 071801 (2023); CODATA 2018.
- K. S. Yee, *IEEE Trans. Antennas Propag.* 14, 302 (1966); A. Taflove, S. C. Hagness, *Computational Electrodynamics*, 3rd ed. (2005).
- R. W. Hellwarth, P. Nouchi, "Focused one-cycle electromagnetic pulses," *Phys. Rev. E* 54, 889 (1996).
- H. Pothier, P. Lafarge, C. Urbina, D. Esteve, M. H. Devoret, "Single-electron pump based on charging effects," *Europhys. Lett.* 17, 249 (1992); D. V. Averin, K. K. Likharev (1986); H. Grabert, M. H. Devoret (eds.), *Single Charge Tunneling* (Plenum, 1992), ch. 2.
- P. T. Farnsworth, "Multipactor oscillator and amplifier," US 2,091,439 (1937); R. A. Kishek, Y. Y. Lau, "Interaction of multipactor discharge and rf circuit," *Phys. Rev. Lett.* 80, 193 (1998); V. Baglin et al., LHC Project Report 472 (2002); A. J. Hatch, H. B. Williams (1954) and J. R. M. Vaughan, *IEEE Trans. Electron Devices* 35, 1172 (1988): the gap model of `em/multipactor.go`; S. Ramo's induced-current theorem.
- W. D. Kilpatrick, "Criterion for vacuum sparking designed to include both rf and dc," *Rev. Sci. Instrum.* 28, 824 (1957).
- J. D. Cockcroft, E. T. S. Walton, *Proc. R. Soc. A* 136, 619 (1932).
- LG 2M214 magnetron specification sheet; Stellant Systems CWM-100kW sell sheet (May 2023).
