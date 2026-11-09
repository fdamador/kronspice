# Harmonic cavitation: an efficient cavitation heater with dimple drums

The machine:
- two dimpled metal drums, one at each outlet, on one shaft;
- in the user's loop, water enters at each drum's outer top; the drums
  raise steam; the steam drives a Tesla turbine on the same shaft; it is
  condensed at the outlet (the home's heat) and pumped back; a motor
  supplies what the turbine does not;
- each drum spins inside a fixed housing, and its rows of hemispherical
  dimples face a matching ring of dimples in the housing across a small
  gap;
- water is cavitated in the gap and leaves through the end caps as hot
  water or steam, to heat homes.

The goal is the most cavitation per joule of shaft work. The layout (a
dimple drum at each outlet) is the user's choice of 2026-11-08. The user's Mandelbulb cavity (H46) is evolved
over its power n as a challenger.

- Plan and rules: [CLAUDE.md](CLAUDE.md).
- Papers and the numbers taken from them: [LITERATURE.md](LITERATURE.md).
- Questions left for the user: [QUESTIONS.md](QUESTIONS.md).
- Report: [harmonic_cavitation_report.pdf](harmonic_cavitation_report.pdf).

```
go run ./examples/harmonic_cavitation           # all ten parts, about 10 s
go run ./examples/harmonic_cavitation -part 8   # one part
go test ./examples/harmonic_cavitation
go run ./examples/harmonic_cavitation -viz examples/harmonic_cavitation/viz   # report data
python examples/harmonic_cavitation/viz/plot_all.py
python examples/harmonic_cavitation/viz/make_report.py
```

Full output: [run_output.txt](run_output.txt).

## What "efficient" means here

η = collapse energy per joule of shaft work. Collapse energy is
(p − p_v)·V_vapour for each cavity cycle, summed over cavities and cycles.

Heat per joule cannot be the measure. Every joule of shaft work ends as
heat, sound and steam whatever the cavity, so heat per joule does not
separate good cavitators from bad ones.

## The model (reduced, calibrated, not CFD)

**Calibration on Xue et al. 2024** (rotor Ø 80 mm, 20 grooves; 1 atm
declared at the rotor):
- cavitation starts between σ = 0.968 and 1.243 (σ_i = 1.105 is taken);
- the pulsation is strongest at σ = 0.775;
- the pulsation Strouhal number fW/U is 0.446.

**Vapour law.** Vapour per cavity is κ_v·g(σ)·V, with g = x·e^(1−x) and
x = (σ_i − σ)/Δ. Δ = 0.272 is fitted so that Xue's collapse rate peaks at
its best speed.

**Collapse energy per event.** (p_local − p_v + housing pulse)·V_vapour.
The housing pulse is C_pulse·½ρU²·(1.5 mm/gap)^0.60. The exponent is
fitted to Petkovšek's 0.7 bar at 3.5 mm and 1.7 bar at 0.8 mm.

**Event rate.**
- On a plain housing, it is the cloud rate St·U/W.
- With housing dimples, it is the passing rate of the other ring, capped
  at the cloud rate.
- In both cases it is never faster than the Rayleigh collapse allows.
- Every design is scored at the worse of St = 0.28 and 0.446.

**Shaft power.**
- Wall friction: the larger of laminar Couette μU/gap and a turbulent
  flat-plate c_f, on the drum surface and both end faces (hub 10 mm).
- Cavity drag: ½ρU²·C_d·A per cavity, with C_d = 0.25 declared.

**Energy consistency.** Collapse energy is part of the cavity drag, and
any design with more collapse energy than drag is rejected.

**The shaft circuit.** In kronspice the shaft is a speed source ω driving
one conductance per loss. The GA ranks designs by the power the operating
point books into the collapse branch over the power the source gives. The
ledger closes to 0.

## Results

**1–2. Calibration and drum sizing.** A dimple works only within a window
of σ = (p − p_v)/(½ρU²). At 1 atm cavitation starts near U = 13.3 m/s
and is strongest near 15.9 m/s.

| Drum radius | Cavitation from | Strongest at |
|---|---|---|
| 50 mm | 2,550 rpm | 3,040 rpm |
| 84 mm | 1,520 rpm | 1,810 rpm |

Every dimple on a drum runs at the same speed, so one shaft speed suits
them all. The resonance rule (N = 2πr·St/W) gives 18.7 grooves on Xue's
rotor against the 20 built.

**3. Reference drums.** R 50 mm, 10 rows of Ø 10 mm, 3050 rpm.
- η = 0.155 at 13.7 kW.
- A plain housing gives 0.110.
- The kronspice shaft circuit gives the same power, with ledger residual 0.

**4. NSGA-II (80 × 80, shaft ≤ 3 kW).** The front is flat in η
(0.1886–0.1891) from 1.1 to 3.0 kW, so the drums can be sized to the load.

| | Heater design (most collapse power within 0.5 % of the best η) |
|---|---|
| Drum radius | 84.3 mm |
| Rows | 2 |
| Dimples | 26 per row, Ø 6.0 mm hemispheres |
| Housing | 73 dimples per ring |
| Gap | 0.8 mm |
| Speed | 1,828 rpm (16.1 m/s) |
| Outlet | 1.18 bar |
| σ | 0.881 |
| Result | η = 0.1886: 565 W of collapse energy from 2,998 W of shaft work |
| Same drums, plain housing | η = 0.117 |

The rules that hold in every rerun:
- full hemispheres;
- dimples on the housing;
- the smallest gap in the data's range (0.8 mm);
- σ ≈ 0.88, a slower speed than the most collapse per second (σ 0.775).

The absolute η depends on the declared constants:

| Changed constant | Best η |
|---|---|
| C_d 0.10 / 0.40 | 0.47 / 0.12 |
| C_pulse 0.15 / 0.50 | 0.16 / 0.23 |
| κ_v 0.10 | 0.09 |
| Xue's St only | 0.31 |

**5. The gap as an acoustic ring (kronspice AC).**
- The modes equal the lattice's exact frequencies to 1e-4.
- The 26-dimple pattern passes at 792 Hz; its own mode 26 is at 73.7 kHz
  in clear water.
- The two would meet only if the bubbly mixture's sound speed fell to
  16.1 m/s. Wood's minimum at 118 kPa is 21.7 m/s, so they never meet.

**6. One collapse and the axis.**
- Vapour 11.1 mm³ (R0 1.38 mm) under 171 kPa.
- Collapse time 96.66 µs against Rayleigh's 96.69 µs; ledger 2.1e-6.
- With equal rotor and housing counts all 26 collapse together, and the
  axis peak rises by up to ×26 with the same energy.
- These are upper bounds. In a solid drum the axis is metal.

**7. Energy books** (declared: motor 0.90, pump 0.70, Tesla turbine 0.30).

| Drive | Electric | Heat to water | Collapse per electric J | Steam |
|---|---|---|---|---|
| Motor on the shaft | 3.3 kW | 3.0 kW | 0.170 | 4.1 kg/h |
| Pump + Tesla turbine | 14.3 kW | 12.8 kW | 0.040 | 17.7 kg/h |

A heat pump on the motor drive's electricity would give 11.7 kW of heat.

**8. Mandelbulb challenger (H46).**
- **Geometry checks** (tested): (n − 1)-fold symmetry in both
  conventions; in the latitude convention the equator is the Multibrot set;
  n = 8 shows 7 lobes.
- **Evolved over n and the machine:** best η = 0.187 (n = 14, 31 mm
  drums, 6,740 rpm).
- **Each n evolved separately:** η 0.158–0.186.
- **Fine detail removed:** 0.191, as good as the dimples. A smooth lobed
  ring matches dimples, and the fractal detail lowers η.
- **n = 8 on the heater's drums** without re-optimizing: 0.099. The bulb
  needs its own size and speed.

**10. Closed steam loop (the user's layout).**

The cycle:
1. Water enters at each drum's outer top.
2. The drums raise steam.
3. The steam drives the Tesla turbine on the same shaft.
4. The exhaust is condensed at the outlet, which is the home's heat.
5. The condensate is pumped back.

The model uses IAPWS-IF97 saturation and IAPWS-95 steam tables. Declared:
turbine 0.30 of ideal, condenser 0.2 bar (60 °C).

| Cycle (3 kW shaft) | Steam | Turbine returns | Electric | Home heat | η |
|---|---|---|---|---|---|
| A: boil in the drums (1.5 bar) | 4.4 kg/h | 115 W | 3.21 kW | 2.89 kW | 0.011 |
| B: subcooled drums (3 bar) + flash tank | 4.5 kg/h | 93 W | 3.24 kW | 2.91 kW | 0.156 |

- **The loop cannot run itself.** Carnot allows 13 % between 111 °C and
  60 °C; the Tesla turbine returns 3–4 % and the motor supplies the rest.
- **The home gets about the electricity**, as with no turbine.
- **Boiling wastes the cavitation.** Steam that leaves has not collapsed,
  and p − p_v → 0 as the water nears boiling.
- **Pressurized drums plus a flash tank keep it.** The books close, and
  kronspice's shaft circuit (motor = speed source, turbine = torque source)
  agrees.

### 11. A 3,000 ft² home at −10 °F, 24/7 (`house.go`)

Loop B (subcooled drums at 3 bar, flash tank, Tesla turbine, condenser).
The drum speed is matched to the demand: a variable-frequency drive sets
the rpm so that the shaft work equals the heat the house needs. Rows are
added until the design day runs at the loop's best speed.

Declared inputs (replace with the real house):
- indoor 70 °F, outdoor −10 °F;
- 3,000 ft², with loads of 15 / 25 / 40 BTU/h per ft² (tight / average /
  leaky). A Manual J calculation replaces these;
- condenser U = 2,500 W/m²·K;
- heat-pump COP 1.5 at −10 °F.

Panel radiators, 70/50 °C water, condenser at 1 atm (100 °C):

| House | Load | Shaft | Rows (Ø 6 mm) | Drum band | rpm | Steam | Electric | kWh/day | Heat pump kWh/day | Condenser | Radiators (rated ΔT50) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| tight | 13.2 kW (45,000 BTU/h) | 13.3 kW | 7 | 67 mm | 1985 | 21 kg/h | 14.7 kW | 352 | 211 | 0.13 m² | 18.3 kW |
| average | 22.0 kW (75,000 BTU/h) | 22.1 kW | 12 | 115 mm | 1967 | 35 kg/h | 24.5 kW | 587 | 352 | 0.22 m² | 30.5 kW |
| leaky | 35.2 kW (120,000 BTU/h) | 35.3 kW | 18 | 173 mm | 2011 | 56 kg/h | 39.1 kW | 939 | 563 | 0.36 m² | 48.8 kW |

**Radiant floor** (45/35 °C water, condenser at 0.2 bar): the floor needs
152 / 254 / 407 m² at its 29 °C surface limit. The house has 279 m², so
the leaky house cannot be heated by floor alone.

**Over the winter** (average house, 12 rows):
- −20 °F: 2046 rpm, slightly over the design speed.
- −10 °F: 1967 rpm.
- 20 °F: 1682 rpm, η 0.05.
- From about 30 °F up, the matched speed falls below cavitation inception,
  so the drive must cycle on and off at the best speed instead of slowing
  further.

**Hole size and the Tesla turbine:**
- **The dimple size does not move the best speed.** Ø 3–25 mm all give
  2000 rpm at R 84 mm. The rim speed sets it (17.6 m/s at 3 bar, through
  σ). A smaller drum turns faster: 4200 rpm at R 40 mm.
- **The turbine wants about 10× the drums' speed.** Its tip should be near
  half the steam-jet speed (365–708 m/s), i.e. 20,700–40,100 rpm at R 84 mm.
  On a shared shaft it runs far below that.
- **Turbine efficiency does not change the electricity a given home heat
  needs.** At η 0.3, 0.7 or 1.0, home heat per unit of electricity stays
  0.899. That is the motor efficiency; the rest is motor heat in the
  basement. Turbine work only circulates round the loop and ends as heat
  (tested).
- **So the machine is a resistance heater (COP ≈ 1).** A cold-climate heat
  pump needs about 1/1.5 of the electricity at −10 °F, and less in milder
  weather.

## What is declared, not measured

- 1 atm at Xue's rotor.
- κ_v = 0.2, C_d = 0.25 and C_pulse = 0.3.
- Gas content 1e-3.
- Motor, pump and Tesla turbine efficiencies.
- Hub radius 10 mm.

Next: CFD outside kronspice (OpenFOAM, Schnerr–Sauer, as Sun et al. 2020),
then a bench with KI dosimetry and a power analyser.
