# CLAUDE.md — examples/harmonic_cavitation (dimple-drum cavitation heater, closed steam loop, Mandelbulb challenger)

Notes for Claude sessions working in this example. The repo's root CLAUDE.md
governs everything here: standard prediction beside every hypothesis, the
energy ledger closes (§3c), every element carries the full port set (§3.1),
cite source and equation for every formula, declared inputs stay labelled.
New example ⇒ a row in `examples/README.md`, its own `README.md`, and a
citation in a `book/` chapter (memory: examples-index-row).

- Papers and the numbers taken from them: [LITERATURE.md](LITERATURE.md).
- Open questions for the user: [QUESTIONS.md](QUESTIONS.md).
- Report: `harmonic_cavitation_report.pdf`, built by `viz/make_report.py`.
- Concept drawing (not in the repo):
  https://claude.ai/artifact/SMdn1Lxob8T9vnKcZtiNhc

## The machine (user's decisions, 2026-11-08)

```
 WATER IN → [dimple drum] ─ TESLA TURBINE (steam) ─ [dimple drum] ← WATER IN
            STEAM OUT ← along the shaft → condenser (home heat) → pump back
```

- **Two dimple drums**, one at each outlet, on one shaft, each turning in a
  fixed housing.
- The drums carry rows of hemispherical dimples facing a matching ring of
  dimples in the housing across a small gap.
- **Closed steam loop:**
  1. Water enters at each drum's outer top.
  2. The drums' cavitation and shaft work raise steam.
  3. The steam drives the Tesla turbine on the same shaft.
  4. The exhaust is condensed at the outlet (the home's heat).
  5. The condensate is pumped back.
  6. A motor on the shaft supplies what the turbine does not.
- **The goal is the most efficient cavitation generator.** η = collapse
  energy, (p − p_v)·V_vapour per cavity cycle, per joule of shaft work. Heat
  per joule is no measure: all shaft work ends as heat whatever the cavity.
- **The Mandelbulb (H46)** is the user's hypothesis: cavities shaped as the
  White–Nylander power-n fractal, whose harmonics and implosions gather in
  the middle. It is evolved over n as a challenger against the best dimple
  drum.
- **Collapses gathering in the middle.** Focusing raises the peak pressure
  on the axis by up to the number of simultaneous collapses (gcd of rotor
  and housing counts). It never raises the energy beyond the collapses' sum,
  which the shaft pays.

The rotor's cavities must move past a fixed surface or fixed cavities
across a small gap (every reactor in LITERATURE.md §1). Water that turns
with its walls makes no cavitation and no heat.

## Scope rule: CFD outside, reduced model inside

Root §1 non-goals keep CFD and bubble dynamics out of the core. Inside
kronspice:
- the cavitation fitness is a **reduced model calibrated on the published
  reactors** (Xue 2024, Sun 2020/2021, Petkovšek 2013);
- the bubble is Rayleigh–Plesset with a ledger;
- the gap acoustics is a kronspice network;
- the shaft is a kronspice circuit (mobility analogy: speed across, torque
  through) whose ledger books every watt.

The best designs are to be checked by **external CFD** (OpenFOAM
`interPhaseChangeFoam`, Schnerr–Sauer, SST k-ω, as Sun 2020/2021), imported
like bench data, never run inside the engine.

## Model constants (fitted or declared; see README and the report)

- **Fitted on Xue 2024** (1 atm declared at the rotor):
  - inception σ 0.968–1.243 (σ_i 1.105);
  - best σ 0.775;
  - pulsation St 0.446;
  - vapour law g = x·e^(1−x), x = (σ_i − σ)/Δ, Δ 0.272.
- **Declared, run at both ends in Sensitivity:**
  - κ_v 0.2;
  - C_d 0.25;
  - housing pulse 0.3·½ρU²·(1.5 mm/gap)^0.60, the exponent fitted to
    Petkovšek;
  - St also 0.28; designs are scored at the worse St.
- **Friction:** max(laminar Couette, turbulent flat plate) on the drum
  surface and both end faces (hub 10 mm).
- **Steam:** IAPWS-IF97 saturation (tested) and IAPWS-95 table states.
- **Declared drive efficiencies:** Tesla turbine 0.30 of ideal, motor 0.90,
  pump 0.70.

## Files

- `sizing.go`: water properties, Xue calibration, Wood's mixture sound
  speed.
- `model.go`: the reduced cavitation model of the dimple drums, `Hemicap`,
  `DimpleMachine`, the kronspice shaft circuit.
- `ga.go`: NSGA-II through `evo` (objective read from the shaft circuit's
  operating point), `Heater`, `Sensitivity`, `PerN`.
- `acoustics.go`: the gap as an acoustic ring (kronspice AC).
- `bubble.go`: Rayleigh–Plesset with a ledger, and focusing.
- `heat.go`: motor drive against pump + turbine.
- `bulb.go`: the Mandelbulb, both angle conventions, valleys as cavities.
- `loop.go`: the closed steam loop (Part 10).
- `house.go`: the home at −10 °F (Part 11).
- `viz.go`: `-viz` writes the report's CSVs.
- `main.go`: parts 1–11.
- `viz/plot_all.py` and `viz/make_report.py`.
- Tests: `main_test.go`, `model_test.go`, `loop_test.go`, `house_test.go`.

## Results (2026-11-08; numbers in README.md and run_output.txt)

- **Dimple drums:** η = 0.189, flat from 1.1 to 3.0 kW.
  - Heater design: R 84 mm, 2 rows of Ø 6 mm hemispheres, a housing
    dimple ring, 0.8 mm gap, 1,830 rpm, σ 0.88.
  - A plain housing gives 0.117.
  - The rules hold in every sensitivity rerun: full hemispheres, housing
    dimples, the smallest gap, σ ≈ 0.88.
  - The absolute η ranges 0.09–0.47 with the declared constants.
- **Gap ring:** the kronspice modes equal the lattice's. The dimple pattern
  cannot reach them (Wood's minimum sound speed is above the drum speed).
- **Collapse:** Rayleigh time to 3×10⁻⁴; ledger 2×10⁻⁶. Focusing gives up
  to ×26 peak with equal counts.
- **Drive:** a motor on the shaft gives 4× the cavitation per electrical
  joule of pump + Tesla turbine.
- **Closed steam loop (Part 10):**
  - The turbine returns 3–4 % of the shaft work (Carnot 13 % between
    111 °C and 60 °C), so the motor supplies the rest.
  - Boiling in the drums drops η to 0.011.
  - Subcooled drums at 3 bar plus a flash tank keep 0.156.
  - Thermal suppression of collapse in hot water is not modelled (upper
    bounds).
- **Mandelbulb (H46):** best 0.187 (n = 14); per n 0.158–0.186; without
  its fine detail 0.191. A smooth lobed ring matches dimples, and the
  fractal detail does not help.

Not done:
- OpenFOAM case writer;
- STL/DXF export;
- `go test -race` (needs cgo; run it in Docker).

## Home heating decided (user, 2026-11-08): 3,000 ft² at −10 °F, 24/7

- **The drum speed is matched to the demand** (user's choice). A VFD sets
  the rpm so that the shaft work equals the house's heat need. Rows are
  added until the design day runs at loop B's best speed.
- **Built:** `house.go`, Part 11.
  - Declared: 15 / 25 / 40 BTU/h·ft² loads (replace with Manual J);
    condenser U 2,500; EN 1264 floor and EN 442 radiators; heat pump COP
    1.5 at −10 °F.
  - Average house (22.0 kW): 12 rows of Ø 6 mm on each R 84 mm drum, a
    115 mm band, 1967 rpm, 24.5 kW electric.
  - Above about 30 °F the matched speed falls below inception, so the
    drive cycles on and off.
- **User question: can the hole size be changed so the Tesla turbine runs
  faster?** No.
  - The drums' best speed is set by the rim speed through σ, not by hole
    size: Ø 3–25 mm all give 2000 rpm at R 84 mm. A smaller drum turns
    faster (4200 rpm at R 40 mm).
  - The turbine wants its tip near half the jet speed: 20,700–40,100 rpm
    at R 84 mm.
  - Turbine efficiency cannot change the electricity per unit of home heat
    (0.899 at η 0.3–1.0, tested). Its work circulates and ends as heat.
  - A faster turbine therefore matters only if it drives something other
    than this heater's own shaft.
- **Not yet in the PDF report.** It still covers parts 1–10.
