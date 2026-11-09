# Questions for the user (saved 2026-11-08 while all stages ran)

Each question says what it changes. Answer in any order.

## Design

1. **Drive.** The pump + Tesla turbine drive gives a quarter of the
   cavitation per electrical joule of a motor on the shaft, because a
   Tesla turbine in water is inefficient (declared 30 %, to be checked).
   Its losses still heat the water, so heat per joule is about the same.
   Options:
   - keep the turbine as the drive;
   - drive the shaft with a motor and use the Tesla disks as the pump that
     moves water through the gap;
   - use both.
2. **Drums (answered 2026-11-08): a dimple drum at each outlet.** Still
   open:
   - One drum at a single outlet, or two (one per outlet, as modelled)?
     The efficiency is the same; only the size per drum changes.
   - The heater design wants 84 mm radius drums at 1,830 rpm for 3 kW.
     Smaller drums turn faster: 39 mm at 3,930 rpm gives the same
     efficiency at 1.1 kW.
3. **Hollow drum.** Collapses lined up on a ring focus on the axis. In a
   solid drum that is metal, which means fatigue rather than useful water
   pressure. Should the drum be hollow with a water core, so the focus
   is in the water, with the steam outlet along the axis?
4. **Housing dimples.** The model says they help (η 0.12 → 0.19). Can
   the housing be machined or printed with a matching dimple ring?
5. **Limits.**
   - Power: 3 kW of shaft power was assumed for a home heater.
   - Size: drums of 31–91 mm radius came out of the searches.
   - Speed: up to 15,000 rpm was allowed, and the designs ran at
     1,600–6,700 rpm.
   - Gap: 0.8 mm was the smallest in the data and the search always went
     there. Can you hold 0.8 mm, or smaller, at speed?
   - Materials and erosion: stainless? bronze? replaceable drum shell?
6. **Output.** Steam at 1 atm, or pressurized hot water? A steam loop is
   a pressure vessel (relief valve and rating needed).
7. **Steam loop (Part 10).**
   - **Boiling or flash?** Should the drums boil (your sketch; cavitation
     then drops to η 0.011), or run pressurized at about 3 bar with a
     flash tank (η 0.156)? The flash version needs a throttle valve and a
     small tank between the drums and the turbine.
   - **Turbine efficiency.** The Tesla turbine on steam was declared at
     30 % of ideal. A measured figure would replace it.
   - **Condenser.** It was set at 60 °C (0.2 bar) so it can feed
     radiators. A colder condenser gives the turbine a little more work but
     cooler heat for the home.
   - **The turbine returns only 3–4 % of the shaft work.** Keep it, or
     drive the shaft with the motor alone and use the steam only for heat?
     The home heat is the same either way.

## Mandelbulb (H46)

7. **Where the bulb goes.** In the model the bulb is the profile of the
   spinning drum, its lobes facing the housing. A closed cavity inside a
   spinning part has no moving surface across a gap, so it would not
   cavitate by rotation. Is the drum-surface reading what you meant, or
   should an inside version be built with a fixed inner core?
8. **Convention.** Polar (arccos, the usual pictures) was used. The
   latitude (arcsin) form is also built; its equator is the Multibrot
   set. Any preference?
9. **Detail.** Detail smaller than 0.3 mm (metal printing) was removed.
   What resolution can your printer or machine hold?

## Data I could not get

10. Full texts I could not read (only abstracts or summaries):
    - Sun et al. 2021 (CGU height definition, measured vapour);
    - Kwon & Yoon 2012 (heater geometry and calorimetry);
    - Liu & Katz 2008;
    - the Fluids 2024 paper (403);
    - the Griggs patent text.
    If you can get the PDFs into `reference/`, the declared constants
    (C_d, κ_v, C_pulse) can be fitted instead of declared. They are the
    main uncertainty: η moves 0.12–0.46 with C_d alone.
11. **Tesla turbine efficiency in water.** A measured figure (paper or
    your own bench) would replace the declared 30 %.

## Repo

12. **Hypothesis register.** H46 is added to the root CLAUDE.md register
    as *bounded* (this model). Change the wording if you want.
13. **Book test.** `book/` fails its test for nine other examples that no
    chapter cites. That was the case before this work: AIAS,
    einstein_teleparallel, ether_cavitation, fdamador_neutrino,
    forward_casimir, hoffmann_motor, kaluza_kron, shoulders_ev,
    zpf_hydrogen. harmonic_cavitation is now cited in chapter 7.
14. **Race test.** `go test -race` needs cgo, which is off on this
    machine. Run it in the Docker tools image (CLAUDE.md §7).
