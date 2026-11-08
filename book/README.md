# The kronspice book

A guided tour of the engine through the examples built on it. Each chapter
takes one part of kronspice, states the question it answers, shows a
netlist you can run with the command-line tool and what it prints, then
the Go program that goes beyond what a netlist can say, the results
against measurement or the standard prediction, and the sources. The
chapters follow the order of CLAUDE.md, so the book is also the long-form
reading of that design document.

How to read it: the netlists in every chapter run from the repository
root with `go run ./cmd/kronspice file.cir` (or `kronspice` once
installed); the Go examples run with `go run ./examples/<name>`. Numbers
quoted in the text are the programs' own output at the time of writing,
and the per-example READMEs under `examples/` carry the full tables.

**Rule.** Every example folder is cited by at least one chapter, and every
cited example exists; `go test ./book` enforces both. When an example is
added, it gets a row in `examples/README.md` and a place in a chapter.

## Contents

| # | Chapter | Examples | Status |
|---|---|---|---|
| 1 | [Circuits and the ledger](01-circuits-and-the-ledger.md) | rc_lowpass_values, negative_resistance_gain, rectifier_filter | written |
| 2 | [Evolution](02-evolution.md) | amplifier_gain, filter_topology, mos_amplifier_pareto, torn_ga | written |
| 3 | [Kron's machines](03-krons-machines.md) | kron_machine, parametric_pump, parametric_lcr, homopolar, homopolar_bearing, bearing_motor, kron_nmr | written |
| 4 | [Power and polyphase](04-power-and-polyphase.md) | power_filter | written |
| 5 | [Fields from Maxwell](05-fields-from-maxwell.md) | pcb_layout, rotating_efield, rate_of_change | written |
| 6 | [Coils, pulses and gaps](06-coils-pulses-and-gaps.md) | bial_bifilar, pulsed_bifilar, water_lamp, sparkgap_diode, sparkgap_diodebridge | written |
| 7 | [Networks beyond circuits](07-networks-beyond-circuits.md) | acoustic_drive, water_vibration, dna_vibration, rna_vibration, kron_schrodinger | written |
| 8 | [Craft and space](08-craft-and-space.md) | tether, coiled_tether, ring_towers, double_hull, machine_spacetime, earth_machine | written |
| 9 | [The electron test bench](09-the-electron-test-bench.md) | kelvin_ring, vortex_ring, space_vortex, kron_electron, photon_electron, light_hydrogen | written |
| 10 | [Hypotheses under test](10-hypotheses-under-test.md) | eckardt_homopolar, parametric_lcr, pais_cavity, space_power_generator, octave_magnet, toroidal_coherence, barrett_osc, magnetostrictive_transformer, al_cell_pump, rodin_resonance | written |

## The shape of a chapter

1. **The question.** One paragraph: what this part of the engine is for
   and what the chapter's examples ask.
2. **Run it.** A netlist (or the example's exported one) and the report
   it prints, read line by line.
3. **Beyond the netlist.** What the Go program adds: elements the netlist
   cannot express, the Kron frame, the field solve, the GA loop, the
   standard prediction beside a hypothesis.
4. **Results.** Tables against analytic values, measurements or ngspice,
   with the ledger's residual.
5. **What to keep.** The rules the chapter illustrates.
6. **Sources.** Only works the code or its documentation cite.
