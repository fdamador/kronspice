# Examples

Every example is a self-contained program (or an `evolve` config) that uses
the engine on one question, with the energy ledger on and the standard
prediction stated beside any hypothesis. Each folder has its own README with
the construction, the results and the sources. This file is the index.

**Rule.** When a new example is added, add a row here in the right group,
with one sentence of scope, the CLAUDE.md section or register hypothesis it
serves, and how to run it.

Run commands are given from the repository root. Config-driven experiments
run through `go run ./cmd/evolve -config examples/<name>/config.json`;
the rest through `go run ./examples/<name>`, and all have
`go test ./examples/<name>`.

## Genetic algorithm experiments

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [rc_lowpass_values](rc_lowpass_values) | Fixed RC topology; evolve two E-series values for a 1 kHz cutoff. The first GA example (M5). | §4, M5 | config |
| [rectifier_filter](rectifier_filter) | Half-wave rectifier; evolve the filter for mean output voltage and ripple (M4 diode, Newton). | §4, M4 | config |
| [negative_resistance_gain](negative_resistance_gain) | A resistive divider asked for gain: evolution with signed values finds negative resistance; its power is booked as model-supplied, never as inflow. | §3c, §4 | config |
| [amplifier_gain](amplifier_gain) | 2N3904 common-emitter stage; eight standard values evolved for 30 dB gain, 100 Hz to 200 kHz, bias targets (M7). | §4, M7 | config |
| [filter_topology](filter_topology) | Low-pass filter grown from an empty embryo by the topology genome; mask objective and parsimony (M8). | §4, M8 | config |
| [mos_amplifier_pareto](mos_amplifier_pareto) | NMOS common-source amplifier; NSGA-II gain / bandwidth / power front (M9). | §4, M9 | config |
| [torn_ga](torn_ga) | The same GA run by MNA, untorn Kron and Kron diakoptics torn at one resistor with cached pieces; same best design, timings (M9b). | §3a, M9b | run |
| [power_filter](power_filter) | Harmonic filters for a six-pulse rectifier, scored by geometric-algebra power (IEEE 1459 twin); PF 0.903 to 0.957 (M9f). | §3b, M9f | run |
| [pcb_layout](pcb_layout) | Microstrip pair whose L, C and K come from `em/extract` field solves; evolved for 50 Ω and low crosstalk (M9i). | §3c, M9i | run |

## Kron's tensor analysis and machines

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [kron_machine](kron_machine) | Kron's moving inductors: induction-motor start-up in frames I and II; a parametric inductance against Mathieu (M9c.1). | §3a, M9c.1 | run |
| [hoffmann_motor](hoffmann_motor) | Hoffmann 1949 on Kron's non-Riemannian electrodynamics: a single-phase squirrel-cage motor on rotating and on stationary rotor coils, the torsion term, three solutions agreeing. | §3a, M9c.1 | run |
| [kaluza_kron](kaluza_kron) | Kaluza 1921: charge as motion along a fifth dimension (exact 5-D geodesic = Lorentz force), his electron difficulty, Klein's circle; the same construction as a Kron machine with cyclic winding charges. | §3a, §3f | run |
| [einstein_teleparallel](einstein_teleparallel) | Einstein 1928 distant parallelism beside Kron: the motor's frame as a 4-bein, one torsion in three notations, φ and curvature, Einstein's linear field equations checked, and a kronspice motor built from the torsion alone. | §3a, M9c.1 | run |
| [parametric_pump](parametric_pump) | Netlist of the LC tank pumped by a capacitor modulated at twice resonance; envelope growth against Mathieu's rate, pump work on the ledger. | §3a, M9c.1 | netlist |
| [parametric_lcr](parametric_lcr) | Amador & Eckardt 2012 series R, L, C(t): every case of the paper with the capacitor pump's power on the ledger. | H16 | run |
| [homopolar](homopolar) | Tesla US 406,968 homopolar dynamo as machine elements from the Lorentz force; EMF, torque, shaft balance; GA on geometry (M9n). | §3f, M9n | run |
| [homopolar_bearing](homopolar_bearing) | Tesla two-disk homopolar with load current through its ball bearings: bore versus magnetic area, bearing drop and damage limit. | §3f, H13, H14 | run |
| [eckardt_homopolar](eckardt_homopolar) | Eckardt 2012 self-excited twin-disk homopolar: the paper's equations beside the machine element with its back torque. | H15 | run |
| [bearing_motor](bearing_motor) | Huber ball-bearing motor: measurement registry, κ and ω0 fits, whether bearing torque can explain the SPG. | H13, H14 | run |
| [space_power_generator](space_power_generator) | Tewari's SPG test (UPSM pp. 159–160) through the energy ledger: the unmeasured loaded-EMF assumption. | H13 | run |
| [kron_nmr](kron_nmr) | NMR and ESR as Kron's frame II(b): Bloch's equations as a three-winding machine, rotating frame, T1/T2, pickup coil. | §3b, H21 | run |
| [machine_spacetime](machine_spacetime) | One rotor in three gravitational fields: clocks, precession, Sagnac, burst and tidal margins; its own curvature and frame dragging. | §3l | run |

## The electron test bench and the hypothesis register

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [kelvin_ring](kelvin_ring) | Kelvin's 1899 charged-ring molecule: induction (eΦ/2π, Larmor) and his Zeeman test, gyrostat against the full Lorentz force. | §6 (Kelvin) | run |
| [vortex_ring](vortex_ring) | Kelvin's 1867 vortex atoms through his ring-electromagnet analogy: ring speed from inductance, leapfrogging rings. | §6 (Kelvin) | run |
| [kron_electron](kron_electron) | A ring of circulating charge as a Kron machine: Kelvin's ring twin test, frame-I torque, electron rings against Dirac (g ≤ 1). | H21 | run |
| [space_vortex](space_vortex) | Tewari's space-vortex electron: charge, moment, g = 15/8, hydrogen orbit, against measurement. | H13 | run |
| [photon_electron](photon_electron) | Williamson & van der Mark toroidal-photon electron; light on the double loop in FDTD; three working electron pumps with ledgers; magnetron-driven multipactor arrays. | H23, §3i | run |
| [fdamador_neutrino](fdamador_neutrino) | The toroidal field inverted (E is the ring, B circles it) as a neutrino candidate: neutral, oscillating; E-ring resonator in kronspice (Q ∝ εr^1.5, ledger), SN1987A confinement bound, inverted (CL) plain and crossed ring lines; zero-point inflow balancing the radiation outflow; resonance chart and sound. | H28 | run |
| [knotted_light](knotted_light) | Rañada's knotted light (hopfion) in the field solver against the exact solution: helicity and linking kept, energy not held. Particles as knots, step 1. Implosion then explosion fed through a Huygens surface, and a pulsing knot with inflow = outflow on the ledger. | H29, H30 | run |
| [zpf_hydrogen](zpf_hydrogen) | Puthoff 1987: hydrogen's ground state as a balance of zero-point absorption and radiation; the oscillator in a simulated zero-point field and as a kronspice circuit in zero-point resistor noise. | H24 | run |
| [forward_casimir](forward_casimir) | Forward 1984 vacuum-fluctuation battery: Casimir cohesion of charged leaves under a servo supply, as a Kron (q, a) circuit with ledger and kronspice twin; one-shot gain, cycle is a net loss. | H25 | run |
| [plasmoid_rings](plasmoid_rings) | Nachamkin 1992 force-free plasmoids as toroidal current rings: his field and numbers checked, radiation lifetimes, a spheromak ring decaying in kronspice, and a collapsing electron ring driven by a kronspice bank. | H26 | run |
| [shoulders_ev](shoulders_ev) | Shoulders 1987 EVs (10¹¹-electron beads) checked against Coulomb, Gauss, the Casimir force and his sources' charge and energy; his source and an EV strike as kronspice circuits with ledgers. | H27 | run |
| [rotating_efield](rotating_efield) | Can a rotating E field make a B field? Displacement current, Rowland disk, Röntgen–Eichenwald disk; the magnetic-charge prediction labelled hypothetical. | §3f, H2 | run |
| [rate_of_change](rate_of_change) | Where higher time derivatives enter: water's Debye response as a circuit twin, and Bopp–Podolsky electrodynamics beside Maxwell. | H19 | run |
| [quantumtricity](quantumtricity) | Work in progress. Circuit theory as U(1) gauge theory on Kron's network (KVL zero holonomy, KCL the Gauss generator, Yee grid = FDTD to round-off), then SU(3) on the same graph: no circuit quantity composes non-abelianly so far; the matter sector (quark ↔ electron) is not yet modeled. | H22 | run |

## Coils, pulses and discharges

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [bial_bifilar](bial_bifilar) | First experiment: Bi–Al capacitor element (thermal, galvanic, thermopile) in parallel with a bifilar coil, evolved for Q against stored energy. | §3e, M5a | run |
| [pulsed_bifilar](pulsed_bifilar) | Coils as 20 coupled sections driven by one-sign pulses: ordinary, Tesla bifilar, conical and flat-spiral; build-up only on a harmonic; tuned pickup. | §3e | run |
| [coiled_tether](coiled_tether) | The 5-mile tether coiled (ordinary and bifilar) versus straight: motional EMF along the wire, solar-driven resonance, every joule from the Sun. | §3e, §3f | run |
| [water_lamp](water_lamp) | The neon flashing-lamp oscillator with water instead of neon: glow tube, liquid spark gap, misted air gap; tuned into the bifilar coil. | §3g, §3j, H5 | run |
| [sparkgap_diode](sparkgap_diode) | A water gap as a diode from its electrodes: aluminium valve cells (Pollak/Nodon) and polarity spark gaps; rectifiers, a 10 J pulse test, amps per volt vs frequency. | §3g, §3h | run |
| [sparkgap_diodebridge](sparkgap_diodebridge) | Water cells and spark gaps combined for six jobs: DC supply, surge protection, spark pulser, voltage doubler, oxide-film varactor with memory, frequency doubler. | §3g, §3h | run |
| [acoustic_drive](acoustic_drive) | Electrical to mechanical to acoustic: a piezo ring driving the breathing mode of a liquid sphere in H₂O and D₂O, with the power path. | §3.1, §3h | run |

## Fields, craft and space

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [tether](tether) | Electrodynamic tether in low Earth orbit: push and pull against a planetary field, momentum ledger. | §3f | run |
| [ring_towers](ring_towers) | Towers on the Earth–ionosphere ring as an LC ladder: traveling, standing and cancelling waves, synchronized pairs, DC bias. | §3e | run |
| [earth_machine](earth_machine) | The Earth as three Kron networks: Rikitake two-disk geodynamo (field reversals), the global atmospheric circuit as leaky shells (250 kV, 1 kA, Earth's charge), the Earth–ionosphere cavity as a line in colatitude (Schumann modes vs Legendre and measurement, two-height ionosphere); then one network for circuit and cavity, the dynamo's cosmic-ray shielding, paper107 drive studies. | §3a, §3f | run |
| [double_hull](double_hull) | A craft's double wall as a ring waveguide: gap voltage per section, DC bias, vacuum multipactor and Kilpatrick checks, GA for a target voltage. | §3i | run |
| [pais_cavity](pais_cavity) | Pais US 10,144,532 B2 annular cavity as a ring line by Happ two-level diakoptics; the patent's figures against physical limits. | H18 | run |

## Chemistry

| Example | Scope | CLAUDE.md | Run |
|---|---|---|---|
| [water_vibration](water_vibration) | H₂O bonds as a circuit: atoms as mass capacitors, bonds as spring inductors; AC peaks are the IR fundamentals; predicts D₂O and HDO. | doc/chemistry-plan.md | run |
| [kron_schrodinger](kron_schrodinger) | Kron's 1945 Schrödinger network on the (ρ, z) plane: H atom and the H₂⁺ bond; LC-circuit twin in kronspice AC. | §3a.1 | run |
| [kron_polyhedron](kron_polyhedron) | Kron's field network and 1962 polyhedron against his printed numbers: 1945 cavity, waveguide, Tables I–VI; frames as integer vectors recovered by the GA; the oscillatory frame in closed form. | §3a.1, H20 | run |
| [dna_vibration](dna_vibration) | The four DNA bases and the Watson–Crick pairs as circuits (AMBER OL15): absorption, pairing modes, stacks, sound renderings. | doc/chemistry-plan.md | run |

## Common conventions

- Every run keeps the energy ledger (CLAUDE.md §3c) and reports its
  residual; a surplus is a bug to find, not a result.
- Hypothesis examples (H-numbers) print the standard prediction beside
  the hypothesis, and hypothetical modes are labelled in output.
- `export/` folders hold each circuit as a SPICE netlist, a KiCad schematic
  and a results sheet comparing MNA with Kron's torn solution, written by
  `TestExport` with `KRONSPICE_EXPORT=1`.
- `run_output.txt`, where present, is the output quoted in the README.
